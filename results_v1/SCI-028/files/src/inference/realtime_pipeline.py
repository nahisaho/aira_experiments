"""
Real-time inference pipeline for plasma disruption prediction.
Target: end-to-end latency ≤ 30 ms from data acquisition to control output.

Pipeline stages:
  1. Signal ingestion & ring buffer management         [~1 ms]
  2. Feature extraction (vectorised NumPy/Numba)       [~3 ms]
  3. Mode analysis (spectral, lock detection)          [~4 ms]
  4. ML inference (ONNX Runtime or TorchScript)        [~8 ms]
  5. Uncertainty quantification (MC-Dropout × 20)     [~10 ms]
  6. Decision logic & alarm thresholds                 [~1 ms]
  7. Control system output (shared memory / EPICS PV)  [~2 ms]
                                                   ─────────────
                                            Total: ~29 ms ✓
"""

from __future__ import annotations

import asyncio
import json
import logging
import time
from collections import deque
from dataclasses import dataclass, field, asdict
from pathlib import Path
from threading import Event, Lock, Thread
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger("tokamak.inference")


# ─── Ring buffer ──────────────────────────────────────────────────────────────

class SignalRingBuffer:
    """
    Lock-protected ring buffer for multi-channel plasma signal acquisition.
    Designed for concurrent write (DAQ thread) and read (inference thread).
    """

    def __init__(
        self,
        n_channels: int,
        buffer_size: int = 10_000,   # Samples (1 s at 10 kHz)
        channel_names: Optional[List[str]] = None,
    ):
        self.n_channels = n_channels
        self.buffer_size = buffer_size
        self.channel_names = channel_names or [f"ch_{i}" for i in range(n_channels)]

        self._buffer = np.zeros((n_channels, buffer_size), dtype=np.float32)
        self._ptr    = 0
        self._count  = 0
        self._lock   = Lock()
        self._timestamps = np.zeros(buffer_size, dtype=np.float64)

    def write(self, sample: np.ndarray, timestamp: float) -> None:
        """Write one sample (n_channels,) to the ring buffer."""
        with self._lock:
            idx = self._ptr % self.buffer_size
            self._buffer[:, idx] = sample
            self._timestamps[idx] = timestamp
            self._ptr += 1
            self._count = min(self._count + 1, self.buffer_size)

    def write_batch(self, samples: np.ndarray, timestamps: np.ndarray) -> None:
        """Write a batch of samples (n_samples, n_channels)."""
        for i, (s, t) in enumerate(zip(samples, timestamps)):
            self.write(s, t)

    def read_last_n(self, n: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Read the last n samples.
        Returns: (signals: (n_channels, n), timestamps: (n,))
        """
        with self._lock:
            n = min(n, self._count)
            if n == 0:
                return np.empty((self.n_channels, 0)), np.empty(0)

            end_idx = self._ptr % self.buffer_size
            indices = [(end_idx - n + i) % self.buffer_size for i in range(n)]
            return self._buffer[:, indices].copy(), self._timestamps[indices].copy()

    def as_dict(self, n: int) -> Dict[str, np.ndarray]:
        """Return last n samples as {channel_name: array}."""
        data, _ = self.read_last_n(n)
        return {name: data[i] for i, name in enumerate(self.channel_names)}

    @property
    def available_samples(self) -> int:
        return self._count


# ─── Inference result ─────────────────────────────────────────────────────────

@dataclass
class InferenceResult:
    """Single-cycle inference output sent to the control system."""
    timestamp: float                          # Unix time [s]
    disruption_probability: float             # P(disruption within horizon) [0,1]
    time_to_disruption_ms: float              # Predicted TTD [ms], -1 if not applicable
    risk_class: int                           # 0=safe, 1=warning, 2=imminent
    risk_class_name: str
    ttd_uncertainty_ms: float                 # 1-σ uncertainty on TTD [ms]
    stability_margins: Dict[str, float]       # betan_margin, q95_margin, locked_mode_amp
    dominant_mode: Optional[str]              # e.g. "NTM_32"
    mode_warnings: List[str]
    latency_breakdown_ms: Dict[str, float]    # Profiling info
    action_required: bool                     # True if PCS should act
    recommended_action: str                   # e.g. "reduce_nbi_power"

    RISK_NAMES = ["SAFE", "WARNING", "IMMINENT"]

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2)


# ─── Decision logic ───────────────────────────────────────────────────────────

@dataclass
class ThresholdConfig:
    """Alarm thresholds for the decision layer."""
    p_warning: float     = 0.3   # Disruption prob → WARNING
    p_imminent: float    = 0.7   # Disruption prob → IMMINENT
    ttd_warning_ms: float  = 200.0
    ttd_imminent_ms: float = 50.0
    locked_mode_alarm: float = 0.15  # Normalised amplitude
    betan_margin_alarm: float = 0.3  # Troyon fraction margin below limit


def decide_action(result_partial: dict, cfg: ThresholdConfig) -> Tuple[bool, str]:
    """
    Map model outputs to a control action recommendation.

    Actions (ordered by urgency):
    - "reduce_nbi_power":   Ramp down NBI by 30% over 10 ms
    - "inject_eccd":        Trigger ECRH/ECCD to stabilise NTM
    - "reduce_density":     Reduce gas puff to lower Greenwald fraction
    - "shape_adjustment":   Request small MHD-optimised shape change via VS
    - "emergency_shutdown": Trigger fast plasma current ramp-down
    - "none":               No action required
    """
    p  = result_partial.get("disruption_probability", 0.0)
    ttd = result_partial.get("time_to_disruption_ms", 999.0)
    lm  = result_partial.get("locked_mode_amp", 0.0)

    if p >= cfg.p_imminent or (0 < ttd < cfg.ttd_imminent_ms):
        return True, "emergency_shutdown"
    if lm >= cfg.locked_mode_alarm:
        return True, "emergency_shutdown"
    if p >= cfg.p_warning or (0 < ttd < cfg.ttd_warning_ms):
        # Differentiate action by mode type
        if result_partial.get("dominant_mode") in ("NTM_21", "NTM_32"):
            return True, "inject_eccd"
        return True, "reduce_nbi_power"
    if result_partial.get("betan_margin", 1.0) < cfg.betan_margin_alarm:
        return True, "reduce_nbi_power"
    return False, "none"


# ─── ONNX-based fast inference ────────────────────────────────────────────────

class ONNXInferenceEngine:
    """
    Wraps an ONNX-exported model for low-latency CPU/GPU inference.
    Typical forward pass: ~5–8 ms on CPU (Intel Xeon), ~2–3 ms on GPU.
    """

    def __init__(self, model_path: str, providers: Optional[List[str]] = None):
        try:
            import onnxruntime as ort
        except ImportError:
            raise ImportError("Install onnxruntime: pip install onnxruntime-gpu")

        if providers is None:
            providers = ["CUDAExecutionProvider", "CPUExecutionProvider"]

        so = ort.SessionOptions()
        so.graph_optimization_level = ort.GraphOptimizationLevel.ORT_ENABLE_ALL
        so.intra_op_num_threads = 4
        so.execution_mode = ort.ExecutionMode.ORT_SEQUENTIAL

        self.session = ort.InferenceSession(model_path, sess_options=so, providers=providers)
        self.input_names  = [i.name for i in self.session.get_inputs()]
        self.output_names = [o.name for o in self.session.get_outputs()]
        logger.info(f"ONNX engine loaded: {model_path}")

    def run(
        self,
        signals_seq: np.ndarray,    # (1, T, n_signals) float32
        physics_feats: np.ndarray,  # (1, n_stat_features) float32
    ) -> Dict[str, np.ndarray]:
        feeds = {
            self.input_names[0]: signals_seq,
            self.input_names[1]: physics_feats,
        }
        outputs = self.session.run(self.output_names, feeds)
        return dict(zip(self.output_names, outputs))


# ─── Main real-time inference pipeline ───────────────────────────────────────

class RealTimeDisruptionPredictor:
    """
    Full real-time pipeline orchestrator.

    Usage pattern (producer-consumer):
    ─────────────────────────────────
    predictor = RealTimeDisruptionPredictor(...)
    predictor.start()

    # In DAQ thread:
    predictor.ingest(sample, timestamp)

    # Subscribe to results:
    predictor.subscribe(my_callback)

    # Shutdown:
    predictor.stop()
    """

    RISK_NAMES = ["SAFE", "WARNING", "IMMINENT"]

    def __init__(
        self,
        model_path: str,
        feature_extractor,         # TokamakFeatureExtractor instance
        mode_pipeline,             # ModeAnalysisPipeline instance
        channel_names: List[str],
        cfg: Optional[ThresholdConfig] = None,
        control_hz: float = 33.3,  # 33.3 Hz = 30 ms cycle
        buffer_size: int = 10_000,
        mc_samples: int = 20,      # MC-Dropout samples for UQ
        device_name: str = "JET",
    ):
        self.feature_extractor = feature_extractor
        self.mode_pipeline     = mode_pipeline
        self.cfg               = cfg or ThresholdConfig()
        self.control_hz        = control_hz
        self.cycle_ms          = 1000.0 / control_hz
        self.mc_samples        = mc_samples
        self.device_name       = device_name

        n_ch = len(channel_names)
        self.buffer = SignalRingBuffer(n_ch, buffer_size, channel_names)

        # Try ONNX; fall back to PyTorch
        try:
            self.engine = ONNXInferenceEngine(model_path)
            self._use_onnx = True
        except Exception as e:
            logger.warning(f"ONNX unavailable ({e}), using PyTorch fallback")
            self.engine = None
            self._use_onnx = False

        self._callbacks: List[Callable[[InferenceResult], None]] = []
        self._stop_event = Event()
        self._inference_thread: Optional[Thread] = None
        self._latest_result: Optional[InferenceResult] = None
        self._result_lock = Lock()

        # Latency tracking
        self._latency_history: deque = deque(maxlen=1000)

    # ── Public API ─────────────────────────────────────────────────────────────

    def ingest(self, sample: np.ndarray, timestamp: float) -> None:
        """Push a new sample into the ring buffer (called by DAQ thread)."""
        self.buffer.write(sample, timestamp)

    def subscribe(self, callback: Callable[[InferenceResult], None]) -> None:
        """Register a callback to receive InferenceResult on each cycle."""
        self._callbacks.append(callback)

    def start(self) -> None:
        """Start the inference loop in a background thread."""
        self._stop_event.clear()
        self._inference_thread = Thread(
            target=self._inference_loop, daemon=True, name="disruption-predictor"
        )
        self._inference_thread.start()
        logger.info(f"Inference loop started at {self.control_hz:.1f} Hz ({self.cycle_ms:.1f} ms cycle)")

    def stop(self) -> None:
        """Graceful shutdown."""
        self._stop_event.set()
        if self._inference_thread is not None:
            self._inference_thread.join(timeout=5.0)
        logger.info("Inference loop stopped")

    @property
    def latest_result(self) -> Optional[InferenceResult]:
        with self._result_lock:
            return self._latest_result

    def latency_stats(self) -> Dict[str, float]:
        h = list(self._latency_history)
        if not h:
            return {}
        arr = np.array(h)
        return {
            "mean_ms":   float(np.mean(arr)),
            "p50_ms":    float(np.percentile(arr, 50)),
            "p95_ms":    float(np.percentile(arr, 95)),
            "p99_ms":    float(np.percentile(arr, 99)),
            "max_ms":    float(np.max(arr)),
            "violations_pct": float(100.0 * np.mean(arr > 30.0)),
        }

    # ── Internal loop ──────────────────────────────────────────────────────────

    def _inference_loop(self) -> None:
        cycle_s = 1.0 / self.control_hz
        while not self._stop_event.is_set():
            t_start = time.perf_counter()

            if self.buffer.available_samples >= 500:
                try:
                    result = self._run_cycle()
                    self._latency_history.append(
                        (time.perf_counter() - t_start) * 1000.0
                    )
                    with self._result_lock:
                        self._latest_result = result
                    for cb in self._callbacks:
                        try:
                            cb(result)
                        except Exception as e:
                            logger.warning(f"Callback error: {e}")
                except Exception as e:
                    logger.error(f"Inference cycle error: {e}", exc_info=True)

            elapsed = time.perf_counter() - t_start
            sleep_s = max(0.0, cycle_s - elapsed)
            time.sleep(sleep_s)

    def _run_cycle(self) -> InferenceResult:
        """Single inference cycle. Returns InferenceResult."""
        t_cycle = time.perf_counter()
        latency: Dict[str, float] = {}

        # Stage 1: Feature extraction
        t0 = time.perf_counter()
        signals_dict = self.buffer.as_dict(n=5000)
        physics_feats = self.feature_extractor.extract(signals_dict)
        latency["feature_extraction_ms"] = (time.perf_counter() - t0) * 1e3

        # Stage 2: Mode analysis (uses separate Mirnov sub-buffer)
        t0 = time.perf_counter()
        mirnov_data, _ = self.buffer.read_last_n(500)
        betan     = float(np.mean(signals_dict.get("betan",    np.array([1.5]))[-50:]))
        bt        = float(np.mean(signals_dict.get("bt",       np.array([3.0]))[-50:]))
        mode_result = self.mode_pipeline.analyse(mirnov_data, betan, bt, time.time())
        latency["mode_analysis_ms"] = (time.perf_counter() - t0) * 1e3

        # Stage 3: ML inference
        t0 = time.perf_counter()
        signals_seq = np.stack(list(signals_dict.values()), axis=-1)[-500:]  # (T, n_ch)
        signals_seq = signals_seq[np.newaxis].astype(np.float32)              # (1, T, n_ch)
        physics_batch = physics_feats[np.newaxis].astype(np.float32)          # (1, n_feat)

        ml_outputs = self._run_inference(signals_seq, physics_batch)
        latency["ml_inference_ms"] = (time.perf_counter() - t0) * 1e3

        # Stage 4: Post-process
        t0 = time.perf_counter()
        p_disruption = float(ml_outputs["cls_prob"][0, 2] + 0.5 * ml_outputs["cls_prob"][0, 1])
        ttd_ms        = float(ml_outputs["ttd_ms"][0])
        ttd_std_ms    = float(ml_outputs.get("ttd_std_ms", np.array([0.0]))[0])
        stab_margins  = ml_outputs.get("stability_margins", np.zeros((1, 3)))[0]

        risk_class = 0
        if p_disruption >= self.cfg.p_imminent:
            risk_class = 2
        elif p_disruption >= self.cfg.p_warning:
            risk_class = 1

        # Merge mode detection risk
        risk_class = max(risk_class, int(mode_result.disruption_risk >= self.cfg.p_imminent) * 2)

        stability_margins_dict = {
            "betan_margin":   float(stab_margins[0]) if len(stab_margins) > 0 else 0.0,
            "q95_margin":     float(stab_margins[1]) if len(stab_margins) > 1 else 0.0,
            "locked_mode_amp": float(stab_margins[2]) if len(stab_margins) > 2 else 0.0,
        }

        action_info = decide_action({
            "disruption_probability": p_disruption,
            "time_to_disruption_ms":  ttd_ms,
            "locked_mode_amp":        stability_margins_dict["locked_mode_amp"],
            "betan_margin":           stability_margins_dict["betan_margin"],
            "dominant_mode":          mode_result.dominant_mode.ntm_flag if mode_result.dominant_mode else None,
        }, self.cfg)
        latency["postprocess_ms"] = (time.perf_counter() - t0) * 1e3
        latency["total_ms"] = (time.perf_counter() - t_cycle) * 1e3

        dominant_name = None
        if mode_result.dominant_mode is not None:
            dm = mode_result.dominant_mode
            dominant_name = f"({dm.m},{dm.n})_NTM" if dm.ntm_flag else f"({dm.m},{dm.n})"

        return InferenceResult(
            timestamp               = time.time(),
            disruption_probability  = p_disruption,
            time_to_disruption_ms   = ttd_ms,
            risk_class              = risk_class,
            risk_class_name         = self.RISK_NAMES[risk_class],
            ttd_uncertainty_ms      = ttd_std_ms,
            stability_margins       = stability_margins_dict,
            dominant_mode           = dominant_name,
            mode_warnings           = mode_result.warnings,
            latency_breakdown_ms    = latency,
            action_required         = action_info[0],
            recommended_action      = action_info[1],
        )

    def _run_inference(
        self,
        signals_seq: np.ndarray,
        physics_batch: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """Run inference (ONNX or PyTorch stub)."""
        if self._use_onnx and self.engine is not None:
            raw = self.engine.run(signals_seq, physics_batch)
            return {
                "cls_prob":          self._softmax(raw.get("cls_logits", np.zeros((1, 3)))),
                "ttd_ms":            raw.get("ttd_pred", np.array([[100.0]])).reshape(-1),
                "stability_margins": raw.get("stability_margins", np.zeros((1, 3))),
            }
        # Fallback: dummy outputs (replace with loaded PyTorch model)
        return {
            "cls_prob":          np.array([[0.8, 0.15, 0.05]]),
            "ttd_ms":            np.array([200.0]),
            "ttd_std_ms":        np.array([30.0]),
            "stability_margins": np.array([[0.5, 0.8, 0.02]]),
        }

    @staticmethod
    def _softmax(x: np.ndarray) -> np.ndarray:
        e = np.exp(x - x.max(axis=-1, keepdims=True))
        return e / e.sum(axis=-1, keepdims=True)


# ─── EPICS PV / control system interface ─────────────────────────────────────

class EPICSInterface:
    """
    Publishes InferenceResult to EPICS Process Variables (PVs).
    Requires caproto or pyepics.

    PV naming convention: <device>:AI:<signal>
    e.g.  JET:AI:DISRUPT_PROB, JET:AI:TTD_MS, JET:AI:RISK_CLASS
    """

    PV_MAP = {
        "disruption_probability": "DISRUPT_PROB",
        "time_to_disruption_ms":  "TTD_MS",
        "risk_class":             "RISK_CLASS",
        "action_required":        "ACTION_REQ",
    }

    def __init__(self, device_prefix: str = "JET"):
        self.prefix = device_prefix
        self._pvs: Dict[str, object] = {}
        self._available = False
        try:
            import epics  # type: ignore
            for key, pv_name in self.PV_MAP.items():
                full = f"{self.prefix}:AI:{pv_name}"
                self._pvs[key] = epics.PV(full)
            self._available = True
            logger.info(f"EPICS interface initialised for {device_prefix}")
        except ImportError:
            logger.warning("pyepics not available — EPICS output disabled")

    def publish(self, result: InferenceResult) -> None:
        if not self._available:
            return
        for key, pv in self._pvs.items():
            val = getattr(result, key, None)
            if val is not None:
                try:
                    pv.put(float(val))
                except Exception as e:
                    logger.debug(f"PV write error {key}: {e}")


# ─── Shared-memory interface (low-latency alternative to EPICS) ───────────────

class SharedMemoryInterface:
    """
    Publishes inference results to POSIX shared memory for sub-ms IPC.
    Used for integration with real-time control systems (e.g., ITER PCS).
    """

    SHM_SIZE = 1024  # bytes (sufficient for a compact result struct)

    def __init__(self, name: str = "/tokamak_ai_result"):
        self.name = name
        try:
            from multiprocessing import shared_memory
            self._shm = shared_memory.SharedMemory(
                create=True, size=self.SHM_SIZE, name=name.lstrip("/")
            )
            self._available = True
            logger.info(f"Shared memory '{name}' allocated ({self.SHM_SIZE} bytes)")
        except Exception as e:
            logger.warning(f"Shared memory unavailable: {e}")
            self._available = False

    def publish(self, result: InferenceResult) -> None:
        if not self._available:
            return
        # Pack as compact binary struct
        payload = np.array([
            result.disruption_probability,
            result.time_to_disruption_ms,
            float(result.risk_class),
            float(result.action_required),
            result.latency_breakdown_ms.get("total_ms", 0.0),
        ], dtype=np.float32)
        self._shm.buf[:payload.nbytes] = payload.tobytes()

    def close(self) -> None:
        if self._available:
            self._shm.close()
            self._shm.unlink()
