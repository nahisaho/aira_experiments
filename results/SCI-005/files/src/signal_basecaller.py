"""Signal-level basecalling primitives for the DeepSV-LR pipeline.

This module implements a lightweight, NumPy-based representation of a
bidirectional GRU (BiGRU) signal basecaller. The implementation is intended for
algorithmic prototyping, deterministic testing, and documentation of the model
architecture without requiring a heavyweight deep-learning runtime.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np
from numpy.typing import ArrayLike, NDArray


@dataclass(frozen=True)
class BeamState:
    """Internal state used during CTC beam-search decoding."""

    prefix: str
    log_prob_blank: float
    log_prob_non_blank: float

    @property
    def score(self) -> float:
        return _logsumexp([self.log_prob_blank, self.log_prob_non_blank])


class SignalBasecaller:
    r"""RNN-based signal-level basecaller using a 5-layer BiGRU encoder.

    Mathematical formulation
    ------------------------
    Preprocessing uses robust median absolute deviation (MAD) normalization:

    .. math::
        \tilde{x}_t = \frac{x_t - \mathrm{median}(x)}{1.4826 \cdot
        \mathrm{median}(|x - \mathrm{median}(x)|) + \varepsilon}

    Each GRU direction computes gates for time-step :math:`t`:

    .. math::
        z_t = \sigma(W_z x_t + U_z h_{t-1} + b_z)

    .. math::
        r_t = \sigma(W_r x_t + U_r h_{t-1} + b_r)

    .. math::
        \hat{h}_t = \tanh(W_h x_t + U_h (r_t \odot h_{t-1}) + b_h)

    .. math::
        h_t = (1 - z_t) \odot h_{t-1} + z_t \odot \hat{h}_t

    The forward and backward hidden states are concatenated and projected into
    CTC logits. Decoding finds the most probable output string under the CTC
    objective:

    .. math::
        y^* = \arg\max_y \sum_{\pi \in \mathcal{B}^{-1}(y)}
        \prod_t p(\pi_t | x)

    where :math:`\mathcal{B}` collapses repeated labels and removes blanks.
    """

    def __init__(
        self,
        input_size: int = 1,
        hidden_size: int = 256,
        num_layers: int = 5,
        alphabet: Sequence[str] = ("-", "A", "C", "G", "T"),
        seed: int = 7,
    ) -> None:
        self.input_size = input_size
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.alphabet = tuple(alphabet)
        self.blank_index = 0
        self.rng = np.random.default_rng(seed)
        self.layers = self._initialize_layers()
        self.output_weight = self.rng.normal(
            0.0,
            0.02,
            size=(len(self.alphabet), hidden_size * 2),
        )
        self.output_bias = np.zeros(len(self.alphabet), dtype=np.float64)

    def _initialize_layers(self) -> List[Dict[str, Dict[str, NDArray[np.float64]]]]:
        layers: List[Dict[str, Dict[str, NDArray[np.float64]]]] = []
        current_dim = self.input_size
        for _ in range(self.num_layers):
            layer: Dict[str, Dict[str, NDArray[np.float64]]] = {}
            for direction in ("forward", "backward"):
                scale = 1.0 / np.sqrt(max(current_dim, 1))
                recurrent_scale = 1.0 / np.sqrt(self.hidden_size)
                layer[direction] = {
                    "W_z": self.rng.normal(0.0, scale, size=(self.hidden_size, current_dim)),
                    "U_z": self.rng.normal(0.0, recurrent_scale, size=(self.hidden_size, self.hidden_size)),
                    "b_z": np.zeros(self.hidden_size, dtype=np.float64),
                    "W_r": self.rng.normal(0.0, scale, size=(self.hidden_size, current_dim)),
                    "U_r": self.rng.normal(0.0, recurrent_scale, size=(self.hidden_size, self.hidden_size)),
                    "b_r": np.zeros(self.hidden_size, dtype=np.float64),
                    "W_h": self.rng.normal(0.0, scale, size=(self.hidden_size, current_dim)),
                    "U_h": self.rng.normal(0.0, recurrent_scale, size=(self.hidden_size, self.hidden_size)),
                    "b_h": np.zeros(self.hidden_size, dtype=np.float64),
                }
            layers.append(layer)
            current_dim = self.hidden_size * 2
        return layers

    def preprocess_signal(
        self,
        signal: ArrayLike,
        smoothing_window: int = 5,
        epsilon: float = 1e-6,
    ) -> NDArray[np.float64]:
        r"""Normalize and smooth raw ionic current.

        Parameters
        ----------
        signal:
            Raw nanopore current samples.
        smoothing_window:
            Width of the moving-average denoiser.
        epsilon:
            Stabilizer added to the MAD term.

        Returns
        -------
        numpy.ndarray
            A two-dimensional array of shape ``(timesteps, 1)``.

        Notes
        -----
        The robust scale estimate is:

        .. math::
            s = 1.4826 \cdot \mathrm{median}(|x - \mathrm{median}(x)|) + \varepsilon

        and the normalized trace is smoothed with a uniform kernel to suppress
        high-frequency pore noise while preserving event boundaries.
        """

        values = np.asarray(signal, dtype=np.float64).reshape(-1)
        if values.size == 0:
            return np.empty((0, self.input_size), dtype=np.float64)

        median = float(np.median(values))
        mad = float(np.median(np.abs(values - median)))
        scale = 1.4826 * mad
        if scale < epsilon:
            scale = float(np.std(values)) + epsilon
        normalized = (values - median) / (scale + epsilon)

        if smoothing_window > 1:
            kernel = np.ones(smoothing_window, dtype=np.float64) / smoothing_window
            padded = np.pad(normalized, (smoothing_window // 2,), mode="edge")
            normalized = np.convolve(padded, kernel, mode="valid")[: values.size]

        return normalized.reshape(-1, self.input_size)

    def forward_pass(self, features: ArrayLike) -> NDArray[np.float64]:
        """Compute CTC logits for a preprocessed signal chunk.

        Parameters
        ----------
        features:
            Array-like input with shape ``(timesteps, input_dim)``.

        Returns
        -------
        numpy.ndarray
            Logits with shape ``(timesteps, alphabet_size)``.
        """

        x = np.asarray(features, dtype=np.float64)
        if x.ndim == 1:
            x = x.reshape(-1, self.input_size)
        if x.size == 0:
            return np.empty((0, len(self.alphabet)), dtype=np.float64)

        current = x
        for layer in self.layers:
            forward_hidden = self._run_gru_direction(current, layer["forward"], reverse=False)
            backward_hidden = self._run_gru_direction(current, layer["backward"], reverse=True)
            current = np.concatenate([forward_hidden, backward_hidden], axis=1)

        return current @ self.output_weight.T + self.output_bias

    def ctc_decode(self, logits: ArrayLike, beam_width: int = 5) -> str:
        r"""Decode CTC logits with prefix beam search.

        Parameters
        ----------
        logits:
            Frame-level network output before softmax.
        beam_width:
            Number of active prefixes retained after each frame.

        Returns
        -------
        str
            Decoded nucleotide sequence.

        Notes
        -----
        The decoder tracks prefix probabilities for paths ending in blank and
        non-blank states. After each step it keeps the top ``beam_width``
        prefixes ranked by

        .. math::
            \log p(\ell) = \log\left(p_b(\ell) + p_{nb}(\ell)\right).
        """

        scores = np.asarray(logits, dtype=np.float64)
        if scores.size == 0:
            return ""

        log_probs = scores - _logsumexp(scores, axis=1, keepdims=True)
        beams: Dict[str, BeamState] = {"": BeamState(prefix="", log_prob_blank=0.0, log_prob_non_blank=-np.inf)}

        for frame in log_probs:
            next_beams: Dict[str, BeamState] = {}
            for prefix, state in beams.items():
                total = state.score
                blank_log_prob = total + frame[self.blank_index]
                self._update_beam(next_beams, prefix, blank_log_prob, is_blank=True)

                for token_index, token in enumerate(self.alphabet[1:], start=1):
                    prev_char = prefix[-1] if prefix else None
                    if token == prev_char:
                        emit = state.log_prob_blank + frame[token_index]
                        self._update_beam(next_beams, prefix, emit, is_blank=False)
                    else:
                        emit = total + frame[token_index]
                        self._update_beam(next_beams, prefix + token, emit, is_blank=False)

            top_prefixes = sorted(next_beams.values(), key=lambda beam: beam.score, reverse=True)[:beam_width]
            beams = {beam.prefix: beam for beam in top_prefixes}

        best = max(beams.values(), key=lambda beam: beam.score)
        return best.prefix

    def quality_consensus(
        self,
        sequences: Sequence[str],
        phred_scores: Sequence[Sequence[float]],
    ) -> Tuple[str, List[float]]:
        """Build a simple quality-weighted consensus across multiple basecalls.

        Parameters
        ----------
        sequences:
            Candidate basecalled sequences.
        phred_scores:
            Per-base Phred quality scores aligned to ``sequences``.

        Returns
        -------
        tuple[str, list[float]]
            Consensus sequence and consensus Phred qualities.
        """

        if not sequences:
            return "", []
        if len(sequences) != len(phred_scores):
            raise ValueError("sequences and phred_scores must have the same length")

        max_len = max(len(sequence) for sequence in sequences)
        consensus_bases: List[str] = []
        consensus_qualities: List[float] = []

        for index in range(max_len):
            weights: Dict[str, float] = {}
            total_weight = 0.0
            for sequence, qualities in zip(sequences, phred_scores):
                if index >= len(sequence):
                    continue
                base = sequence[index]
                quality = float(qualities[index] if index < len(qualities) else 10.0)
                support = 1.0 - 10.0 ** (-quality / 10.0)
                weights[base] = weights.get(base, 0.0) + support
                total_weight += support

            if not weights:
                continue
            base, weight = max(weights.items(), key=lambda item: item[1])
            posterior = weight / max(total_weight, 1e-8)
            posterior = min(max(posterior, 1e-6), 1.0 - 1e-6)
            consensus_bases.append(base)
            consensus_qualities.append(-10.0 * np.log10(1.0 - posterior))

        return "".join(consensus_bases), consensus_qualities

    def _run_gru_direction(
        self,
        inputs: NDArray[np.float64],
        weights: Dict[str, NDArray[np.float64]],
        reverse: bool,
    ) -> NDArray[np.float64]:
        ordered = inputs[::-1] if reverse else inputs
        hidden = np.zeros(self.hidden_size, dtype=np.float64)
        outputs: List[NDArray[np.float64]] = []

        for timestep in ordered:
            z_t = _sigmoid(weights["W_z"] @ timestep + weights["U_z"] @ hidden + weights["b_z"])
            r_t = _sigmoid(weights["W_r"] @ timestep + weights["U_r"] @ hidden + weights["b_r"])
            proposal = np.tanh(weights["W_h"] @ timestep + weights["U_h"] @ (r_t * hidden) + weights["b_h"])
            hidden = (1.0 - z_t) * hidden + z_t * proposal
            outputs.append(hidden.copy())

        stacked = np.vstack(outputs)
        return stacked[::-1] if reverse else stacked

    def _update_beam(
        self,
        beams: Dict[str, BeamState],
        prefix: str,
        log_prob: float,
        is_blank: bool,
    ) -> None:
        current = beams.get(prefix, BeamState(prefix=prefix, log_prob_blank=-np.inf, log_prob_non_blank=-np.inf))
        if is_blank:
            beams[prefix] = BeamState(
                prefix=prefix,
                log_prob_blank=_logsumexp([current.log_prob_blank, log_prob]),
                log_prob_non_blank=current.log_prob_non_blank,
            )
        else:
            beams[prefix] = BeamState(
                prefix=prefix,
                log_prob_blank=current.log_prob_blank,
                log_prob_non_blank=_logsumexp([current.log_prob_non_blank, log_prob]),
            )


def _sigmoid(values: NDArray[np.float64]) -> NDArray[np.float64]:
    return 1.0 / (1.0 + np.exp(-np.clip(values, -60.0, 60.0)))


def _logsumexp(
    values: ArrayLike,
    axis: int | None = None,
    keepdims: bool = False,
) -> NDArray[np.float64] | float:
    array = np.asarray(values)
    if array.dtype == object:
        flattened = [np.asarray(value, dtype=np.float64) for value in values]  # type: ignore[arg-type]
        array = np.stack(flattened, axis=0)
    array = np.asarray(array, dtype=np.float64)
    if array.size == 0:
        return -np.inf
    maximum = np.max(array, axis=axis, keepdims=True)
    finite_maximum = np.where(np.isfinite(maximum), maximum, 0.0)
    stable = np.exp(array - finite_maximum)
    stable = np.where(np.isfinite(array), stable, 0.0)
    summed = np.sum(stable, axis=axis, keepdims=True)
    result = np.where(np.isfinite(maximum), maximum + np.log(np.maximum(summed, 1e-300)), -np.inf)
    if not keepdims and axis is not None:
        result = np.squeeze(result, axis=axis)
    if not keepdims and axis is None:
        result = np.squeeze(result)
    if np.isscalar(result) or getattr(result, "ndim", 0) == 0:
        return float(np.asarray(result))
    return np.asarray(result, dtype=np.float64)


__all__ = ["BeamState", "SignalBasecaller"]
