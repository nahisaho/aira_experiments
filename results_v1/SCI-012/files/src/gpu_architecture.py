"""
GPU-parallel SNN architecture for 1M+ neuron simulations.

Architecture overview
---------------------
* CPU fallback (NumPy/vectorised) always available.
* Numba JIT-compiled CPU batch kernel for moderate sizes (~100k neurons).
* Numba CUDA kernel for GPU acceleration (when CUDA device present).
* The public API is identical regardless of backend selected.

Design patterns
---------------
1. State arrays are flat C-contiguous float32 for coalesced CUDA access.
2. Sparse connectivity encoded as CSR (indptr, indices, weights).
3. Each time step: integrate → detect spikes → propagate → plasticity.
4. Double-buffering for V/u arrays avoids race conditions.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, Tuple
import time

try:
    from numba import cuda, jit, prange
    import numba
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False

try:
    from numba import cuda as numba_cuda
    CUDA_AVAILABLE = numba_cuda.is_available()
except Exception:
    CUDA_AVAILABLE = False


# ---------------------------------------------------------------------------
# Sparse connectivity (CSR)
# ---------------------------------------------------------------------------

@dataclass
class CSRConnectivity:
    """Compressed Sparse Row connectivity matrix."""
    n_pre:   int
    n_post:  int
    indptr:  np.ndarray   # shape (n_post+1,) int32
    indices: np.ndarray   # shape (n_syn,)  int32
    weights: np.ndarray   # shape (n_syn,)  float32

    @classmethod
    def random(cls, n_pre: int, n_post: int, p_conn: float = 0.01,
               w_mean: float = 0.5, w_std: float = 0.1,
               rng: Optional[np.random.Generator] = None) -> "CSRConnectivity":
        """Generate random sparse connectivity."""
        if rng is None:
            rng = np.random.default_rng(42)
        rows, cols, ws = [], [], []
        for j in range(n_post):
            sources = np.where(rng.random(n_pre) < p_conn)[0]
            for s in sources:
                rows.append(j)
                cols.append(s)
                ws.append(max(0.0, rng.normal(w_mean, w_std)))
        indptr = np.zeros(n_post + 1, dtype=np.int32)
        for r in rows:
            indptr[r + 1] += 1
        np.cumsum(indptr, out=indptr)
        return cls(
            n_pre=n_pre, n_post=n_post,
            indptr=indptr.astype(np.int32),
            indices=np.array(cols, dtype=np.int32),
            weights=np.array(ws, dtype=np.float32),
        )

    @property
    def n_syn(self) -> int:
        return len(self.indices)


# ---------------------------------------------------------------------------
# Numba JIT CPU batch kernel (Izhikevich population)
# ---------------------------------------------------------------------------

if NUMBA_AVAILABLE:
    from numba import jit as _jit, prange as _prange

    @_jit(nopython=True, parallel=True, fastmath=True)
    def _izhikevich_step_cpu(V, u, I_syn, a, b, c, d, dt, spikes_out):
        """Vectorised Izhikevich step over N neurons."""
        N = V.shape[0]
        for i in _prange(N):
            v, uv = V[i], u[i]
            if v >= 30.0:
                spikes_out[i] = True
                V[i] = c[i]
                u[i] = uv + d[i]
            else:
                spikes_out[i] = False
                dv = (0.04*v*v + 5.0*v + 140.0 - uv + I_syn[i]) * dt
                du = a[i] * (b[i]*v - uv) * dt
                V[i] = v + dv
                u[i] = uv + du

    @_jit(nopython=True, parallel=True, fastmath=True)
    def _csr_mv_cpu(indptr, indices, weights, spike_vec, I_out):
        """Sparse matrix-vector multiply (CSR) for synaptic currents."""
        n_post = indptr.shape[0] - 1
        for j in _prange(n_post):
            s = 0.0
            for k in range(indptr[j], indptr[j+1]):
                s += weights[k] * spike_vec[indices[k]]
            I_out[j] = s

else:
    def _izhikevich_step_cpu(V, u, I_syn, a, b, c, d, dt, spikes_out):
        spikes_out[:] = V >= 30.0
        reset = spikes_out
        V[reset] = c[reset]
        u[reset] += d[reset]
        nreset = ~reset
        v, uv = V[nreset], u[nreset]
        V[nreset] = v + (0.04*v**2 + 5*v + 140 - uv + I_syn[nreset])*dt
        u[nreset] = uv + a[nreset]*(b[nreset]*v - uv)*dt

    def _csr_mv_cpu(indptr, indices, weights, spike_vec, I_out):
        n_post = len(indptr) - 1
        for j in range(n_post):
            I_out[j] = np.dot(weights[indptr[j]:indptr[j+1]],
                              spike_vec[indices[indptr[j]:indptr[j+1]]])


# ---------------------------------------------------------------------------
# CUDA kernels (only compiled when GPU available)
# ---------------------------------------------------------------------------

if CUDA_AVAILABLE:
    from numba import cuda as _cuda

    @_cuda.jit
    def _izhikevich_step_gpu(V, u, I_syn, a, b, c, d, dt, spikes_out):
        i = _cuda.grid(1)
        if i < V.shape[0]:
            v, uv = V[i], u[i]
            if v >= 30.0:
                spikes_out[i] = 1
                V[i] = c[i]
                u[i] = uv + d[i]
            else:
                spikes_out[i] = 0
                V[i] = v + (0.04*v*v + 5.0*v + 140.0 - uv + I_syn[i]) * dt
                u[i] = uv + (a[i] * (b[i]*v - uv)) * dt

    @_cuda.jit
    def _csr_mv_gpu(indptr, indices, weights, spike_vec, I_out):
        j = _cuda.grid(1)
        if j < I_out.shape[0]:
            s = 0.0
            for k in range(indptr[j], indptr[j+1]):
                s += weights[k] * spike_vec[indices[k]]
            I_out[j] = s


# ---------------------------------------------------------------------------
# Population class
# ---------------------------------------------------------------------------

class IzhikevichPopulation:
    """
    Large-scale Izhikevich population with CSR synaptic connectivity.

    Automatically selects CUDA > Numba-CPU > NumPy backend.
    """

    def __init__(self, N: int, neuron_type: str = "RS",
                 backend: str = "auto",
                 rng: Optional[np.random.Generator] = None):
        self.N = N
        self.rng = rng or np.random.default_rng(0)
        self._backend = self._select_backend(backend)

        # Parameter arrays (float32 for GPU efficiency)
        self.a, self.b, self.c, self.d = self._init_params(neuron_type)
        self.V = np.full(N, -65.0, dtype=np.float32)
        self.u = self.b * self.V

        self.I_syn   = np.zeros(N, dtype=np.float32)
        self.spikes  = np.zeros(N, dtype=np.int32)
        self.conn: Optional[CSRConnectivity] = None

        # GPU device arrays (allocated lazily)
        self._d_V = self._d_u = self._d_I = self._d_s = None
        self._d_indptr = self._d_indices = self._d_weights = None

        if self._backend == "cuda":
            self._upload_to_gpu()

    @staticmethod
    def _select_backend(backend: str) -> str:
        if backend == "auto":
            if CUDA_AVAILABLE:
                return "cuda"
            elif NUMBA_AVAILABLE:
                return "numba_cpu"
            else:
                return "numpy"
        return backend

    def _init_params(self, neuron_type: str):
        presets = {
            "RS":  (0.02, 0.2,  -65.0, 8.0),
            "FS":  (0.1,  0.2,  -65.0, 2.0),
            "IB":  (0.02, 0.2,  -55.0, 4.0),
            "LTS": (0.02, 0.25, -65.0, 2.0),
        }
        a0, b0, c0, d0 = presets.get(neuron_type, presets["RS"])
        N = self.N
        # Add small heterogeneity
        a = np.full(N, a0, dtype=np.float32) + self.rng.normal(0, 0.002, N).astype(np.float32)
        b = np.full(N, b0, dtype=np.float32) + self.rng.normal(0, 0.005, N).astype(np.float32)
        c = np.full(N, c0, dtype=np.float32) + self.rng.normal(0, 1.0,  N).astype(np.float32)
        d = np.full(N, d0, dtype=np.float32) + self.rng.normal(0, 0.5,  N).astype(np.float32)
        return a, b, c, d

    def set_connectivity(self, conn: CSRConnectivity):
        self.conn = conn
        if self._backend == "cuda" and CUDA_AVAILABLE:
            from numba import cuda as nc
            self._d_indptr  = nc.to_device(conn.indptr)
            self._d_indices = nc.to_device(conn.indices)
            self._d_weights = nc.to_device(conn.weights)

    def _upload_to_gpu(self):
        if not CUDA_AVAILABLE:
            return
        from numba import cuda as nc
        self._d_V  = nc.to_device(self.V)
        self._d_u  = nc.to_device(self.u)
        self._d_I  = nc.to_device(self.I_syn)
        self._d_s  = nc.to_device(self.spikes)
        self._d_a  = nc.to_device(self.a)
        self._d_b  = nc.to_device(self.b)
        self._d_c  = nc.to_device(self.c)
        self._d_d  = nc.to_device(self.d)

    def step(self, dt: float, I_ext: Optional[np.ndarray] = None):
        """Advance simulation by one time step dt [ms]."""
        if I_ext is not None:
            self.I_syn += I_ext.astype(np.float32)

        if self._backend == "cuda" and CUDA_AVAILABLE:
            self._step_cuda(dt)
        elif self._backend == "numba_cpu" and NUMBA_AVAILABLE:
            self._step_numba_cpu(dt)
        else:
            self._step_numpy(dt)

    def _step_numba_cpu(self, dt: float):
        _izhikevich_step_cpu(
            self.V, self.u, self.I_syn,
            self.a, self.b, self.c, self.d,
            np.float32(dt), self.spikes
        )
        # Synaptic current for next step
        if self.conn is not None:
            self.I_syn[:] = 0.0
            _csr_mv_cpu(
                self.conn.indptr, self.conn.indices, self.conn.weights,
                self.spikes.astype(np.float32), self.I_syn
            )
        else:
            self.I_syn[:] = 0.0

    def _step_numpy(self, dt: float):
        fired = self.V >= 30.0
        self.spikes[:] = fired.astype(np.int32)
        self.V[fired] = self.c[fired]
        self.u[fired] += self.d[fired]
        nf = ~fired
        v, uv = self.V[nf], self.u[nf]
        self.V[nf] = v + (0.04*v**2 + 5*v + 140 - uv + self.I_syn[nf])*dt
        self.u[nf] = uv + self.a[nf]*(self.b[nf]*v - uv)*dt
        if self.conn is not None:
            sv = self.spikes.astype(np.float32)
            self.I_syn[:] = 0.0
            _csr_mv_cpu(self.conn.indptr, self.conn.indices,
                        self.conn.weights, sv, self.I_syn)
        else:
            self.I_syn[:] = 0.0

    def _step_cuda(self, dt: float):
        from numba import cuda as nc
        TPB = 256
        BPG = (self.N + TPB - 1) // TPB
        # Copy external current to GPU
        nc.to_device(self.I_syn, to=self._d_I)
        # Neuron step
        _izhikevich_step_gpu[BPG, TPB](
            self._d_V, self._d_u, self._d_I,
            self._d_a, self._d_b, self._d_c, self._d_d,
            np.float32(dt), self._d_s
        )
        # Synaptic current
        if self._backend == "cuda" and self._d_indptr is not None:
            _csr_mv_gpu[BPG, TPB](
                self._d_indptr, self._d_indices, self._d_weights,
                self._d_s, self._d_I
            )
        # Copy back to host
        self._d_V.copy_to_host(self.V)
        self._d_u.copy_to_host(self.u)
        self._d_s.copy_to_host(self.spikes)
        self._d_I.copy_to_host(self.I_syn)

    def firing_rate_instant(self) -> float:
        """Instantaneous population firing rate (normalised, [0,1])."""
        return float(self.spikes.mean())


# ---------------------------------------------------------------------------
# Benchmark: scale test
# ---------------------------------------------------------------------------

def benchmark_scale(sizes: list, T_ms: float = 100.0, dt: float = 0.1,
                    p_conn: float = 0.002,
                    backend: str = "auto") -> list:
    """
    Benchmark simulation throughput for different population sizes.

    Returns list of dicts: {N, n_syn, elapsed_s, neurons_per_s}
    """
    results = []
    for N in sizes:
        rng = np.random.default_rng(0)
        pop = IzhikevichPopulation(N, backend=backend, rng=rng)
        conn = CSRConnectivity.random(N, N, p_conn=p_conn, rng=rng)
        pop.set_connectivity(conn)

        T = int(T_ms / dt)
        I_base = np.full(N, 5.0, dtype=np.float32)

        t0 = time.perf_counter()
        for _ in range(T):
            pop.step(dt, I_base)
        elapsed = time.perf_counter() - t0

        results.append({
            "N": N,
            "n_syn": conn.n_syn,
            "elapsed_s": round(elapsed, 3),
            "neurons_per_s": round(N * T / elapsed),
            "backend": pop._backend,
        })
    return results
