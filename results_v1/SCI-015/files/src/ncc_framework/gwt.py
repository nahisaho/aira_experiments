"""
Global Workspace Theory (GWT) — Network-based consciousness analysis.

Implements:
  - Functional connectivity estimation (correlation, MI, coherence)
  - Graph-theoretic metrics (efficiency, clustering, betweenness)
  - Global workspace ignition detection
  - Information broadcast capacity
  - Thalamo-cortical integration index
"""
import numpy as np
from itertools import combinations
from scipy import signal as scipy_signal
from .utils import mutual_information, transfer_entropy


def phase_coherence(x: np.ndarray, y: np.ndarray,
                    fs: float = 256.0, band: tuple = (8, 13)) -> float:
    """
    Phase-locking value (PLV) between two signals in a frequency band.
    """
    from scipy.signal import butter, filtfilt, hilbert

    lo, hi = band
    nyq = fs / 2
    b, a = butter(4, [lo / nyq, hi / nyq], btype='band')

    try:
        xf = filtfilt(b, a, x)
        yf = filtfilt(b, a, y)
        phi_x = np.angle(hilbert(xf))
        phi_y = np.angle(hilbert(yf))
        return float(np.abs(np.mean(np.exp(1j * (phi_x - phi_y)))))
    except Exception:
        return 0.0


def build_connectivity_matrix(
    data: np.ndarray,
    method: str = "correlation",
    fs: float = 256.0,
    band: tuple = (8, 13),
) -> np.ndarray:
    """
    Build n×n functional connectivity matrix.

    methods: "correlation", "mutual_information", "phase_coherence", "transfer_entropy"
    """
    n = data.shape[0]
    mat = np.zeros((n, n))

    for i, j in combinations(range(n), 2):
        if method == "correlation":
            val = float(np.corrcoef(data[i], data[j])[0, 1])
            val = abs(val)
        elif method == "mutual_information":
            val = mutual_information(data[i], data[j])
        elif method == "phase_coherence":
            val = phase_coherence(data[i], data[j], fs=fs, band=band)
        elif method == "transfer_entropy":
            val = transfer_entropy(data[i], data[j])
        else:
            raise ValueError(f"Unknown method: {method}")
        mat[i, j] = mat[j, i] = val

    return mat


def threshold_matrix(mat: np.ndarray, threshold: float = 0.2) -> np.ndarray:
    """Threshold connectivity matrix to get binary adjacency matrix."""
    adj = (mat > threshold).astype(float)
    np.fill_diagonal(adj, 0)
    return adj


def global_efficiency(adj: np.ndarray) -> float:
    """
    Global efficiency = mean inverse shortest path length.
    Uses Floyd-Warshall algorithm.
    """
    n = adj.shape[0]
    if n == 0:
        return 0.0

    # Initialize distance matrix
    dist = np.full((n, n), np.inf)
    np.fill_diagonal(dist, 0)
    dist[adj > 0] = 1.0

    # Floyd-Warshall
    for k in range(n):
        for i in range(n):
            for j in range(n):
                if dist[i, k] + dist[k, j] < dist[i, j]:
                    dist[i, j] = dist[i, k] + dist[k, j]

    inv_dist = 1.0 / dist
    inv_dist[np.isinf(inv_dist)] = 0
    np.fill_diagonal(inv_dist, 0)

    if n > 1:
        return float(inv_dist.sum() / (n * (n - 1)))
    return 0.0


def local_efficiency(adj: np.ndarray) -> float:
    """Local efficiency = mean global efficiency of node neighborhoods."""
    n = adj.shape[0]
    local_effs = []
    for i in range(n):
        neighbors = np.where(adj[i] > 0)[0]
        if len(neighbors) >= 2:
            sub_adj = adj[np.ix_(neighbors, neighbors)]
            local_effs.append(global_efficiency(sub_adj))
        else:
            local_effs.append(0.0)
    return float(np.mean(local_effs))


def clustering_coefficient(adj: np.ndarray) -> float:
    """Mean clustering coefficient of the network."""
    n = adj.shape[0]
    cc = []
    for i in range(n):
        neighbors = np.where(adj[i] > 0)[0]
        k = len(neighbors)
        if k < 2:
            cc.append(0.0)
            continue
        # Edges among neighbors
        sub = adj[np.ix_(neighbors, neighbors)]
        e = sub.sum() / 2
        cc.append(2 * e / (k * (k - 1)))
    return float(np.mean(cc))


def small_world_index(adj: np.ndarray) -> float:
    """
    Small-world index: σ = (C/C_rand) / (L/L_rand)
    Approximates C_rand and L_rand using random graph estimates.
    """
    n = adj.shape[0]
    k = adj.sum() / n  # mean degree

    c = clustering_coefficient(adj)
    e = global_efficiency(adj)

    if k < 1:
        return 0.0

    # Random graph approximation
    p = k / (n - 1)
    c_rand = max(p, 1e-10)
    l_rand = max(np.log(n) / np.log(max(k, 2)), 1e-10)

    l = 1.0 / max(e, 1e-10)  # approximate mean path length
    sigma = (c / c_rand) / (l / l_rand)
    return float(sigma)


def information_broadcast_capacity(
    data: np.ndarray, threshold: float = 0.3
) -> float:
    """
    GWT information broadcast capacity:
    Fraction of neural population receiving information above threshold connectivity.
    """
    mat = build_connectivity_matrix(data, method="correlation")
    adj = threshold_matrix(mat, threshold)
    degree = adj.sum(axis=1)
    return float(degree.mean() / (data.shape[0] - 1))


def ignition_index(data: np.ndarray, threshold: float = 0.3) -> float:
    """
    GWT ignition index: measures sudden recruitment of a global workspace.
    Approximated by the maximum eigenvalue of the connectivity matrix.
    (High eigenvalue → all-or-none ignition dynamics)
    """
    mat = build_connectivity_matrix(data, method="correlation")
    try:
        eigenvalues = np.linalg.eigvalsh(mat)
        max_eig = eigenvalues[-1]
        # Normalize by network size
        return float(max_eig / data.shape[0])
    except np.linalg.LinAlgError:
        return 0.0


class GlobalWorkspaceAnalyzer:
    """
    Global Workspace Theory (GWT) network analyzer.

    Computes graph-theoretic and information-theoretic indices
    that operationalize GWT predictions for consciousness.

    The GWT predicts:
    - Conscious access → global broadcast (high efficiency, ignition)
    - Unconscious processing → local, clustered activity (low efficiency)
    - Consciousness level correlates with global efficiency and ignition index

    Parameters
    ----------
    threshold : float
        Connectivity threshold for binarizing adjacency matrix
    fs : float
        Sampling rate of input data
    """

    def __init__(self, threshold: float = 0.25, fs: float = 256.0):
        self.threshold = threshold
        self.fs = fs

    def analyze(self, data: np.ndarray) -> dict:
        """
        Full GWT analysis of multi-channel data.

        Returns dict with all graph-theoretic and information metrics.
        """
        mat = build_connectivity_matrix(data, method="correlation")
        adj = threshold_matrix(mat, self.threshold)

        g_eff = global_efficiency(adj)
        l_eff = local_efficiency(adj)
        cc = clustering_coefficient(adj)
        sw = small_world_index(adj)
        broadcast = information_broadcast_capacity(data, self.threshold)
        ignition = ignition_index(data, self.threshold)

        # Transfer entropy network (directed)
        te_mat = np.zeros((data.shape[0], data.shape[0]))
        n = data.shape[0]
        for i in range(min(n, 8)):  # limit for tractability
            for j in range(min(n, 8)):
                if i != j:
                    te_mat[i, j] = transfer_entropy(data[i], data[j])

        te_asymmetry = float(np.mean(np.abs(te_mat - te_mat.T)))

        return {
            "global_efficiency": g_eff,
            "local_efficiency": l_eff,
            "clustering_coefficient": cc,
            "small_world_index": sw,
            "information_broadcast_capacity": broadcast,
            "ignition_index": ignition,
            "te_asymmetry": te_asymmetry,
            "connectivity_matrix": mat,
            "adjacency_matrix": adj,
        }

    def gwt_consciousness_index(self, data: np.ndarray) -> float:
        """
        Composite GWT-based consciousness index (GCI):
        Weighted combination of global efficiency, ignition, and broadcast.
        """
        metrics = self.analyze(data)
        gci = (
            0.4 * metrics["global_efficiency"]
            + 0.4 * metrics["ignition_index"]
            + 0.2 * metrics["information_broadcast_capacity"]
        )
        return float(np.clip(gci, 0, 1))
