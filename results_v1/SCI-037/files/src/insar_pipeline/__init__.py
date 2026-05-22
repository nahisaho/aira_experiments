"""InSAR time-series analysis package for crustal deformation monitoring."""

from .ps_insar import build_ps_network, estimate_ps_coherence, estimate_velocity, select_psc, unwrap_phase
from .sbas import construct_network, estimate_deformation, select_pairs, svd_inversion
from .atmospheric_correction import era5_correction, estimate_aps, gacos_correction, spatial_filter
from .trend_decomposition import extract_transient, fit_linear_trend, fit_seasonal, kalman_filter_decompose
from .precursor_detection import classify_alert_level, detect_acceleration, detect_strain_anomaly, spatial_clustering
from .displacement_3d import decompose_enu, integrate_gps, los_to_3d, propagate_errors
from .nankai_monitoring import detect_sse, estimate_coupling, generate_alert, monitor_stress_accumulation

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "select_psc",
    "estimate_ps_coherence",
    "unwrap_phase",
    "estimate_velocity",
    "build_ps_network",
    "select_pairs",
    "construct_network",
    "svd_inversion",
    "estimate_deformation",
    "era5_correction",
    "gacos_correction",
    "spatial_filter",
    "estimate_aps",
    "fit_linear_trend",
    "fit_seasonal",
    "extract_transient",
    "kalman_filter_decompose",
    "detect_strain_anomaly",
    "detect_acceleration",
    "spatial_clustering",
    "classify_alert_level",
    "los_to_3d",
    "decompose_enu",
    "propagate_errors",
    "integrate_gps",
    "estimate_coupling",
    "detect_sse",
    "monitor_stress_accumulation",
    "generate_alert",
]
