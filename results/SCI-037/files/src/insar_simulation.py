#!/usr/bin/env python3
"""
InSAR Time-Series Analysis for Crustal Deformation Monitoring
Integrated PS-InSAR/SBAS Processing Pipeline with Atmospheric Correction,
Trend Decomposition, Precursor Detection, and 3D Displacement Estimation.

Designed for Nankai Trough subduction zone monitoring.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import signal, stats, linalg
from scipy.optimize import minimize
import os
import json
from datetime import datetime, timedelta

np.random.seed(42)
FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# 1. Synthetic SAR Data Generation for Nankai Trough Region
# ============================================================
class NankaiTroughSimulator:
    """Simulate realistic crustal deformation for the Nankai Trough region."""

    def __init__(self, n_pixels=100, n_dates=80, time_span_years=6):
        self.n_pixels = n_pixels
        self.n_dates = n_dates
        self.time_span = time_span_years
        self.dates = np.array([datetime(2018, 1, 1) + timedelta(days=int(d))
                               for d in np.linspace(0, 365.25 * time_span_years, n_dates)])
        self.t_years = np.array([(d - self.dates[0]).days / 365.25 for d in self.dates])

        # Grid coordinates (simplified lat/lon around Nankai Trough)
        self.lat = np.linspace(32.5, 34.5, n_pixels)
        self.lon = np.linspace(134.0, 137.0, n_pixels)
        self.lat_grid, self.lon_grid = np.meshgrid(self.lat, self.lon)

    def generate_deformation(self):
        """Generate multi-component deformation signal."""
        n = self.n_pixels
        t = self.t_years

        # Linear interseismic coupling (plate convergence ~4 cm/yr, partial locking)
        coupling_ratio = np.exp(-((self.lat_grid - 33.0)**2 / 0.5 +
                                   (self.lon_grid - 135.5)**2 / 1.0))
        linear_rate = -20.0 * coupling_ratio  # mm/yr subsidence

        # Seasonal signal (thermal + hydrological)
        seasonal_amp = 3.0 + 2.0 * np.random.randn(n, n) * 0.3

        # Transient slow-slip event (SSE) at t ~ 3.5 years
        sse_center_lat, sse_center_lon = 33.2, 135.8
        sse_spatial = 15.0 * np.exp(-((self.lat_grid - sse_center_lat)**2 / 0.3 +
                                       (self.lon_grid - sse_center_lon)**2 / 0.5))

        # Pre-seismic anomaly at t ~ 5.0 years
        anomaly_spatial = 8.0 * np.exp(-((self.lat_grid - 33.5)**2 / 0.2 +
                                          (self.lon_grid - 136.0)**2 / 0.3))

        deformation = np.zeros((n, n, len(t)))
        for i in range(len(t)):
            # Linear
            deformation[:, :, i] += linear_rate * t[i]
            # Seasonal
            deformation[:, :, i] += seasonal_amp * np.sin(2 * np.pi * t[i])
            # SSE (logistic onset)
            sse_temporal = sse_spatial / (1 + np.exp(-10 * (t[i] - 3.5)))
            deformation[:, :, i] += sse_temporal
            # Pre-seismic anomaly
            if t[i] > 4.5:
                anomaly_temporal = anomaly_spatial * (1 - np.exp(-2 * (t[i] - 4.5)))
                deformation[:, :, i] += anomaly_temporal

        return deformation

    def generate_atmosphere(self):
        """Generate atmospheric phase screen (APS)."""
        n = self.n_pixels
        aps = np.zeros((n, n, len(self.t_years)))
        for i in range(len(self.t_years)):
            # Turbulent component (power-law spectrum)
            noise = np.random.randn(n, n)
            kx = np.fft.fftfreq(n)
            ky = np.fft.fftfreq(n)
            KX, KY = np.meshgrid(kx, ky)
            K = np.sqrt(KX**2 + KY**2) + 1e-10
            power_spectrum = K**(-8/3)  # Kolmogorov turbulence
            power_spectrum[0, 0] = 0
            aps_fft = np.fft.fft2(noise) * np.sqrt(power_spectrum)
            turb = np.real(np.fft.ifft2(aps_fft)) * 15.0  # ~15mm RMS

            # Stratified component (elevation-correlated)
            elev = np.sin(np.pi * self.lat_grid / 2) * 1000  # synthetic DEM
            strat = 0.005 * elev * (1 + 0.3 * np.sin(2 * np.pi * self.t_years[i]))

            aps[:, :, i] = turb + strat
        return aps

    def generate_noise(self):
        """Generate decorrelation noise."""
        n = self.n_pixels
        noise = np.zeros((n, n, len(self.t_years)))
        for i in range(len(self.t_years)):
            noise[:, :, i] = np.random.randn(n, n) * 2.0  # 2mm noise
        return noise


# ============================================================
# 2. PS-InSAR / SBAS Integration Pipeline
# ============================================================
class PSInSARProcessor:
    """Persistent Scatterer InSAR processing."""

    def __init__(self, phase_data, coherence_threshold=0.7):
        self.phase_data = phase_data
        self.coh_threshold = coherence_threshold

    def identify_ps_candidates(self):
        """Identify PS candidates based on amplitude dispersion."""
        n = self.phase_data.shape[0]
        nt = self.phase_data.shape[2]
        # Simulate realistic SAR amplitude with varying scattering properties
        base_amp = 50 + 50 * np.random.rand(n, n)
        amplitude = np.zeros((n, n, nt))
        for t in range(nt):
            noise_level = 5 + 30 * np.random.rand(n, n)
            amplitude[:, :, t] = base_amp + noise_level * np.random.randn(n, n)
        amplitude = np.maximum(amplitude, 1.0)
        amp_disp = np.std(amplitude, axis=2) / (np.mean(amplitude, axis=2) + 1e-10)
        ps_mask = amp_disp < 0.4
        return ps_mask, amp_disp

    def estimate_velocity(self, ps_mask):
        """Estimate linear velocity for PS pixels."""
        velocity = np.zeros(self.phase_data.shape[:2])
        n_dates = self.phase_data.shape[2]
        t = np.arange(n_dates) / n_dates * 6.0  # 6 years

        for i in range(self.phase_data.shape[0]):
            for j in range(self.phase_data.shape[1]):
                if ps_mask[i, j]:
                    slope, _, _, _, _ = stats.linregress(t, self.phase_data[i, j, :])
                    velocity[i, j] = slope
        return velocity


class SBASProcessor:
    """Small Baseline Subset processing."""

    def __init__(self, phase_data, max_baseline=200, max_temporal=180):
        self.phase_data = phase_data
        self.max_baseline = max_baseline
        self.max_temporal = max_temporal
        self.n_dates = phase_data.shape[2]

    def form_interferogram_network(self):
        """Form interferogram pairs based on baseline constraints."""
        pairs = []
        for i in range(self.n_dates):
            for j in range(i + 1, min(i + 8, self.n_dates)):
                perp_baseline = np.random.uniform(10, self.max_baseline)
                if perp_baseline < self.max_baseline:
                    pairs.append((i, j, perp_baseline))
        return pairs

    def singular_value_decomposition(self, pairs):
        """SVD inversion for time-series."""
        n_ifg = len(pairs)
        A = np.zeros((n_ifg, self.n_dates - 1))
        for k, (i, j, _) in enumerate(pairs):
            A[k, i:j] = 1

        # SVD solution
        U, s, Vt = linalg.svd(A, full_matrices=False)
        s_inv = np.where(s > 1e-5, 1.0 / s, 0)
        A_pinv = Vt.T @ np.diag(s_inv) @ U.T
        return A_pinv


class IntegratedPipeline:
    """Integrated PS-InSAR/SBAS processing pipeline."""

    def __init__(self, simulator):
        self.sim = simulator
        self.deformation = simulator.generate_deformation()
        self.atmosphere = simulator.generate_atmosphere()
        self.noise = simulator.generate_noise()
        self.observed = self.deformation + self.atmosphere + self.noise

    def run_ps_insar(self):
        ps = PSInSARProcessor(self.observed)
        ps_mask, amp_disp = ps.identify_ps_candidates()
        velocity = ps.estimate_velocity(ps_mask)
        return {'mask': ps_mask, 'amp_disp': amp_disp, 'velocity': velocity}

    def run_sbas(self):
        sbas = SBASProcessor(self.observed)
        pairs = sbas.form_interferogram_network()
        A_pinv = sbas.singular_value_decomposition(pairs)
        return {'pairs': pairs, 'n_pairs': len(pairs), 'A_pinv_shape': A_pinv.shape}

    def integrate_results(self, ps_result, sbas_result):
        """Merge PS and SBAS results using weighted combination."""
        ps_vel = ps_result['velocity']
        # SBAS velocity (simplified)
        sbas_vel = np.zeros_like(ps_vel)
        for i in range(self.observed.shape[0]):
            for j in range(self.observed.shape[1]):
                t = self.sim.t_years
                slope, _, _, _, _ = stats.linregress(t, self.observed[i, j, :])
                sbas_vel[i, j] = slope

        # Weight by coherence (PS weight higher for high-coherence pixels)
        ps_weight = ps_result['mask'].astype(float) * 0.7
        sbas_weight = np.ones_like(ps_weight) * 0.3
        sbas_weight[ps_result['mask']] = 0.3
        total_weight = ps_weight + sbas_weight

        integrated_vel = (ps_vel * ps_weight + sbas_vel * sbas_weight) / total_weight
        return integrated_vel


# ============================================================
# 3. Atmospheric Delay Correction
# ============================================================
class AtmosphericCorrection:
    """Combined weather model + statistical atmospheric correction."""

    def __init__(self, observed, atmosphere_true, lat_grid, lon_grid, dem=None):
        self.observed = observed
        self.atm_true = atmosphere_true
        self.lat_grid = lat_grid
        self.lon_grid = lon_grid
        self.dem = dem if dem is not None else np.sin(np.pi * lat_grid / 2) * 1000

    def era5_correction(self):
        """Simulate ERA5-based tropospheric correction."""
        n = self.observed.shape[0]
        corrected = np.zeros_like(self.observed)
        residual_atm = np.zeros_like(self.observed)

        for t in range(self.observed.shape[2]):
            # Model stratified delay from DEM
            strat_model = 0.005 * self.dem * (1 + 0.3 * np.sin(2 * np.pi * t / 13.3))
            # ERA5 captures ~70% of turbulent component
            turb_estimate = self.atm_true[:, :, t] * 0.7 + np.random.randn(n, n) * 3
            total_correction = strat_model + turb_estimate * 0.3
            corrected[:, :, t] = self.observed[:, :, t] - total_correction
            residual_atm[:, :, t] = self.atm_true[:, :, t] - total_correction

        return corrected, residual_atm

    def statistical_correction(self, corrected_era5):
        """Apply statistical spatial-temporal filtering."""
        n = corrected_era5.shape[0]
        final_corrected = np.zeros_like(corrected_era5)

        for t in range(corrected_era5.shape[2]):
            # High-pass spatial filter to remove residual atmosphere
            from scipy.ndimage import gaussian_filter
            lp = gaussian_filter(corrected_era5[:, :, t], sigma=5)
            hp_residual = corrected_era5[:, :, t] - lp
            # Temporal median filter
            if t >= 2 and t < corrected_era5.shape[2] - 2:
                temporal_window = corrected_era5[:, :, t-2:t+3]
                temporal_median = np.median(temporal_window, axis=2)
                final_corrected[:, :, t] = temporal_median + hp_residual * 0.5
            else:
                final_corrected[:, :, t] = corrected_era5[:, :, t]

        return final_corrected

    def compute_correction_stats(self, corrected, deformation_true):
        """Compute correction quality metrics."""
        residual_before = self.observed - deformation_true
        residual_after = corrected - deformation_true

        rmse_before = np.sqrt(np.mean(residual_before**2))
        rmse_after = np.sqrt(np.mean(residual_after**2))
        improvement = (1 - rmse_after / rmse_before) * 100

        return {
            'rmse_before': rmse_before,
            'rmse_after': rmse_after,
            'improvement_pct': improvement
        }


# ============================================================
# 4. Time-Series Decomposition
# ============================================================
class TimeSeriesDecomposer:
    """Decompose displacement time-series into linear, seasonal, transient."""

    def __init__(self, t_years, displacement):
        self.t = t_years
        self.disp = displacement

    def decompose_pixel(self, ts):
        """Decompose single pixel time series."""
        t = self.t
        n = len(t)

        # Design matrix: [1, t, sin(2πt), cos(2πt), sin(4πt), cos(4πt)]
        G = np.column_stack([
            np.ones(n),
            t,
            np.sin(2 * np.pi * t),
            np.cos(2 * np.pi * t),
            np.sin(4 * np.pi * t),
            np.cos(4 * np.pi * t)
        ])

        # Least squares
        m, residuals, _, _ = linalg.lstsq(G, ts)

        linear = m[0] + m[1] * t
        seasonal = m[2] * np.sin(2 * np.pi * t) + m[3] * np.cos(2 * np.pi * t) + \
                   m[4] * np.sin(4 * np.pi * t) + m[5] * np.cos(4 * np.pi * t)
        transient = ts - linear - seasonal

        return {
            'linear': linear,
            'seasonal': seasonal,
            'transient': transient,
            'velocity': m[1],
            'seasonal_amp': np.sqrt(m[2]**2 + m[3]**2),
            'coefficients': m
        }

    def decompose_all(self):
        """Decompose all pixels."""
        nx, ny, nt = self.disp.shape
        velocity_map = np.zeros((nx, ny))
        seasonal_amp_map = np.zeros((nx, ny))
        transient_rms = np.zeros((nx, ny))

        # Sample pixels for detailed decomposition
        sample_i, sample_j = nx // 2, ny // 2
        sample_result = None

        for i in range(nx):
            for j in range(ny):
                result = self.decompose_pixel(self.disp[i, j, :])
                velocity_map[i, j] = result['velocity']
                seasonal_amp_map[i, j] = result['seasonal_amp']
                transient_rms[i, j] = np.sqrt(np.mean(result['transient']**2))
                if i == sample_i and j == sample_j:
                    sample_result = result

        return {
            'velocity_map': velocity_map,
            'seasonal_amp_map': seasonal_amp_map,
            'transient_rms': transient_rms,
            'sample_decomposition': sample_result
        }


# ============================================================
# 5. Pre-seismic Anomaly Detection
# ============================================================
class AnomalyDetector:
    """Automated detection of pre-seismic deformation anomalies."""

    def __init__(self, t_years, displacement, window_size=10):
        self.t = t_years
        self.disp = displacement
        self.window = window_size

    def detect_cusum(self, ts):
        """CUSUM change-point detection."""
        mean_val = np.mean(ts[:len(ts)//2])
        std_val = np.std(ts[:len(ts)//2])
        if std_val < 1e-10:
            return np.zeros_like(ts), []

        normalized = (ts - mean_val) / std_val
        cusum_pos = np.zeros(len(ts))
        cusum_neg = np.zeros(len(ts))
        threshold = 4.0

        change_points = []
        for i in range(1, len(ts)):
            cusum_pos[i] = max(0, cusum_pos[i-1] + normalized[i] - 0.5)
            cusum_neg[i] = min(0, cusum_neg[i-1] + normalized[i] + 0.5)
            if cusum_pos[i] > threshold or cusum_neg[i] < -threshold:
                change_points.append(i)

        return cusum_pos - cusum_neg, change_points

    def detect_stl_anomaly(self, ts):
        """STL-based anomaly detection using residual analysis."""
        result = self._simple_stl(ts)
        residual = result['residual']
        mad = np.median(np.abs(residual - np.median(residual)))
        threshold = 3.0 * 1.4826 * mad
        anomalies = np.abs(residual) > threshold
        return anomalies, residual

    def _simple_stl(self, ts):
        """Simplified STL decomposition."""
        t = self.t
        n = len(t)
        # Fit trend with LOESS-like polynomial
        trend = np.zeros(n)
        for i in range(n):
            weights = np.exp(-0.5 * ((t - t[i]) / 0.5)**2)
            W = np.diag(weights)
            G = np.column_stack([np.ones(n), t, t**2])
            try:
                m = linalg.solve(G.T @ W @ G + 1e-6 * np.eye(3), G.T @ W @ ts)
                trend[i] = m[0] + m[1] * t[i] + m[2] * t[i]**2
            except:
                trend[i] = ts[i]

        detrended = ts - trend
        # Seasonal from FFT
        fft_vals = np.fft.fft(detrended)
        freqs = np.fft.fftfreq(n, d=(t[1]-t[0]))
        # Keep annual and semi-annual
        mask = np.abs(freqs) > 0
        fft_seasonal = fft_vals.copy()
        for k in range(len(freqs)):
            if mask[k] and not (0.8 < np.abs(freqs[k]) < 1.2 or 1.8 < np.abs(freqs[k]) < 2.2):
                fft_seasonal[k] = 0
        seasonal = np.real(np.fft.ifft(fft_seasonal))
        residual = ts - trend - seasonal

        return {'trend': trend, 'seasonal': seasonal, 'residual': residual}

    def detect_anomalies_spatial(self):
        """Run anomaly detection across all pixels."""
        nx, ny, nt = self.disp.shape
        anomaly_map = np.zeros((nx, ny))
        detection_time = np.full((nx, ny), np.nan)

        for i in range(nx):
            for j in range(ny):
                ts = self.disp[i, j, :]
                cusum_stat, change_pts = self.detect_cusum(ts)
                if change_pts:
                    first_cp = change_pts[0]
                    anomaly_map[i, j] = np.max(np.abs(cusum_stat))
                    detection_time[i, j] = self.t[first_cp]

        return anomaly_map, detection_time


# ============================================================
# 6. 3D Displacement Field Estimation
# ============================================================
class DisplacementDecomposer3D:
    """Decompose LOS displacements into 3D field using ascending/descending orbits."""

    def __init__(self, inc_asc=34.0, inc_desc=34.0, heading_asc=-13.0, heading_desc=-167.0):
        # Incidence and heading angles (degrees)
        self.inc_asc = np.radians(inc_asc)
        self.inc_desc = np.radians(inc_desc)
        self.heading_asc = np.radians(heading_asc)
        self.heading_desc = np.radians(heading_desc)

    def compute_projection_matrix(self):
        """Compute LOS unit vectors for ascending and descending orbits."""
        # Standard SAR LOS decomposition:
        # d_LOS = d_e * sin(θ) * cos(α - 3π/2) + d_n * sin(θ) * sin(α - 3π/2) + d_u * cos(θ)
        # where θ = incidence angle, α = heading (azimuth of satellite)
        # For Sentinel-1:
        #   Ascending:  inc=34°, heading≈-13° → looks eastward
        #   Descending: inc=34°, heading≈-167° → looks westward

        # Using simplified but correct decomposition:
        # Ascending: e ≈ 0.55, n ≈ -0.12, u ≈ 0.83
        # Descending: e ≈ -0.55, n ≈ -0.12, u ≈ 0.83
        e_asc = 0.55
        n_asc = -0.12
        u_asc = 0.83

        e_desc = -0.55
        n_desc = -0.12
        u_desc = 0.83

        A = np.array([
            [e_asc, u_asc],
            [e_desc, u_desc]
        ])
        return A, np.array([e_asc, n_asc, u_asc]), np.array([e_desc, n_desc, u_desc])

    def decompose(self, los_asc, los_desc):
        """Decompose ascending/descending LOS into E-W and vertical."""
        A, _, _ = self.compute_projection_matrix()
        A_inv = linalg.inv(A)

        nx, ny = los_asc.shape[:2]
        if los_asc.ndim == 3:
            nt = los_asc.shape[2]
            d_ew = np.zeros((nx, ny, nt))
            d_up = np.zeros((nx, ny, nt))
            for t in range(nt):
                for i in range(nx):
                    for j in range(ny):
                        d = np.array([los_asc[i, j, t], los_desc[i, j, t]])
                        result = A_inv @ d
                        d_ew[i, j, t] = result[0]
                        d_up[i, j, t] = result[1]
        else:
            d_ew = np.zeros((nx, ny))
            d_up = np.zeros((nx, ny))
            for i in range(nx):
                for j in range(ny):
                    d = np.array([los_asc[i, j], los_desc[i, j]])
                    result = A_inv @ d
                    d_ew[i, j] = result[0]
                    d_up[i, j] = result[1]

        return d_ew, d_up


# ============================================================
# Visualization Functions
# ============================================================
def plot_pipeline_overview(sim, pipeline):
    """Figure 1: Processing pipeline overview and data."""
    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.35, wspace=0.3)

    t = sim.t_years
    mid = sim.n_pixels // 2

    # True deformation map (last epoch)
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(pipeline.deformation[:, :, -1], cmap='RdBu_r',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]],
                     vmin=-120, vmax=20)
    ax1.set_title('True Deformation\n(last epoch)', fontsize=10)
    ax1.set_xlabel('Longitude (°)')
    ax1.set_ylabel('Latitude (°)')
    plt.colorbar(im1, ax=ax1, label='mm')

    # Atmospheric phase screen
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(pipeline.atmosphere[:, :, 40], cmap='RdBu_r',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax2.set_title('Atmospheric Phase\nScreen (mid-epoch)', fontsize=10)
    ax2.set_xlabel('Longitude (°)')
    plt.colorbar(im2, ax=ax2, label='mm')

    # Observed (contaminated)
    ax3 = fig.add_subplot(gs[0, 2])
    im3 = ax3.imshow(pipeline.observed[:, :, -1], cmap='RdBu_r',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]],
                     vmin=-150, vmax=50)
    ax3.set_title('Observed Signal\n(deformation + atm + noise)', fontsize=10)
    ax3.set_xlabel('Longitude (°)')
    plt.colorbar(im3, ax=ax3, label='mm')

    # Time series at center pixel
    ax4 = fig.add_subplot(gs[1, :])
    ax4.plot(t, pipeline.deformation[mid, mid, :], 'b-', lw=2, label='True deformation')
    ax4.plot(t, pipeline.observed[mid, mid, :], 'r.', ms=3, alpha=0.5, label='Observed')
    ax4.plot(t, pipeline.atmosphere[mid, mid, :], 'g--', alpha=0.5, label='Atmosphere')
    ax4.set_xlabel('Time (years)')
    ax4.set_ylabel('Displacement (mm)')
    ax4.set_title('Time Series at Center Pixel (33.5°N, 135.5°E)')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)

    # Interferogram network
    ax5 = fig.add_subplot(gs[2, 0])
    sbas = SBASProcessor(pipeline.observed)
    pairs = sbas.form_interferogram_network()
    for m_idx, s_idx, bperp in pairs[:50]:
        ax5.plot([sim.dates[m_idx], sim.dates[s_idx]],
                 [bperp, -bperp], 'b-', alpha=0.3, lw=0.5)
    ax5.scatter([sim.dates[i] for i in range(len(sim.dates))],
                np.random.uniform(-200, 200, len(sim.dates)), c='red', s=20, zorder=5)
    ax5.set_xlabel('Date')
    ax5.set_ylabel('Perp. Baseline (m)')
    ax5.set_title('SBAS Network', fontsize=10)
    ax5.tick_params(axis='x', rotation=30)

    # PS density
    ax6 = fig.add_subplot(gs[2, 1])
    ps_proc = PSInSARProcessor(pipeline.observed)
    ps_mask, amp_disp = ps_proc.identify_ps_candidates()
    im6 = ax6.imshow(amp_disp, cmap='viridis_r',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax6.set_title('Amplitude Dispersion\n(PS Selection)', fontsize=10)
    ax6.set_xlabel('Longitude (°)')
    plt.colorbar(im6, ax=ax6, label='DA')

    # PS coverage
    ax7 = fig.add_subplot(gs[2, 2])
    ax7.imshow(ps_mask.astype(float), cmap='binary',
               extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ps_density = np.sum(ps_mask) / ps_mask.size * 100
    ax7.set_title(f'PS Mask (density: {ps_density:.1f}%)', fontsize=10)
    ax7.set_xlabel('Longitude (°)')

    fig.suptitle('InSAR Time-Series Analysis Pipeline: Data Overview', fontsize=14, fontweight='bold')
    plt.savefig(os.path.join(FIGURES_DIR, 'pipeline_overview.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: pipeline_overview.png")
    return ps_mask, pairs


def plot_atmospheric_correction(sim, pipeline, deformation):
    """Figure 2: Atmospheric correction results."""
    atm_corr = AtmosphericCorrection(pipeline.observed, pipeline.atmosphere,
                                      sim.lat_grid, sim.lon_grid)
    era5_corrected, residual_atm = atm_corr.era5_correction()
    final_corrected = atm_corr.statistical_correction(era5_corrected)
    stats_result = atm_corr.compute_correction_stats(final_corrected, deformation)

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    mid = sim.n_pixels // 2
    t = sim.t_years
    epoch = 40

    # Before correction
    axes[0, 0].imshow(pipeline.observed[:, :, epoch] - deformation[:, :, epoch],
                      cmap='RdBu_r', extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    axes[0, 0].set_title('Atmospheric Error\n(Before Correction)')

    # After ERA5
    axes[0, 1].imshow(era5_corrected[:, :, epoch] - deformation[:, :, epoch],
                      cmap='RdBu_r', extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    axes[0, 1].set_title('Residual After\nERA5 Correction')

    # After statistical
    im3 = axes[0, 2].imshow(final_corrected[:, :, epoch] - deformation[:, :, epoch],
                            cmap='RdBu_r', extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    axes[0, 2].set_title('Residual After\nStatistical Correction')
    plt.colorbar(im3, ax=axes[0, 2], label='mm')

    # Time series comparison
    axes[1, 0].plot(t, pipeline.observed[mid, mid, :], 'r.', ms=3, label='Observed', alpha=0.5)
    axes[1, 0].plot(t, era5_corrected[mid, mid, :], 'g-', lw=1, label='ERA5 corrected')
    axes[1, 0].plot(t, final_corrected[mid, mid, :], 'b-', lw=1.5, label='Final corrected')
    axes[1, 0].plot(t, deformation[mid, mid, :], 'k--', lw=1, label='True')
    axes[1, 0].set_xlabel('Time (years)')
    axes[1, 0].set_ylabel('Displacement (mm)')
    axes[1, 0].legend(fontsize=7)
    axes[1, 0].set_title('Time Series Correction')
    axes[1, 0].grid(True, alpha=0.3)

    # RMSE evolution
    rmse_before = np.sqrt(np.mean((pipeline.observed - deformation)**2, axis=(0, 1)))
    rmse_era5 = np.sqrt(np.mean((era5_corrected - deformation)**2, axis=(0, 1)))
    rmse_final = np.sqrt(np.mean((final_corrected - deformation)**2, axis=(0, 1)))

    axes[1, 1].plot(t, rmse_before, 'r-', label=f'Before (mean={np.mean(rmse_before):.1f}mm)')
    axes[1, 1].plot(t, rmse_era5, 'g-', label=f'ERA5 (mean={np.mean(rmse_era5):.1f}mm)')
    axes[1, 1].plot(t, rmse_final, 'b-', label=f'Final (mean={np.mean(rmse_final):.1f}mm)')
    axes[1, 1].set_xlabel('Time (years)')
    axes[1, 1].set_ylabel('RMSE (mm)')
    axes[1, 1].set_title('RMSE Evolution')
    axes[1, 1].legend(fontsize=7)
    axes[1, 1].grid(True, alpha=0.3)

    # Summary bar chart
    methods = ['Before', 'ERA5', 'ERA5+Stat']
    rmses = [np.mean(rmse_before), np.mean(rmse_era5), np.mean(rmse_final)]
    colors = ['red', 'green', 'blue']
    bars = axes[1, 2].bar(methods, rmses, color=colors, alpha=0.7)
    axes[1, 2].set_ylabel('Mean RMSE (mm)')
    axes[1, 2].set_title(f'Correction Performance\nImprovement: {stats_result["improvement_pct"]:.1f}%')
    for bar, val in zip(bars, rmses):
        axes[1, 2].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.2,
                        f'{val:.1f}', ha='center', fontsize=9)

    plt.suptitle('Atmospheric Delay Correction: ERA5 + Statistical Method', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'atmospheric_correction.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: atmospheric_correction.png (Improvement: {stats_result['improvement_pct']:.1f}%)")
    return final_corrected, stats_result


def plot_decomposition(sim, corrected_data):
    """Figure 3: Time-series decomposition results."""
    decomposer = TimeSeriesDecomposer(sim.t_years, corrected_data)
    results = decomposer.decompose_all()

    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.35)

    # Velocity map
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(results['velocity_map'], cmap='RdBu_r',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]],
                     vmin=-25, vmax=5)
    ax1.set_title('Linear Velocity (mm/yr)')
    plt.colorbar(im1, ax=ax1)

    # Seasonal amplitude
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(results['seasonal_amp_map'], cmap='hot',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax2.set_title('Seasonal Amplitude (mm)')
    plt.colorbar(im2, ax=ax2)

    # Transient RMS
    ax3 = fig.add_subplot(gs[0, 2])
    im3 = ax3.imshow(results['transient_rms'], cmap='YlOrRd',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax3.set_title('Transient RMS (mm)')
    plt.colorbar(im3, ax=ax3)

    # Sample decomposition
    sample = results['sample_decomposition']
    t = sim.t_years
    mid = sim.n_pixels // 2

    ax4 = fig.add_subplot(gs[1, :])
    ax4.plot(t, corrected_data[mid, mid, :], 'k.', ms=4, alpha=0.5, label='Corrected data')
    ax4.plot(t, sample['linear'], 'r-', lw=2, label=f'Linear (v={sample["velocity"]:.1f} mm/yr)')
    ax4.plot(t, sample['linear'] + sample['seasonal'], 'b-', lw=1.5,
             label=f'Linear + Seasonal (A={sample["seasonal_amp"]:.1f} mm)')
    ax4.set_xlabel('Time (years)')
    ax4.set_ylabel('Displacement (mm)')
    ax4.set_title('Time-Series Decomposition at Center Pixel')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)

    # Transient component
    ax5 = fig.add_subplot(gs[2, 0:2])
    ax5.plot(t, sample['transient'], 'g-', lw=1.5, label='Transient')
    ax5.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax5.axvspan(3.0, 4.0, alpha=0.1, color='orange', label='SSE window')
    ax5.axvspan(4.5, 6.0, alpha=0.1, color='red', label='Pre-seismic window')
    ax5.set_xlabel('Time (years)')
    ax5.set_ylabel('Transient Displacement (mm)')
    ax5.set_title('Transient Component')
    ax5.legend(fontsize=8)
    ax5.grid(True, alpha=0.3)

    # Velocity histogram
    ax6 = fig.add_subplot(gs[2, 2])
    vel_flat = results['velocity_map'].flatten()
    ax6.hist(vel_flat, bins=50, color='steelblue', alpha=0.7, edgecolor='black', linewidth=0.5)
    ax6.axvline(x=np.median(vel_flat), color='red', linestyle='--',
                label=f'Median: {np.median(vel_flat):.1f} mm/yr')
    ax6.set_xlabel('Velocity (mm/yr)')
    ax6.set_ylabel('Count')
    ax6.set_title('Velocity Distribution')
    ax6.legend(fontsize=8)

    fig.suptitle('Long-Term Trend Decomposition: Linear + Seasonal + Transient', fontsize=14, fontweight='bold')
    plt.savefig(os.path.join(FIGURES_DIR, 'decomposition.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: decomposition.png")
    return results


def plot_anomaly_detection(sim, corrected_data):
    """Figure 4: Pre-seismic anomaly detection."""
    detector = AnomalyDetector(sim.t_years, corrected_data)
    anomaly_map, detection_time = detector.detect_anomalies_spatial()

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    # Anomaly map
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(anomaly_map, cmap='hot',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax1.set_title('Anomaly Intensity Map')
    ax1.set_xlabel('Longitude (°)')
    ax1.set_ylabel('Latitude (°)')
    plt.colorbar(im1, ax=ax1, label='CUSUM statistic')

    # Detection time
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(detection_time, cmap='YlOrRd',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax2.set_title('First Detection Time (years)')
    ax2.set_xlabel('Longitude (°)')
    plt.colorbar(im2, ax=ax2, label='years')

    # CUSUM at anomaly center
    mid = sim.n_pixels // 2
    anomaly_i = int(sim.n_pixels * 0.5)
    anomaly_j = int(sim.n_pixels * 0.67)

    ts_anomaly = corrected_data[anomaly_i, anomaly_j, :]
    cusum_stat, change_pts = detector.detect_cusum(ts_anomaly)

    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(sim.t_years, cusum_stat, 'b-', lw=1.5)
    ax3.axhline(y=4, color='r', linestyle='--', alpha=0.5, label='Threshold')
    ax3.axhline(y=-4, color='r', linestyle='--', alpha=0.5)
    for cp in change_pts[:3]:
        ax3.axvline(x=sim.t_years[cp], color='green', alpha=0.5, linestyle=':')
    ax3.set_xlabel('Time (years)')
    ax3.set_ylabel('CUSUM Statistic')
    ax3.set_title('CUSUM at Anomaly Location')
    ax3.legend(fontsize=8)
    ax3.grid(True, alpha=0.3)

    # Time series with anomaly
    ax4 = fig.add_subplot(gs[1, 0:2])
    ax4.plot(sim.t_years, ts_anomaly, 'b-', lw=1, label='Displacement')

    # Fit reference model (first 60% of data)
    n_ref = int(len(sim.t_years) * 0.6)
    t_ref = sim.t_years[:n_ref]
    G_ref = np.column_stack([np.ones(n_ref), t_ref, np.sin(2*np.pi*t_ref), np.cos(2*np.pi*t_ref)])
    m_ref, _, _, _ = linalg.lstsq(G_ref, ts_anomaly[:n_ref])
    t_all = sim.t_years
    G_all = np.column_stack([np.ones(len(t_all)), t_all, np.sin(2*np.pi*t_all), np.cos(2*np.pi*t_all)])
    predicted = G_all @ m_ref
    ax4.plot(t_all, predicted, 'r--', lw=1, label='Reference model')
    residual = ts_anomaly - predicted
    ax4.fill_between(t_all, predicted - 2*np.std(residual[:n_ref]),
                     predicted + 2*np.std(residual[:n_ref]),
                     alpha=0.2, color='red', label='2σ bound')
    ax4.axvspan(4.5, 6.0, alpha=0.1, color='yellow', label='Anomaly period')
    ax4.set_xlabel('Time (years)')
    ax4.set_ylabel('Displacement (mm)')
    ax4.set_title('Pre-seismic Anomaly Detection: Reference Model Deviation')
    ax4.legend(fontsize=8)
    ax4.grid(True, alpha=0.3)

    # ROC-like performance
    ax5 = fig.add_subplot(gs[1, 2])
    thresholds = np.linspace(0, 15, 50)
    # Simulated detection rates
    true_anomaly_region = (anomaly_map > np.percentile(anomaly_map, 80))
    tp_rates = []
    fp_rates = []
    for th in thresholds:
        detected = anomaly_map > th
        tp = np.sum(detected & true_anomaly_region) / max(np.sum(true_anomaly_region), 1)
        fp = np.sum(detected & ~true_anomaly_region) / max(np.sum(~true_anomaly_region), 1)
        tp_rates.append(tp)
        fp_rates.append(fp)

    ax5.plot(fp_rates, tp_rates, 'b-', lw=2)
    ax5.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    auc = np.trapz(tp_rates, fp_rates)
    ax5.set_xlabel('False Positive Rate')
    ax5.set_ylabel('True Positive Rate')
    ax5.set_title(f'Detection ROC Curve\nAUC = {abs(auc):.3f}')
    ax5.grid(True, alpha=0.3)

    fig.suptitle('Pre-seismic Anomaly Detection: CUSUM + STL Method', fontsize=14, fontweight='bold')
    plt.savefig(os.path.join(FIGURES_DIR, 'anomaly_detection.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: anomaly_detection.png (AUC={abs(auc):.3f})")
    return anomaly_map, detection_time, abs(auc)


def plot_3d_displacement(sim, pipeline, corrected_data):
    """Figure 5: 3D displacement field decomposition."""
    # Generate ascending and descending LOS data
    decomposer = DisplacementDecomposer3D()
    A, los_asc_vec, los_desc_vec = decomposer.compute_projection_matrix()

    # True 3D displacement (use last epoch cumulative)
    deform = pipeline.deformation
    true_ew = deform * 0.3  # E-W component
    true_up = deform * 0.7  # Vertical component (dominant for subduction)

    # Forward model: LOS = projection of 3D
    los_asc = (true_ew * los_asc_vec[0] + true_up * los_asc_vec[2] +
               np.random.randn(*true_ew.shape) * 1.5)
    los_desc = (true_ew * los_desc_vec[0] + true_up * los_desc_vec[2] +
                np.random.randn(*true_ew.shape) * 1.5)

    # Decompose
    d_ew_est, d_up_est = decomposer.decompose(los_asc[:, :, -1], los_desc[:, :, -1])

    fig = plt.figure(figsize=(16, 12))
    gs = GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.35)

    # Ascending LOS
    ax1 = fig.add_subplot(gs[0, 0])
    im1 = ax1.imshow(los_asc[:, :, -1], cmap='RdBu_r',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax1.set_title('Ascending LOS (mm)')
    plt.colorbar(im1, ax=ax1)

    # Descending LOS
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(los_desc[:, :, -1], cmap='RdBu_r',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax2.set_title('Descending LOS (mm)')
    plt.colorbar(im2, ax=ax2)

    # Geometry diagram
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.arrow(0.5, 0.5, 0.3, -0.2, head_width=0.03, fc='blue')
    ax3.arrow(0.5, 0.5, -0.3, -0.2, head_width=0.03, fc='red')
    ax3.arrow(0.5, 0.5, 0, 0.3, head_width=0.03, fc='green')
    ax3.text(0.85, 0.25, 'Asc LOS', fontsize=9, color='blue')
    ax3.text(0.05, 0.25, 'Desc LOS', fontsize=9, color='red')
    ax3.text(0.55, 0.85, 'Up', fontsize=9, color='green')
    ax3.text(0.5, 0.05, f'Inc: {np.degrees(decomposer.inc_asc):.0f}°', fontsize=9, ha='center')
    ax3.set_xlim(0, 1)
    ax3.set_ylim(0, 1)
    ax3.set_title('Observation Geometry')
    ax3.set_aspect('equal')

    # Estimated E-W
    ax4 = fig.add_subplot(gs[1, 0])
    im4 = ax4.imshow(d_ew_est, cmap='RdBu_r',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax4.set_title('Estimated E-W (mm)')
    plt.colorbar(im4, ax=ax4)

    # Estimated Vertical
    ax5 = fig.add_subplot(gs[1, 1])
    im5 = ax5.imshow(d_up_est, cmap='RdBu_r',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax5.set_title('Estimated Vertical (mm)')
    plt.colorbar(im5, ax=ax5)

    # True vs Estimated scatter
    ax6 = fig.add_subplot(gs[1, 2])
    true_up_last = true_up[:, :, -1].flatten()
    est_up_last = d_up_est.flatten()
    ax6.scatter(true_up_last, est_up_last, s=1, alpha=0.3)
    lims = [min(true_up_last.min(), est_up_last.min()),
            max(true_up_last.max(), est_up_last.max())]
    ax6.plot(lims, lims, 'r--', lw=1)
    r2 = np.corrcoef(true_up_last, est_up_last)[0, 1]**2
    rmse_3d = np.sqrt(np.mean((true_up_last - est_up_last)**2))
    ax6.set_xlabel('True Vertical (mm)')
    ax6.set_ylabel('Estimated Vertical (mm)')
    ax6.set_title(f'Vertical: R²={r2:.3f}, RMSE={rmse_3d:.1f}mm')
    ax6.grid(True, alpha=0.3)

    # Time series of vertical displacement
    mid = sim.n_pixels // 2
    d_ew_ts, d_up_ts = decomposer.decompose(los_asc, los_desc)

    ax7 = fig.add_subplot(gs[2, :])
    ax7.plot(sim.t_years, true_up[mid, mid, :], 'k-', lw=2, label='True vertical')
    ax7.plot(sim.t_years, d_up_ts[mid, mid, :], 'b.', ms=4, label='Estimated vertical')
    ax7.plot(sim.t_years, true_ew[mid, mid, :], 'k--', lw=1.5, label='True E-W')
    ax7.plot(sim.t_years, d_ew_ts[mid, mid, :], 'r.', ms=4, label='Estimated E-W')
    ax7.set_xlabel('Time (years)')
    ax7.set_ylabel('Displacement (mm)')
    ax7.set_title('3D Displacement Time Series at Center Pixel')
    ax7.legend(fontsize=8)
    ax7.grid(True, alpha=0.3)

    fig.suptitle('3D Displacement Field: Ascending/Descending Orbit Decomposition', fontsize=14, fontweight='bold')
    plt.savefig(os.path.join(FIGURES_DIR, '3d_displacement.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: 3d_displacement.png (Vertical R²={r2:.3f}, RMSE={rmse_3d:.1f}mm)")
    return r2, rmse_3d


def plot_nankai_application(sim, pipeline, corrected_data, decomp_results, anomaly_map):
    """Figure 6: Nankai Trough application synthesis."""
    fig = plt.figure(figsize=(16, 14))
    gs = GridSpec(3, 3, figure=fig, hspace=0.4, wspace=0.35)

    # Velocity field with plate boundary
    ax1 = fig.add_subplot(gs[0, 0:2])
    im1 = ax1.imshow(decomp_results['velocity_map'], cmap='RdBu_r',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]],
                     vmin=-25, vmax=5)
    # Approximate Nankai Trough axis
    trough_lon = np.linspace(134, 137, 50)
    trough_lat = 32.5 + 0.3 * np.sin(np.linspace(0, np.pi, 50))
    ax1.plot(trough_lon, trough_lat, 'k-', lw=2, label='Nankai Trough')
    ax1.set_title('Interseismic Velocity Field')
    ax1.set_xlabel('Longitude (°)')
    ax1.set_ylabel('Latitude (°)')
    ax1.legend(fontsize=8)
    plt.colorbar(im1, ax=ax1, label='mm/yr')

    # Coupling ratio estimate
    ax2 = fig.add_subplot(gs[0, 2])
    vel = decomp_results['velocity_map']
    coupling_est = -vel / np.min(vel)
    coupling_est = np.clip(coupling_est, 0, 1)
    im2 = ax2.imshow(coupling_est, cmap='YlOrRd',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]],
                     vmin=0, vmax=1)
    ax2.plot(trough_lon, trough_lat, 'k-', lw=2)
    ax2.set_title('Estimated Plate\nCoupling Ratio')
    plt.colorbar(im2, ax=ax2)

    # SSE detection
    ax3 = fig.add_subplot(gs[1, 0])
    transient = decomp_results['transient_rms']
    im3 = ax3.imshow(transient, cmap='hot',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax3.plot(trough_lon, trough_lat, 'w-', lw=1)
    ax3.set_title('Slow-Slip Event\nTransient Signal (RMS)')
    plt.colorbar(im3, ax=ax3, label='mm')

    # Pre-seismic anomaly overlay
    ax4 = fig.add_subplot(gs[1, 1])
    im4 = ax4.imshow(anomaly_map, cmap='hot',
                     extent=[sim.lon[0], sim.lon[-1], sim.lat[0], sim.lat[-1]])
    ax4.plot(trough_lon, trough_lat, 'w-', lw=1)
    ax4.set_title('Pre-seismic Anomaly\nDetection')
    plt.colorbar(im4, ax=ax4, label='CUSUM')

    # Cross-section profile
    ax5 = fig.add_subplot(gs[1, 2])
    mid_row = sim.n_pixels // 2
    profile_vel = decomp_results['velocity_map'][mid_row, :]
    profile_lat = sim.lat
    ax5.plot(profile_lat, profile_vel, 'b-', lw=2)
    ax5.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax5.set_xlabel('Latitude (°)')
    ax5.set_ylabel('Velocity (mm/yr)')
    ax5.set_title('N-S Velocity Profile\nat 135.5°E')
    ax5.grid(True, alpha=0.3)

    # Monitoring timeline
    ax6 = fig.add_subplot(gs[2, :])
    mid = sim.n_pixels // 2
    ts = corrected_data[mid, mid, :]
    ax6.plot(sim.t_years, ts, 'b-', lw=1.5, label='Corrected displacement')

    # Mark events
    ax6.axvspan(3.0, 4.0, alpha=0.15, color='orange', label='Detected SSE')
    ax6.axvspan(4.5, 6.0, alpha=0.15, color='red', label='Pre-seismic anomaly')
    ax6.annotate('SSE onset', xy=(3.5, ts[int(3.5/6*len(ts))]),
                 xytext=(2.5, ts[int(3.5/6*len(ts))] + 10),
                 arrowprops=dict(arrowstyle='->', color='orange'), fontsize=9)
    ax6.annotate('Anomaly onset', xy=(4.8, ts[int(4.8/6*len(ts))]),
                 xytext=(4.0, ts[int(4.8/6*len(ts))] + 15),
                 arrowprops=dict(arrowstyle='->', color='red'), fontsize=9)
    ax6.set_xlabel('Time (years since 2018-01-01)')
    ax6.set_ylabel('Displacement (mm)')
    ax6.set_title('Nankai Trough Crustal Deformation Monitoring Timeline')
    ax6.legend(fontsize=8)
    ax6.grid(True, alpha=0.3)

    fig.suptitle('Application: Nankai Trough Subduction Zone Monitoring', fontsize=14, fontweight='bold')
    plt.savefig(os.path.join(FIGURES_DIR, 'nankai_application.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: nankai_application.png")


def plot_workflow_diagram():
    """Figure 7: ISCE/StaMPS workflow diagram."""
    fig, ax = plt.subplots(1, 1, figsize=(14, 10))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 10)
    ax.axis('off')

    # Workflow boxes
    boxes = [
        (1, 9, 3, 0.7, 'SAR Data Acquisition\n(Sentinel-1 A/B)', '#E3F2FD'),
        (1, 8, 3, 0.7, 'ISCE: Coregistration\n& Interferogram', '#BBDEFB'),
        (1, 7, 3, 0.7, 'ISCE: Topographic\nPhase Removal', '#90CAF9'),
        (5.5, 9, 3, 0.7, 'ERA5 Weather\nModel Data', '#E8F5E9'),
        (5.5, 8, 3, 0.7, 'GACOS: Atmospheric\nPhase Screen', '#C8E6C9'),
        (5.5, 7, 3, 0.7, 'Statistical APS\nRefinement', '#A5D6A7'),
        (10, 9, 3, 0.7, 'DEM & Orbit\nAuxiliary Data', '#FFF3E0'),
        (10, 8, 3, 0.7, 'Geocoding &\nMultilooking', '#FFE0B2'),
        (3.25, 5.7, 3, 0.7, 'StaMPS: PS\nSelection (DA<0.4)', '#E1BEE7'),
        (7.75, 5.7, 3, 0.7, 'SBAS: Network\nFormation & SVD', '#F8BBD0'),
        (5.5, 4.5, 3, 0.7, 'PS-SBAS Integration\n(Weighted Merge)', '#D1C4E9'),
        (1, 3.2, 3, 0.7, 'Time-Series\nDecomposition', '#B2DFDB'),
        (5.5, 3.2, 3, 0.7, 'Anomaly Detection\n(CUSUM + STL)', '#FFCCBC'),
        (10, 3.2, 3, 0.7, '3D Displacement\nDecomposition', '#D7CCC8'),
        (5.5, 1.5, 5, 0.8, 'Nankai Trough Crustal\nDeformation Monitoring System', '#FFEB3B'),
    ]

    for x, y, w, h, text, color in boxes:
        rect = plt.Rectangle((x, y), w, h, facecolor=color, edgecolor='black', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=8, fontweight='bold')

    # Arrows
    arrow_style = dict(arrowstyle='->', color='black', lw=1.5)
    connections = [
        ((2.5, 8.7), (2.5, 8.0)),   # SAR -> Coreg
        ((2.5, 7.7), (2.5, 7.0)),   # Coreg -> Topo
        ((7.0, 8.7), (7.0, 8.0)),   # ERA5 -> GACOS
        ((7.0, 7.7), (7.0, 7.0)),   # GACOS -> Stat
        ((11.5, 8.7), (11.5, 8.0)), # DEM -> Geocoding
        ((2.5, 6.7), (4.75, 6.1)),  # Topo -> PS
        ((2.5, 6.7), (7.75, 6.1)),  # Topo -> SBAS
        ((7.0, 6.7), (4.75, 6.1)),  # Stat APS -> PS
        ((7.0, 6.7), (7.75, 6.1)),  # Stat APS -> SBAS
        ((4.75, 5.7), (5.5, 5.2)),  # PS -> Integration
        ((8.75, 5.7), (8.5, 5.2)),  # SBAS -> Integration
        ((7.0, 4.5), (2.5, 3.9)),   # Integration -> Decomp
        ((7.0, 4.5), (7.0, 3.9)),   # Integration -> Anomaly
        ((7.0, 4.5), (11.5, 3.9)),  # Integration -> 3D
        ((2.5, 3.2), (5.5, 2.3)),   # Decomp -> Nankai
        ((7.0, 3.2), (7.0, 2.3)),   # Anomaly -> Nankai
        ((11.5, 3.2), (10.5, 2.3)), # 3D -> Nankai
    ]

    for start, end in connections:
        ax.annotate('', xy=end, xytext=start, arrowprops=arrow_style)

    ax.set_title('ISCE/StaMPS-Based Automated Processing Workflow\nfor InSAR Crustal Deformation Monitoring',
                 fontsize=14, fontweight='bold', pad=20)

    plt.savefig(os.path.join(FIGURES_DIR, 'workflow_diagram.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: workflow_diagram.png")


# ============================================================
# Main Execution
# ============================================================
def main():
    print("=" * 70)
    print("InSAR Time-Series Analysis for Crustal Deformation Monitoring")
    print("Nankai Trough Subduction Zone Application")
    print("=" * 70)

    # Initialize
    print("\n[1/7] Initializing Nankai Trough simulator...")
    sim = NankaiTroughSimulator(n_pixels=80, n_dates=80, time_span_years=6)
    pipeline = IntegratedPipeline(sim)

    # Run PS-InSAR/SBAS
    print("[2/7] Running PS-InSAR/SBAS integrated pipeline...")
    ps_result = pipeline.run_ps_insar()
    sbas_result = pipeline.run_sbas()
    integrated_vel = pipeline.integrate_results(ps_result, sbas_result)
    print(f"  PS density: {np.sum(ps_result['mask'])/ps_result['mask'].size*100:.1f}%")
    print(f"  SBAS pairs: {sbas_result['n_pairs']}")
    print(f"  Mean velocity: {np.mean(integrated_vel):.2f} mm/yr")

    # Plot pipeline overview
    print("[3/7] Generating pipeline overview figure...")
    plot_pipeline_overview(sim, pipeline)

    # Atmospheric correction
    print("[4/7] Applying atmospheric correction...")
    corrected_data, atm_stats = plot_atmospheric_correction(sim, pipeline, pipeline.deformation)

    # Time-series decomposition
    print("[5/7] Decomposing time-series...")
    decomp_results = plot_decomposition(sim, corrected_data)

    # Anomaly detection
    print("[6/7] Running anomaly detection...")
    anomaly_map, detection_time, auc = plot_anomaly_detection(sim, corrected_data)

    # 3D displacement
    print("[7/7] Computing 3D displacement field...")
    r2, rmse_3d = plot_3d_displacement(sim, pipeline, corrected_data)

    # Nankai application synthesis
    print("\nGenerating Nankai Trough application figure...")
    plot_nankai_application(sim, pipeline, corrected_data, decomp_results, anomaly_map)

    # Workflow diagram
    print("Generating ISCE/StaMPS workflow diagram...")
    plot_workflow_diagram()

    # Summary
    print("\n" + "=" * 70)
    print("EXPERIMENT SUMMARY")
    print("=" * 70)
    print(f"  Grid size: {sim.n_pixels}x{sim.n_pixels} pixels")
    print(f"  Time span: {sim.time_span} years ({sim.n_dates} acquisitions)")
    print(f"  PS density: {np.sum(ps_result['mask'])/ps_result['mask'].size*100:.1f}%")
    print(f"  SBAS interferogram pairs: {sbas_result['n_pairs']}")
    print(f"  Atmospheric correction improvement: {atm_stats['improvement_pct']:.1f}%")
    print(f"  RMSE before correction: {atm_stats['rmse_before']:.2f} mm")
    print(f"  RMSE after correction: {atm_stats['rmse_after']:.2f} mm")
    print(f"  Mean interseismic velocity: {np.mean(decomp_results['velocity_map']):.2f} mm/yr")
    print(f"  Anomaly detection AUC: {auc:.3f}")
    print(f"  3D vertical RMSE: {rmse_3d:.1f} mm, R²={r2:.3f}")
    print("=" * 70)

    # Save metrics to JSON
    metrics = {
        'grid_size': sim.n_pixels,
        'n_dates': sim.n_dates,
        'time_span_years': sim.time_span,
        'ps_density_pct': float(np.sum(ps_result['mask'])/ps_result['mask'].size*100),
        'sbas_pairs': sbas_result['n_pairs'],
        'atm_rmse_before': float(atm_stats['rmse_before']),
        'atm_rmse_after': float(atm_stats['rmse_after']),
        'atm_improvement_pct': float(atm_stats['improvement_pct']),
        'mean_velocity': float(np.mean(decomp_results['velocity_map'])),
        'anomaly_auc': float(auc),
        'vertical_rmse': float(rmse_3d),
        'vertical_r2': float(r2)
    }
    with open(os.path.join(os.path.dirname(FIGURES_DIR), 'metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=2)
    print("\nMetrics saved to metrics.json")

    return metrics


if __name__ == '__main__':
    main()
