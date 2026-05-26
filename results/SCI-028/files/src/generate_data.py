"""
Synthetic Tokamak Plasma Data Generator
Generates realistic time-series data mimicking JET and KSTAR diagnostics
for disruption prediction experiments.
"""
import numpy as np
import os

np.random.seed(42)

def generate_plasma_shot(shot_id, device='JET', disruption=False, n_timesteps=500, dt=0.002):
    """Generate a single plasma shot with multiple diagnostic channels."""
    t = np.arange(n_timesteps) * dt  # time in seconds

    # Base plasma parameters (device-dependent scaling)
    if device == 'JET':
        Ip_base, ne_base, Te_base, beta_base = 3.0, 4.0, 3.5, 1.8
        B_tor = 3.45
        R0, a = 2.96, 1.25
    elif device == 'KSTAR':
        Ip_base, ne_base, Te_base, beta_base = 0.6, 3.0, 2.0, 1.2
        B_tor = 3.5
        R0, a = 1.8, 0.5
    else:  # ITER-like
        Ip_base, ne_base, Te_base, beta_base = 15.0, 10.0, 8.0, 2.5
        B_tor = 5.3
        R0, a = 6.2, 2.0

    # Ramp-up, flat-top, ramp-down phases
    ramp_up = np.minimum(t / 0.2, 1.0)
    ramp_down = np.where(t > 0.8, np.maximum(1.0 - (t - 0.8) / 0.2, 0.0), 1.0)
    envelope = ramp_up * ramp_down

    # Plasma current Ip (MA)
    Ip = Ip_base * envelope + np.random.normal(0, 0.01 * Ip_base, n_timesteps)

    # Electron density ne (1e19 m^-3)
    ne = ne_base * envelope + np.random.normal(0, 0.05 * ne_base, n_timesteps)

    # Electron temperature Te (keV)
    Te = Te_base * envelope + np.random.normal(0, 0.03 * Te_base, n_timesteps)

    # Normalized beta
    beta_N = beta_base * envelope + np.random.normal(0, 0.05, n_timesteps)

    # Internal inductance li
    li = 0.9 + 0.2 * np.sin(2 * np.pi * t / 0.5) + np.random.normal(0, 0.02, n_timesteps)

    # Radiated power fraction
    P_rad_frac = 0.3 + 0.1 * np.sin(2 * np.pi * t / 0.3) + np.random.normal(0, 0.02, n_timesteps)

    # Safety factor q95
    q95 = 3.0 + 0.5 * np.sin(2 * np.pi * t / 0.4) + np.random.normal(0, 0.1, n_timesteps)

    # Greenwald fraction
    n_G = ne / (Ip_base * 1e6 / (np.pi * a**2) * 1e-20)

    # Locked mode amplitude (mT)
    lm_amp = np.abs(np.random.normal(0, 0.05, n_timesteps))

    # MHD activity (Mirnov signal RMS)
    mirnov_rms = np.abs(np.random.normal(0, 0.1, n_timesteps))

    # Stored energy (MJ)
    W_mhd = 0.5 * beta_N * Ip * B_tor * a**2 / R0

    if disruption:
        # Add disruption precursors
        t_disrupt = 0.6 + np.random.uniform(-0.1, 0.1)
        precursor_start = t_disrupt - 0.15  # ~150ms before disruption

        for i, ti in enumerate(t):
            if ti > precursor_start:
                decay = np.exp(3.0 * (ti - precursor_start))
                # Radiation spike
                P_rad_frac[i] += 0.3 * decay / np.exp(3.0 * 0.15)
                # Density increase (Greenwald limit approach)
                ne[i] *= 1.0 + 0.5 * decay / np.exp(3.0 * 0.15)
                # Temperature drop
                Te[i] *= max(0.1, 1.0 - 0.8 * decay / np.exp(3.0 * 0.15))
                # Locked mode growth
                lm_amp[i] += 2.0 * decay / np.exp(3.0 * 0.15)
                # MHD activity increase
                mirnov_rms[i] += 1.5 * decay / np.exp(3.0 * 0.15)
                # Beta collapse
                beta_N[i] *= max(0.1, 1.0 - 0.7 * decay / np.exp(3.0 * 0.15))
                # q95 drop
                q95[i] -= 1.0 * decay / np.exp(3.0 * 0.15)

            if ti > t_disrupt:
                # Thermal quench
                Te[i] *= 0.01
                Ip[i] *= max(0.0, 1.0 - 10.0 * (ti - t_disrupt))
                W_mhd[i] *= 0.01

    # Labels: 1 if within 300ms of disruption, 0 otherwise
    labels = np.zeros(n_timesteps)
    if disruption:
        t_disrupt_idx = int(t_disrupt / dt)
        warning_window = int(0.3 / dt)  # 300ms
        start_idx = max(0, t_disrupt_idx - warning_window)
        labels[start_idx:t_disrupt_idx] = 1.0

    features = np.column_stack([Ip, ne, Te, beta_N, li, P_rad_frac, q95, n_G, lm_amp, mirnov_rms, W_mhd])
    feature_names = ['Ip', 'ne', 'Te', 'beta_N', 'li', 'P_rad_frac', 'q95', 'n_G', 'lm_amp', 'mirnov_rms', 'W_mhd']

    return {
        'shot_id': shot_id,
        'device': device,
        'disruption': disruption,
        't': t,
        'features': features,
        'feature_names': feature_names,
        'labels': labels,
    }


def add_tearing_mode(shot_data, mode='NTM', onset_time=0.4):
    """Add tearing mode / NTM signatures to a plasma shot."""
    t = shot_data['t']
    features = shot_data['features'].copy()
    n = len(t)
    dt = t[1] - t[0]
    onset_idx = int(onset_time / dt)

    # Create TM/NTM labels
    tm_labels = np.zeros(n)

    for i in range(onset_idx, n):
        progress = (t[i] - onset_time) / (t[-1] - onset_time + 1e-6)
        if mode == 'NTM':
            # NTM: slower growth, island width evolution
            island_width = 0.05 * (1 - np.exp(-5 * progress))
            features[i, 9] += 0.8 * island_width * 20  # mirnov
            features[i, 3] -= 0.3 * island_width  # beta_N drop
            features[i, 2] *= (1.0 - 0.2 * island_width)  # Te reduction
            freq = 8000 * (1 - 0.3 * progress)  # frequency chirp
            features[i, 9] += 0.5 * np.sin(2 * np.pi * freq * t[i]) * island_width
        else:
            # Classical TM
            growth_rate = 0.03 * np.exp(2 * progress)
            features[i, 9] += growth_rate * 5
            features[i, 3] -= growth_rate * 0.5

        tm_labels[i] = 1.0

    shot_data['features'] = features
    shot_data['tm_labels'] = tm_labels
    shot_data['tm_type'] = mode
    return shot_data


def generate_dataset(n_shots=200, device='JET', disrupt_frac=0.3, tm_frac=0.2):
    """Generate a full dataset of plasma shots."""
    shots = []
    n_disrupt = int(n_shots * disrupt_frac)
    n_tm = int(n_shots * tm_frac)

    for i in range(n_shots):
        disruption = i < n_disrupt
        shot = generate_plasma_shot(f"{device}_{i:05d}", device=device, disruption=disruption)

        if i >= n_disrupt and i < n_disrupt + n_tm:
            mode = 'NTM' if np.random.random() > 0.4 else 'TM'
            shot = add_tearing_mode(shot, mode=mode)
        elif disruption and np.random.random() > 0.5:
            shot = add_tearing_mode(shot, mode='NTM')

        shots.append(shot)

    return shots


def prepare_sequences(shots, seq_len=50, stride=10):
    """Convert shots into sliding window sequences for ML training."""
    X_list, y_list, y_tm_list = [], [], []

    for shot in shots:
        features = shot['features']
        labels = shot['labels']
        tm_labels = shot.get('tm_labels', np.zeros(len(labels)))
        n = len(labels)

        # Normalize features per-shot
        mu = features.mean(axis=0, keepdims=True)
        std = features.std(axis=0, keepdims=True) + 1e-8
        features_norm = (features - mu) / std

        for start in range(0, n - seq_len, stride):
            end = start + seq_len
            X_list.append(features_norm[start:end])
            y_list.append(labels[end - 1])
            y_tm_list.append(tm_labels[end - 1])

    X = np.array(X_list, dtype=np.float32)
    y = np.array(y_list, dtype=np.float32)
    y_tm = np.array(y_tm_list, dtype=np.float32)
    return X, y, y_tm


if __name__ == '__main__':
    os.makedirs('data', exist_ok=True)

    for device in ['JET', 'KSTAR']:
        print(f"Generating {device} dataset...")
        shots = generate_dataset(n_shots=200, device=device)
        X, y, y_tm = prepare_sequences(shots)
        np.savez(f'data/{device}_data.npz', X=X, y=y, y_tm=y_tm)
        print(f"  {device}: X={X.shape}, y={y.shape} (disrupt={y.sum():.0f}), y_tm={y_tm.sum():.0f} TM samples")

    # ITER synthetic (small)
    print("Generating ITER synthetic dataset...")
    shots_iter = generate_dataset(n_shots=30, device='ITER', disrupt_frac=0.2)
    X_iter, y_iter, y_tm_iter = prepare_sequences(shots_iter)
    np.savez('data/ITER_data.npz', X=X_iter, y=y_iter, y_tm=y_tm_iter)
    print(f"  ITER: X={X_iter.shape}")

    print("Data generation complete.")
