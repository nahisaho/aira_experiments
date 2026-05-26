"""
Synthetic data generation for mHealth neurodegenerative disease biomarker study.
Generates realistic sensor data for gait, voice, touch, and longitudinal monitoring.
"""

import numpy as np
import pandas as pd
from pathlib import Path

np.random.seed(42)

DATA_DIR = Path(__file__).parent.parent / "data"
DATA_DIR.mkdir(exist_ok=True)


def generate_gait_data(n_subjects=200, n_steps=500, sampling_rate=100):
    """Generate synthetic accelerometer and gyroscope gait data.
    
    Healthy subjects have regular gait patterns; PD subjects show
    increased variability, asymmetry, and reduced stride length.
    """
    records = []
    labels = []
    
    n_pd = n_subjects // 2
    n_healthy = n_subjects - n_pd
    
    for i in range(n_subjects):
        is_pd = i < n_pd
        severity = np.random.uniform(0.3, 1.0) if is_pd else 0.0
        t = np.linspace(0, n_steps / sampling_rate, n_steps)
        
        # Base gait frequency (~2 Hz for healthy)
        freq = 2.0 - 0.5 * severity + np.random.normal(0, 0.1)
        
        # Accelerometer (x=forward, y=lateral, z=vertical)
        ax = np.sin(2 * np.pi * freq * t) * (1.0 - 0.3 * severity)
        ay = np.sin(2 * np.pi * freq * t + np.pi/3) * (0.5 + 0.3 * severity)
        az = np.abs(np.sin(2 * np.pi * freq * t)) * (9.8 + 0.5 * severity)
        
        # Add PD-specific irregularity
        if is_pd:
            jitter = np.random.normal(0, 0.2 * severity, n_steps)
            ax += jitter
            ay += np.random.normal(0, 0.15 * severity, n_steps)
            az += np.random.normal(0, 0.3 * severity, n_steps)
            # Freezing episodes
            n_freezes = np.random.poisson(2 * severity)
            for _ in range(n_freezes):
                start = np.random.randint(0, n_steps - 50)
                ax[start:start+50] *= 0.1
                ay[start:start+50] *= 0.1
        else:
            ax += np.random.normal(0, 0.05, n_steps)
            ay += np.random.normal(0, 0.05, n_steps)
            az += np.random.normal(0, 0.1, n_steps)
        
        # Gyroscope
        gx = np.gradient(ax) * sampling_rate + np.random.normal(0, 0.1 * (1 + severity), n_steps)
        gy = np.gradient(ay) * sampling_rate + np.random.normal(0, 0.1 * (1 + severity), n_steps)
        gz = np.gradient(az) * sampling_rate + np.random.normal(0, 0.05 * (1 + severity), n_steps)
        
        # Extract features
        features = extract_gait_features(ax, ay, az, gx, gy, gz, sampling_rate)
        features['subject_id'] = i
        features['severity'] = severity
        records.append(features)
        labels.append(1 if is_pd else 0)
    
    df = pd.DataFrame(records)
    df['label'] = labels
    df.to_csv(DATA_DIR / "gait_features.csv", index=False)
    return df


def extract_gait_features(ax, ay, az, gx, gy, gz, sr):
    """Extract gait features from raw IMU signals."""
    mag_acc = np.sqrt(ax**2 + ay**2 + az**2)
    mag_gyro = np.sqrt(gx**2 + gy**2 + gz**2)
    
    features = {
        # Acceleration features
        'acc_mean': np.mean(mag_acc),
        'acc_std': np.std(mag_acc),
        'acc_range': np.ptp(mag_acc),
        'acc_rms': np.sqrt(np.mean(mag_acc**2)),
        'acc_skew': float(pd.Series(mag_acc).skew()),
        'acc_kurt': float(pd.Series(mag_acc).kurtosis()),
        
        # Step regularity (autocorrelation-based)
        'step_regularity': compute_step_regularity(mag_acc, sr),
        'stride_regularity': compute_stride_regularity(mag_acc, sr),
        
        # Asymmetry
        'lateral_asymmetry': np.abs(np.mean(ay[:len(ay)//2]) - np.mean(ay[len(ay)//2:])),
        
        # Gyroscope features
        'gyro_mean': np.mean(mag_gyro),
        'gyro_std': np.std(mag_gyro),
        'gyro_range': np.ptp(mag_gyro),
        'gyro_rms': np.sqrt(np.mean(mag_gyro**2)),
        
        # Frequency domain
        'dominant_freq': get_dominant_freq(mag_acc, sr),
        'spectral_entropy': compute_spectral_entropy(mag_acc, sr),
        
        # Jerk (derivative of acceleration)
        'jerk_mean': np.mean(np.abs(np.diff(mag_acc) * sr)),
        'jerk_std': np.std(np.abs(np.diff(mag_acc) * sr)),
    }
    return features


def compute_step_regularity(signal, sr):
    ac = np.correlate(signal - np.mean(signal), signal - np.mean(signal), mode='full')
    ac = ac[len(ac)//2:]
    ac /= ac[0] + 1e-10
    # Look for first peak after ~0.4s (step time)
    start = int(0.3 * sr)
    end = int(0.7 * sr)
    if end > len(ac):
        end = len(ac)
    if start >= end:
        return 0.0
    return float(np.max(ac[start:end]))


def compute_stride_regularity(signal, sr):
    ac = np.correlate(signal - np.mean(signal), signal - np.mean(signal), mode='full')
    ac = ac[len(ac)//2:]
    ac /= ac[0] + 1e-10
    start = int(0.8 * sr)
    end = int(1.4 * sr)
    if end > len(ac):
        end = len(ac)
    if start >= end:
        return 0.0
    return float(np.max(ac[start:end]))


def get_dominant_freq(signal, sr):
    fft_vals = np.abs(np.fft.rfft(signal - np.mean(signal)))
    freqs = np.fft.rfftfreq(len(signal), 1.0/sr)
    return float(freqs[np.argmax(fft_vals[1:]) + 1])


def compute_spectral_entropy(signal, sr):
    fft_vals = np.abs(np.fft.rfft(signal - np.mean(signal)))
    psd = fft_vals**2
    psd_norm = psd / (np.sum(psd) + 1e-10)
    psd_norm = psd_norm[psd_norm > 0]
    return float(-np.sum(psd_norm * np.log2(psd_norm + 1e-10)))


def generate_voice_data(n_subjects=150, n_sessions=10):
    """Generate synthetic voice feature data for ALS progression monitoring.
    
    Features: jitter, shimmer, MFCC coefficients, HNR, F0.
    ALS subjects show progressive deterioration over sessions.
    """
    records = []
    n_als = n_subjects // 2
    
    for i in range(n_subjects):
        is_als = i < n_als
        base_severity = np.random.uniform(0.1, 0.5) if is_als else 0.0
        progression_rate = np.random.uniform(0.02, 0.08) if is_als else 0.0
        
        for session in range(n_sessions):
            severity = min(base_severity + progression_rate * session, 1.0)
            
            # Fundamental frequency
            f0 = 150 - 30 * severity + np.random.normal(0, 10)
            
            # Jitter (frequency perturbation)
            jitter = 0.5 + 2.5 * severity + np.random.normal(0, 0.3)
            jitter = max(jitter, 0.1)
            
            # Shimmer (amplitude perturbation)
            shimmer = 1.0 + 4.0 * severity + np.random.normal(0, 0.5)
            shimmer = max(shimmer, 0.2)
            
            # HNR (harmonic-to-noise ratio, decreases with severity)
            hnr = 25 - 10 * severity + np.random.normal(0, 2)
            
            # MFCC features (13 coefficients)
            mfcc = np.zeros(13)
            mfcc[0] = 12.0 - 3.0 * severity + np.random.normal(0, 1)
            for k in range(1, 13):
                mfcc[k] = np.random.normal(0, 1) * (1 + 0.5 * severity)
            
            record = {
                'subject_id': i,
                'session': session,
                'is_als': int(is_als),
                'severity': severity,
                'f0': f0,
                'jitter': jitter,
                'shimmer': shimmer,
                'hnr': hnr,
            }
            for k in range(13):
                record[f'mfcc_{k}'] = mfcc[k]
            
            records.append(record)
    
    df = pd.DataFrame(records)
    df.to_csv(DATA_DIR / "voice_features.csv", index=False)
    return df


def generate_touch_data(n_subjects=180, n_sessions=8):
    """Generate synthetic touchscreen interaction data for cognitive decline detection.
    
    Features: tap accuracy, reaction time, swipe velocity, typing patterns.
    Cognitive decline leads to increased latency and reduced accuracy.
    """
    records = []
    n_impaired = n_subjects // 3
    n_mci = n_subjects // 3
    n_healthy = n_subjects - n_impaired - n_mci
    
    for i in range(n_subjects):
        if i < n_impaired:
            group = 'impaired'
            cog_score = np.random.uniform(0.5, 1.0)
        elif i < n_impaired + n_mci:
            group = 'mci'
            cog_score = np.random.uniform(0.2, 0.5)
        else:
            group = 'healthy'
            cog_score = np.random.uniform(0.0, 0.15)
        
        for session in range(n_sessions):
            decline = cog_score * (1 + 0.02 * session)
            
            # Reaction time (ms)
            reaction_time = 300 + 400 * decline + np.random.normal(0, 50)
            
            # Tap accuracy (proportion)
            tap_accuracy = max(0.5, 0.98 - 0.3 * decline + np.random.normal(0, 0.03))
            
            # Swipe velocity (pixels/s)
            swipe_velocity = 800 - 300 * decline + np.random.normal(0, 80)
            
            # Double-tap interval variability
            dt_variability = 20 + 60 * decline + np.random.normal(0, 10)
            
            # Typing speed (chars/min)
            typing_speed = max(10, 60 - 30 * decline + np.random.normal(0, 5))
            
            # Error rate
            error_rate = 0.02 + 0.15 * decline + np.random.normal(0, 0.02)
            error_rate = max(0, min(1, error_rate))
            
            # Pressure variability
            pressure_var = 0.05 + 0.2 * decline + np.random.normal(0, 0.03)
            
            # Trail-making time (simulated)
            trail_time = 30 + 80 * decline + np.random.normal(0, 10)
            
            records.append({
                'subject_id': i,
                'session': session,
                'group': group,
                'cog_score': cog_score,
                'reaction_time': reaction_time,
                'tap_accuracy': tap_accuracy,
                'swipe_velocity': swipe_velocity,
                'dt_variability': dt_variability,
                'typing_speed': typing_speed,
                'error_rate': error_rate,
                'pressure_var': pressure_var,
                'trail_time': trail_time,
            })
    
    df = pd.DataFrame(records)
    df.to_csv(DATA_DIR / "touch_features.csv", index=False)
    return df


def generate_longitudinal_data(n_subjects=100, n_timepoints=52):
    """Generate longitudinal multimodal data for change point detection.
    
    Simulates weekly measurements over 1 year with disease onset at random points.
    """
    records = []
    
    for i in range(n_subjects):
        has_onset = np.random.random() < 0.6
        onset_week = np.random.randint(10, 40) if has_onset else n_timepoints + 10
        
        for t in range(n_timepoints):
            is_post_onset = t >= onset_week
            time_since_onset = max(0, t - onset_week)
            decline = 0.05 * time_since_onset if is_post_onset else 0.0
            
            # Composite signals
            gait_score = max(0, 1.0 - decline + np.random.normal(0, 0.05))
            voice_score = max(0, 1.0 - 0.8 * decline + np.random.normal(0, 0.04))
            touch_score = max(0, 1.0 - 0.6 * decline + np.random.normal(0, 0.06))
            clinical_score = max(0, 30 - 2.0 * time_since_onset + np.random.normal(0, 1)) if is_post_onset else 30 + np.random.normal(0, 1)
            
            records.append({
                'subject_id': i,
                'week': t,
                'has_onset': int(has_onset),
                'true_onset_week': onset_week if has_onset else -1,
                'gait_score': gait_score,
                'voice_score': voice_score,
                'touch_score': touch_score,
                'clinical_score': clinical_score,
            })
    
    df = pd.DataFrame(records)
    df.to_csv(DATA_DIR / "longitudinal_data.csv", index=False)
    return df


if __name__ == "__main__":
    print("Generating gait data...")
    gait_df = generate_gait_data()
    print(f"  Gait: {len(gait_df)} subjects, {gait_df.shape[1]} features")
    
    print("Generating voice data...")
    voice_df = generate_voice_data()
    print(f"  Voice: {len(voice_df)} records")
    
    print("Generating touch data...")
    touch_df = generate_touch_data()
    print(f"  Touch: {len(touch_df)} records")
    
    print("Generating longitudinal data...")
    long_df = generate_longitudinal_data()
    print(f"  Longitudinal: {len(long_df)} records")
    
    print("Data generation complete.")
