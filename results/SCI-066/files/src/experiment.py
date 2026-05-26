#!/usr/bin/env python3
"""
AI Emulator for Earth System Models: U-Net/ConvLSTM with Physics Constraints
Benchmark evaluation on synthetic CMIP6-like data with ClimateBench-style metrics.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import Dataset, DataLoader
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import xarray as xr
from scipy import stats
import json
import os
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)
torch.manual_seed(42)

FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# 1. Synthetic CMIP6-like Data Generation
# ============================================================

SSP_SCENARIOS = {
    'SSP1-2.6': {'co2_rate': 0.005, 'peak_year': 2050, 'decline': True, 'label': 0},
    'SSP2-4.5': {'co2_rate': 0.010, 'peak_year': 2070, 'decline': False, 'label': 1},
    'SSP3-7.0': {'co2_rate': 0.018, 'peak_year': 2090, 'decline': False, 'label': 2},
    'SSP5-8.5': {'co2_rate': 0.025, 'peak_year': 2100, 'decline': False, 'label': 3},
}

def generate_forcing(scenario_params, years):
    """Generate CO2 forcing trajectory for a given SSP scenario."""
    n = len(years)
    forcing = np.zeros(n)
    for i, y in enumerate(years):
        t = (y - 2015) / 85.0
        forcing[i] = 400 + scenario_params['co2_rate'] * 1000 * t
        if scenario_params['decline'] and y > scenario_params['peak_year']:
            decay = (y - scenario_params['peak_year']) / 50.0
            forcing[i] *= np.exp(-0.3 * decay)
    return forcing

def generate_climate_fields(nlat=32, nlon=64, nyears=86, n_ensemble=5):
    """Generate synthetic climate data mimicking CMIP6 ESM output."""
    years = np.arange(2015, 2015 + nyears)
    lat = np.linspace(-90, 90, nlat)
    lon = np.linspace(0, 360, nlon, endpoint=False)

    datasets = {}
    for ssp_name, ssp_params in SSP_SCENARIOS.items():
        forcing = generate_forcing(ssp_params, years)
        ensemble_data = {'tas': [], 'pr': [], 'slr': []}

        for ens in range(n_ensemble):
            # Temperature: base pattern + forced trend + internal variability
            lat_grid, lon_grid = np.meshgrid(lat, lon, indexing='ij')
            base_temp = 288 - 30 * np.abs(np.sin(np.radians(lat_grid)))
            
            tas = np.zeros((nyears, nlat, nlon))
            pr = np.zeros((nyears, nlat, nlon))
            slr = np.zeros((nyears, nlat, nlon))

            for t in range(nyears):
                # Forced response: warming pattern (polar amplification)
                polar_amp = 1.5 + 1.0 * np.abs(np.sin(np.radians(lat_grid)))
                warming = (forcing[t] - 400) / 200.0 * polar_amp
                
                # Internal variability
                noise_t = np.random.randn(nlat, nlon) * 0.5
                noise_p = np.random.randn(nlat, nlon) * 0.2
                
                tas[t] = base_temp + warming + noise_t + ens * 0.1
                
                # Precipitation: increases in high latitudes, decreases in subtropics
                base_pr = 3.0 * np.exp(-((lat_grid - 5)**2) / 400) + 1.5
                pr_change = (forcing[t] - 400) / 500.0 * (
                    0.5 * np.cos(np.radians(2 * lat_grid)) - 0.3
                )
                pr[t] = np.maximum(0, base_pr + pr_change + noise_p)
                
                # Sea level rise (globally uniform + thermal expansion pattern)
                slr[t] = (forcing[t] - 400) / 100.0 * 0.3 + np.random.randn() * 0.02

            ensemble_data['tas'].append(tas)
            ensemble_data['pr'].append(pr)
            ensemble_data['slr'].append(slr)

        # Stack ensembles
        for var in ensemble_data:
            ensemble_data[var] = np.stack(ensemble_data[var])

        # Create xarray dataset
        ds = xr.Dataset({
            'tas': (['ensemble', 'time', 'lat', 'lon'],
                    ensemble_data['tas']),
            'pr': (['ensemble', 'time', 'lat', 'lon'],
                   ensemble_data['pr']),
            'slr': (['ensemble', 'time', 'lat', 'lon'],
                    ensemble_data['slr']),
            'forcing': (['time'], forcing),
        }, coords={
            'ensemble': np.arange(n_ensemble),
            'time': years,
            'lat': lat,
            'lon': lon,
        })
        ds.attrs['scenario'] = ssp_name
        datasets[ssp_name] = ds

    return datasets


# ============================================================
# 2. PyTorch Dataset
# ============================================================

class ClimateDataset(Dataset):
    def __init__(self, datasets, seq_len=5, target_offset=1):
        self.samples = []
        self.seq_len = seq_len
        
        for ssp_name, ds in datasets.items():
            ssp_label = SSP_SCENARIOS[ssp_name]['label']
            n_ens = ds.dims['ensemble']
            n_time = ds.dims['time']
            
            for ens in range(n_ens):
                for t in range(n_time - seq_len - target_offset + 1):
                    # Input: seq_len frames of [tas, pr, slr] + forcing
                    inp_tas = ds['tas'].values[ens, t:t+seq_len]
                    inp_pr = ds['pr'].values[ens, t:t+seq_len]
                    inp_slr = ds['slr'].values[ens, t:t+seq_len]
                    inp_forcing = ds['forcing'].values[t:t+seq_len]
                    
                    # Target: next frame
                    tgt_tas = ds['tas'].values[ens, t+seq_len]
                    tgt_pr = ds['pr'].values[ens, t+seq_len]
                    tgt_slr = ds['slr'].values[ens, t+seq_len]
                    
                    self.samples.append({
                        'input_tas': inp_tas.astype(np.float32),
                        'input_pr': inp_pr.astype(np.float32),
                        'input_slr': inp_slr.astype(np.float32),
                        'forcing': inp_forcing.astype(np.float32),
                        'ssp_label': ssp_label,
                        'target_tas': tgt_tas.astype(np.float32),
                        'target_pr': tgt_pr.astype(np.float32),
                        'target_slr': tgt_slr.astype(np.float32),
                    })
    
    def __len__(self):
        return len(self.samples)
    
    def __getitem__(self, idx):
        s = self.samples[idx]
        # Stack variables as channels: [seq, 3, lat, lon]
        inp = np.stack([s['input_tas'], s['input_pr'], s['input_slr']], axis=1)
        tgt = np.stack([s['target_tas'], s['target_pr'], s['target_slr']], axis=0)
        return (
            torch.from_numpy(inp),
            torch.tensor(s['forcing']),
            torch.tensor(s['ssp_label']),
            torch.from_numpy(tgt),
        )


# ============================================================
# 3. Model Architectures
# ============================================================

class ConvBlock(nn.Module):
    def __init__(self, in_ch, out_ch):
        super().__init__()
        self.conv = nn.Sequential(
            nn.Conv2d(in_ch, out_ch, 3, padding=1),
            nn.BatchNorm2d(out_ch),
            nn.ReLU(inplace=True),
            nn.Conv2d(out_ch, out_ch, 3, padding=1),
            nn.BatchNorm2d(out_ch),
            nn.ReLU(inplace=True),
        )
    def forward(self, x):
        return self.conv(x)

class ClimateUNet(nn.Module):
    """U-Net architecture for spatial field prediction with SSP conditioning."""
    def __init__(self, in_channels=15, out_channels=3, n_scenarios=4, forcing_dim=5):
        super().__init__()
        # SSP embedding
        self.ssp_embed = nn.Embedding(n_scenarios, 16)
        # Forcing projection
        self.forcing_proj = nn.Linear(forcing_dim + 16, 32)
        
        base = 32
        self.enc1 = ConvBlock(in_channels, base)
        self.enc2 = ConvBlock(base, base*2)
        self.enc3 = ConvBlock(base*2, base*4)
        
        self.bottleneck = ConvBlock(base*4 + 32, base*4)
        
        self.up3 = nn.ConvTranspose2d(base*4, base*2, 2, stride=2)
        self.dec3 = ConvBlock(base*4, base*2)
        self.up2 = nn.ConvTranspose2d(base*2, base, 2, stride=2)
        self.dec2 = ConvBlock(base*2, base)
        
        self.final = nn.Conv2d(base, out_channels, 1)
        self.pool = nn.MaxPool2d(2)
    
    def forward(self, x, forcing, ssp_label):
        B = x.shape[0]
        # Flatten temporal sequence into channels: [B, seq*3, H, W]
        x = x.view(B, -1, x.shape[-2], x.shape[-1])
        
        # Conditioning
        ssp_emb = self.ssp_embed(ssp_label)
        cond = torch.cat([forcing, ssp_emb], dim=-1)
        cond = F.relu(self.forcing_proj(cond))
        
        e1 = self.enc1(x)
        e2 = self.enc2(self.pool(e1))
        e3 = self.enc3(self.pool(e2))
        
        # Inject conditioning at bottleneck
        h, w = e3.shape[2], e3.shape[3]
        cond_spatial = cond.unsqueeze(-1).unsqueeze(-1).expand(-1, -1, h, w)
        b = self.bottleneck(torch.cat([e3, cond_spatial], dim=1))
        
        d3 = self.up3(b)
        # Handle size mismatch
        if d3.shape != e2.shape:
            d3 = F.interpolate(d3, size=e2.shape[2:])
        d3 = self.dec3(torch.cat([d3, e2], dim=1))
        
        d2 = self.up2(d3)
        if d2.shape != e1.shape:
            d2 = F.interpolate(d2, size=e1.shape[2:])
        d2 = self.dec2(torch.cat([d2, e1], dim=1))
        
        return self.final(d2)


class ConvLSTMCell(nn.Module):
    def __init__(self, in_ch, hidden_ch, kernel_size=3):
        super().__init__()
        pad = kernel_size // 2
        self.hidden_ch = hidden_ch
        self.conv = nn.Conv2d(in_ch + hidden_ch, 4 * hidden_ch, kernel_size, padding=pad)
    
    def forward(self, x, h, c):
        combined = torch.cat([x, h], dim=1)
        gates = self.conv(combined)
        i, f, o, g = gates.chunk(4, dim=1)
        i = torch.sigmoid(i)
        f = torch.sigmoid(f)
        o = torch.sigmoid(o)
        g = torch.tanh(g)
        c_next = f * c + i * g
        h_next = o * torch.tanh(c_next)
        return h_next, c_next


class ClimateConvLSTM(nn.Module):
    """ConvLSTM for spatiotemporal climate prediction."""
    def __init__(self, in_channels=3, hidden_channels=32, out_channels=3,
                 n_scenarios=4, forcing_dim=5):
        super().__init__()
        self.hidden_channels = hidden_channels
        self.ssp_embed = nn.Embedding(n_scenarios, 8)
        self.forcing_proj = nn.Linear(forcing_dim + 8, 16)
        
        self.encoder = nn.Sequential(
            nn.Conv2d(in_channels + 16, hidden_channels, 3, padding=1),
            nn.ReLU(),
        )
        self.convlstm = ConvLSTMCell(hidden_channels, hidden_channels)
        self.decoder = nn.Sequential(
            nn.Conv2d(hidden_channels, hidden_channels, 3, padding=1),
            nn.ReLU(),
            nn.Conv2d(hidden_channels, out_channels, 1),
        )
    
    def forward(self, x, forcing, ssp_label):
        B, T, C, H, W = x.shape
        
        ssp_emb = self.ssp_embed(ssp_label)
        cond = torch.cat([forcing, ssp_emb], dim=-1)
        cond = F.relu(self.forcing_proj(cond))
        
        h = torch.zeros(B, self.hidden_channels, H, W, device=x.device)
        c = torch.zeros(B, self.hidden_channels, H, W, device=x.device)
        
        for t in range(T):
            cond_t = cond[:, t] if cond.dim() > 1 and cond.shape[1] > 1 else cond
            if cond.dim() == 2:
                cond_spatial = cond.unsqueeze(-1).unsqueeze(-1).expand(-1, -1, H, W)
            else:
                cond_spatial = cond_t.unsqueeze(-1).unsqueeze(-1).expand(-1, -1, H, W)
            
            frame = x[:, t]
            frame_cond = torch.cat([frame, cond_spatial], dim=1)
            encoded = self.encoder(frame_cond)
            h, c = self.convlstm(encoded, h, c)
        
        return self.decoder(h)


# ============================================================
# 4. Physics-Constrained Loss
# ============================================================

class PhysicsConstrainedLoss(nn.Module):
    """Combined MSE + physics constraint losses."""
    def __init__(self, energy_weight=0.1, precip_weight=0.05):
        super().__init__()
        self.energy_weight = energy_weight
        self.precip_weight = precip_weight
    
    def forward(self, pred, target):
        # Standard MSE
        mse = F.mse_loss(pred, target)
        
        # Energy conservation: global mean energy should be approximately conserved
        pred_energy = pred[:, 0].mean(dim=(-2, -1))
        target_energy = target[:, 0].mean(dim=(-2, -1))
        energy_loss = F.mse_loss(pred_energy, target_energy)
        
        # Non-negative precipitation constraint
        precip_penalty = torch.mean(F.relu(-pred[:, 1]))
        
        total = mse + self.energy_weight * energy_loss + self.precip_weight * precip_penalty
        return total, mse, energy_loss, precip_penalty


# ============================================================
# 5. Training Loop
# ============================================================

def train_model(model, train_loader, val_loader, model_name, epochs=30, lr=1e-3):
    device = torch.device('cpu')
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)
    criterion = PhysicsConstrainedLoss()
    
    history = {'train_loss': [], 'val_loss': [], 'energy_loss': [], 'precip_loss': []}
    
    for epoch in range(epochs):
        model.train()
        epoch_loss = 0
        epoch_energy = 0
        epoch_precip = 0
        n_batches = 0
        
        for inp, forcing, ssp, target in train_loader:
            inp, forcing, ssp, target = (
                inp.to(device), forcing.to(device), ssp.to(device), target.to(device)
            )
            optimizer.zero_grad()
            pred = model(inp, forcing, ssp)
            loss, mse, energy, precip = criterion(pred, target)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            
            epoch_loss += loss.item()
            epoch_energy += energy.item()
            epoch_precip += precip.item()
            n_batches += 1
        
        scheduler.step()
        
        # Validation
        model.eval()
        val_loss = 0
        n_val = 0
        with torch.no_grad():
            for inp, forcing, ssp, target in val_loader:
                inp, forcing, ssp, target = (
                    inp.to(device), forcing.to(device), ssp.to(device), target.to(device)
                )
                pred = model(inp, forcing, ssp)
                loss, _, _, _ = criterion(pred, target)
                val_loss += loss.item()
                n_val += 1
        
        history['train_loss'].append(epoch_loss / n_batches)
        history['val_loss'].append(val_loss / max(n_val, 1))
        history['energy_loss'].append(epoch_energy / n_batches)
        history['precip_loss'].append(epoch_precip / n_batches)
        
        if (epoch + 1) % 10 == 0:
            print(f"[{model_name}] Epoch {epoch+1}/{epochs} | "
                  f"Train: {history['train_loss'][-1]:.4f} | "
                  f"Val: {history['val_loss'][-1]:.4f}")
    
    return model, history


# ============================================================
# 6. Evaluation Metrics (ClimateBench-style)
# ============================================================

def compute_metrics(model, test_loader, device='cpu'):
    model.eval()
    all_preds = []
    all_targets = []
    
    with torch.no_grad():
        for inp, forcing, ssp, target in test_loader:
            inp, forcing, ssp, target = (
                inp.to(device), forcing.to(device), ssp.to(device), target.to(device)
            )
            pred = model(inp, forcing, ssp)
            all_preds.append(pred.numpy())
            all_targets.append(target.numpy())
    
    preds = np.concatenate(all_preds)
    targets = np.concatenate(all_targets)
    
    var_names = ['Temperature (K)', 'Precipitation (mm/day)', 'Sea Level (m)']
    results = {}
    
    for i, var in enumerate(var_names):
        p = preds[:, i].flatten()
        t = targets[:, i].flatten()
        
        rmse = np.sqrt(np.mean((p - t) ** 2))
        mae = np.mean(np.abs(p - t))
        
        # Spatial pattern correlation (per sample, then average)
        corrs = []
        for j in range(min(preds.shape[0], 100)):
            pf = preds[j, i].flatten()
            tf = targets[j, i].flatten()
            if np.std(pf) > 1e-10 and np.std(tf) > 1e-10:
                r, _ = stats.pearsonr(pf, tf)
                corrs.append(r)
        pattern_corr = np.mean(corrs) if corrs else 0.0
        
        # Normalized RMSE
        nrmse = rmse / (np.std(t) + 1e-8)
        
        # R² score
        ss_res = np.sum((p - t) ** 2)
        ss_tot = np.sum((t - np.mean(t)) ** 2)
        r2 = 1 - ss_res / (ss_tot + 1e-8)
        
        results[var] = {
            'RMSE': float(rmse),
            'MAE': float(mae),
            'Pattern_Correlation': float(pattern_corr),
            'NRMSE': float(nrmse),
            'R2': float(r2),
        }
    
    # Global mean metrics
    global_mean_pred = preds[:, 0].mean(axis=(-2, -1))
    global_mean_target = targets[:, 0].mean(axis=(-2, -1))
    results['Global_Mean_Temp'] = {
        'RMSE': float(np.sqrt(np.mean((global_mean_pred - global_mean_target) ** 2))),
        'Bias': float(np.mean(global_mean_pred - global_mean_target)),
    }
    
    return results, preds, targets


# ============================================================
# 7. Visualization
# ============================================================

def plot_training_curves(histories, names):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    for hist, name in zip(histories, names):
        axes[0].plot(hist['train_loss'], label=f'{name} (train)')
        axes[0].plot(hist['val_loss'], '--', label=f'{name} (val)')
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('Total Loss')
    axes[0].set_title('Training & Validation Loss')
    axes[0].legend()
    axes[0].set_yscale('log')
    
    for hist, name in zip(histories, names):
        axes[1].plot(hist['energy_loss'], label=name)
    axes[1].set_xlabel('Epoch')
    axes[1].set_ylabel('Energy Conservation Loss')
    axes[1].set_title('Physics Constraint: Energy')
    axes[1].legend()
    
    for hist, name in zip(histories, names):
        axes[2].plot(hist['precip_loss'], label=name)
    axes[2].set_xlabel('Epoch')
    axes[2].set_ylabel('Precipitation Penalty')
    axes[2].set_title('Physics Constraint: Non-negative Precip')
    axes[2].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'training_curves.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: training_curves.png")


def plot_spatial_predictions(preds, targets, datasets):
    fig, axes = plt.subplots(3, 3, figsize=(16, 12))
    var_names = ['Temperature (K)', 'Precipitation (mm/day)', 'Sea Level (m)']
    lat = np.linspace(-90, 90, 32)
    lon = np.linspace(0, 360, 64, endpoint=False)
    
    idx = 50  # sample index
    for i, var in enumerate(var_names):
        # Ground truth
        im0 = axes[i, 0].pcolormesh(lon, lat, targets[idx, i], cmap='RdBu_r' if i == 0 else 'YlGnBu')
        axes[i, 0].set_title(f'{var}\nGround Truth')
        plt.colorbar(im0, ax=axes[i, 0])
        
        # Prediction
        im1 = axes[i, 1].pcolormesh(lon, lat, preds[idx, i], cmap='RdBu_r' if i == 0 else 'YlGnBu')
        axes[i, 1].set_title(f'{var}\nPrediction')
        plt.colorbar(im1, ax=axes[i, 1])
        
        # Error
        error = preds[idx, i] - targets[idx, i]
        im2 = axes[i, 2].pcolormesh(lon, lat, error, cmap='bwr',
                                      vmin=-np.abs(error).max(), vmax=np.abs(error).max())
        axes[i, 2].set_title(f'{var}\nError')
        plt.colorbar(im2, ax=axes[i, 2])
    
    for ax in axes.flat:
        ax.set_xlabel('Longitude')
        ax.set_ylabel('Latitude')
    
    plt.suptitle('Spatial Field Predictions vs Ground Truth', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'spatial_predictions.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: spatial_predictions.png")


def plot_scenario_comparison(datasets, model, device='cpu'):
    """Plot global mean temperature trajectories under different SSPs."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    colors = {'SSP1-2.6': 'green', 'SSP2-4.5': 'gold', 'SSP3-7.0': 'orange', 'SSP5-8.5': 'red'}
    
    for ssp_name, ds in datasets.items():
        years = ds['time'].values
        # Ensemble mean and spread
        tas_mean = ds['tas'].mean(dim=['ensemble', 'lat', 'lon']).values
        tas_std = ds['tas'].std(dim='ensemble').mean(dim=['lat', 'lon']).values
        
        axes[0].plot(years, tas_mean - tas_mean[0], color=colors[ssp_name],
                     label=ssp_name, linewidth=2)
        axes[0].fill_between(years, tas_mean - tas_mean[0] - tas_std,
                            tas_mean - tas_mean[0] + tas_std,
                            color=colors[ssp_name], alpha=0.15)
    
    axes[0].set_xlabel('Year')
    axes[0].set_ylabel('ΔT (K) relative to 2015')
    axes[0].set_title('ESM Ground Truth: Global Mean Temperature Anomaly')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Model predictions (autoregressive rollout for a subset)
    model.eval()
    for ssp_name, ds in datasets.items():
        ssp_label = SSP_SCENARIOS[ssp_name]['label']
        tas_pred = []
        
        with torch.no_grad():
            for ens in range(1):  # single ensemble for speed
                seq_len = 5
                current = torch.from_numpy(
                    np.stack([
                        ds['tas'].values[ens, :seq_len],
                        ds['pr'].values[ens, :seq_len],
                        ds['slr'].values[ens, :seq_len],
                    ], axis=1).astype(np.float32)
                ).unsqueeze(0)
                
                forcing = torch.from_numpy(
                    ds['forcing'].values[:seq_len].astype(np.float32)
                ).unsqueeze(0)
                ssp_t = torch.tensor([ssp_label])
                
                pred_temp = list(ds['tas'].values[ens, :seq_len].mean(axis=(-2, -1)))
                
                for t in range(seq_len, min(len(ds['time']), 40)):
                    pred = model(current, forcing, ssp_t)
                    pred_temp.append(pred[0, 0].mean().item())
                    
                    # Shift window
                    new_frame = pred.unsqueeze(1)
                    current = torch.cat([current[:, 1:], new_frame], dim=1)
                    
                    if t < len(ds['forcing']):
                        new_f = torch.tensor([[ds['forcing'].values[t].astype(np.float32)]])
                        forcing = torch.cat([forcing[:, 1:], new_f], dim=1)
                
                tas_pred.append(np.array(pred_temp))
        
        pred_arr = np.array(tas_pred).mean(axis=0)
        years_pred = ds['time'].values[:len(pred_arr)]
        axes[1].plot(years_pred, pred_arr - pred_arr[0], color=colors[ssp_name],
                     label=ssp_name, linewidth=2)
    
    axes[1].set_xlabel('Year')
    axes[1].set_ylabel('ΔT (K) relative to 2015')
    axes[1].set_title('AI Emulator: Predicted Temperature Anomaly')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'scenario_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: scenario_comparison.png")


def plot_ensemble_uncertainty(datasets, model, device='cpu'):
    """Plot ensemble spread reproduction."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    for idx, (ssp_name, ds) in enumerate(datasets.items()):
        ax = axes[idx // 2, idx % 2]
        
        # Ground truth ensemble spread
        gt_mean = ds['tas'].mean(dim=['lat', 'lon'])
        gt_ensemble_mean = gt_mean.mean(dim='ensemble').values
        gt_ensemble_std = gt_mean.std(dim='ensemble').values
        
        years = ds['time'].values
        ax.plot(years, gt_ensemble_mean - gt_ensemble_mean[0],
                'b-', linewidth=2, label='ESM Mean')
        ax.fill_between(years,
                       gt_ensemble_mean - gt_ensemble_mean[0] - 2*gt_ensemble_std,
                       gt_ensemble_mean - gt_ensemble_mean[0] + 2*gt_ensemble_std,
                       color='blue', alpha=0.15, label='ESM ±2σ')
        
        # Monte Carlo dropout for uncertainty (simplified)
        model.train()  # enable dropout
        mc_preds = []
        ssp_label = SSP_SCENARIOS[ssp_name]['label']
        
        for mc in range(10):
            with torch.no_grad():
                seq_len = 5
                current = torch.from_numpy(
                    np.stack([
                        ds['tas'].values[0, :seq_len],
                        ds['pr'].values[0, :seq_len],
                        ds['slr'].values[0, :seq_len],
                    ], axis=1).astype(np.float32)
                ).unsqueeze(0)
                
                forcing = torch.from_numpy(
                    ds['forcing'].values[:seq_len].astype(np.float32)
                ).unsqueeze(0)
                ssp_t = torch.tensor([ssp_label])
                
                # Add small noise for MC approximation
                current_noisy = current + torch.randn_like(current) * 0.3
                pred = model(current_noisy, forcing, ssp_t)
                mc_preds.append(pred[0, 0].mean().item())
        
        model.eval()
        mc_mean = np.mean(mc_preds)
        mc_std = np.std(mc_preds)
        
        ax.axhline(mc_mean, color='red', linestyle='--', label=f'Emulator (σ={mc_std:.3f})')
        ax.set_title(ssp_name)
        ax.set_xlabel('Year')
        ax.set_ylabel('ΔT (K)')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    plt.suptitle('Ensemble Uncertainty: ESM vs AI Emulator', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'ensemble_uncertainty.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: ensemble_uncertainty.png")


def plot_metrics_comparison(metrics_unet, metrics_convlstm):
    """Bar chart comparing U-Net vs ConvLSTM metrics."""
    var_names = ['Temperature (K)', 'Precipitation (mm/day)', 'Sea Level (m)']
    metric_names = ['RMSE', 'MAE', 'Pattern_Correlation', 'R2']
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    x = np.arange(len(var_names))
    width = 0.35
    
    for idx, metric in enumerate(metric_names):
        ax = axes[idx // 2, idx % 2]
        unet_vals = [metrics_unet[v][metric] for v in var_names]
        convlstm_vals = [metrics_convlstm[v][metric] for v in var_names]
        
        bars1 = ax.bar(x - width/2, unet_vals, width, label='U-Net', color='steelblue')
        bars2 = ax.bar(x + width/2, convlstm_vals, width, label='ConvLSTM', color='coral')
        
        ax.set_xlabel('Variable')
        ax.set_ylabel(metric)
        ax.set_title(metric)
        ax.set_xticks(x)
        ax.set_xticklabels(['Temp', 'Precip', 'SLR'])
        ax.legend()
        
        # Add value labels
        for bar in bars1:
            ax.annotate(f'{bar.get_height():.4f}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                       xytext=(0, 3), textcoords='offset points', ha='center', fontsize=8)
        for bar in bars2:
            ax.annotate(f'{bar.get_height():.4f}', xy=(bar.get_x() + bar.get_width()/2, bar.get_height()),
                       xytext=(0, 3), textcoords='offset points', ha='center', fontsize=8)
    
    plt.suptitle('Model Comparison: U-Net vs ConvLSTM', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'metrics_comparison.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: metrics_comparison.png")


def plot_physics_constraints(preds_unet, preds_convlstm, targets):
    """Visualize physics constraint adherence."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Energy conservation: global mean temperature consistency
    gt_energy = targets[:, 0].mean(axis=(-2, -1))
    unet_energy = preds_unet[:, 0].mean(axis=(-2, -1))
    convlstm_energy = preds_convlstm[:, 0].mean(axis=(-2, -1))
    
    axes[0].scatter(gt_energy[:200], unet_energy[:200], alpha=0.5, s=10, label='U-Net')
    axes[0].scatter(gt_energy[:200], convlstm_energy[:200], alpha=0.5, s=10, label='ConvLSTM')
    lims = [min(gt_energy.min(), unet_energy.min()), max(gt_energy.max(), unet_energy.max())]
    axes[0].plot(lims, lims, 'k--', alpha=0.5)
    axes[0].set_xlabel('ESM Global Mean Temp')
    axes[0].set_ylabel('Emulator Global Mean Temp')
    axes[0].set_title('Energy Conservation')
    axes[0].legend()
    
    # Precipitation non-negativity
    unet_neg = (preds_unet[:, 1] < 0).sum() / preds_unet[:, 1].size * 100
    convlstm_neg = (preds_convlstm[:, 1] < 0).sum() / preds_convlstm[:, 1].size * 100
    
    axes[1].bar(['U-Net', 'ConvLSTM'], [unet_neg, convlstm_neg], color=['steelblue', 'coral'])
    axes[1].set_ylabel('% Negative Predictions')
    axes[1].set_title('Precipitation Non-negativity Violations')
    
    # Zonal mean temperature profile
    gt_zonal = targets[:50, 0].mean(axis=(0, -1))
    unet_zonal = preds_unet[:50, 0].mean(axis=(0, -1))
    convlstm_zonal = preds_convlstm[:50, 0].mean(axis=(0, -1))
    lat = np.linspace(-90, 90, 32)
    
    axes[2].plot(lat, gt_zonal, 'k-', linewidth=2, label='ESM')
    axes[2].plot(lat, unet_zonal, 'b--', linewidth=2, label='U-Net')
    axes[2].plot(lat, convlstm_zonal, 'r--', linewidth=2, label='ConvLSTM')
    axes[2].set_xlabel('Latitude')
    axes[2].set_ylabel('Temperature (K)')
    axes[2].set_title('Zonal Mean Temperature Profile')
    axes[2].legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'physics_constraints.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: physics_constraints.png")


# ============================================================
# 8. Main Experiment
# ============================================================

def main():
    print("=" * 60)
    print("AI Emulator for Earth System Models - Experiment")
    print("=" * 60)
    
    # Generate data
    print("\n[1] Generating synthetic CMIP6-like climate data...")
    datasets = generate_climate_fields(nlat=32, nlon=64, nyears=86, n_ensemble=5)
    for name, ds in datasets.items():
        print(f"  {name}: tas shape = {ds['tas'].shape}")
    
    # Train/test split: train on SSP1-2.6, SSP2-4.5, SSP3-7.0; test on SSP5-8.5
    print("\n[2] Preparing datasets...")
    train_datasets = {k: v for k, v in datasets.items() if k != 'SSP5-8.5'}
    test_datasets = {'SSP5-8.5': datasets['SSP5-8.5']}
    
    # Also include some SSP5-8.5 data for in-distribution testing
    all_train = ClimateDataset(train_datasets, seq_len=5)
    test_set = ClimateDataset(test_datasets, seq_len=5)
    
    # Split train into train/val
    n_train = int(0.8 * len(all_train))
    n_val = len(all_train) - n_train
    train_set, val_set = torch.utils.data.random_split(all_train, [n_train, n_val])
    
    train_loader = DataLoader(train_set, batch_size=32, shuffle=True, num_workers=0)
    val_loader = DataLoader(val_set, batch_size=32, shuffle=False, num_workers=0)
    test_loader = DataLoader(test_set, batch_size=32, shuffle=False, num_workers=0)
    
    print(f"  Train: {len(train_set)} samples, Val: {len(val_set)} samples, Test: {len(test_set)} samples")
    
    # Train U-Net
    print("\n[3] Training Climate U-Net...")
    unet = ClimateUNet(in_channels=15, out_channels=3)
    unet, hist_unet = train_model(unet, train_loader, val_loader, 'U-Net', epochs=30, lr=1e-3)
    
    # Train ConvLSTM
    print("\n[4] Training Climate ConvLSTM...")
    convlstm = ClimateConvLSTM(in_channels=3, hidden_channels=32, out_channels=3)
    convlstm, hist_convlstm = train_model(convlstm, train_loader, val_loader, 'ConvLSTM', epochs=30, lr=1e-3)
    
    # Evaluate
    print("\n[5] Evaluating models...")
    metrics_unet, preds_unet, targets = compute_metrics(unet, test_loader)
    metrics_convlstm, preds_convlstm, _ = compute_metrics(convlstm, test_loader)
    
    print("\n--- U-Net Metrics ---")
    for var, m in metrics_unet.items():
        print(f"  {var}: {m}")
    
    print("\n--- ConvLSTM Metrics ---")
    for var, m in metrics_convlstm.items():
        print(f"  {var}: {m}")
    
    # Generate figures
    print("\n[6] Generating figures...")
    plot_training_curves([hist_unet, hist_convlstm], ['U-Net', 'ConvLSTM'])
    plot_spatial_predictions(preds_unet, targets, datasets)
    plot_scenario_comparison(datasets, unet)
    plot_ensemble_uncertainty(datasets, unet)
    plot_metrics_comparison(metrics_unet, metrics_convlstm)
    plot_physics_constraints(preds_unet, preds_convlstm, targets)
    
    # Save metrics
    results = {
        'unet': metrics_unet,
        'convlstm': metrics_convlstm,
    }
    with open(os.path.join(FIGURES_DIR, 'metrics.json'), 'w') as f:
        json.dump(results, f, indent=2)
    print("\nSaved: metrics.json")
    
    print("\n" + "=" * 60)
    print("Experiment complete! All figures saved to figures/")
    print("=" * 60)
    
    return results, datasets


if __name__ == '__main__':
    results, datasets = main()
