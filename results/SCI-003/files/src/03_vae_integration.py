import json
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, TensorDataset

warnings.filterwarnings('ignore')
np.random.seed(42)
torch.manual_seed(42)

try:
    import scanpy as sc
except Exception:
    sc = None

BASE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = BASE_DIR / 'data'
RESULTS_DIR = BASE_DIR / 'results'
FIG_DIR = BASE_DIR / 'figures'
LOG_PATH = BASE_DIR / 'logs' / 'process-log.jsonl'
SKILL_NAME = 'co-scientist-multi-omics'
DEVICE = torch.device('cuda' if torch.cuda.is_available() else 'cpu')


def log_event(phase, event_type, status='ok', files_written=None, extra=None):
    entry = {
        'timestamp': datetime.utcnow().isoformat(),
        'phase': phase,
        'event_type': event_type,
        'actor': 'co-scientist',
        'skill': SKILL_NAME,
        'status': status,
        'files_written': files_written or [],
    }
    if extra:
        entry.update(extra)
    with open(LOG_PATH, 'a', encoding='utf-8') as fh:
        fh.write(json.dumps(entry, ensure_ascii=False) + '\n')


class ModalityEncoder(nn.Module):
    def __init__(self, input_dim, h1, h2):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, h1),
            nn.ReLU(),
            nn.Linear(h1, h2),
            nn.ReLU(),
        )

    def forward(self, x):
        return self.net(x)


class ModalityDecoder(nn.Module):
    def __init__(self, latent_dim, h1, h2, output_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(latent_dim, h2),
            nn.ReLU(),
            nn.Linear(h2, h1),
            nn.ReLU(),
            nn.Linear(h1, output_dim),
        )

    def forward(self, z):
        return self.net(z)


class MultimodalVAE(nn.Module):
    def __init__(self):
        super().__init__()
        self.rna_encoder = ModalityEncoder(2000, 512, 256)
        self.atac_encoder = ModalityEncoder(5000, 512, 256)
        self.meth_encoder = ModalityEncoder(1000, 256, 128)
        self.fc_mu = nn.Linear(256 + 256 + 128, 32)
        self.fc_logvar = nn.Linear(256 + 256 + 128, 32)
        self.rna_decoder = ModalityDecoder(32, 512, 256, 2000)
        self.atac_decoder = ModalityDecoder(32, 512, 256, 5000)
        self.meth_decoder = ModalityDecoder(32, 256, 128, 1000)

    def encode(self, rna_x, atac_x, meth_x):
        h = torch.cat([
            self.rna_encoder(rna_x),
            self.atac_encoder(atac_x),
            self.meth_encoder(meth_x),
        ], dim=1)
        return self.fc_mu(h), self.fc_logvar(h)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z):
        return self.rna_decoder(z), self.atac_decoder(z), self.meth_decoder(z)

    def forward(self, rna_x, atac_x, meth_x):
        mu, logvar = self.encode(rna_x, atac_x, meth_x)
        z = self.reparameterize(mu, logvar)
        recon = self.decode(z)
        return recon, mu, logvar


def main():
    RESULTS_DIR.mkdir(exist_ok=True)
    FIG_DIR.mkdir(exist_ok=True)
    log_event('vae_integration', 'module_started')

    if sc is None:
        raise ImportError('scanpy is required to load AnnData objects')

    rna = sc.read_h5ad(DATA_DIR / 'rna_processed.h5ad')
    atac = sc.read_h5ad(DATA_DIR / 'atac_processed.h5ad')
    meth = sc.read_h5ad(DATA_DIR / 'methylation_processed.h5ad')

    rna_x = np.asarray(rna.X, dtype=np.float32)
    atac_x = np.log1p(np.asarray(atac.layers['counts'], dtype=np.float32))
    meth_x = np.asarray(meth.layers['beta'], dtype=np.float32)

    rna_x = StandardScaler().fit_transform(rna_x).astype(np.float32)
    atac_x = StandardScaler().fit_transform(atac_x).astype(np.float32)
    meth_x = StandardScaler().fit_transform(meth_x).astype(np.float32)

    dataset = TensorDataset(
        torch.tensor(rna_x),
        torch.tensor(atac_x),
        torch.tensor(meth_x),
    )
    loader = DataLoader(dataset, batch_size=64, shuffle=True)

    model = MultimodalVAE().to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

    history = []
    for epoch in range(50):
        model.train()
        total_loss = 0.0
        for batch_rna, batch_atac, batch_meth in loader:
            batch_rna = batch_rna.to(DEVICE)
            batch_atac = batch_atac.to(DEVICE)
            batch_meth = batch_meth.to(DEVICE)
            optimizer.zero_grad()
            (recon_rna, recon_atac, recon_meth), mu, logvar = model(batch_rna, batch_atac, batch_meth)
            recon_loss = (
                F.mse_loss(recon_rna, batch_rna, reduction='mean') +
                F.mse_loss(recon_atac, batch_atac, reduction='mean') +
                F.mse_loss(recon_meth, batch_meth, reduction='mean')
            )
            kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
            loss = recon_loss + 1e-3 * kl
            loss.backward()
            optimizer.step()
            total_loss += float(loss.item())
        history.append(total_loss / max(len(loader), 1))

    model.eval()
    with torch.no_grad():
        mu, _ = model.encode(
            torch.tensor(rna_x, device=DEVICE),
            torch.tensor(atac_x, device=DEVICE),
            torch.tensor(meth_x, device=DEVICE),
        )
    latent = mu.cpu().numpy()

    latent_df = pd.DataFrame(latent, index=rna.obs_names, columns=[f'latent_{i+1}' for i in range(latent.shape[1])])
    latent_df['cell_type'] = rna.obs['cell_type'].values
    latent_path = RESULTS_DIR / 'vae_latent.csv'
    latent_df.to_csv(latent_path)

    model_path = RESULTS_DIR / 'vae_model.pt'
    torch.save(model.state_dict(), model_path)

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(np.arange(1, 51), history, color='royalblue', lw=2)
    ax.set_title('Multimodal VAE training loss')
    ax.set_xlabel('Epoch')
    ax.set_ylabel('ELBO loss')
    loss_fig = FIG_DIR / 'vae_training_loss.png'
    plt.savefig(loss_fig, dpi=150, bbox_inches='tight')
    plt.close(fig)

    files_written = [str(model_path), str(latent_path), str(loss_fig)]

    if sc is not None:
        latent_adata = sc.AnnData(latent, obs=rna.obs.copy())
        sc.pp.neighbors(latent_adata, n_neighbors=20)
        sc.tl.umap(latent_adata)
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.scatter(latent_adata.obsm['X_umap'][:, 0], latent_adata.obsm['X_umap'][:, 1], c=pd.Categorical(latent_adata.obs['cell_type']).codes, cmap='tab10', s=12)
        ax.set_title('VAE latent UMAP')
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        umap_fig = FIG_DIR / 'vae_latent_umap.png'
        plt.savefig(umap_fig, dpi=150, bbox_inches='tight')
        plt.close(fig)
        latent_adata.write_h5ad(DATA_DIR / 'vae_latent.h5ad')
        files_written.extend([str(umap_fig), str(DATA_DIR / 'vae_latent.h5ad')])

    log_event('vae_integration', 'module_completed', files_written=files_written, extra={'summary': {'final_loss': float(history[-1]), 'latent_dim': 32}})
    return {
        'final_loss': float(history[-1]),
        'files_written': files_written,
    }


if __name__ == '__main__':
    main()
