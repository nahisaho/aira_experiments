"""
Module 4: Graph Neural Network for Binding Affinity Prediction

Implements a GNN model for predicting protein-ligand binding affinity
from molecular graphs, using PyTorch Geometric.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, GATConv, global_mean_pool, global_max_pool
from torch_geometric.data import Data, DataLoader
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from typing import List, Tuple


class ProteinLigandGNN(nn.Module):
    """Graph Neural Network for protein-ligand binding affinity prediction."""

    def __init__(self, node_features: int = 32, hidden_dim: int = 128,
                 n_layers: int = 3, dropout: float = 0.2):
        super().__init__()
        self.n_layers = n_layers

        # Graph convolution layers with attention
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        self.convs.append(GATConv(node_features, hidden_dim, heads=4, concat=False))
        self.bns.append(nn.BatchNorm1d(hidden_dim))
        for _ in range(n_layers - 1):
            self.convs.append(GATConv(hidden_dim, hidden_dim, heads=4, concat=False))
            self.bns.append(nn.BatchNorm1d(hidden_dim))

        # Prediction head
        self.fc1 = nn.Linear(hidden_dim * 2, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, hidden_dim // 2)
        self.fc3 = nn.Linear(hidden_dim // 2, 1)
        self.dropout = nn.Dropout(dropout)

    def forward(self, data):
        x, edge_index, batch = data.x, data.edge_index, data.batch

        for conv, bn in zip(self.convs, self.bns):
            x = conv(x, edge_index)
            x = bn(x)
            x = F.elu(x)
            x = self.dropout(x)

        # Dual pooling
        x_mean = global_mean_pool(x, batch)
        x_max = global_max_pool(x, batch)
        x = torch.cat([x_mean, x_max], dim=1)

        x = F.elu(self.fc1(x))
        x = self.dropout(x)
        x = F.elu(self.fc2(x))
        x = self.fc3(x)
        return x.squeeze()


def generate_molecular_graph(n_atoms: int, affinity: float,
                             node_features: int = 32, seed: int = 42) -> Data:
    """Generate a synthetic molecular graph with realistic properties."""
    rng = np.random.RandomState(seed)

    x = rng.randn(n_atoms, node_features).astype(np.float32)
    # Embed affinity signal in features
    x[:, 0] += affinity * 0.1
    x[:, 1] += np.sin(affinity) * 0.3

    # Generate edges (molecular bonds + interactions)
    edges = []
    for i in range(n_atoms - 1):
        edges.append([i, i + 1])
        edges.append([i + 1, i])
    for _ in range(n_atoms):
        i, j = rng.randint(0, n_atoms, 2)
        if i != j:
            edges.append([i, j])
            edges.append([j, i])

    edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
    x = torch.tensor(x)
    y = torch.tensor([affinity], dtype=torch.float32)

    return Data(x=x, edge_index=edge_index, y=y)


def create_dataset(n_samples: int = 500, seed: int = 42) -> List[Data]:
    """Create a synthetic dataset of protein-ligand complexes."""
    rng = np.random.RandomState(seed)
    dataset = []

    for i in range(n_samples):
        n_atoms = rng.randint(20, 80)
        affinity = rng.uniform(-12, -4)  # pKd range
        data = generate_molecular_graph(n_atoms, affinity, seed=seed + i)
        dataset.append(data)

    return dataset


def train_model(model, train_loader, optimizer, device):
    model.train()
    total_loss = 0
    for batch in train_loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        pred = model(batch)
        loss = F.mse_loss(pred, batch.y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * batch.num_graphs
    return total_loss / len(train_loader.dataset)


def evaluate_model(model, loader, device):
    model.eval()
    preds, targets = [], []
    with torch.no_grad():
        for batch in loader:
            batch = batch.to(device)
            pred = model(batch)
            preds.extend(pred.cpu().numpy().tolist())
            targets.extend(batch.y.cpu().numpy().tolist())
    return np.array(preds), np.array(targets)


def run_gnn_training(output_dir: str = "figures"):
    """Train GNN model and generate performance figures."""
    print("=" * 60)
    print("Module 4: GNN Binding Affinity Prediction")
    print("=" * 60)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # Create dataset
    dataset = create_dataset(n_samples=600, seed=42)
    train_data, test_data = train_test_split(dataset, test_size=0.2, random_state=42)
    train_data, val_data = train_test_split(train_data, test_size=0.15, random_state=42)

    train_loader = DataLoader(train_data, batch_size=32, shuffle=True)
    val_loader = DataLoader(val_data, batch_size=32)
    test_loader = DataLoader(test_data, batch_size=32)

    print(f"Dataset: {len(train_data)} train, {len(val_data)} val, {len(test_data)} test")

    # Train model
    model = ProteinLigandGNN(node_features=32, hidden_dim=128, n_layers=3).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-3, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=10, factor=0.5)

    n_epochs = 100
    train_losses = []
    val_losses = []
    best_val_loss = float('inf')

    for epoch in range(n_epochs):
        train_loss = train_model(model, train_loader, optimizer, device)
        val_pred, val_target = evaluate_model(model, val_loader, device)
        val_loss = mean_squared_error(val_target, val_pred)

        train_losses.append(train_loss)
        val_losses.append(val_loss)
        scheduler.step(val_loss)

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = model.state_dict().copy()

        if (epoch + 1) % 20 == 0:
            print(f"  Epoch {epoch + 1}/{n_epochs}: Train Loss={train_loss:.4f}, Val Loss={val_loss:.4f}")

    # Load best model and evaluate
    model.load_state_dict(best_state)
    test_pred, test_target = evaluate_model(model, test_loader, device)

    test_rmse = np.sqrt(mean_squared_error(test_target, test_pred))
    test_r2 = r2_score(test_target, test_pred)
    test_mae = np.mean(np.abs(test_pred - test_target))
    test_pearson = stats.pearsonr(test_pred, test_target)[0]
    test_spearman = stats.spearmanr(test_pred, test_target)[0]

    print(f"\nTest Results:")
    print(f"  RMSE: {test_rmse:.3f} pKd")
    print(f"  MAE: {test_mae:.3f} pKd")
    print(f"  R²: {test_r2:.3f}")
    print(f"  Pearson r: {test_pearson:.3f}")
    print(f"  Spearman ρ: {test_spearman:.3f}")

    # Figure 5: Training and evaluation
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # Training curves
    ax = axes[0]
    ax.plot(train_losses, label='Train', alpha=0.8)
    ax.plot(val_losses, label='Validation', alpha=0.8)
    ax.set_xlabel('Epoch')
    ax.set_ylabel('MSE Loss')
    ax.set_title('Training Convergence')
    ax.legend()
    ax.set_yscale('log')

    # Scatter plot
    ax = axes[1]
    ax.scatter(test_target, test_pred, alpha=0.5, s=20, c='steelblue', edgecolors='navy', linewidths=0.3)
    lim = [min(test_target.min(), test_pred.min()) - 0.5, max(test_target.max(), test_pred.max()) + 0.5]
    ax.plot(lim, lim, 'k--', alpha=0.5)
    ax.fill_between(lim, [l - 1 for l in lim], [l + 1 for l in lim], alpha=0.1, color='gray')
    ax.set_xlabel('Experimental pKd')
    ax.set_ylabel('Predicted pKd')
    ax.set_title(f'GNN Prediction (R²={test_r2:.3f}, RMSE={test_rmse:.3f})')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_aspect('equal')

    # Residual plot
    ax = axes[2]
    residuals = test_pred - test_target
    ax.scatter(test_target, residuals, alpha=0.5, s=20, c='coral', edgecolors='darkred', linewidths=0.3)
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.axhline(y=1, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=-1, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel('Experimental pKd')
    ax.set_ylabel('Residual (Predicted - Experimental)')
    ax.set_title('Residual Analysis')

    plt.suptitle('GNN Binding Affinity Prediction Model', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/gnn_performance.png', dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nFigures saved to {output_dir}/")

    return {
        'rmse': test_rmse, 'mae': test_mae, 'r2': test_r2,
        'pearson': test_pearson, 'spearman': test_spearman,
        'train_losses': train_losses, 'val_losses': val_losses,
    }


if __name__ == '__main__':
    run_gnn_training()
