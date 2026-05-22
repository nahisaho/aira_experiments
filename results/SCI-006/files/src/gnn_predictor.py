"""
Module 4: Graph Neural Network (GNN) for Binding Affinity Prediction

Implements a heterogeneous GNN architecture for predicting protein-ligand
binding affinities from molecular graphs, incorporating:
- Molecular graph construction from SMILES/3D structures
- Protein-ligand interaction graph with distance-based edges
- Message-passing neural network with attention
- Multi-task learning for affinity + selectivity
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional


@dataclass
class AtomFeatures:
    """Node features for molecular graph atoms."""
    atomic_number: int
    degree: int
    formal_charge: int
    hybridization: int  # sp=1, sp2=2, sp3=3
    aromaticity: bool
    num_hydrogens: int
    in_ring: bool
    ring_size: int = 0
    gasteiger_charge: float = 0.0
    
    def to_vector(self) -> List[float]:
        return [
            self.atomic_number / 53.0,  # Normalize by I
            self.degree / 6.0,
            (self.formal_charge + 2) / 4.0,
            self.hybridization / 3.0,
            float(self.aromaticity),
            self.num_hydrogens / 4.0,
            float(self.in_ring),
            self.ring_size / 8.0,
            self.gasteiger_charge,
        ]


@dataclass
class BondFeatures:
    """Edge features for molecular graph bonds."""
    bond_type: int  # single=1, double=2, triple=3, aromatic=4
    is_conjugated: bool
    is_in_ring: bool
    stereo: int = 0
    
    def to_vector(self) -> List[float]:
        return [
            self.bond_type / 4.0,
            float(self.is_conjugated),
            float(self.is_in_ring),
            self.stereo / 6.0,
        ]


@dataclass
class InteractionEdge:
    """Edge representing protein-ligand interaction."""
    protein_atom_idx: int
    ligand_atom_idx: int
    distance: float
    interaction_type: str  # hydrogen_bond, hydrophobic, pi_stacking, etc.
    
    def to_vector(self) -> List[float]:
        type_encoding = {
            "hydrogen_bond": [1, 0, 0, 0, 0],
            "hydrophobic": [0, 1, 0, 0, 0],
            "pi_stacking": [0, 0, 1, 0, 0],
            "salt_bridge": [0, 0, 0, 1, 0],
            "halogen_bond": [0, 0, 0, 0, 1],
        }
        enc = type_encoding.get(self.interaction_type, [0] * 5)
        return [self.distance / 10.0] + enc


@dataclass
class MolecularGraph:
    """Graph representation of a molecule."""
    node_features: np.ndarray  # (n_atoms, n_node_features)
    edge_index: np.ndarray     # (2, n_edges)
    edge_features: np.ndarray  # (n_edges, n_edge_features)
    n_atoms: int = 0
    n_bonds: int = 0


@dataclass
class ProteinLigandGraph:
    """Heterogeneous graph for protein-ligand complex."""
    protein_graph: MolecularGraph
    ligand_graph: MolecularGraph
    interaction_edges: np.ndarray   # (2, n_interactions)
    interaction_features: np.ndarray  # (n_interactions, n_interaction_features)


@dataclass
class GNNArchitecture:
    """GNN model architecture specification."""
    # Encoder
    node_feature_dim: int = 9
    edge_feature_dim: int = 4
    interaction_feature_dim: int = 6
    
    # Message passing
    hidden_dim: int = 128
    n_message_passing_steps: int = 6
    attention_heads: int = 4
    dropout: float = 0.1
    
    # Readout
    readout_type: str = "attention_weighted"  # sum, mean, attention_weighted
    readout_hidden_dim: int = 256
    
    # Output
    n_output_tasks: int = 1  # 1 for single affinity, >1 for multi-task
    output_dim: int = 1


@dataclass
class TrainingConfig:
    """Training configuration for GNN model."""
    batch_size: int = 32
    learning_rate: float = 1e-4
    weight_decay: float = 1e-5
    epochs: int = 200
    patience: int = 20
    lr_schedule: str = "cosine_warmup"
    warmup_epochs: int = 10
    
    # Data
    train_split: float = 0.8
    val_split: float = 0.1
    test_split: float = 0.1
    
    # Augmentation
    use_conformer_augmentation: bool = True
    n_conformers: int = 5
    
    # Loss
    loss_function: str = "huber"  # mse, huber, smooth_l1
    huber_delta: float = 1.0


@dataclass
class GNNPrediction:
    """Single prediction result."""
    ligand_id: str
    predicted_pki: float
    experimental_pki: Optional[float] = None
    uncertainty: float = 0.0
    attention_weights: Optional[Dict[int, float]] = None


@dataclass
class GNNEvaluation:
    """Model evaluation metrics."""
    rmse: float = 0.0
    mae: float = 0.0
    r_squared: float = 0.0
    pearson_r: float = 0.0
    spearman_rho: float = 0.0
    kendall_tau: float = 0.0
    
    # Per-target metrics
    per_target_rmse: Dict[str, float] = field(default_factory=dict)
    
    # Calibration
    mean_uncertainty: float = 0.0
    calibration_error: float = 0.0


def generate_pytorch_geometric_model() -> str:
    """Generate PyTorch Geometric GNN model code."""
    return '''
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import (
    GATv2Conv, TransformerConv, global_attention_pool,
    HeteroConv, Linear
)
from torch_geometric.data import HeteroData


class ProteinLigandGNN(nn.Module):
    """
    Heterogeneous Graph Neural Network for binding affinity prediction.
    
    Architecture:
    1. Node embedding layers for protein and ligand atoms
    2. Intra-molecular message passing (protein-protein, ligand-ligand)
    3. Inter-molecular message passing (protein-ligand interactions)
    4. Attention-weighted readout
    5. Multi-layer prediction head
    """
    
    def __init__(self, config: dict):
        super().__init__()
        hidden = config["hidden_dim"]
        heads = config["attention_heads"]
        dropout = config["dropout"]
        n_layers = config["n_message_passing_steps"]
        
        # Node embeddings
        self.protein_embed = nn.Sequential(
            nn.Linear(config["protein_node_dim"], hidden),
            nn.LayerNorm(hidden),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        self.ligand_embed = nn.Sequential(
            nn.Linear(config["ligand_node_dim"], hidden),
            nn.LayerNorm(hidden),
            nn.ReLU(),
            nn.Dropout(dropout)
        )
        
        # Message passing layers
        self.conv_layers = nn.ModuleList()
        self.norms = nn.ModuleList()
        
        for _ in range(n_layers):
            conv_dict = {
                # Intra-molecular
                ("protein", "bond", "protein"): GATv2Conv(
                    hidden, hidden // heads, heads=heads,
                    edge_dim=config["edge_dim"], dropout=dropout
                ),
                ("ligand", "bond", "ligand"): GATv2Conv(
                    hidden, hidden // heads, heads=heads,
                    edge_dim=config["edge_dim"], dropout=dropout
                ),
                # Inter-molecular
                ("protein", "interacts", "ligand"): TransformerConv(
                    hidden, hidden // heads, heads=heads,
                    edge_dim=config["interaction_dim"], dropout=dropout
                ),
                ("ligand", "interacts", "protein"): TransformerConv(
                    hidden, hidden // heads, heads=heads,
                    edge_dim=config["interaction_dim"], dropout=dropout
                ),
            }
            self.conv_layers.append(HeteroConv(conv_dict, aggr="sum"))
            self.norms.append(nn.ModuleDict({
                "protein": nn.LayerNorm(hidden),
                "ligand": nn.LayerNorm(hidden),
            }))
        
        # Attention pooling
        gate_protein = nn.Sequential(nn.Linear(hidden, 1))
        gate_ligand = nn.Sequential(nn.Linear(hidden, 1))
        self.pool_protein = lambda x, batch: global_attention_pool(x, batch, gate_protein)
        self.pool_ligand = lambda x, batch: global_attention_pool(x, batch, gate_ligand)
        
        # Prediction head
        self.predictor = nn.Sequential(
            nn.Linear(hidden * 2, config["readout_hidden_dim"]),
            nn.LayerNorm(config["readout_hidden_dim"]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(config["readout_hidden_dim"], config["readout_hidden_dim"] // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(config["readout_hidden_dim"] // 2, 1)
        )
        
        # Uncertainty estimation (MC Dropout or evidential)
        self.uncertainty_head = nn.Sequential(
            nn.Linear(hidden * 2, 64),
            nn.ReLU(),
            nn.Linear(64, 1),
            nn.Softplus()
        )
    
    def forward(self, data: HeteroData):
        # Embed nodes
        x_prot = self.protein_embed(data["protein"].x)
        x_lig = self.ligand_embed(data["ligand"].x)
        
        x_dict = {"protein": x_prot, "ligand": x_lig}
        
        # Message passing with residual connections
        for conv, norm in zip(self.conv_layers, self.norms):
            x_new = conv(x_dict, data.edge_index_dict, data.edge_attr_dict)
            x_dict = {
                key: norm[key](x_new[key] + x_dict[key])
                for key in x_dict
            }
        
        # Pool
        prot_repr = self.pool_protein(
            x_dict["protein"], data["protein"].batch
        )
        lig_repr = self.pool_ligand(
            x_dict["ligand"], data["ligand"].batch
        )
        
        # Concatenate and predict
        complex_repr = torch.cat([prot_repr, lig_repr], dim=-1)
        affinity = self.predictor(complex_repr)
        uncertainty = self.uncertainty_head(complex_repr)
        
        return affinity, uncertainty


class AffinityLoss(nn.Module):
    """Evidential deep learning loss for uncertainty-aware regression."""
    
    def __init__(self, coeff=0.01):
        super().__init__()
        self.coeff = coeff
    
    def forward(self, pred, target, uncertainty):
        # Huber loss
        diff = pred - target
        huber = torch.where(
            diff.abs() < 1.0,
            0.5 * diff ** 2,
            diff.abs() - 0.5
        )
        # NLL with learned uncertainty
        nll = 0.5 * torch.log(uncertainty) + 0.5 * diff ** 2 / uncertainty
        return huber.mean() + self.coeff * nll.mean()
'''


def simulate_gnn_training(
    n_epochs: int = 200,
    seed: int = 42
) -> Dict[str, List[float]]:
    """Simulate GNN training curves for demonstration."""
    rng = np.random.RandomState(seed)
    
    train_losses = []
    val_losses = []
    val_rmses = []
    val_r2s = []
    
    for epoch in range(n_epochs):
        # Exponential decay with noise
        t = epoch / n_epochs
        base_train = 2.5 * np.exp(-4 * t) + 0.3
        base_val = 2.5 * np.exp(-3.5 * t) + 0.45
        
        train_loss = base_train + rng.normal(0, 0.05)
        val_loss = base_val + rng.normal(0, 0.08)
        
        # Add slight overfitting after epoch 150
        if epoch > 150:
            val_loss += 0.002 * (epoch - 150)
        
        val_rmse = np.sqrt(max(0.01, val_loss)) * 0.7
        val_r2 = min(0.95, max(0, 1 - val_loss / 2.5 + rng.normal(0, 0.02)))
        
        train_losses.append(max(0.01, train_loss))
        val_losses.append(max(0.01, val_loss))
        val_rmses.append(max(0.1, val_rmse))
        val_r2s.append(val_r2)
    
    return {
        "epoch": list(range(n_epochs)),
        "train_loss": train_losses,
        "val_loss": val_losses,
        "val_rmse": val_rmses,
        "val_r2": val_r2s,
    }


def simulate_gnn_predictions(
    n_compounds: int = 200,
    seed: int = 42
) -> List[GNNPrediction]:
    """Generate simulated GNN predictions for demonstration."""
    rng = np.random.RandomState(seed)
    
    predictions = []
    for i in range(n_compounds):
        exp_pki = rng.uniform(4, 10)
        noise = rng.normal(0, 0.5)
        pred_pki = exp_pki + noise
        unc = abs(noise) * 0.5 + rng.uniform(0.1, 0.3)
        
        predictions.append(GNNPrediction(
            ligand_id=f"CMPD_{i+1:04d}",
            predicted_pki=float(pred_pki),
            experimental_pki=float(exp_pki),
            uncertainty=float(unc),
        ))
    
    return predictions


def evaluate_predictions(predictions: List[GNNPrediction]) -> GNNEvaluation:
    """Evaluate GNN predictions against experimental values."""
    pred = np.array([p.predicted_pki for p in predictions])
    exp = np.array([p.experimental_pki for p in predictions
                    if p.experimental_pki is not None])
    pred_matched = np.array([p.predicted_pki for p in predictions
                             if p.experimental_pki is not None])
    
    errors = pred_matched - exp
    rmse = float(np.sqrt(np.mean(errors ** 2)))
    mae = float(np.mean(np.abs(errors)))
    
    corr = np.corrcoef(pred_matched, exp)[0, 1]
    r2 = float(corr ** 2)
    pearson_r = float(corr)
    
    # Spearman rank correlation
    from scipy.stats import spearmanr, kendalltau
    try:
        spearman_rho = float(spearmanr(pred_matched, exp).correlation)
        kendall_tau_val = float(kendalltau(pred_matched, exp).correlation)
    except Exception:
        spearman_rho = float(np.corrcoef(
            np.argsort(np.argsort(pred_matched)),
            np.argsort(np.argsort(exp))
        )[0, 1])
        kendall_tau_val = 0.0
    
    uncertainties = np.array([p.uncertainty for p in predictions])
    
    return GNNEvaluation(
        rmse=rmse,
        mae=mae,
        r_squared=r2,
        pearson_r=pearson_r,
        spearman_rho=spearman_rho,
        kendall_tau=kendall_tau_val,
        mean_uncertainty=float(np.mean(uncertainties)),
        calibration_error=float(abs(rmse - np.mean(uncertainties))),
    )
