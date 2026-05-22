"""
gnn_processor.py — Graph Neural Network processor for atmospheric message passing.

Implements the core GNN processor with multi-scale message passing,
inspired by GraphCast's encode-process-decode architecture.
"""

import torch
import torch.nn as nn
from torch_geometric.nn import MessagePassing
from torch_geometric.data import Data
from typing import Optional, Tuple


class AtmosphericMessagePassing(MessagePassing):
    """
    Custom message passing layer for atmospheric fields.
    Uses edge features (distance, relative position) and
    multi-head attention for directional information flow.
    """

    def __init__(self, d_model: int = 256, n_heads: int = 4, dropout: float = 0.1):
        super().__init__(aggr='mean')
        self.d_model = d_model
        self.n_heads = n_heads
        head_dim = d_model // n_heads

        self.q_proj = nn.Linear(d_model, d_model)
        self.k_proj = nn.Linear(d_model, d_model)
        self.v_proj = nn.Linear(d_model, d_model)
        self.o_proj = nn.Linear(d_model, d_model)

        self.edge_mlp = nn.Sequential(
            nn.Linear(3, d_model),  # edge features: dx, dy, dz
            nn.SiLU(),
            nn.Linear(d_model, d_model),
        )

        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)

        self.ffn = nn.Sequential(
            nn.Linear(d_model, d_model * 4),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(d_model * 4, d_model),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor,
                edge_attr: Optional[torch.Tensor] = None) -> torch.Tensor:
        # Pre-norm residual connection for message passing
        h = self.norm1(x)
        h = self.propagate(edge_index, x=h, edge_attr=edge_attr)
        x = x + h

        # Pre-norm residual connection for FFN
        h = self.norm2(x)
        h = self.ffn(h)
        x = x + h

        return x

    def message(self, x_i: torch.Tensor, x_j: torch.Tensor,
                edge_attr: Optional[torch.Tensor] = None) -> torch.Tensor:
        """Compute attention-weighted messages."""
        q = self.q_proj(x_i)
        k = self.k_proj(x_j)
        v = self.v_proj(x_j)

        # Scale dot-product attention
        d_k = q.shape[-1] / self.n_heads
        attn = (q * k).sum(dim=-1, keepdim=True) / (d_k ** 0.5)
        attn = torch.sigmoid(attn)

        msg = attn * v

        # Incorporate edge features if available
        if edge_attr is not None:
            edge_emb = self.edge_mlp(edge_attr)
            msg = msg + edge_emb

        return self.o_proj(msg)


class GNNProcessor(nn.Module):
    """
    Multi-layer GNN processor with residual connections.

    This is the core "processor" in the encode-process-decode architecture.
    It performs multiple rounds of message passing on the multi-scale mesh
    to propagate information across spatial scales and distances.
    """

    def __init__(
        self,
        d_model: int = 256,
        n_layers: int = 8,
        n_heads: int = 4,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.layers = nn.ModuleList([
            AtmosphericMessagePassing(d_model, n_heads, dropout)
            for _ in range(n_layers)
        ])
        self.final_norm = nn.LayerNorm(d_model)

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor,
                edge_attr: Optional[torch.Tensor] = None) -> torch.Tensor:
        """
        Args:
            x: (N, d_model) node features
            edge_index: (2, E) edge indices
            edge_attr: (E, 3) optional edge features

        Returns:
            (N, d_model) processed node features
        """
        for layer in self.layers:
            x = layer(x, edge_index, edge_attr)
        return self.final_norm(x)
