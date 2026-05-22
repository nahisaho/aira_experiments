"""
ConvLSTM architecture for temporal evolution of climate fields.

Captures spatiotemporal dynamics by combining convolutional operations
with LSTM gating for sequential climate field prediction.
"""

import torch
import torch.nn as nn
from typing import Optional


class ConvLSTMCell(nn.Module):
    """
    Convolutional LSTM cell.

    Replaces fully-connected gates with convolutions to preserve
    spatial structure in hidden/cell states.
    """

    def __init__(self, input_dim: int, hidden_dim: int, kernel_size: int = 3):
        super().__init__()
        self.hidden_dim = hidden_dim
        padding = kernel_size // 2

        self.gates = nn.Conv2d(
            input_dim + hidden_dim, 4 * hidden_dim,
            kernel_size, padding=padding, bias=True
        )
        self.layer_norm = nn.GroupNorm(4, 4 * hidden_dim)

    def forward(self, x: torch.Tensor,
                h: Optional[torch.Tensor] = None,
                c: Optional[torch.Tensor] = None):
        B, _, H, W = x.shape
        if h is None:
            h = torch.zeros(B, self.hidden_dim, H, W, device=x.device)
        if c is None:
            c = torch.zeros(B, self.hidden_dim, H, W, device=x.device)

        combined = torch.cat([x, h], dim=1)
        gates = self.layer_norm(self.gates(combined))

        i, f, o, g = gates.chunk(4, dim=1)
        i = torch.sigmoid(i)
        f = torch.sigmoid(f)
        o = torch.sigmoid(o)
        g = torch.tanh(g)

        c_next = f * c + i * g
        h_next = o * torch.tanh(c_next)

        return h_next, c_next


class ConvLSTMPredictor(nn.Module):
    """
    Multi-layer ConvLSTM for temporal climate field prediction.

    Processes a sequence of climate fields and predicts the next
    timestep(s), maintaining spatial coherence through convolutional
    recurrence.
    """

    def __init__(self, input_dim: int = 3, hidden_dims: list = None,
                 n_layers: int = 3, kernel_size: int = 3,
                 out_channels: int = 3):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [64, 64, 64]
        assert len(hidden_dims) == n_layers

        self.n_layers = n_layers
        self.hidden_dims = hidden_dims

        cells = []
        for i in range(n_layers):
            in_dim = input_dim if i == 0 else hidden_dims[i - 1]
            cells.append(ConvLSTMCell(in_dim, hidden_dims[i], kernel_size))
        self.cells = nn.ModuleList(cells)

        self.output_conv = nn.Sequential(
            nn.Conv2d(hidden_dims[-1], hidden_dims[-1] // 2, 3, padding=1),
            nn.GELU(),
            nn.Conv2d(hidden_dims[-1] // 2, out_channels, 1),
        )

    def forward(self, x_seq: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x_seq: (B, T, C, H, W) - sequence of climate fields

        Returns:
            output: (B, C, H, W) - predicted next timestep
        """
        B, T, C, H, W = x_seq.shape

        h_states = [None] * self.n_layers
        c_states = [None] * self.n_layers

        for t in range(T):
            x_t = x_seq[:, t]
            for layer_idx in range(self.n_layers):
                h, c = self.cells[layer_idx](
                    x_t if layer_idx == 0 else h_states[layer_idx - 1],
                    h_states[layer_idx],
                    c_states[layer_idx],
                )
                h_states[layer_idx] = h
                c_states[layer_idx] = c

        return self.output_conv(h_states[-1])

    def predict_sequence(self, x_seq: torch.Tensor,
                         n_steps: int = 10) -> torch.Tensor:
        """Autoregressive multi-step prediction."""
        predictions = []
        current_seq = x_seq.clone()

        for _ in range(n_steps):
            pred = self.forward(current_seq)
            predictions.append(pred.unsqueeze(1))
            current_seq = torch.cat([current_seq[:, 1:], pred.unsqueeze(1)], dim=1)

        return torch.cat(predictions, dim=1)
