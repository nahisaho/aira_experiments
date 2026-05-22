"""
Module 1: Multi-scale Fourier Feature Embedding for PINNs.

Addresses spectral bias in neural networks by mapping inputs through
random Fourier features, enabling learning of high-frequency solutions.

References:
- Tancik et al., "Fourier Features Let Networks Learn High Frequency
  Functions in Low Dimensional Domains" (NeurIPS 2020)
- Wang et al., "On the eigenvector bias of Fourier feature networks" (2021)
"""

import numpy as np
import torch
import torch.nn as nn
from typing import Optional, List, Tuple


class FourierFeatureEmbedding(nn.Module):
    """Random Fourier Feature mapping for multi-scale PINN inputs."""

    def __init__(
        self,
        input_dim: int,
        num_frequencies: int = 128,
        sigma_values: Optional[List[float]] = None,
        trainable: bool = False,
    ):
        super().__init__()
        self.input_dim = input_dim
        self.num_frequencies = num_frequencies

        if sigma_values is None:
            sigma_values = [1.0, 10.0, 100.0]
        self.sigma_values = sigma_values
        self.num_scales = len(sigma_values)
        self.output_dim = 2 * num_frequencies * self.num_scales

        B_list = []
        for sigma in sigma_values:
            B = torch.randn(input_dim, num_frequencies) * sigma
            B_list.append(B)

        B_all = torch.cat(B_list, dim=1)
        if trainable:
            self.B = nn.Parameter(B_all)
        else:
            self.register_buffer("B", B_all)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        projection = 2.0 * np.pi * x @ self.B
        return torch.cat([torch.sin(projection), torch.cos(projection)], dim=-1)


class MultiScalePINN(nn.Module):
    """PINN with multi-scale Fourier feature embedding."""

    def __init__(
        self,
        input_dim: int = 2,
        output_dim: int = 1,
        hidden_dims: List[int] = [256, 256, 256, 256],
        num_frequencies: int = 128,
        sigma_values: List[float] = [1.0, 10.0, 100.0],
        activation: str = "tanh",
    ):
        super().__init__()
        self.embedding = FourierFeatureEmbedding(
            input_dim, num_frequencies, sigma_values
        )

        act_fn = {"tanh": nn.Tanh, "gelu": nn.GELU, "silu": nn.SiLU}[activation]

        layers = []
        in_dim = self.embedding.output_dim
        for h in hidden_dims:
            layers.append(nn.Linear(in_dim, h))
            layers.append(act_fn())
            in_dim = h
        layers.append(nn.Linear(in_dim, output_dim))
        self.net = nn.Sequential(*layers)

        self._init_weights()

    def _init_weights(self):
        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        features = self.embedding(x)
        return self.net(features)


class AdaptiveFourierPINN(nn.Module):
    """PINN with learnable Fourier feature scales (adaptive multi-scale)."""

    def __init__(
        self,
        input_dim: int = 2,
        output_dim: int = 1,
        hidden_dims: List[int] = [256, 256, 256, 256],
        num_frequencies: int = 64,
        num_scales: int = 4,
    ):
        super().__init__()
        self.log_sigmas = nn.Parameter(torch.linspace(-1, 3, num_scales))

        self.B_base = nn.Parameter(
            torch.randn(input_dim, num_frequencies), requires_grad=False
        )
        self.output_dim_ff = 2 * num_frequencies * num_scales

        layers = []
        in_dim = self.output_dim_ff
        for h in hidden_dims:
            layers.append(nn.Linear(in_dim, h))
            layers.append(nn.Tanh())
            in_dim = h
        layers.append(nn.Linear(in_dim, output_dim))
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        sigmas = torch.exp(self.log_sigmas)
        projections = []
        for s in sigmas:
            proj = 2.0 * np.pi * x @ (self.B_base * s)
            projections.extend([torch.sin(proj), torch.cos(proj)])
        features = torch.cat(projections, dim=-1)
        return self.net(features)


def run_multiscale_benchmark():
    """Benchmark: solve u_xx + u_yy = f with multi-scale source term."""
    torch.manual_seed(42)
    device = torch.device("cpu")

    def exact_solution(x, y):
        return (
            torch.sin(2 * np.pi * x) * torch.sin(2 * np.pi * y)
            + 0.1 * torch.sin(20 * np.pi * x) * torch.sin(20 * np.pi * y)
            + 0.01 * torch.sin(50 * np.pi * x) * torch.sin(50 * np.pi * y)
        )

    def source_term(x, y):
        return (
            -2 * (2 * np.pi) ** 2 * torch.sin(2 * np.pi * x) * torch.sin(2 * np.pi * y)
            - 0.1 * 2 * (20 * np.pi) ** 2 * torch.sin(20 * np.pi * x) * torch.sin(20 * np.pi * y)
            - 0.01 * 2 * (50 * np.pi) ** 2 * torch.sin(50 * np.pi * x) * torch.sin(50 * np.pi * y)
        )

    configs = {
        "Standard MLP": {"sigma_values": [1.0], "num_frequencies": 128},
        "Multi-scale (sigma=1,10,50)": {"sigma_values": [1.0, 10.0, 50.0], "num_frequencies": 64},
        "Multi-scale (sigma=1,10,50,100)": {"sigma_values": [1.0, 10.0, 50.0, 100.0], "num_frequencies": 48},
    }

    results = {}
    n_interior = 2000
    n_boundary = 400
    n_epochs = 2000
    lr = 1e-3

    for name, cfg in configs.items():
        model = MultiScalePINN(
            input_dim=2,
            output_dim=1,
            hidden_dims=[256, 256, 256],
            num_frequencies=cfg["num_frequencies"],
            sigma_values=cfg["sigma_values"],
        ).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, n_epochs)

        loss_history = []
        for epoch in range(n_epochs):
            x_int = torch.rand(n_interior, 1, requires_grad=True, device=device)
            y_int = torch.rand(n_interior, 1, requires_grad=True, device=device)
            xy_int = torch.cat([x_int, y_int], dim=1)

            u_pred = model(xy_int)

            u_x = torch.autograd.grad(u_pred, x_int, torch.ones_like(u_pred), create_graph=True)[0]
            u_xx = torch.autograd.grad(u_x, x_int, torch.ones_like(u_x), create_graph=True)[0]
            u_y = torch.autograd.grad(u_pred, y_int, torch.ones_like(u_pred), create_graph=True)[0]
            u_yy = torch.autograd.grad(u_y, y_int, torch.ones_like(u_y), create_graph=True)[0]

            f_pred = u_xx + u_yy
            f_exact = source_term(x_int, y_int)
            loss_pde = torch.mean((f_pred - f_exact) ** 2)

            x_bc = torch.cat([
                torch.zeros(n_boundary // 4, 1, device=device),
                torch.ones(n_boundary // 4, 1, device=device),
                torch.rand(n_boundary // 4, 1, device=device),
                torch.rand(n_boundary // 4, 1, device=device),
            ])
            y_bc = torch.cat([
                torch.rand(n_boundary // 4, 1, device=device),
                torch.rand(n_boundary // 4, 1, device=device),
                torch.zeros(n_boundary // 4, 1, device=device),
                torch.ones(n_boundary // 4, 1, device=device),
            ])
            xy_bc = torch.cat([x_bc, y_bc], dim=1)
            u_bc = model(xy_bc)
            loss_bc = torch.mean(u_bc ** 2)

            loss = loss_pde + 100.0 * loss_bc

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            scheduler.step()

            if epoch % 100 == 0:
                loss_history.append({"epoch": epoch, "loss": loss.item(), "pde": loss_pde.item(), "bc": loss_bc.item()})

        x_test = torch.linspace(0, 1, 100, device=device)
        y_test = torch.linspace(0, 1, 100, device=device)
        X, Y = torch.meshgrid(x_test, y_test, indexing="ij")
        xy_test = torch.stack([X.flatten(), Y.flatten()], dim=1)

        with torch.no_grad():
            u_pred_test = model(xy_test).reshape(100, 100)
            u_exact_test = exact_solution(X, Y)
            l2_error = torch.norm(u_pred_test - u_exact_test) / torch.norm(u_exact_test)

        results[name] = {
            "l2_relative_error": l2_error.item(),
            "final_loss": loss_history[-1]["loss"],
            "loss_history": loss_history,
            "params": sum(p.numel() for p in model.parameters()),
        }
        print(f"{name}: L2 error = {l2_error.item():.6f}, params = {results[name]['params']}")

    return results


if __name__ == "__main__":
    results = run_multiscale_benchmark()
