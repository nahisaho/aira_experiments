"""
Module 4: Adaptive Collocation Point Placement.

Implements residual-based adaptive refinement (RAR) and error-indicator
driven strategies to concentrate collocation points where PDE residuals
are large.

References:
- Lu et al., "DeepXDE" (SIAM Review, 2021)
- Wu et al., "Comprehensive study of non-adaptive and residual-based
  adaptive sampling for PINNs" (CMAME, 2023)
"""

import numpy as np
import torch
import torch.nn as nn
from typing import List, Tuple, Optional, Dict


class AdaptiveCollocationPINN(nn.Module):
    """PINN with adaptive collocation point strategy."""

    def __init__(
        self,
        input_dim: int = 2,
        output_dim: int = 1,
        hidden_dims: List[int] = [128, 128, 128, 128],
    ):
        super().__init__()
        layers = []
        in_dim = input_dim
        for h in hidden_dims:
            layers.append(nn.Linear(in_dim, h))
            layers.append(nn.Tanh())
            in_dim = h
        layers.append(nn.Linear(in_dim, output_dim))
        self.net = nn.Sequential(*layers)

        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class ResidualAdaptiveRefinement:
    """Residual-based Adaptive Refinement (RAR) for collocation points."""

    def __init__(
        self,
        domain_bounds: List[Tuple[float, float]],
        n_initial: int = 1000,
        n_candidate: int = 10000,
        n_add: int = 200,
        device: torch.device = torch.device("cpu"),
    ):
        self.domain_bounds = domain_bounds
        self.n_initial = n_initial
        self.n_candidate = n_candidate
        self.n_add = n_add
        self.device = device
        self.dim = len(domain_bounds)

        self.points = self._sample_uniform(n_initial)

    def _sample_uniform(self, n: int) -> torch.Tensor:
        points = []
        for lo, hi in self.domain_bounds:
            points.append(torch.rand(n, 1, device=self.device) * (hi - lo) + lo)
        return torch.cat(points, dim=1)

    def refine(
        self, model: nn.Module, residual_fn, k: float = 2.0
    ) -> torch.Tensor:
        """Add points where residual is largest."""
        candidates = self._sample_uniform(self.n_candidate)

        cols = [candidates[:, i:i+1].clone().requires_grad_(True) for i in range(self.dim)]

        residuals = residual_fn(model, *cols)
        res_magnitude = torch.abs(residuals.detach()).squeeze()

        probs = res_magnitude ** k
        probs = probs / probs.sum()

        indices = torch.multinomial(probs, self.n_add, replacement=False)
        new_points = candidates[indices].detach()

        self.points = torch.cat([self.points, new_points], dim=0)
        return self.points

    def get_points(self) -> torch.Tensor:
        return self.points


def run_adaptive_collocation_benchmark():
    """
    Benchmark: Poisson equation with sharp internal layer.
    u_xx + u_yy = f(x,y) on [0,1]^2
    """
    torch.manual_seed(42)
    device = torch.device("cpu")

    def exact_solution(x, y):
        return torch.tanh(30 * (x - 0.5)) * torch.sin(np.pi * y)

    def source_term(x, y):
        t = torch.tanh(30 * (x - 0.5))
        s = 1 - t ** 2
        sin_y = torch.sin(np.pi * y)
        u_xx = -1800.0 * t * s * sin_y
        u_yy = -np.pi ** 2 * t * sin_y
        return u_xx + u_yy

    strategies = {
        "uniform": {"use_rar": False, "n_points": 2000},
        "RAR (k=1)": {"use_rar": True, "k": 1.0, "n_points": 1000},
        "RAR (k=2)": {"use_rar": True, "k": 2.0, "n_points": 1000},
    }

    results = {}

    for name, cfg in strategies.items():
        model = AdaptiveCollocationPINN(
            input_dim=2, output_dim=1, hidden_dims=[128, 128, 128, 128]
        ).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

        if cfg["use_rar"]:
            rar = ResidualAdaptiveRefinement(
                domain_bounds=[(0, 1), (0, 1)],
                n_initial=cfg["n_points"],
                n_candidate=5000,
                n_add=100,
                device=device,
            )

        n_epochs = 2000
        n_bc = 400
        refine_interval = 500
        loss_history = []

        for epoch in range(n_epochs):
            if cfg["use_rar"] and epoch > 0 and epoch % refine_interval == 0:
                def residual_fn(m, x, y):
                    xy = torch.cat([x, y], dim=1)
                    u = m(xy)
                    u_x = torch.autograd.grad(u, x, torch.ones_like(u), create_graph=True)[0]
                    u_xx = torch.autograd.grad(u_x, x, torch.ones_like(u_x), create_graph=True)[0]
                    u_y = torch.autograd.grad(u, y, torch.ones_like(u), create_graph=True)[0]
                    u_yy = torch.autograd.grad(u_y, y, torch.ones_like(u_y), create_graph=True)[0]
                    f = source_term(x.detach(), y.detach())
                    return u_xx + u_yy - f
                rar.refine(model, residual_fn, k=cfg.get("k", 2.0))

            if cfg["use_rar"]:
                pts = rar.get_points()
                x_pde = pts[:, 0:1].requires_grad_(True)
                y_pde = pts[:, 1:2].requires_grad_(True)
            else:
                x_pde = torch.rand(cfg["n_points"], 1, requires_grad=True, device=device)
                y_pde = torch.rand(cfg["n_points"], 1, requires_grad=True, device=device)

            xy_pde = torch.cat([x_pde, y_pde], dim=1)
            u = model(xy_pde)
            u_x = torch.autograd.grad(u, x_pde, torch.ones_like(u), create_graph=True)[0]
            u_xx = torch.autograd.grad(u_x, x_pde, torch.ones_like(u_x), create_graph=True)[0]
            u_y = torch.autograd.grad(u, y_pde, torch.ones_like(u), create_graph=True)[0]
            u_yy = torch.autograd.grad(u_y, y_pde, torch.ones_like(u_y), create_graph=True)[0]

            f = source_term(x_pde.detach(), y_pde.detach())
            loss_pde = torch.mean((u_xx + u_yy - f) ** 2)

            x_bc = torch.cat([
                torch.zeros(n_bc // 4, 1), torch.ones(n_bc // 4, 1),
                torch.rand(n_bc // 4, 1), torch.rand(n_bc // 4, 1),
            ]).to(device)
            y_bc = torch.cat([
                torch.rand(n_bc // 4, 1), torch.rand(n_bc // 4, 1),
                torch.zeros(n_bc // 4, 1), torch.ones(n_bc // 4, 1),
            ]).to(device)
            xy_bc = torch.cat([x_bc, y_bc], dim=1)
            u_bc_pred = model(xy_bc)
            u_bc_exact = exact_solution(x_bc, y_bc)
            loss_bc = torch.mean((u_bc_pred - u_bc_exact) ** 2)

            loss = loss_pde + 100.0 * loss_bc

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            if epoch % 100 == 0:
                loss_history.append({"epoch": epoch, "loss": loss.item()})

        x_test = torch.linspace(0, 1, 100, device=device)
        y_test = torch.linspace(0, 1, 100, device=device)
        X, Y = torch.meshgrid(x_test, y_test, indexing="ij")
        xy_test = torch.stack([X.flatten(), Y.flatten()], dim=1)

        with torch.no_grad():
            u_pred = model(xy_test).reshape(100, 100)
            u_exact = exact_solution(X, Y)
            l2_error = torch.norm(u_pred - u_exact) / torch.norm(u_exact)

        n_final_pts = rar.get_points().shape[0] if cfg["use_rar"] else cfg["n_points"]

        results[name] = {
            "l2_relative_error": l2_error.item(),
            "final_loss": loss_history[-1]["loss"],
            "n_collocation_points": n_final_pts,
            "loss_history": loss_history,
        }
        print(f"{name}: L2 error = {l2_error.item():.6f}, points = {n_final_pts}")

    return results


if __name__ == "__main__":
    results = run_adaptive_collocation_benchmark()
