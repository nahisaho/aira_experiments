"""
Module 3: Causal Training (Time-Discrete PINN).

Enforces temporal causality during training by weighting PDE residual
losses so that earlier time steps are learned before later ones.

References:
- Wang et al., "Respecting Causality is All You Need for Training
  Physics-Informed Neural Networks" (2022)
"""

import numpy as np
import torch
import torch.nn as nn
from typing import List, Dict, Optional, Tuple


class CausalPINN(nn.Module):
    """PINN with causal training weighting."""

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


class CausalTrainer:
    """
    Causal training strategy for time-dependent PDEs.

    Divides the time domain into N_t segments. Each segment's loss is
    weighted by exp(-epsilon * cumulative_loss_of_previous_segments).
    """

    def __init__(
        self,
        model: CausalPINN,
        t_min: float = 0.0,
        t_max: float = 1.0,
        n_time_segments: int = 20,
        epsilon: float = 1.0,
        device: torch.device = torch.device("cpu"),
    ):
        self.model = model.to(device)
        self.device = device
        self.t_min = t_min
        self.t_max = t_max
        self.n_time_segments = n_time_segments
        self.epsilon = epsilon

        dt = (t_max - t_min) / n_time_segments
        self.t_boundaries = [t_min + i * dt for i in range(n_time_segments + 1)]

    def compute_causal_weights(
        self, segment_losses: List[float]
    ) -> torch.Tensor:
        weights = torch.ones(len(segment_losses), device=self.device)
        cumulative = 0.0
        for i in range(len(segment_losses)):
            weights[i] = torch.exp(torch.tensor(-self.epsilon * cumulative, device=self.device))
            cumulative += segment_losses[i]
        return weights

    def train_step(
        self,
        pde_residual_fn,
        n_points_per_segment: int = 100,
        x_range: Tuple[float, float] = (0.0, 1.0),
    ) -> Tuple[torch.Tensor, List[float], torch.Tensor]:
        segment_losses = []
        segment_residuals = []

        for i in range(self.n_time_segments):
            t_lo = self.t_boundaries[i]
            t_hi = self.t_boundaries[i + 1]

            x = torch.rand(n_points_per_segment, 1, requires_grad=True, device=self.device)
            x = x * (x_range[1] - x_range[0]) + x_range[0]
            t = torch.rand(n_points_per_segment, 1, requires_grad=True, device=self.device)
            t = t * (t_hi - t_lo) + t_lo

            residual = pde_residual_fn(self.model, x, t)
            seg_loss = torch.mean(residual ** 2)
            segment_losses.append(seg_loss.item())
            segment_residuals.append(seg_loss)

        weights = self.compute_causal_weights(segment_losses)

        total_loss = sum(w * r for w, r in zip(weights, segment_residuals))
        return total_loss, segment_losses, weights


def run_causal_training_benchmark():
    """
    Benchmark: Convection equation u_t + c*u_x = 0
    Compare causal vs standard training.
    """
    torch.manual_seed(42)
    device = torch.device("cpu")

    c = 1.0
    T_final = 1.0

    def exact_solution(x, t):
        return torch.sin(2 * np.pi * (x - c * t))

    def pde_residual(model, x, t):
        xt = torch.cat([x, t], dim=1)
        u = model(xt)
        u_t = torch.autograd.grad(u, t, torch.ones_like(u), create_graph=True)[0]
        u_x = torch.autograd.grad(u, x, torch.ones_like(u), create_graph=True)[0]
        return u_t + c * u_x

    results = {}

    for mode in ["standard", "causal"]:
        model = CausalPINN(input_dim=2, output_dim=1, hidden_dims=[128, 128, 128, 128])
        model = model.to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

        if mode == "causal":
            trainer = CausalTrainer(
                model, t_min=0.0, t_max=T_final,
                n_time_segments=20, epsilon=10.0, device=device,
            )

        n_epochs = 3000
        loss_history = []
        n_pde = 2000
        n_ic = 500
        n_bc = 200

        for epoch in range(n_epochs):
            if mode == "causal":
                loss_pde, seg_losses, weights = trainer.train_step(
                    pde_residual, n_points_per_segment=100
                )
            else:
                x_pde = torch.rand(n_pde, 1, requires_grad=True, device=device)
                t_pde = torch.rand(n_pde, 1, requires_grad=True, device=device) * T_final
                residual = pde_residual(model, x_pde, t_pde)
                loss_pde = torch.mean(residual ** 2)

            x_ic = torch.rand(n_ic, 1, device=device)
            t_ic = torch.zeros(n_ic, 1, device=device)
            xt_ic = torch.cat([x_ic, t_ic], dim=1)
            u_ic = model(xt_ic)
            u_ic_exact = torch.sin(2 * np.pi * x_ic)
            loss_ic = torch.mean((u_ic - u_ic_exact) ** 2)

            t_bc = torch.rand(n_bc, 1, device=device) * T_final
            xt_0 = torch.cat([torch.zeros(n_bc, 1, device=device), t_bc], dim=1)
            xt_1 = torch.cat([torch.ones(n_bc, 1, device=device), t_bc], dim=1)
            loss_bc = torch.mean((model(xt_0) - model(xt_1)) ** 2)

            loss = loss_pde + 10.0 * loss_ic + 10.0 * loss_bc

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            if epoch % 100 == 0:
                loss_history.append({"epoch": epoch, "loss": loss.item()})

        errors_by_time = []
        time_slices = [0.0, 0.25, 0.5, 0.75, 1.0]
        for t_val in time_slices:
            x_test = torch.linspace(0, 1, 200, device=device).unsqueeze(1)
            t_test = torch.full_like(x_test, t_val)
            xt_test = torch.cat([x_test, t_test], dim=1)
            with torch.no_grad():
                u_pred = model(xt_test)
            u_exact = exact_solution(x_test, t_test)
            l2_err = torch.norm(u_pred - u_exact) / torch.norm(u_exact)
            errors_by_time.append({"t": t_val, "l2_error": l2_err.item()})

        results[mode] = {
            "loss_history": loss_history,
            "errors_by_time": errors_by_time,
            "final_loss": loss_history[-1]["loss"],
            "avg_l2_error": np.mean([e["l2_error"] for e in errors_by_time]),
        }

        print(f"\n{mode.upper()} training:")
        for e in errors_by_time:
            print(f"  t={e['t']:.2f}: L2 error = {e['l2_error']:.6f}")
        print(f"  Average L2 error = {results[mode]['avg_l2_error']:.6f}")

    improvement = (
        (results["standard"]["avg_l2_error"] - results["causal"]["avg_l2_error"])
        / results["standard"]["avg_l2_error"]
    ) * 100
    print(f"\nCausal training improvement: {improvement:.1f}%")
    results["improvement_pct"] = improvement

    return results


if __name__ == "__main__":
    results = run_causal_training_benchmark()
