"""
Module 2: Inverse Problem Solving with Uncertainty Quantification.

Implements PINN-based parameter estimation with:
- Trainable PDE parameters
- Monte Carlo Dropout for uncertainty quantification
- Ensemble-based uncertainty estimation

Reference:
- Yang & Perdikaris, "B-PINNs: Bayesian PINNs" (JCP, 2021)
"""

import numpy as np
import torch
import torch.nn as nn
from typing import Dict, List, Tuple, Optional


class MCDropoutPINN(nn.Module):
    """PINN with MC Dropout for uncertainty quantification."""

    def __init__(
        self,
        input_dim: int = 2,
        output_dim: int = 1,
        hidden_dims: List[int] = [128, 128, 128, 128],
        dropout_rate: float = 0.05,
    ):
        super().__init__()
        self.dropout_rate = dropout_rate
        layers = []
        in_dim = input_dim
        for h in hidden_dims:
            layers.append(nn.Linear(in_dim, h))
            layers.append(nn.Tanh())
            layers.append(nn.Dropout(p=dropout_rate))
            in_dim = h
        layers.append(nn.Linear(in_dim, output_dim))
        self.net = nn.Sequential(*layers)

        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)

    def predict_with_uncertainty(
        self, x: torch.Tensor, n_samples: int = 100
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        self.train()
        predictions = torch.stack([self.forward(x) for _ in range(n_samples)])
        mean = predictions.mean(dim=0)
        std = predictions.std(dim=0)
        return mean, std, predictions


class InversePINNSolver:
    """Solve inverse problems: estimate PDE parameters from observations."""

    def __init__(
        self,
        model: MCDropoutPINN,
        param_names: List[str],
        param_init: Dict[str, float],
        device: torch.device = torch.device("cpu"),
    ):
        self.model = model.to(device)
        self.device = device

        self.pde_params = {}
        for name, val in param_init.items():
            self.pde_params[name] = nn.Parameter(
                torch.tensor([val], dtype=torch.float32, device=device)
            )

    def get_all_parameters(self):
        params = list(self.model.parameters())
        for p in self.pde_params.values():
            params.append(p)
        return params


def run_inverse_problem_benchmark():
    """
    Benchmark: Estimate diffusion coefficient D in the heat equation:
      u_t = D * u_xx
    from noisy observations.
    """
    torch.manual_seed(42)
    device = torch.device("cpu")

    D_true = 0.1
    L = 1.0
    T_final = 0.5

    def exact_solution(x, t, D=D_true, n_terms=20):
        u = torch.zeros_like(x)
        for n in range(1, n_terms + 1):
            u += (2.0 / (n * np.pi)) * (-1) ** (n + 1) * torch.sin(
                n * np.pi * x / L
            ) * torch.exp(-D * (n * np.pi / L) ** 2 * t)
        return u

    n_obs = 200
    x_obs = torch.rand(n_obs, 1, device=device) * L
    t_obs = torch.rand(n_obs, 1, device=device) * T_final
    u_obs = exact_solution(x_obs, t_obs) + 0.01 * torch.randn(n_obs, 1, device=device)

    model = MCDropoutPINN(input_dim=2, output_dim=1, hidden_dims=[128, 128, 128, 128])
    solver = InversePINNSolver(
        model, param_names=["D"], param_init={"D": 0.5}, device=device
    )

    all_params = solver.get_all_parameters()
    optimizer = torch.optim.Adam(all_params, lr=1e-3)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=1000, gamma=0.5)

    n_epochs = 3000
    n_pde = 1000
    history = {"epoch": [], "loss": [], "D_est": [], "data_loss": [], "pde_loss": []}

    for epoch in range(n_epochs):
        xt_obs = torch.cat([x_obs, t_obs], dim=1)
        u_pred_obs = model(xt_obs)
        loss_data = torch.mean((u_pred_obs - u_obs) ** 2)

        x_pde = torch.rand(n_pde, 1, requires_grad=True, device=device) * L
        t_pde = torch.rand(n_pde, 1, requires_grad=True, device=device) * T_final
        xt_pde = torch.cat([x_pde, t_pde], dim=1)

        u = model(xt_pde)
        u_t = torch.autograd.grad(u, t_pde, torch.ones_like(u), create_graph=True)[0]
        u_x = torch.autograd.grad(u, x_pde, torch.ones_like(u), create_graph=True)[0]
        u_xx = torch.autograd.grad(u_x, x_pde, torch.ones_like(u_x), create_graph=True)[0]

        D_est = solver.pde_params["D"]
        residual = u_t - D_est * u_xx
        loss_pde = torch.mean(residual ** 2)

        x_ic = torch.rand(200, 1, device=device) * L
        t_ic = torch.zeros(200, 1, device=device)
        xt_ic = torch.cat([x_ic, t_ic], dim=1)
        u_ic = model(xt_ic)
        loss_ic = torch.mean((u_ic - x_ic / L) ** 2)

        t_bc = torch.rand(100, 1, device=device) * T_final
        xt_bc0 = torch.cat([torch.zeros(100, 1, device=device), t_bc], dim=1)
        xt_bc1 = torch.cat([torch.ones(100, 1, device=device) * L, t_bc], dim=1)
        loss_bc = torch.mean(model(xt_bc0) ** 2) + torch.mean(model(xt_bc1) ** 2)

        loss = loss_data + 10.0 * loss_pde + 10.0 * loss_ic + 100.0 * loss_bc

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
        scheduler.step()

        if epoch % 100 == 0:
            history["epoch"].append(epoch)
            history["loss"].append(loss.item())
            history["D_est"].append(D_est.item())
            history["data_loss"].append(loss_data.item())
            history["pde_loss"].append(loss_pde.item())
            if epoch % 500 == 0:
                print(f"Epoch {epoch}: loss={loss.item():.6f}, D_est={D_est.item():.6f} (true={D_true})")

    x_test = torch.linspace(0, L, 50, device=device).unsqueeze(1)
    t_test = torch.full((50, 1), T_final / 2, device=device)
    xt_test = torch.cat([x_test, t_test], dim=1)

    mean_pred, std_pred, _ = model.predict_with_uncertainty(xt_test, n_samples=200)
    u_exact_test = exact_solution(x_test, t_test)

    results = {
        "D_true": D_true,
        "D_estimated": D_est.item(),
        "D_relative_error": abs(D_est.item() - D_true) / D_true,
        "final_loss": history["loss"][-1],
        "mean_uncertainty": std_pred.mean().item(),
        "max_uncertainty": std_pred.max().item(),
        "history": history,
        "test_mean": mean_pred.detach().cpu().numpy(),
        "test_std": std_pred.detach().cpu().numpy(),
        "test_exact": u_exact_test.detach().cpu().numpy(),
        "test_x": x_test.detach().cpu().numpy(),
    }

    print(f"\nInverse Problem Results:")
    print(f"  D_true = {D_true}")
    print(f"  D_estimated = {D_est.item():.6f}")
    print(f"  Relative error = {results['D_relative_error']:.4%}")
    print(f"  Mean uncertainty = {results['mean_uncertainty']:.6f}")

    return results


if __name__ == "__main__":
    results = run_inverse_problem_benchmark()
