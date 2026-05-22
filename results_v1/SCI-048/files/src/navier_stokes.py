"""
Module 6: Navier-Stokes Turbulence Case Study.

PINN for 2D incompressible Navier-Stokes (lid-driven cavity).

References:
- Raissi et al., "PINNs" (JCP, 2019)
"""

import numpy as np
import torch
import torch.nn as nn
from typing import List, Dict


class NavierStokesPINN(nn.Module):
    """PINN for 2D incompressible Navier-Stokes."""

    def __init__(
        self,
        hidden_dims: List[int] = [128, 128, 128, 128, 128],
        use_fourier: bool = False,
        n_fourier: int = 64,
        sigma: float = 1.0,
    ):
        super().__init__()
        self.use_fourier = use_fourier

        if use_fourier:
            self.B = nn.Parameter(
                torch.randn(3, n_fourier) * sigma, requires_grad=False
            )
            input_dim = 2 * n_fourier
        else:
            input_dim = 3

        layers = []
        in_dim = input_dim
        for h in hidden_dims:
            layers.append(nn.Linear(in_dim, h))
            layers.append(nn.Tanh())
            in_dim = h
        layers.append(nn.Linear(in_dim, 2))  # (psi, p)
        self.net = nn.Sequential(*layers)

        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x, y, t):
        xyt = torch.cat([x, y, t], dim=1)
        if self.use_fourier:
            proj = 2.0 * np.pi * xyt @ self.B
            xyt = torch.cat([torch.sin(proj), torch.cos(proj)], dim=-1)
        out = self.net(xyt)
        psi = out[:, 0:1]
        p = out[:, 1:2]
        u = torch.autograd.grad(psi, y, torch.ones_like(psi), create_graph=True)[0]
        v = -torch.autograd.grad(psi, x, torch.ones_like(psi), create_graph=True)[0]
        return u, v, p


def ns_residuals(model, x, y, t, nu):
    u, v, p = model(x, y, t)
    u_t = torch.autograd.grad(u, t, torch.ones_like(u), create_graph=True)[0]
    u_x = torch.autograd.grad(u, x, torch.ones_like(u), create_graph=True)[0]
    u_y = torch.autograd.grad(u, y, torch.ones_like(u), create_graph=True)[0]
    v_t = torch.autograd.grad(v, t, torch.ones_like(v), create_graph=True)[0]
    v_x = torch.autograd.grad(v, x, torch.ones_like(v), create_graph=True)[0]
    v_y = torch.autograd.grad(v, y, torch.ones_like(v), create_graph=True)[0]
    p_x = torch.autograd.grad(p, x, torch.ones_like(p), create_graph=True)[0]
    p_y = torch.autograd.grad(p, y, torch.ones_like(p), create_graph=True)[0]
    u_xx = torch.autograd.grad(u_x, x, torch.ones_like(u_x), create_graph=True)[0]
    u_yy = torch.autograd.grad(u_y, y, torch.ones_like(u_y), create_graph=True)[0]
    v_xx = torch.autograd.grad(v_x, x, torch.ones_like(v_x), create_graph=True)[0]
    v_yy = torch.autograd.grad(v_y, y, torch.ones_like(v_y), create_graph=True)[0]

    res_u = u_t + u * u_x + v * u_y + p_x - nu * (u_xx + u_yy)
    res_v = v_t + u * v_x + v * v_y + p_y - nu * (v_xx + v_yy)
    res_cont = u_x + v_y
    return res_u, res_v, res_cont


def run_navier_stokes_benchmark():
    """Benchmark: 2D lid-driven cavity at Re=100 and Re=400."""
    torch.manual_seed(42)
    device = torch.device("cpu")

    reynolds_numbers = [100, 400]
    results = {}

    for Re in reynolds_numbers:
        nu = 1.0 / Re
        print(f"\n=== Re = {Re} (nu = {nu:.4f}) ===")

        model = NavierStokesPINN(
            hidden_dims=[128, 128, 128, 128, 128],
            use_fourier=True, n_fourier=64, sigma=2.0,
        ).to(device)

        optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=3000)

        n_epochs = 3000
        n_pde = 2000
        n_bc = 200
        loss_history = []

        for epoch in range(n_epochs):
            x_int = torch.rand(n_pde, 1, requires_grad=True, device=device)
            y_int = torch.rand(n_pde, 1, requires_grad=True, device=device)
            t_int = torch.zeros(n_pde, 1, requires_grad=True, device=device)

            res_u, res_v, res_cont = ns_residuals(model, x_int, y_int, t_int, nu)
            loss_pde = (
                torch.mean(res_u ** 2) + torch.mean(res_v ** 2)
                + 10.0 * torch.mean(res_cont ** 2)
            )

            # Bottom: u=0, v=0
            x_b = torch.rand(n_bc, 1, requires_grad=True, device=device)
            y_b = torch.zeros(n_bc, 1, requires_grad=True, device=device)
            t_b = torch.zeros(n_bc, 1, requires_grad=True, device=device)
            u_b, v_b, _ = model(x_b, y_b, t_b)
            loss_bc = torch.mean(u_b ** 2) + torch.mean(v_b ** 2)

            # Top (lid): u=1, v=0
            x_t = torch.rand(n_bc, 1, requires_grad=True, device=device)
            y_t = torch.ones(n_bc, 1, requires_grad=True, device=device)
            t_t = torch.zeros(n_bc, 1, requires_grad=True, device=device)
            u_t_bc, v_t_bc, _ = model(x_t, y_t, t_t)
            loss_bc += torch.mean((u_t_bc - 1.0) ** 2) + torch.mean(v_t_bc ** 2)

            # Left: u=0, v=0
            x_l = torch.zeros(n_bc, 1, requires_grad=True, device=device)
            y_l = torch.rand(n_bc, 1, requires_grad=True, device=device)
            t_l = torch.zeros(n_bc, 1, requires_grad=True, device=device)
            u_l, v_l, _ = model(x_l, y_l, t_l)
            loss_bc += torch.mean(u_l ** 2) + torch.mean(v_l ** 2)

            # Right: u=0, v=0
            x_r = torch.ones(n_bc, 1, requires_grad=True, device=device)
            y_r = torch.rand(n_bc, 1, requires_grad=True, device=device)
            t_r = torch.zeros(n_bc, 1, requires_grad=True, device=device)
            u_r, v_r, _ = model(x_r, y_r, t_r)
            loss_bc += torch.mean(u_r ** 2) + torch.mean(v_r ** 2)

            loss = loss_pde + 100.0 * loss_bc

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            scheduler.step()

            if epoch % 200 == 0:
                loss_history.append({
                    "epoch": epoch, "loss": loss.item(),
                    "pde": loss_pde.item(), "bc": loss_bc.item(),
                })
                if epoch % 1000 == 0:
                    print(f"  Epoch {epoch}: loss={loss.item():.6f}")

        # Evaluate centerlines
        n_eval = 100
        x_vc = torch.full((n_eval, 1), 0.5, requires_grad=True, device=device)
        y_vc = torch.linspace(0, 1, n_eval, device=device).unsqueeze(1).clone().requires_grad_(True)
        t_vc = torch.zeros(n_eval, 1, requires_grad=True, device=device)
        u_vc, v_vc, p_vc = model(x_vc, y_vc, t_vc)

        x_hc = torch.linspace(0, 1, n_eval, device=device).unsqueeze(1).clone().requires_grad_(True)
        y_hc = torch.full((n_eval, 1), 0.5, requires_grad=True, device=device)
        t_hc = torch.zeros(n_eval, 1, requires_grad=True, device=device)
        u_hc, v_hc, p_hc = model(x_hc, y_hc, t_hc)

        # Continuity check
        x_check = torch.rand(1000, 1, requires_grad=True, device=device)
        y_check = torch.rand(1000, 1, requires_grad=True, device=device)
        t_check = torch.zeros(1000, 1, requires_grad=True, device=device)
        _, _, res_c = ns_residuals(model, x_check, y_check, t_check, nu)
        cont_error = torch.mean(res_c ** 2).item()

        results[f"Re{Re}"] = {
            "Re": Re, "nu": nu,
            "final_loss": loss_history[-1]["loss"],
            "final_pde_loss": loss_history[-1]["pde"],
            "final_bc_loss": loss_history[-1]["bc"],
            "continuity_mse": cont_error,
            "loss_history": loss_history,
            "u_centerline_vert": u_vc.detach().cpu().numpy(),
            "v_centerline_horiz": v_hc.detach().cpu().numpy(),
            "y_coords": y_vc.detach().cpu().numpy(),
            "x_coords": x_hc.detach().cpu().numpy(),
            "params": sum(p.numel() for p in model.parameters()),
        }
        print(f"  Re={Re}: final_loss={loss_history[-1]['loss']:.6f}, continuity_mse={cont_error:.8f}")

    return results


if __name__ == "__main__":
    results = run_navier_stokes_benchmark()
