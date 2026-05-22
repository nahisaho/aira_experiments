"""
Module 5: Operator Learning — DeepONet vs FNO Comparison.

References:
- Lu et al., "Learning nonlinear operators via DeepONet" (Nature MI, 2021)
- Li et al., "Fourier Neural Operator for Parametric PDEs" (ICLR 2021)
"""

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import List, Tuple, Dict


class DeepONet(nn.Module):
    """Deep Operator Network."""

    def __init__(
        self,
        branch_input_dim: int = 100,
        trunk_input_dim: int = 1,
        hidden_dim: int = 128,
        n_basis: int = 64,
        n_layers: int = 4,
    ):
        super().__init__()
        branch_layers = [nn.Linear(branch_input_dim, hidden_dim), nn.Tanh()]
        for _ in range(n_layers - 2):
            branch_layers.extend([nn.Linear(hidden_dim, hidden_dim), nn.Tanh()])
        branch_layers.append(nn.Linear(hidden_dim, n_basis))
        self.branch = nn.Sequential(*branch_layers)

        trunk_layers = [nn.Linear(trunk_input_dim, hidden_dim), nn.Tanh()]
        for _ in range(n_layers - 2):
            trunk_layers.extend([nn.Linear(hidden_dim, hidden_dim), nn.Tanh()])
        trunk_layers.append(nn.Linear(hidden_dim, n_basis))
        self.trunk = nn.Sequential(*trunk_layers)

        self.bias = nn.Parameter(torch.zeros(1))

    def forward(self, u_input, y_query):
        branch_out = self.branch(u_input)
        trunk_out = self.trunk(y_query)
        if branch_out.dim() == 2 and trunk_out.dim() == 2:
            output = torch.einsum("bp,qp->bq", branch_out, trunk_out) + self.bias
        else:
            output = (branch_out * trunk_out).sum(dim=-1, keepdim=True) + self.bias
        return output


class SpectralConv1d(nn.Module):
    """1D Fourier layer for FNO."""

    def __init__(self, in_channels: int, out_channels: int, modes: int):
        super().__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.modes = modes
        scale = 1.0 / (in_channels * out_channels)
        self.weights = nn.Parameter(
            scale * torch.rand(in_channels, out_channels, modes, dtype=torch.cfloat)
        )

    def forward(self, x):
        batch_size = x.shape[0]
        x_ft = torch.fft.rfft(x, dim=-1)
        out_ft = torch.zeros(
            batch_size, self.out_channels, x_ft.shape[-1],
            dtype=torch.cfloat, device=x.device,
        )
        out_ft[:, :, :self.modes] = torch.einsum(
            "bix,iox->box", x_ft[:, :, :self.modes], self.weights
        )
        return torch.fft.irfft(out_ft, n=x.shape[-1])


class FNO1d(nn.Module):
    """1D Fourier Neural Operator."""

    def __init__(self, modes=16, width=64, n_layers=4, input_dim=2, output_dim=1):
        super().__init__()
        self.modes = modes
        self.width = width
        self.n_layers = n_layers

        self.fc0 = nn.Linear(input_dim, width)
        self.spectral_convs = nn.ModuleList()
        self.w_convs = nn.ModuleList()
        for _ in range(n_layers):
            self.spectral_convs.append(SpectralConv1d(width, width, modes))
            self.w_convs.append(nn.Conv1d(width, width, 1))

        self.fc1 = nn.Linear(width, 128)
        self.fc2 = nn.Linear(128, output_dim)

    def forward(self, x):
        x = self.fc0(x)
        x = x.permute(0, 2, 1)
        for i in range(self.n_layers):
            x1 = self.spectral_convs[i](x)
            x2 = self.w_convs[i](x)
            x = x1 + x2
            if i < self.n_layers - 1:
                x = F.gelu(x)
        x = x.permute(0, 2, 1)
        x = F.gelu(self.fc1(x))
        x = self.fc2(x)
        return x


def generate_parametric_poisson_data(n_samples=500, n_grid=100, device=torch.device("cpu")):
    """Generate data: -u_xx = f(x) on [0,1], u(0)=u(1)=0."""
    x = torch.linspace(0, 1, n_grid, device=device)
    dx = x[1] - x[0]

    f_all, u_all = [], []

    for _ in range(n_samples):
        n_modes = np.random.randint(1, 6)
        f = torch.zeros(n_grid, device=device)
        for _ in range(n_modes):
            k = np.random.uniform(1, 10)
            a = np.random.uniform(-2, 2)
            phi = np.random.uniform(0, 2 * np.pi)
            f += a * torch.sin(torch.tensor(k * np.pi) * x + phi)

        n_inner = n_grid - 2
        A = torch.zeros(n_inner, n_inner, device=device)
        for i in range(n_inner):
            A[i, i] = 2.0
            if i > 0:
                A[i, i - 1] = -1.0
            if i < n_inner - 1:
                A[i, i + 1] = -1.0
        A = A / dx ** 2

        rhs = f[1:-1]
        u_inner = torch.linalg.solve(A, rhs)
        u = torch.zeros(n_grid, device=device)
        u[1:-1] = u_inner

        f_all.append(f)
        u_all.append(u)

    return torch.stack(f_all), torch.stack(u_all), x


def run_operator_learning_benchmark():
    """Benchmark DeepONet vs FNO on parametric Poisson equation."""
    torch.manual_seed(42)
    device = torch.device("cpu")
    n_grid = 64

    print("Generating training data...")
    f_train, u_train, x_grid = generate_parametric_poisson_data(
        n_samples=500, n_grid=n_grid, device=device
    )
    f_test, u_test, _ = generate_parametric_poisson_data(
        n_samples=100, n_grid=n_grid, device=device
    )

    results = {}

    # DeepONet
    print("\nTraining DeepONet...")
    deeponet = DeepONet(
        branch_input_dim=n_grid, trunk_input_dim=1,
        hidden_dim=128, n_basis=64, n_layers=4,
    ).to(device)
    optimizer_don = torch.optim.Adam(deeponet.parameters(), lr=1e-3)
    x_trunk = x_grid.unsqueeze(1)
    loss_history_don = []

    for epoch in range(2000):
        idx = torch.randperm(f_train.shape[0])[:64]
        u_pred = deeponet(f_train[idx], x_trunk)
        loss = F.mse_loss(u_pred, u_train[idx])
        optimizer_don.zero_grad()
        loss.backward()
        optimizer_don.step()
        if epoch % 200 == 0:
            loss_history_don.append({"epoch": epoch, "loss": loss.item()})

    with torch.no_grad():
        u_pred_test = deeponet(f_test, x_trunk)
        l2_don = (torch.norm(u_pred_test - u_test) / torch.norm(u_test)).item()
        mse_don = F.mse_loss(u_pred_test, u_test).item()

    results["DeepONet"] = {
        "l2_relative_error": l2_don, "mse": mse_don,
        "params": sum(p.numel() for p in deeponet.parameters()),
        "loss_history": loss_history_don,
    }
    print(f"DeepONet: L2 error = {l2_don:.6f}")

    # FNO
    print("\nTraining FNO...")
    fno = FNO1d(modes=16, width=64, n_layers=4, input_dim=2, output_dim=1).to(device)
    optimizer_fno = torch.optim.Adam(fno.parameters(), lr=1e-3)
    x_grid_expanded = x_grid.unsqueeze(0).unsqueeze(-1).expand(f_train.shape[0], -1, -1)
    loss_history_fno = []

    for epoch in range(2000):
        idx = torch.randperm(f_train.shape[0])[:32]
        f_batch = f_train[idx].unsqueeze(-1)
        input_batch = torch.cat([x_grid_expanded[idx], f_batch], dim=-1)
        u_batch = u_train[idx].unsqueeze(-1)
        u_pred = fno(input_batch)
        loss = F.mse_loss(u_pred, u_batch)
        optimizer_fno.zero_grad()
        loss.backward()
        optimizer_fno.step()
        if epoch % 200 == 0:
            loss_history_fno.append({"epoch": epoch, "loss": loss.item()})

    with torch.no_grad():
        f_test_in = f_test.unsqueeze(-1)
        x_test_in = x_grid.unsqueeze(0).unsqueeze(-1).expand(f_test.shape[0], -1, -1)
        input_test = torch.cat([x_test_in, f_test_in], dim=-1)
        u_pred_fno = fno(input_test).squeeze(-1)
        l2_fno = (torch.norm(u_pred_fno - u_test) / torch.norm(u_test)).item()
        mse_fno = F.mse_loss(u_pred_fno, u_test).item()

    results["FNO"] = {
        "l2_relative_error": l2_fno, "mse": mse_fno,
        "params": sum(p.numel() for p in fno.parameters()),
        "loss_history": loss_history_fno,
    }
    print(f"FNO: L2 error = {l2_fno:.6f}")

    print("\n=== Operator Learning Comparison ===")
    for name, r in results.items():
        print(f"{name}: L2={r['l2_relative_error']:.6f}, MSE={r['mse']:.8f}, Params={r['params']}")

    return results


if __name__ == "__main__":
    results = run_operator_learning_benchmark()
