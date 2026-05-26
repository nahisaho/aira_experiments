"""
PINN Extension Experiments: Multi-scale, Inverse, Causal, Adaptive, Operator, and Navier-Stokes
JAX-based implementation framework for Physics-Informed Neural Networks.
"""

import jax
import jax.numpy as jnp
from jax import random, grad, jit, vmap
import optax
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from functools import partial
import time
import os

FIGURES_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGURES_DIR, exist_ok=True)

# ============================================================
# Core PINN Building Blocks
# ============================================================

def init_mlp(key, layer_sizes):
    """Initialize MLP parameters with Xavier initialization."""
    params = []
    for i in range(len(layer_sizes) - 1):
        key, subkey = random.split(key)
        fan_in, fan_out = layer_sizes[i], layer_sizes[i+1]
        scale = jnp.sqrt(2.0 / (fan_in + fan_out))
        w = scale * random.normal(subkey, (fan_in, fan_out))
        b = jnp.zeros(fan_out)
        params.append((w, b))
    return params

def mlp_forward(params, x):
    """Forward pass through MLP with tanh activation."""
    for w, b in params[:-1]:
        x = jnp.tanh(x @ w + b)
    w, b = params[-1]
    return x @ w + b

def fourier_feature_mapping(x, B):
    """Map inputs through random Fourier features: [sin(2πBx), cos(2πBx)]."""
    proj = 2.0 * jnp.pi * x @ B.T
    return jnp.concatenate([jnp.sin(proj), jnp.cos(proj)], axis=-1)

def init_fourier_pinn(key, input_dim, hidden_sizes, output_dim, n_fourier=64, sigma=1.0):
    """Initialize a Fourier-feature PINN."""
    key1, key2 = random.split(key)
    B = random.normal(key1, (n_fourier, input_dim)) * sigma
    layer_sizes = [2 * n_fourier] + list(hidden_sizes) + [output_dim]
    mlp_params = init_mlp(key2, layer_sizes)
    return {'B': B, 'mlp': mlp_params}

def fourier_pinn_forward(params, x):
    """Forward pass: Fourier features -> MLP."""
    feat = fourier_feature_mapping(x, params['B'])
    return mlp_forward(params['mlp'], feat)

# ============================================================
# Experiment 1: Multi-scale Problem (Helmholtz Equation)
# ============================================================

def experiment_multiscale():
    """
    Multi-scale Helmholtz: -Δu - k²u = f on [0,1]
    Compare standard PINN vs Fourier-feature PINN.
    """
    print("=" * 60)
    print("Experiment 1: Multi-scale Helmholtz (Fourier Features)")
    print("=" * 60)

    k_val = 20.0  # high frequency
    exact_fn = lambda x: jnp.sin(k_val * jnp.pi * x)
    source_fn = lambda x: (k_val**2 * jnp.pi**2 - k_val**2) * jnp.sin(k_val * jnp.pi * x)

    key = random.PRNGKey(42)
    n_colloc = 500
    n_bc = 2

    # Standard PINN
    key, k1 = random.split(key)
    std_params = init_mlp(k1, [1, 128, 128, 128, 1])

    # Fourier PINN
    key, k2 = random.split(key)
    ff_params = init_fourier_pinn(k2, 1, [128, 128, 128], 1, n_fourier=64, sigma=10.0)

    def helmholtz_loss_std(params, key):
        x_c = random.uniform(key, (n_colloc, 1))
        u = lambda x: mlp_forward(params, x.reshape(1, 1)).squeeze()
        u_xx = vmap(grad(grad(u)))(x_c.squeeze())
        f_val = source_fn(x_c.squeeze())
        pde_res = -u_xx - k_val**2 * vmap(u)(x_c.squeeze()) - f_val
        pde_loss = jnp.mean(pde_res**2)
        # BC
        u0 = mlp_forward(params, jnp.array([[0.0]])).squeeze()
        u1 = mlp_forward(params, jnp.array([[1.0]])).squeeze()
        bc_loss = u0**2 + u1**2
        return pde_loss + 100.0 * bc_loss

    def helmholtz_loss_ff(params, key):
        x_c = random.uniform(key, (n_colloc, 1))
        u = lambda x: fourier_pinn_forward(params, x.reshape(1, 1)).squeeze()
        u_xx = vmap(grad(grad(u)))(x_c.squeeze())
        f_val = source_fn(x_c.squeeze())
        pde_res = -u_xx - k_val**2 * vmap(u)(x_c.squeeze()) - f_val
        pde_loss = jnp.mean(pde_res**2)
        u0 = fourier_pinn_forward(params, jnp.array([[0.0]])).squeeze()
        u1 = fourier_pinn_forward(params, jnp.array([[1.0]])).squeeze()
        bc_loss = u0**2 + u1**2
        return pde_loss + 100.0 * bc_loss

    # Train both
    n_epochs = 3000
    lr = 1e-3
    optimizer_std = optax.adam(lr)
    optimizer_ff = optax.adam(lr)

    opt_state_std = optimizer_std.init(std_params)
    opt_state_ff = optimizer_ff.init(ff_params)

    losses_std, losses_ff = [], []

    @jit
    def step_std(params, opt_state, key):
        loss, grads = jax.value_and_grad(helmholtz_loss_std)(params, key)
        updates, opt_state_new = optimizer_std.update(grads, opt_state)
        params_new = optax.apply_updates(params, updates)
        return params_new, opt_state_new, loss

    @jit
    def step_ff(params, opt_state, key):
        loss, grads = jax.value_and_grad(helmholtz_loss_ff)(params, key)
        updates, opt_state_new = optimizer_ff.update(grads, opt_state)
        params_new = optax.apply_updates(params, updates)
        return params_new, opt_state_new, loss

    t0 = time.time()
    for epoch in range(n_epochs):
        key, k1, k2 = random.split(key, 3)
        std_params, opt_state_std, loss_s = step_std(std_params, opt_state_std, k1)
        ff_params, opt_state_ff, loss_f = step_ff(ff_params, opt_state_ff, k2)
        losses_std.append(float(loss_s))
        losses_ff.append(float(loss_f))
        if (epoch + 1) % 500 == 0:
            print(f"  Epoch {epoch+1}: Std Loss={loss_s:.6f}, FF Loss={loss_f:.6f}")
    t1 = time.time()
    print(f"  Training time: {t1-t0:.1f}s")

    # Evaluate and plot
    x_test = jnp.linspace(0, 1, 200).reshape(-1, 1)
    u_exact = exact_fn(x_test.squeeze())
    u_std = vmap(lambda x: mlp_forward(std_params, x.reshape(1, 1)).squeeze())(x_test.squeeze())
    u_ff = vmap(lambda x: fourier_pinn_forward(ff_params, x.reshape(1, 1)).squeeze())(x_test.squeeze())

    err_std = float(jnp.sqrt(jnp.mean((u_std - u_exact)**2)))
    err_ff = float(jnp.sqrt(jnp.mean((u_ff - u_exact)**2)))
    print(f"  RMSE Standard: {err_std:.6f}")
    print(f"  RMSE Fourier:  {err_ff:.6f}")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].plot(np.array(x_test), np.array(u_exact), 'k-', label='Exact', lw=2)
    axes[0].plot(np.array(x_test), np.array(u_std), 'r--', label=f'Std PINN (RMSE={err_std:.4f})')
    axes[0].plot(np.array(x_test), np.array(u_ff), 'b-.', label=f'FF-PINN (RMSE={err_ff:.4f})')
    axes[0].set_xlabel('x'); axes[0].set_ylabel('u(x)')
    axes[0].set_title('Multi-scale Helmholtz Solution'); axes[0].legend()
    axes[1].semilogy(losses_std, 'r-', alpha=0.7, label='Standard PINN')
    axes[1].semilogy(losses_ff, 'b-', alpha=0.7, label='Fourier PINN')
    axes[1].set_xlabel('Epoch'); axes[1].set_ylabel('Loss')
    axes[1].set_title('Training Convergence'); axes[1].legend()
    axes[2].plot(np.array(x_test), np.abs(np.array(u_std - u_exact)), 'r-', label='Std PINN')
    axes[2].plot(np.array(x_test), np.abs(np.array(u_ff - u_exact)), 'b-', label='FF-PINN')
    axes[2].set_xlabel('x'); axes[2].set_ylabel('|Error|')
    axes[2].set_title('Pointwise Error'); axes[2].legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'exp1_multiscale_helmholtz.png'), dpi=150)
    plt.close()

    return {'rmse_std': err_std, 'rmse_ff': err_ff, 'losses_std': losses_std, 'losses_ff': losses_ff}


# ============================================================
# Experiment 2: Inverse Problem + Uncertainty Quantification
# ============================================================

def experiment_inverse_uq():
    """
    Inverse heat equation: u_t = D * u_xx
    Estimate diffusion coefficient D from noisy observations.
    Ensemble-based uncertainty quantification.
    """
    print("\n" + "=" * 60)
    print("Experiment 2: Inverse Problem + Uncertainty Quantification")
    print("=" * 60)

    D_true = 0.05
    exact_u = lambda x, t: jnp.exp(-D_true * jnp.pi**2 * t) * jnp.sin(jnp.pi * x)

    key = random.PRNGKey(123)
    n_obs = 100
    key, k1, k2, k3 = random.split(key, 4)
    x_obs = random.uniform(k1, (n_obs, 1))
    t_obs = random.uniform(k2, (n_obs, 1)) * 0.5
    u_obs = vmap(lambda xt: exact_u(xt[0], xt[1]))(jnp.hstack([x_obs, t_obs]))
    noise = 0.02 * random.normal(k3, u_obs.shape)
    u_obs_noisy = u_obs + noise

    n_ensemble = 5
    ensemble_D = []
    all_losses = []

    for ens_idx in range(n_ensemble):
        key, k1 = random.split(key)
        params = init_mlp(k1, [2, 64, 64, 64, 1])
        log_D = jnp.array(-3.0)  # Initialize near true value log(0.05)≈-3

        def loss_fn(params, log_D, key):
            D_est = jnp.exp(log_D)
            xt = jnp.hstack([x_obs, t_obs])
            u_pred = vmap(lambda x: mlp_forward(params, x.reshape(1, 2)).squeeze())(xt)
            data_loss = jnp.mean((u_pred - u_obs_noisy)**2)
            # PDE residual
            n_c = 300
            key, k1, k2 = random.split(key, 3)
            xc = random.uniform(k1, (n_c, 1))
            tc = random.uniform(k2, (n_c, 1)) * 0.5

            def u_fn(x, t):
                inp = jnp.array([x, t])
                return mlp_forward(params, inp.reshape(1, 2)).squeeze()

            def pde_residual(x, t):
                u_t = grad(u_fn, argnums=1)(x, t)
                u_xx = grad(grad(u_fn, argnums=0), argnums=0)(x, t)
                return u_t - D_est * u_xx

            res = vmap(pde_residual)(xc.squeeze(), tc.squeeze())
            pde_loss = jnp.mean(res**2)
            # BC: u(0,t)=0, u(1,t)=0
            key, k3 = random.split(key)
            t_bc = random.uniform(k3, (50, 1)) * 0.5
            bc0 = vmap(lambda t: u_fn(0.0, t))(t_bc.squeeze())
            bc1 = vmap(lambda t: u_fn(1.0, t))(t_bc.squeeze())
            bc_loss = jnp.mean(bc0**2) + jnp.mean(bc1**2)
            # IC: u(x,0) = sin(πx)
            key, k4 = random.split(key)
            x_ic = random.uniform(k4, (50, 1))
            ic_pred = vmap(lambda x: u_fn(x, 0.0))(x_ic.squeeze())
            ic_true = jnp.sin(jnp.pi * x_ic.squeeze())
            ic_loss = jnp.mean((ic_pred - ic_true)**2)
            return 10.0 * data_loss + pde_loss + 20.0 * (bc_loss + ic_loss)

        optimizer = optax.adam(5e-4)
        opt_state = optimizer.init((params, log_D))
        losses = []

        @jit
        def step(params, log_D, opt_state, key):
            loss, grads = jax.value_and_grad(loss_fn, argnums=(0, 1))(params, log_D, key)
            updates, opt_state_new = optimizer.update(grads, opt_state)
            params_new, log_D_new = optax.apply_updates((params, log_D), updates)
            return params_new, log_D_new, opt_state_new, loss

        for epoch in range(4000):
            key, k = random.split(key)
            params, log_D, opt_state, loss = step(params, log_D, opt_state, k)
            losses.append(float(loss))
            if (epoch + 1) % 1000 == 0:
                D_est = float(jnp.exp(log_D))
                print(f"  Ensemble {ens_idx+1}, Epoch {epoch+1}: Loss={loss:.6f}, D_est={D_est:.5f}")

        D_est = float(jnp.exp(log_D))
        ensemble_D.append(D_est)
        all_losses.append(losses)
        print(f"  Ensemble {ens_idx+1} final D: {D_est:.5f} (true: {D_true})")

    D_mean = np.mean(ensemble_D)
    D_std = np.std(ensemble_D)
    print(f"\n  Ensemble D: {D_mean:.5f} ± {D_std:.5f} (true: {D_true})")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for i, losses in enumerate(all_losses):
        axes[0].semilogy(losses, alpha=0.7, label=f'Ens {i+1}')
    axes[0].set_xlabel('Epoch'); axes[0].set_ylabel('Loss')
    axes[0].set_title('Training Loss (Ensemble)'); axes[0].legend()

    axes[1].bar(range(n_ensemble), ensemble_D, color='steelblue', alpha=0.7)
    axes[1].axhline(D_true, color='r', linestyle='--', label=f'True D={D_true}')
    axes[1].errorbar(n_ensemble//2, D_mean, yerr=2*D_std, fmt='ko', capsize=5, label=f'Mean±2σ')
    axes[1].set_xlabel('Ensemble Member'); axes[1].set_ylabel('Estimated D')
    axes[1].set_title('Parameter Estimation'); axes[1].legend()

    axes[2].hist(ensemble_D, bins=max(3, n_ensemble//2), color='steelblue', alpha=0.7, edgecolor='black')
    axes[2].axvline(D_true, color='r', linestyle='--', lw=2, label=f'True D={D_true}')
    axes[2].axvline(D_mean, color='k', linestyle='-', lw=2, label=f'Mean={D_mean:.4f}')
    axes[2].set_xlabel('D'); axes[2].set_ylabel('Count')
    axes[2].set_title('UQ Distribution'); axes[2].legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'exp2_inverse_uq.png'), dpi=150)
    plt.close()

    return {'D_true': D_true, 'D_mean': D_mean, 'D_std': D_std, 'ensemble_D': ensemble_D}


# ============================================================
# Experiment 3: Causal Training
# ============================================================

def experiment_causal_training():
    """
    Compare standard vs causal PINN training on advection equation:
    u_t + c * u_x = 0, u(x,0) = sin(2πx)
    Causal weighting: w_i = exp(-ε * Σ_{j<i} L_j)
    """
    print("\n" + "=" * 60)
    print("Experiment 3: Causal Training (Advection Equation)")
    print("=" * 60)

    c = 1.0
    exact_u = lambda x, t: jnp.sin(2 * jnp.pi * (x - c * t))

    key = random.PRNGKey(77)

    def train_pinn(key, causal=False, n_epochs=3000):
        key, k1 = random.split(key)
        params = init_mlp(k1, [2, 64, 64, 64, 1])
        optimizer = optax.adam(1e-3)
        opt_state = optimizer.init(params)

        n_t_steps = 20
        n_x = 30
        eps = 10.0  # causal weighting strength

        def loss_fn(params, key):
            t_vals = jnp.linspace(0.01, 1, n_t_steps)
            key, k1 = random.split(key)
            x_vals = random.uniform(k1, (n_x,))

            def u_fn(x, t):
                inp = jnp.array([x, t])
                return mlp_forward(params, inp.reshape(1, 2)).squeeze()

            # Compute PDE residual at each time level
            def residual_at_t(t):
                def res_at_x(x):
                    u_t = grad(u_fn, argnums=1)(x, t)
                    u_x = grad(u_fn, argnums=0)(x, t)
                    return (u_t + c * u_x)**2
                return jnp.mean(vmap(res_at_x)(x_vals))

            time_losses = vmap(residual_at_t)(t_vals)

            if causal:
                # Causal weighting: weight earlier time steps more
                cumsum = jnp.cumsum(time_losses)
                weights = jnp.exp(-eps * jnp.concatenate([jnp.array([0.0]), cumsum[:-1]]))
                weights = weights / jnp.sum(weights) * n_t_steps
                pde_loss = jnp.mean(weights * time_losses)
            else:
                pde_loss = jnp.mean(time_losses)

            # IC
            x_ic = jnp.linspace(0, 1, 50)
            ic_pred = vmap(lambda x: u_fn(x, 0.0))(x_ic)
            ic_true = jnp.sin(2 * jnp.pi * x_ic)
            ic_loss = jnp.mean((ic_pred - ic_true)**2)

            # Periodic BC
            t_bc = jnp.linspace(0, 1, 30)
            bc_left = vmap(lambda t: u_fn(0.0, t))(t_bc)
            bc_right = vmap(lambda t: u_fn(1.0, t))(t_bc)
            bc_loss = jnp.mean((bc_left - bc_right)**2)

            return pde_loss + 10.0 * ic_loss + 5.0 * bc_loss

        @jit
        def step(params, opt_state, key):
            loss, grads = jax.value_and_grad(loss_fn)(params, key)
            updates, new_state = optimizer.update(grads, opt_state)
            new_params = optax.apply_updates(params, updates)
            return new_params, new_state, loss

        losses = []
        for epoch in range(n_epochs):
            key, k = random.split(key)
            params, opt_state, loss = step(params, opt_state, k)
            losses.append(float(loss))
        return params, losses

    # Train both
    key, k1, k2 = random.split(key, 3)
    t0 = time.time()
    params_std, losses_std = train_pinn(k1, causal=False)
    t_std = time.time() - t0
    t0 = time.time()
    params_causal, losses_causal = train_pinn(k2, causal=True)
    t_causal = time.time() - t0

    print(f"  Standard: final loss = {losses_std[-1]:.6f}, time = {t_std:.1f}s")
    print(f"  Causal:   final loss = {losses_causal[-1]:.6f}, time = {t_causal:.1f}s")

    # Evaluate
    nx, nt = 100, 100
    x_grid = jnp.linspace(0, 1, nx)
    t_grid = jnp.linspace(0, 1, nt)
    X, T = jnp.meshgrid(x_grid, t_grid)
    xt_flat = jnp.stack([X.ravel(), T.ravel()], axis=-1)

    u_exact_grid = vmap(lambda xt: exact_u(xt[0], xt[1]))(xt_flat).reshape(nt, nx)
    u_std_grid = vmap(lambda xt: mlp_forward(params_std, xt.reshape(1, 2)).squeeze())(xt_flat).reshape(nt, nx)
    u_causal_grid = vmap(lambda xt: mlp_forward(params_causal, xt.reshape(1, 2)).squeeze())(xt_flat).reshape(nt, nx)

    err_std = float(jnp.sqrt(jnp.mean((u_std_grid - u_exact_grid)**2)))
    err_causal = float(jnp.sqrt(jnp.mean((u_causal_grid - u_exact_grid)**2)))
    print(f"  RMSE Standard: {err_std:.6f}")
    print(f"  RMSE Causal:   {err_causal:.6f}")

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    X_np, T_np = np.array(X), np.array(T)
    vmin, vmax = -1, 1
    for ax, data, title in zip(axes[0], 
                                [u_exact_grid, u_std_grid, u_causal_grid],
                                ['Exact', f'Standard (RMSE={err_std:.4f})', f'Causal (RMSE={err_causal:.4f})']):
        im = ax.pcolormesh(X_np, T_np, np.array(data), vmin=vmin, vmax=vmax, cmap='RdBu_r', shading='auto')
        ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_title(title)
        plt.colorbar(im, ax=ax)

    err_std_map = np.abs(np.array(u_std_grid - u_exact_grid))
    err_causal_map = np.abs(np.array(u_causal_grid - u_exact_grid))
    for ax, data, title in zip(axes[1, :2], [err_std_map, err_causal_map],
                                ['Standard Error', 'Causal Error']):
        im = ax.pcolormesh(X_np, T_np, data, cmap='hot', shading='auto')
        ax.set_xlabel('x'); ax.set_ylabel('t'); ax.set_title(title)
        plt.colorbar(im, ax=ax)

    axes[1, 2].semilogy(losses_std, 'r-', alpha=0.7, label='Standard')
    axes[1, 2].semilogy(losses_causal, 'b-', alpha=0.7, label='Causal')
    axes[1, 2].set_xlabel('Epoch'); axes[1, 2].set_ylabel('Loss')
    axes[1, 2].set_title('Training Convergence'); axes[1, 2].legend()
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'exp3_causal_training.png'), dpi=150)
    plt.close()

    return {'rmse_std': err_std, 'rmse_causal': err_causal, 'time_std': t_std, 'time_causal': t_causal}


# ============================================================
# Experiment 4: Adaptive Collocation Points
# ============================================================

def experiment_adaptive_collocation():
    """
    Burgers equation: u_t + u*u_x = ν*u_xx
    Compare uniform vs residual-based adaptive collocation.
    """
    print("\n" + "=" * 60)
    print("Experiment 4: Adaptive Collocation (Burgers Equation)")
    print("=" * 60)

    nu = 0.01 / jnp.pi

    key = random.PRNGKey(99)

    def train_burgers(key, adaptive=False, n_epochs=3000):
        key, k1 = random.split(key)
        params = init_mlp(k1, [2, 64, 64, 64, 1])
        optimizer = optax.adam(1e-3)
        opt_state = optimizer.init(params)
        n_colloc = 500

        def u_fn(params, x, t):
            inp = jnp.array([x, t])
            return mlp_forward(params, inp.reshape(1, 2)).squeeze()

        def pde_residual(params, x, t):
            u = u_fn(params, x, t)
            u_t = grad(u_fn, argnums=2)(params, x, t)
            u_x = grad(u_fn, argnums=1)(params, x, t)
            u_xx = grad(grad(u_fn, argnums=1), argnums=1)(params, x, t)
            return u_t + u * u_x - nu * u_xx

        def loss_fn(params, x_c, t_c, key):
            res = vmap(lambda x, t: pde_residual(params, x, t))(x_c, t_c)
            pde_loss = jnp.mean(res**2)
            # IC: u(x,0) = -sin(πx)
            x_ic = jnp.linspace(-1, 1, 100)
            ic_pred = vmap(lambda x: u_fn(params, x, 0.0))(x_ic)
            ic_true = -jnp.sin(jnp.pi * x_ic)
            ic_loss = jnp.mean((ic_pred - ic_true)**2)
            # BC: u(-1,t)=0, u(1,t)=0
            t_bc = jnp.linspace(0, 1, 50)
            bc_left = vmap(lambda t: u_fn(params, -1.0, t))(t_bc)
            bc_right = vmap(lambda t: u_fn(params, 1.0, t))(t_bc)
            bc_loss = jnp.mean(bc_left**2) + jnp.mean(bc_right**2)
            return pde_loss + 10.0 * ic_loss + 10.0 * bc_loss

        # Initial collocation points
        key, k1, k2 = random.split(key, 3)
        x_c = random.uniform(k1, (n_colloc,), minval=-1.0, maxval=1.0)
        t_c = random.uniform(k2, (n_colloc,), minval=0.0, maxval=1.0)

        losses = []
        point_snapshots = [(np.array(x_c), np.array(t_c))]

        for epoch in range(n_epochs):
            key, k = random.split(key)

            # Adaptive resampling every 500 epochs
            if adaptive and (epoch + 1) % 500 == 0 and epoch > 0:
                # Compute residuals at a dense grid
                key, k1, k2 = random.split(key, 3)
                x_dense = random.uniform(k1, (2000,), minval=-1.0, maxval=1.0)
                t_dense = random.uniform(k2, (2000,), minval=0.0, maxval=1.0)
                res_dense = vmap(lambda x, t: pde_residual(params, x, t)**2)(x_dense, t_dense)
                # Sample proportionally to residual magnitude
                probs = jnp.abs(res_dense) + 1e-8
                probs = probs / jnp.sum(probs)
                key, k3 = random.split(key)
                # Keep 70% uniform, 30% adaptive
                n_keep = int(n_colloc * 0.7)
                n_adaptive = n_colloc - n_keep
                key, k4, k5 = random.split(key, 3)
                x_uniform = random.uniform(k4, (n_keep,), minval=-1.0, maxval=1.0)
                t_uniform = random.uniform(k5, (n_keep,), minval=0.0, maxval=1.0)
                idx = random.choice(k3, 2000, shape=(n_adaptive,), p=probs)
                x_c = jnp.concatenate([x_uniform, x_dense[idx]])
                t_c = jnp.concatenate([t_uniform, t_dense[idx]])
                point_snapshots.append((np.array(x_c), np.array(t_c)))

            loss, grads = jax.value_and_grad(loss_fn)(params, x_c, t_c, k)
            updates, opt_state = optimizer.update(grads, opt_state)
            params = optax.apply_updates(params, updates)
            losses.append(float(loss))

            if (epoch + 1) % 1000 == 0:
                print(f"  {'Adaptive' if adaptive else 'Uniform'} Epoch {epoch+1}: Loss={loss:.6f}")

        return params, losses, point_snapshots

    key, k1, k2 = random.split(key, 3)
    params_uni, losses_uni, pts_uni = train_burgers(k1, adaptive=False)
    params_ada, losses_ada, pts_ada = train_burgers(k2, adaptive=True)

    print(f"  Uniform final loss:   {losses_uni[-1]:.6f}")
    print(f"  Adaptive final loss:  {losses_ada[-1]:.6f}")

    # Evaluate
    nx, nt = 100, 100
    x_grid = jnp.linspace(-1, 1, nx)
    t_grid = jnp.linspace(0, 1, nt)
    X, T = jnp.meshgrid(x_grid, t_grid)
    xt_flat = jnp.stack([X.ravel(), T.ravel()], axis=-1)

    u_uni = vmap(lambda xt: mlp_forward(params_uni, xt.reshape(1, 2)).squeeze())(xt_flat).reshape(nt, nx)
    u_ada = vmap(lambda xt: mlp_forward(params_ada, xt.reshape(1, 2)).squeeze())(xt_flat).reshape(nt, nx)

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    X_np, T_np = np.array(X), np.array(T)

    im = axes[0, 0].pcolormesh(X_np, T_np, np.array(u_uni), cmap='RdBu_r', shading='auto')
    axes[0, 0].set_title('Uniform Collocation'); axes[0, 0].set_xlabel('x'); axes[0, 0].set_ylabel('t')
    plt.colorbar(im, ax=axes[0, 0])

    im = axes[0, 1].pcolormesh(X_np, T_np, np.array(u_ada), cmap='RdBu_r', shading='auto')
    axes[0, 1].set_title('Adaptive Collocation'); axes[0, 1].set_xlabel('x'); axes[0, 1].set_ylabel('t')
    plt.colorbar(im, ax=axes[0, 1])

    axes[0, 2].semilogy(losses_uni, 'r-', alpha=0.7, label='Uniform')
    axes[0, 2].semilogy(losses_ada, 'b-', alpha=0.7, label='Adaptive')
    axes[0, 2].set_xlabel('Epoch'); axes[0, 2].set_ylabel('Loss')
    axes[0, 2].set_title('Training Convergence'); axes[0, 2].legend()

    # Show collocation point distribution
    for i, (pts, title) in enumerate([(pts_uni[-1], 'Uniform Points'), (pts_ada[-1], 'Adaptive Points (Final)')]):
        axes[1, i].scatter(pts[0], pts[1], s=2, alpha=0.5)
        axes[1, i].set_xlim(-1, 1); axes[1, i].set_ylim(0, 1)
        axes[1, i].set_xlabel('x'); axes[1, i].set_ylabel('t')
        axes[1, i].set_title(title)

    if len(pts_ada) > 1:
        axes[1, 2].scatter(pts_ada[0][0], pts_ada[0][1], s=2, alpha=0.5, label='Initial', c='gray')
        axes[1, 2].scatter(pts_ada[-1][0], pts_ada[-1][1], s=2, alpha=0.5, label='Final', c='blue')
        axes[1, 2].set_xlim(-1, 1); axes[1, 2].set_ylim(0, 1)
        axes[1, 2].set_xlabel('x'); axes[1, 2].set_ylabel('t')
        axes[1, 2].set_title('Adaptive Points Evolution'); axes[1, 2].legend()
    else:
        axes[1, 2].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'exp4_adaptive_collocation.png'), dpi=150)
    plt.close()

    return {'loss_uniform': losses_uni[-1], 'loss_adaptive': losses_ada[-1]}


# ============================================================
# Experiment 5: Operator Learning Comparison (DeepONet-style vs FNO-style)
# ============================================================

def experiment_operator_comparison():
    """
    Compare PINN, DeepONet-like, and FNO-like approaches on 1D Poisson.
    -u'' = f(x), u(0)=u(1)=0, for multiple source terms f.
    """
    print("\n" + "=" * 60)
    print("Experiment 5: Operator Learning Comparison")
    print("=" * 60)

    key = random.PRNGKey(55)

    # Generate training data: parametric Poisson
    def poisson_exact(a, k, x):
        """Exact solution for f(x) = a*sin(k*pi*x) -> u = a*sin(k*pi*x)/(k*pi)^2"""
        return a * jnp.sin(k * jnp.pi * x) / (k * jnp.pi)**2

    n_functions = 50
    n_x = 64
    x_grid = jnp.linspace(0, 1, n_x)

    key, k1, k2 = random.split(key, 3)
    amplitudes = random.uniform(k1, (n_functions,), minval=0.5, maxval=3.0)
    wavenumbers = random.choice(k2, jnp.array([1, 2, 3, 4, 5]), shape=(n_functions,))

    # Generate input-output pairs
    f_data = vmap(lambda a, k: a * jnp.sin(k * jnp.pi * x_grid))(amplitudes, wavenumbers)
    u_data = vmap(lambda a, k: vmap(lambda x: poisson_exact(a, k, x))(x_grid))(amplitudes, wavenumbers)

    # PINN approach: train per-function (evaluate on a few test functions)
    n_test = 5
    test_idx = list(range(n_functions - n_test, n_functions))
    train_idx = list(range(n_functions - n_test))

    # DeepONet-like: branch (encodes f) + trunk (encodes x)
    key, k1 = random.split(key)
    branch_params = init_mlp(k1, [n_x, 64, 64, 32])
    key, k2 = random.split(key)
    trunk_params = init_mlp(k2, [1, 64, 64, 32])

    def deeponet_forward(branch_p, trunk_p, f_input, x):
        branch_out = mlp_forward(branch_p, f_input.reshape(1, -1))  # (1, 32)
        trunk_out = mlp_forward(trunk_p, x.reshape(1, 1))  # (1, 32)
        return jnp.sum(branch_out * trunk_out)

    def deeponet_loss(bp, tp, f_batch, u_batch):
        def single_loss(f_in, u_true):
            u_pred = vmap(lambda x: deeponet_forward(bp, tp, f_in, x))(x_grid)
            return jnp.mean((u_pred - u_true)**2)
        return jnp.mean(vmap(single_loss)(f_batch, u_batch))

    optimizer_don = optax.adam(1e-3)
    opt_state_don = optimizer_don.init((branch_params, trunk_params))
    f_train = f_data[jnp.array(train_idx)]
    u_train = u_data[jnp.array(train_idx)]

    @jit
    def step_don(bp, tp, opt_state):
        loss, grads = jax.value_and_grad(deeponet_loss, argnums=(0, 1))(bp, tp, f_train, u_train)
        updates, opt_state_new = optimizer_don.update(grads, opt_state)
        (bp_new, tp_new) = optax.apply_updates((bp, tp), updates)
        return bp_new, tp_new, opt_state_new, loss

    losses_don = []
    for epoch in range(2000):
        branch_params, trunk_params, opt_state_don, loss = step_don(branch_params, trunk_params, opt_state_don)
        losses_don.append(float(loss))
        if (epoch + 1) % 500 == 0:
            print(f"  DeepONet Epoch {epoch+1}: Loss={loss:.6f}")

    # FNO-like: spectral convolution layer (simplified)
    key, k1 = random.split(key)
    fno_params = {
        'lift': init_mlp(k1, [1, 32]),
        'spectral_weights': random.normal(random.split(k1)[0], (16, 32, 32)) * 0.01,
        'proj': init_mlp(random.split(k1)[1], [32, 1])
    }

    def fno_forward(params, f_input):
        """Simplified 1D FNO: lift -> spectral conv -> project."""
        # Lift: (n_x,) -> (n_x, 32)
        x = vmap(lambda fi: mlp_forward(params['lift'], fi.reshape(1, 1)).squeeze())(f_input)
        # Spectral convolution on real-valued features
        # Use DCT-like approach: real FFT per feature channel
        n_x, d = x.shape[0], x.shape[1]
        n_modes = min(16, n_x // 2)
        # Real-valued spectral: truncate, multiply by weight matrix, pad back
        x_fft = jnp.fft.rfft(x, axis=0)  # (n_x//2+1, 32) complex
        x_fft_trunc = x_fft[:n_modes]  # (n_modes, 32) complex
        # Apply complex weight: (n_modes, 32, 32) x (n_modes, 32) -> (n_modes, 32)
        sw = params['spectral_weights'][:n_modes]  # (n_modes, 32, 32) real
        x_real = jnp.einsum('kij,kj->ki', sw, x_fft_trunc.real)
        x_imag = jnp.einsum('kij,kj->ki', sw, x_fft_trunc.imag)
        x_filtered = x_real + 1j * x_imag  # (n_modes, 32) complex
        x_padded = jnp.zeros_like(x_fft)
        x_padded = x_padded.at[:n_modes].set(x_filtered)
        x_out = jnp.fft.irfft(x_padded, n=n_x, axis=0)  # (n_x, 32) real
        x = jax.nn.gelu(x + x_out)
        # Project
        u_pred = vmap(lambda xi: mlp_forward(params['proj'], xi.reshape(1, -1)).squeeze())(x)
        return u_pred

    def fno_loss(params, f_batch, u_batch):
        def single_loss(f_in, u_true):
            u_pred = fno_forward(params, f_in)
            return jnp.mean((u_pred - u_true)**2)
        return jnp.mean(vmap(single_loss)(f_batch, u_batch))

    optimizer_fno = optax.adam(1e-3)
    opt_state_fno = optimizer_fno.init(fno_params)

    @jit
    def step_fno(params, opt_state):
        loss, grads = jax.value_and_grad(fno_loss)(params, f_train, u_train)
        updates, opt_state_new = optimizer_fno.update(grads, opt_state)
        params_new = optax.apply_updates(params, updates)
        return params_new, opt_state_new, loss

    losses_fno = []
    for epoch in range(2000):
        fno_params, opt_state_fno, loss = step_fno(fno_params, opt_state_fno)
        losses_fno.append(float(loss))
        if (epoch + 1) % 500 == 0:
            print(f"  FNO Epoch {epoch+1}: Loss={loss:.6f}")

    # Evaluate on test set
    f_test = f_data[jnp.array(test_idx)]
    u_test = u_data[jnp.array(test_idx)]

    errors_don = []
    errors_fno = []
    for i in range(n_test):
        u_pred_don = vmap(lambda x: deeponet_forward(branch_params, trunk_params, f_test[i], x))(x_grid)
        u_pred_fno = fno_forward(fno_params, f_test[i])
        err_don = float(jnp.sqrt(jnp.mean((u_pred_don - u_test[i])**2)))
        err_fno = float(jnp.sqrt(jnp.mean((u_pred_fno - u_test[i])**2)))
        errors_don.append(err_don)
        errors_fno.append(err_fno)

    mean_err_don = np.mean(errors_don)
    mean_err_fno = np.mean(errors_fno)
    print(f"\n  DeepONet test RMSE: {mean_err_don:.6f}")
    print(f"  FNO test RMSE:     {mean_err_fno:.6f}")

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes[0, 0].semilogy(losses_don, 'b-', alpha=0.7, label='DeepONet')
    axes[0, 0].semilogy(losses_fno, 'r-', alpha=0.7, label='FNO')
    axes[0, 0].set_xlabel('Epoch'); axes[0, 0].set_ylabel('Loss')
    axes[0, 0].set_title('Training Convergence'); axes[0, 0].legend()

    x_np = np.array(x_grid)
    for i in range(min(3, n_test)):
        u_pred_don = np.array(vmap(lambda x: deeponet_forward(branch_params, trunk_params, f_test[i], x))(x_grid))
        u_pred_fno = np.array(fno_forward(fno_params, f_test[i]))
        ax = axes[0, 1] if i == 0 else (axes[0, 2] if i == 1 else axes[1, 0])
        ax.plot(x_np, np.array(u_test[i]), 'k-', lw=2, label='Exact')
        ax.plot(x_np, u_pred_don, 'b--', label='DeepONet')
        ax.plot(x_np, u_pred_fno, 'r-.', label='FNO')
        ax.set_xlabel('x'); ax.set_ylabel('u(x)')
        ax.set_title(f'Test Case {i+1}'); ax.legend()

    # Bar chart comparison
    bar_x = np.arange(n_test)
    width = 0.35
    axes[1, 1].bar(bar_x - width/2, errors_don, width, label='DeepONet', color='steelblue')
    axes[1, 1].bar(bar_x + width/2, errors_fno, width, label='FNO', color='coral')
    axes[1, 1].set_xlabel('Test Case'); axes[1, 1].set_ylabel('RMSE')
    axes[1, 1].set_title('Test Error Comparison'); axes[1, 1].legend()

    # Summary
    methods = ['DeepONet', 'FNO']
    mean_errors = [mean_err_don, mean_err_fno]
    axes[1, 2].bar(methods, mean_errors, color=['steelblue', 'coral'])
    axes[1, 2].set_ylabel('Mean RMSE'); axes[1, 2].set_title('Overall Comparison')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'exp5_operator_comparison.png'), dpi=150)
    plt.close()

    return {'deeponet_rmse': mean_err_don, 'fno_rmse': mean_err_fno, 'errors_don': errors_don, 'errors_fno': errors_fno}


# ============================================================
# Experiment 6: Navier-Stokes Case Study (Lid-driven Cavity)
# ============================================================

def experiment_navier_stokes():
    """
    2D Steady Navier-Stokes (lid-driven cavity, simplified).
    u * u_x + v * u_y = -p_x + ν(u_xx + u_yy)
    u * v_x + v * v_y = -p_y + ν(v_xx + v_yy)
    u_x + v_y = 0 (continuity)
    """
    print("\n" + "=" * 60)
    print("Experiment 6: Navier-Stokes (Lid-Driven Cavity)")
    print("=" * 60)

    Re = 100.0
    nu_val = 1.0 / Re
    key = random.PRNGKey(42)

    # Network outputs: [u, v, p]
    key, k1 = random.split(key)
    params = init_fourier_pinn(k1, 2, [128, 128, 128, 128], 3, n_fourier=64, sigma=4.0)

    def uvp_fn(params, x, y):
        inp = jnp.array([x, y])
        out = fourier_pinn_forward(params, inp.reshape(1, 2)).squeeze()
        return out[0], out[1], out[2]

    def ns_residuals(params, x, y):
        u_fn = lambda x, y: uvp_fn(params, x, y)[0]
        v_fn = lambda x, y: uvp_fn(params, x, y)[1]
        p_fn = lambda x, y: uvp_fn(params, x, y)[2]

        u = u_fn(x, y)
        v = v_fn(x, y)

        u_x = grad(u_fn, argnums=0)(x, y)
        u_y = grad(u_fn, argnums=1)(x, y)
        u_xx = grad(grad(u_fn, argnums=0), argnums=0)(x, y)
        u_yy = grad(grad(u_fn, argnums=1), argnums=1)(x, y)

        v_x = grad(v_fn, argnums=0)(x, y)
        v_y = grad(v_fn, argnums=1)(x, y)
        v_xx = grad(grad(v_fn, argnums=0), argnums=0)(x, y)
        v_yy = grad(grad(v_fn, argnums=1), argnums=1)(x, y)

        p_x = grad(p_fn, argnums=0)(x, y)
        p_y = grad(p_fn, argnums=1)(x, y)

        res_x = u * u_x + v * u_y + p_x - nu_val * (u_xx + u_yy)
        res_y = u * v_x + v * v_y + p_y - nu_val * (v_xx + v_yy)
        res_cont = u_x + v_y

        return res_x, res_y, res_cont

    def loss_fn(params, key):
        n_interior = 500
        key, k1, k2 = random.split(key, 3)
        x_int = random.uniform(k1, (n_interior,))
        y_int = random.uniform(k2, (n_interior,))

        res = vmap(lambda x, y: ns_residuals(params, x, y))(x_int, y_int)
        pde_loss = jnp.mean(res[0]**2) + jnp.mean(res[1]**2) + jnp.mean(res[2]**2)

        # BCs: walls no-slip, top lid u=1
        n_bc = 50
        key, k1 = random.split(key)
        s = random.uniform(k1, (n_bc,))

        # Bottom: y=0, u=v=0
        bc_bottom = vmap(lambda x: uvp_fn(params, x, 0.0)[:2])(s)
        # Top: y=1, u=1, v=0
        bc_top = vmap(lambda x: uvp_fn(params, x, 1.0)[:2])(s)
        # Left: x=0, u=v=0
        bc_left = vmap(lambda y: uvp_fn(params, 0.0, y)[:2])(s)
        # Right: x=1, u=v=0
        bc_right = vmap(lambda y: uvp_fn(params, 1.0, y)[:2])(s)

        bc_loss = (jnp.mean(bc_bottom[0]**2) + jnp.mean(bc_bottom[1]**2) +
                   jnp.mean((bc_top[0] - 1.0)**2) + jnp.mean(bc_top[1]**2) +
                   jnp.mean(bc_left[0]**2) + jnp.mean(bc_left[1]**2) +
                   jnp.mean(bc_right[0]**2) + jnp.mean(bc_right[1]**2))

        return pde_loss + 20.0 * bc_loss

    optimizer = optax.adam(1e-3)
    opt_state = optimizer.init(params)

    @jit
    def step(params, opt_state, key):
        loss, grads = jax.value_and_grad(loss_fn)(params, key)
        updates, opt_state_new = optimizer.update(grads, opt_state)
        params_new = optax.apply_updates(params, updates)
        return params_new, opt_state_new, loss

    losses = []
    t0 = time.time()
    n_epochs = 3000
    for epoch in range(n_epochs):
        key, k = random.split(key)
        params, opt_state, loss = step(params, opt_state, k)
        losses.append(float(loss))
        if (epoch + 1) % 1000 == 0:
            print(f"  Epoch {epoch+1}: Loss={loss:.6f}")
    t_train = time.time() - t0
    print(f"  Training time: {t_train:.1f}s")

    # Evaluate on grid
    n_grid = 50
    x_grid = jnp.linspace(0, 1, n_grid)
    y_grid = jnp.linspace(0, 1, n_grid)
    X, Y = jnp.meshgrid(x_grid, y_grid)
    xy_flat = jnp.stack([X.ravel(), Y.ravel()], axis=-1)

    uvp = vmap(lambda xy: jnp.array(uvp_fn(params, xy[0], xy[1])))(xy_flat)
    U = uvp[:, 0].reshape(n_grid, n_grid)
    V = uvp[:, 1].reshape(n_grid, n_grid)
    P = uvp[:, 2].reshape(n_grid, n_grid)
    speed = jnp.sqrt(U**2 + V**2)

    # Compute divergence as quality metric
    res_all = vmap(lambda xy: jnp.array(ns_residuals(params, xy[0], xy[1])))(xy_flat)
    div_field = res_all[:, 2].reshape(n_grid, n_grid)
    mean_div = float(jnp.mean(jnp.abs(div_field)))
    print(f"  Mean |divergence|: {mean_div:.6f}")

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    X_np, Y_np = np.array(X), np.array(Y)

    im = axes[0, 0].pcolormesh(X_np, Y_np, np.array(U), cmap='RdBu_r', shading='auto')
    axes[0, 0].set_title('u-velocity'); axes[0, 0].set_xlabel('x'); axes[0, 0].set_ylabel('y')
    plt.colorbar(im, ax=axes[0, 0])

    im = axes[0, 1].pcolormesh(X_np, Y_np, np.array(V), cmap='RdBu_r', shading='auto')
    axes[0, 1].set_title('v-velocity'); axes[0, 1].set_xlabel('x'); axes[0, 1].set_ylabel('y')
    plt.colorbar(im, ax=axes[0, 1])

    im = axes[0, 2].pcolormesh(X_np, Y_np, np.array(P), cmap='coolwarm', shading='auto')
    axes[0, 2].set_title('Pressure'); axes[0, 2].set_xlabel('x'); axes[0, 2].set_ylabel('y')
    plt.colorbar(im, ax=axes[0, 2])

    im = axes[1, 0].pcolormesh(X_np, Y_np, np.array(speed), cmap='viridis', shading='auto')
    skip = 3
    axes[1, 0].quiver(X_np[::skip, ::skip], Y_np[::skip, ::skip],
                       np.array(U[::skip, ::skip]), np.array(V[::skip, ::skip]),
                       color='white', alpha=0.7, scale=15)
    axes[1, 0].set_title('Velocity Magnitude + Streamlines'); axes[1, 0].set_xlabel('x'); axes[1, 0].set_ylabel('y')
    plt.colorbar(im, ax=axes[1, 0])

    im = axes[1, 1].pcolormesh(X_np, Y_np, np.abs(np.array(div_field)), cmap='hot', shading='auto')
    axes[1, 1].set_title(f'|Divergence| (mean={mean_div:.4f})'); axes[1, 1].set_xlabel('x'); axes[1, 1].set_ylabel('y')
    plt.colorbar(im, ax=axes[1, 1])

    axes[1, 2].semilogy(losses, 'b-', alpha=0.7)
    axes[1, 2].set_xlabel('Epoch'); axes[1, 2].set_ylabel('Loss')
    axes[1, 2].set_title('Training Convergence')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'exp6_navier_stokes.png'), dpi=150)
    plt.close()

    return {'final_loss': losses[-1], 'mean_divergence': mean_div, 'training_time': t_train}


# ============================================================
# Summary Figure
# ============================================================

def create_summary_figure(results):
    """Create a summary comparison figure."""
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # Exp 1: Multi-scale bar
    methods = ['Standard\nPINN', 'Fourier\nPINN']
    rmses = [results['exp1']['rmse_std'], results['exp1']['rmse_ff']]
    colors = ['coral', 'steelblue']
    axes[0, 0].bar(methods, rmses, color=colors)
    axes[0, 0].set_ylabel('RMSE')
    axes[0, 0].set_title('Exp1: Multi-scale Helmholtz')

    # Exp 2: Inverse problem
    axes[0, 1].errorbar(['Estimated D'], [results['exp2']['D_mean']], 
                         yerr=[2*results['exp2']['D_std']], fmt='bo', capsize=10, markersize=10)
    axes[0, 1].axhline(results['exp2']['D_true'], color='r', ls='--', label=f'True D={results["exp2"]["D_true"]}')
    axes[0, 1].set_ylabel('D'); axes[0, 1].set_title('Exp2: Inverse Problem + UQ')
    axes[0, 1].legend()

    # Exp 3: Causal training
    methods = ['Standard', 'Causal']
    rmses = [results['exp3']['rmse_std'], results['exp3']['rmse_causal']]
    axes[0, 2].bar(methods, rmses, color=['coral', 'steelblue'])
    axes[0, 2].set_ylabel('RMSE'); axes[0, 2].set_title('Exp3: Causal Training')

    # Exp 4: Adaptive collocation
    methods = ['Uniform', 'Adaptive']
    losses = [results['exp4']['loss_uniform'], results['exp4']['loss_adaptive']]
    axes[1, 0].bar(methods, losses, color=['coral', 'steelblue'])
    axes[1, 0].set_ylabel('Final Loss'); axes[1, 0].set_title('Exp4: Adaptive Collocation')

    # Exp 5: Operator comparison
    methods = ['DeepONet', 'FNO']
    rmses = [results['exp5']['deeponet_rmse'], results['exp5']['fno_rmse']]
    axes[1, 1].bar(methods, rmses, color=['steelblue', 'coral'])
    axes[1, 1].set_ylabel('Mean RMSE'); axes[1, 1].set_title('Exp5: Operator Learning')

    # Exp 6: NS
    metrics = ['Loss', 'Mean |div|']
    vals = [results['exp6']['final_loss'], results['exp6']['mean_divergence']]
    axes[1, 2].bar(metrics, vals, color=['steelblue', 'coral'])
    axes[1, 2].set_ylabel('Value'); axes[1, 2].set_title('Exp6: Navier-Stokes')

    plt.suptitle('PINN Extension Methods: Summary of Results', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, 'summary_comparison.png'), dpi=150)
    plt.close()


# ============================================================
# Main
# ============================================================

if __name__ == '__main__':
    results = {}
    
    results['exp1'] = experiment_multiscale()
    results['exp2'] = experiment_inverse_uq()
    results['exp3'] = experiment_causal_training()
    results['exp4'] = experiment_adaptive_collocation()
    results['exp5'] = experiment_operator_comparison()
    results['exp6'] = experiment_navier_stokes()
    
    create_summary_figure(results)
    
    print("\n" + "=" * 60)
    print("ALL EXPERIMENTS COMPLETE")
    print("=" * 60)
    print(f"\nResults Summary:")
    print(f"  Exp1 - Multi-scale: Std RMSE={results['exp1']['rmse_std']:.6f}, FF RMSE={results['exp1']['rmse_ff']:.6f}")
    print(f"  Exp2 - Inverse: D={results['exp2']['D_mean']:.5f}±{results['exp2']['D_std']:.5f} (true={results['exp2']['D_true']})")
    print(f"  Exp3 - Causal: Std RMSE={results['exp3']['rmse_std']:.6f}, Causal RMSE={results['exp3']['rmse_causal']:.6f}")
    print(f"  Exp4 - Adaptive: Uniform loss={results['exp4']['loss_uniform']:.6f}, Adaptive loss={results['exp4']['loss_adaptive']:.6f}")
    print(f"  Exp5 - Operators: DeepONet RMSE={results['exp5']['deeponet_rmse']:.6f}, FNO RMSE={results['exp5']['fno_rmse']:.6f}")
    print(f"  Exp6 - NS: Final loss={results['exp6']['final_loss']:.6f}, Mean |div|={results['exp6']['mean_divergence']:.6f}")
    
    # Save results for report generation
    import json
    results_serializable = {}
    for k, v in results.items():
        results_serializable[k] = {}
        for kk, vv in v.items():
            if isinstance(vv, (list, float, int, str)):
                if isinstance(vv, list) and len(vv) > 20:
                    results_serializable[k][kk] = {'first_10': vv[:10], 'last_10': vv[-10:], 'len': len(vv)}
                else:
                    results_serializable[k][kk] = vv
            else:
                results_serializable[k][kk] = str(vv)
    
    with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), 'results.json'), 'w') as f:
        json.dump(results_serializable, f, indent=2)
    print("\nResults saved to results.json")
    print("Figures saved to figures/")
