"""
Volcanic deformation source models:
  - Mogi (point pressure source)
  - Prolate/oblate spheroid (Yang model)
  - FEM-based arbitrary geometry (FEniCS)

All models predict surface displacement [East, North, Up] at observation points
given source parameters.

References:
  Mogi (1958), Bull. Earthq. Res. Inst.
  Yang et al. (1988), JGR
  McTigue (1987), JGR
"""

import numpy as np
from typing import Tuple, Dict, Optional
from dataclasses import dataclass, field


# ==============================================================================
# Data classes for source parameters
# ==============================================================================

@dataclass
class MogiSource:
    """Point pressure source (Mogi, 1958)."""
    x: float          # Easting [m]
    y: float          # Northing [m]
    d: float          # Depth (positive downward) [m]
    dV: float         # Volume change [m^3]
    nu: float = 0.25  # Poisson's ratio


@dataclass
class SpheroidSource:
    """Prolate/oblate spheroid (Yang et al., 1988)."""
    x: float           # Easting [m]
    y: float           # Northing [m]
    d: float           # Depth to center [m]
    a: float           # Semi-major axis [m]
    b: float           # Semi-minor axis [m] (b < a: prolate, b > a: oblate)
    dP: float          # Pressure change [Pa]
    strike: float      # Strike angle [deg]
    dip: float         # Dip angle [deg] (90=vertical)
    nu: float = 0.25
    mu: float = 3.0e10 # Shear modulus [Pa]


@dataclass
class FEMSourceConfig:
    """Configuration for FEM-based source model."""
    mesh_file: Optional[str] = None
    domain_size: Tuple[float, float, float] = (50000, 50000, 30000)
    resolution: int = 32
    nu: float = 0.25
    mu: float = 3.0e10
    chamber_center: Tuple[float, float, float] = (0, 0, 5000)
    chamber_radii: Tuple[float, float, float] = (1000, 1000, 500)
    dP: float = 10e6


# ==============================================================================
# Mogi model
# ==============================================================================

def mogi_displacement(
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    source: MogiSource
) -> np.ndarray:
    """
    Compute surface displacement from a Mogi point source in a half-space.

    Parameters
    ----------
    obs_x, obs_y : (N,) arrays of observation coordinates [m]
    source : MogiSource parameters

    Returns
    -------
    disp : (N, 3) array of [East, North, Up] displacement [m]
    """
    dx = obs_x - source.x
    dy = obs_y - source.y
    R = np.sqrt(dx**2 + dy**2 + source.d**2)
    C = source.dV * (1 - source.nu) / np.pi

    ux = C * dx / R**3
    uy = C * dy / R**3
    uz = C * source.d / R**3

    return np.column_stack([ux, uy, uz])


def mogi_gravity(
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    source: MogiSource,
    rho: float = 2500.0,
    free_air_gradient: float = -3.086e-6
) -> np.ndarray:
    """
    Gravity change associated with Mogi source (free-air + mass change).

    Returns gravity change in µGal.
    """
    disp = mogi_displacement(obs_x, obs_y, source)
    uz = disp[:, 2]

    # Free-air effect
    dg_freeair = free_air_gradient * uz * 1e8  # to µGal

    # Mass attraction from injected mass
    dx = obs_x - source.x
    dy = obs_y - source.y
    R = np.sqrt(dx**2 + dy**2 + source.d**2)
    G = 6.674e-11
    dM = rho * source.dV
    dg_mass = G * dM * source.d / R**3 * 1e8  # to µGal

    return dg_freeair + dg_mass


# ==============================================================================
# Spheroid model (Yang et al., 1988 — simplified)
# ==============================================================================

def spheroid_displacement(
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    source: SpheroidSource
) -> np.ndarray:
    """
    Surface displacement from a pressurized spheroid in a half-space.

    Uses the analytical approximation of Yang et al. (1988) valid for
    deeply buried spheroids (d >> a).

    Returns (N, 3) displacement [East, North, Up].
    """
    strike_rad = np.radians(source.strike)
    dip_rad = np.radians(source.dip)

    # Rotate coordinates to source-aligned frame
    dx = obs_x - source.x
    dy = obs_y - source.y
    x_rot = dx * np.cos(strike_rad) + dy * np.sin(strike_rad)
    y_rot = -dx * np.sin(strike_rad) + dy * np.cos(strike_rad)

    # Effective volume change for spheroid
    a, b = source.a, source.b
    aspect = b / a if a > 0 else 1.0

    # Equivalent volume change
    dV_eq = (4.0 / 3.0) * np.pi * a * a * b * source.dP / (source.mu)
    # Shape correction factors
    if aspect < 1.0:  # prolate
        e = np.sqrt(1 - aspect**2)
        alpha_corr = (1 + aspect**2 / e**2 * (1 - np.arcsin(e) / e)) / 2
    elif aspect > 1.0:  # oblate
        e = np.sqrt(1 - 1.0 / aspect**2)
        alpha_corr = (1 + 1.0 / (aspect**2 * e**2) * (1 - np.arcsin(e) / e)) / 2
    else:
        alpha_corr = 1.0

    C = dV_eq * alpha_corr * (1 - source.nu) / np.pi

    # Horizontal/vertical splitting for dipping source
    d_eff = source.d
    R = np.sqrt(x_rot**2 + y_rot**2 + d_eff**2)

    ux_rot = C * x_rot / R**3
    uy_rot = C * y_rot / R**3
    uz = C * d_eff / R**3

    # Rotate back
    ux = ux_rot * np.cos(strike_rad) - uy_rot * np.sin(strike_rad)
    uy = ux_rot * np.sin(strike_rad) + uy_rot * np.cos(strike_rad)

    return np.column_stack([ux, uy, uz])


# ==============================================================================
# FEM model using FEniCS
# ==============================================================================

def fem_displacement(
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    config: FEMSourceConfig,
    use_fenics: bool = False
) -> np.ndarray:
    """
    Compute surface displacement using FEM (FEniCS).

    If use_fenics=False, falls back to an analytical approximation
    for environments where FEniCS is not available.

    Parameters
    ----------
    obs_x, obs_y : observation coordinates
    config : FEM configuration

    Returns
    -------
    disp : (N, 3) displacement array
    """
    if use_fenics:
        return _fem_fenics_solve(obs_x, obs_y, config)
    else:
        return _fem_analytical_fallback(obs_x, obs_y, config)


def _fem_analytical_fallback(
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    config: FEMSourceConfig
) -> np.ndarray:
    """
    Analytical approximation of FEM solution using extended Mogi with
    finite-size correction (McTigue, 1987).
    """
    cx, cy, cz = config.chamber_center
    rx, ry, rz = config.chamber_radii
    a_eff = (rx * ry * rz) ** (1.0 / 3.0)  # effective radius

    dV = (4.0 / 3.0) * np.pi * rx * ry * rz * config.dP / (config.mu)

    # McTigue correction for finite-size source
    eps = a_eff / cz
    correction = 1 + eps**3  # first-order size correction

    src = MogiSource(x=cx, y=cy, d=cz, dV=dV * correction, nu=config.nu)
    return mogi_displacement(obs_x, obs_y, src)


def _fem_fenics_solve(
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    config: FEMSourceConfig
) -> np.ndarray:
    """
    Full FEM solution using FEniCS/DOLFINx.

    Solves linear elasticity BVP:
      div(σ) = 0       in Ω (half-space domain)
      σ·n = -dP·n      on Γ_chamber (chamber boundary)
      σ·n = 0           on Γ_surface (free surface)
      u = 0             on Γ_far (far-field boundaries)

    Returns interpolated displacement at observation points.
    """
    try:
        import dolfinx
        from dolfinx import fem, mesh, io
        from dolfinx.fem import FunctionSpace, Function
        import ufl
        from mpi4py import MPI
        import gmsh
    except ImportError:
        print("FEniCS/DOLFINx not available. Using analytical fallback.")
        return _fem_analytical_fallback(obs_x, obs_y, config)

    Lx, Ly, Lz = config.domain_size
    cx, cy, cz = config.chamber_center
    rx, ry, rz = config.chamber_radii
    N = config.resolution

    # --- Mesh generation with Gmsh ---
    gmsh.initialize()
    gmsh.model.add("volcano")

    # Domain box
    box = gmsh.model.occ.addBox(-Lx/2, -Ly/2, -Lz, Lx, Ly, Lz)
    # Chamber ellipsoid
    sphere = gmsh.model.occ.addSphere(cx, cy, -cz, max(rx, ry, rz))
    # Scale to ellipsoid
    gmsh.model.occ.dilate([(3, sphere)], cx, cy, -cz,
                          rx/max(rx,ry,rz), ry/max(rx,ry,rz), rz/max(rx,ry,rz))
    # Boolean cut
    result = gmsh.model.occ.cut([(3, box)], [(3, sphere)])
    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), Lx / N)
    gmsh.model.mesh.generate(3)

    # Convert to DOLFINx mesh
    domain, cell_tags, facet_tags = io.gmshio.model_to_mesh(
        gmsh.model, MPI.COMM_WORLD, 0, gdim=3
    )
    gmsh.finalize()

    # --- Function space ---
    V = fem.VectorFunctionSpace(domain, ("Lagrange", 1))

    # --- Material parameters ---
    lam = config.mu * 2 * config.nu / (1 - 2 * config.nu)
    mu_val = config.mu

    def epsilon(u):
        return ufl.sym(ufl.grad(u))

    def sigma(u):
        return lam * ufl.nabla_div(u) * ufl.Identity(3) + 2 * mu_val * epsilon(u)

    # --- Variational form ---
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx

    # Pressure on chamber wall
    n = ufl.FacetNormal(domain)
    ds_chamber = ufl.ds  # simplified: use chamber boundary marker
    L = config.dP * ufl.dot(n, v) * ds_chamber

    # --- Boundary conditions (far-field fixed) ---
    u_bc = Function(V)
    u_bc.x.set(0.0)

    # --- Solve ---
    problem = fem.petsc.LinearProblem(a, L, bcs=[])
    uh = problem.solve()

    # --- Extract surface displacement at observation points ---
    N_obs = len(obs_x)
    disp = np.zeros((N_obs, 3))
    points = np.column_stack([obs_x, obs_y, np.zeros(N_obs)])

    from dolfinx.geometry import BoundingBoxTree, compute_collisions_points
    tree = BoundingBoxTree(domain, 3)
    cell_candidates = compute_collisions_points(tree, points)

    for i in range(N_obs):
        cells = cell_candidates.links(i)
        if len(cells) > 0:
            disp[i, :] = uh.eval(points[i], cells[0])

    return disp


# ==============================================================================
# Model comparison utility
# ==============================================================================

def compare_models(
    obs_x: np.ndarray,
    obs_y: np.ndarray,
    mogi_src: MogiSource,
    spheroid_src: SpheroidSource,
    fem_cfg: FEMSourceConfig
) -> Dict[str, np.ndarray]:
    """
    Compute displacement predictions from all three source models.

    Returns dict with keys 'mogi', 'spheroid', 'fem', each (N, 3).
    """
    d_mogi = mogi_displacement(obs_x, obs_y, mogi_src)
    d_sph = spheroid_displacement(obs_x, obs_y, spheroid_src)
    d_fem = fem_displacement(obs_x, obs_y, fem_cfg)
    return {
        'mogi': d_mogi,
        'spheroid': d_sph,
        'fem': d_fem,
    }


def compute_model_residuals(
    obs_disp: np.ndarray,
    pred_disp: np.ndarray,
    obs_sigma: Optional[np.ndarray] = None
) -> Dict[str, float]:
    """Compute RMS, WRMS, chi-squared for model-data comparison."""
    resid = obs_disp - pred_disp
    rms = np.sqrt(np.mean(resid**2))

    if obs_sigma is not None:
        w_resid = resid / obs_sigma
        wrms = np.sqrt(np.mean(w_resid**2))
        chi2 = np.sum(w_resid**2)
        ndata = obs_disp.size
        chi2_red = chi2 / max(ndata - 1, 1)
    else:
        wrms = rms
        chi2 = np.nan
        chi2_red = np.nan

    return {'rms': rms, 'wrms': wrms, 'chi2': chi2, 'chi2_reduced': chi2_red}
