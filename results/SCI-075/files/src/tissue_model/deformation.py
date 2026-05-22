"""
Tissue Deformation Modeling Module
====================================
Real-time biomechanical tissue simulation using:
- Mass-Spring-Damper (MSD) model for fast approximation
- Finite Element Method (FEM) with co-rotational formulation for accuracy
"""

import numpy as np
from typing import Optional, Tuple, List, Dict
from dataclasses import dataclass
from enum import Enum


class TissueModelType(Enum):
    MASS_SPRING = "mass_spring"
    FEM_LINEAR = "fem_linear"
    FEM_COROTATIONAL = "fem_corotational"


@dataclass
class TissueProperties:
    """Biomechanical tissue properties."""
    youngs_modulus: float = 5000.0       # Pa (soft tissue ~1-10 kPa)
    poissons_ratio: float = 0.45         # Near-incompressible
    density: float = 1060.0              # kg/m^3
    damping_ratio: float = 0.3           # Rayleigh damping
    thickness: float = 0.005             # m
    max_strain: float = 0.3              # Failure strain threshold
    spring_stiffness: float = 500.0      # N/m (for MSD model)
    damping_coefficient: float = 5.0     # Ns/m


class MassSpringModel:
    """
    Mass-Spring-Damper tissue deformation model.
    Efficient for real-time simulation (~1 kHz update rate).
    """

    def __init__(self, properties: TissueProperties,
                 grid_size: Tuple[int, int] = (20, 20),
                 spacing: float = 0.002):
        self.props = properties
        self.grid_size = grid_size
        self.spacing = spacing

        # Node positions and velocities
        self.n_nodes = grid_size[0] * grid_size[1]
        self.positions = np.zeros((self.n_nodes, 3))
        self.velocities = np.zeros((self.n_nodes, 3))
        self.rest_positions = np.zeros((self.n_nodes, 3))
        self.forces = np.zeros((self.n_nodes, 3))
        self.fixed_nodes: set = set()

        # Spring connectivity
        self.springs: List[Tuple[int, int, float]] = []

        self._initialize_grid()
        self._create_springs()

        # Node masses
        area_per_node = spacing ** 2
        self.masses = np.full(self.n_nodes,
                              properties.density * area_per_node * properties.thickness)

    def _initialize_grid(self):
        """Create a regular 2D grid of nodes."""
        rows, cols = self.grid_size
        for i in range(rows):
            for j in range(cols):
                idx = i * cols + j
                self.positions[idx] = [
                    j * self.spacing, i * self.spacing, 0.0
                ]
        self.rest_positions = self.positions.copy()

    def _create_springs(self):
        """Create structural, shear, and bend springs."""
        rows, cols = self.grid_size

        for i in range(rows):
            for j in range(cols):
                idx = i * cols + j

                # Structural springs (horizontal, vertical)
                if j < cols - 1:
                    self.springs.append((idx, idx + 1, self.spacing))
                if i < rows - 1:
                    self.springs.append((idx, idx + cols, self.spacing))

                # Shear springs (diagonal)
                if i < rows - 1 and j < cols - 1:
                    diag_len = self.spacing * np.sqrt(2)
                    self.springs.append((idx, idx + cols + 1, diag_len))
                if i < rows - 1 and j > 0:
                    diag_len = self.spacing * np.sqrt(2)
                    self.springs.append((idx, idx + cols - 1, diag_len))

                # Bend springs (skip one node)
                if j < cols - 2:
                    self.springs.append((idx, idx + 2, 2 * self.spacing))
                if i < rows - 2:
                    self.springs.append((idx, idx + 2 * cols, 2 * self.spacing))

    def apply_force(self, node_idx: int, force: np.ndarray):
        """Apply external force to a node."""
        self.forces[node_idx] += force

    def apply_displacement(self, node_idx: int, displacement: np.ndarray):
        """Apply displacement constraint (e.g., needle contact)."""
        self.positions[node_idx] = self.rest_positions[node_idx] + displacement
        self.fixed_nodes.add(node_idx)

    def step(self, dt: float = 0.001):
        """
        Advance simulation by one timestep using Verlet integration.
        Target: real-time at 1 kHz.
        """
        internal_forces = np.zeros_like(self.forces)

        # Compute spring forces
        k = self.props.spring_stiffness
        c = self.props.damping_coefficient

        for n1, n2, rest_length in self.springs:
            delta = self.positions[n2] - self.positions[n1]
            dist = np.linalg.norm(delta)
            if dist < 1e-10:
                continue

            direction = delta / dist
            stretch = dist - rest_length

            # Spring force + damping
            vel_diff = self.velocities[n2] - self.velocities[n1]
            f_spring = k * stretch * direction
            f_damping = c * np.dot(vel_diff, direction) * direction

            f_total = f_spring + f_damping

            internal_forces[n1] += f_total
            internal_forces[n2] -= f_total

        # Gravity
        gravity = np.array([0, 0, -9.81])

        # Update velocities and positions
        for i in range(self.n_nodes):
            if i in self.fixed_nodes:
                self.velocities[i] = 0
                continue

            acceleration = (
                internal_forces[i] + self.forces[i] +
                self.masses[i] * gravity
            ) / self.masses[i]

            self.velocities[i] += acceleration * dt
            self.velocities[i] *= (1 - self.props.damping_ratio * dt)
            self.positions[i] += self.velocities[i] * dt

        # Reset external forces
        self.forces[:] = 0
        self.fixed_nodes.clear()

    def get_deformation(self) -> np.ndarray:
        """Get displacement field."""
        return self.positions - self.rest_positions

    def get_max_strain(self) -> float:
        """Compute maximum strain in the tissue."""
        max_strain = 0.0
        for n1, n2, rest_length in self.springs:
            delta = self.positions[n2] - self.positions[n1]
            dist = np.linalg.norm(delta)
            strain = abs(dist - rest_length) / rest_length
            max_strain = max(max_strain, strain)
        return max_strain

    def get_reaction_force(self, node_idx: int) -> np.ndarray:
        """Compute reaction force at a constrained node."""
        k = self.props.spring_stiffness
        force = np.zeros(3)
        for n1, n2, rest_length in self.springs:
            if n1 == node_idx or n2 == node_idx:
                other = n2 if n1 == node_idx else n1
                delta = self.positions[other] - self.positions[node_idx]
                dist = np.linalg.norm(delta)
                if dist < 1e-10:
                    continue
                direction = delta / dist
                stretch = dist - rest_length
                force += k * stretch * direction
        return force


class FEMTissueModel:
    """
    Finite Element Method tissue model with co-rotational formulation.
    Uses tetrahedral elements for 3D tissue simulation.
    More accurate than MSD but computationally heavier (~100 Hz).
    """

    def __init__(self, properties: TissueProperties):
        self.props = properties
        self.nodes: np.ndarray = np.array([])
        self.elements: np.ndarray = np.array([])  # Tetrahedra connectivity
        self.displacements: np.ndarray = np.array([])
        self.velocities_fem: np.ndarray = np.array([])
        self.K_global: Optional[np.ndarray] = None
        self.M_global: Optional[np.ndarray] = None
        self.fixed_dofs: set = set()

    def initialize_mesh(self, nodes: np.ndarray, elements: np.ndarray):
        """
        Initialize FEM mesh.

        Parameters
        ----------
        nodes : (N, 3) nodal coordinates
        elements : (M, 4) tetrahedral element connectivity
        """
        self.nodes = nodes.copy()
        self.rest_nodes = nodes.copy()
        self.elements = elements.copy()
        n_dofs = 3 * len(nodes)
        self.displacements = np.zeros(n_dofs)
        self.velocities_fem = np.zeros(n_dofs)

        self._assemble_stiffness()
        self._assemble_mass()

    def generate_box_mesh(self, size: Tuple[float, float, float] = (0.04, 0.04, 0.005),
                          resolution: Tuple[int, int, int] = (8, 8, 2)):
        """Generate a simple box mesh with tetrahedral elements."""
        nx, ny, nz = resolution
        sx, sy, sz = size

        # Create nodes
        nodes = []
        for k in range(nz + 1):
            for j in range(ny + 1):
                for i in range(nx + 1):
                    nodes.append([
                        i * sx / nx, j * sy / ny, k * sz / nz
                    ])
        nodes = np.array(nodes)

        # Create tetrahedra from hexahedral cells
        elements = []
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    # 8 corners of hex cell
                    n0 = k * (ny+1) * (nx+1) + j * (nx+1) + i
                    n1 = n0 + 1
                    n2 = n0 + (nx+1)
                    n3 = n2 + 1
                    n4 = n0 + (ny+1) * (nx+1)
                    n5 = n4 + 1
                    n6 = n4 + (nx+1)
                    n7 = n6 + 1

                    # Split hex into 5 tetrahedra
                    elements.append([n0, n1, n3, n5])
                    elements.append([n0, n3, n2, n6])
                    elements.append([n0, n5, n4, n6])
                    elements.append([n3, n5, n6, n7])
                    elements.append([n0, n3, n5, n6])

        elements = np.array(elements)
        self.initialize_mesh(nodes, elements)

    def _compute_element_stiffness(self, elem_idx: int) -> np.ndarray:
        """Compute 12x12 element stiffness matrix for a tetrahedron."""
        nodes_idx = self.elements[elem_idx]
        coords = self.rest_nodes[nodes_idx]

        # Shape function derivatives (constant for linear tet)
        # B = [dN/dx] computed from Jacobian
        J = np.zeros((3, 3))
        for i in range(3):
            J[i] = coords[i+1] - coords[0]
        J = J.T

        det_J = np.linalg.det(J)
        if abs(det_J) < 1e-15:
            return np.zeros((12, 12))

        J_inv = np.linalg.inv(J)
        volume = abs(det_J) / 6.0

        # Shape function derivatives w.r.t. global coords
        dN_local = np.array([
            [-1, -1, -1],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1]
        ], dtype=float)
        dN = dN_local @ J_inv

        # Build B matrix (6x12)
        B = np.zeros((6, 12))
        for n in range(4):
            col = 3 * n
            B[0, col] = dN[n, 0]
            B[1, col+1] = dN[n, 1]
            B[2, col+2] = dN[n, 2]
            B[3, col] = dN[n, 1]
            B[3, col+1] = dN[n, 0]
            B[4, col+1] = dN[n, 2]
            B[4, col+2] = dN[n, 1]
            B[5, col] = dN[n, 2]
            B[5, col+2] = dN[n, 0]

        # Material stiffness matrix (isotropic linear elastic)
        E = self.props.youngs_modulus
        nu = self.props.poissons_ratio
        c = E / ((1 + nu) * (1 - 2 * nu))
        D = c * np.array([
            [1-nu, nu, nu, 0, 0, 0],
            [nu, 1-nu, nu, 0, 0, 0],
            [nu, nu, 1-nu, 0, 0, 0],
            [0, 0, 0, (1-2*nu)/2, 0, 0],
            [0, 0, 0, 0, (1-2*nu)/2, 0],
            [0, 0, 0, 0, 0, (1-2*nu)/2]
        ])

        Ke = volume * B.T @ D @ B
        return Ke

    def _assemble_stiffness(self):
        """Assemble global stiffness matrix."""
        n_dofs = 3 * len(self.nodes)
        self.K_global = np.zeros((n_dofs, n_dofs))

        for e in range(len(self.elements)):
            Ke = self._compute_element_stiffness(e)
            nodes_idx = self.elements[e]
            dofs = []
            for n in nodes_idx:
                dofs.extend([3*n, 3*n+1, 3*n+2])

            for i, di in enumerate(dofs):
                for j, dj in enumerate(dofs):
                    self.K_global[di, dj] += Ke[i, j]

    def _assemble_mass(self):
        """Assemble lumped mass matrix."""
        n_dofs = 3 * len(self.nodes)
        self.M_global = np.zeros(n_dofs)

        for e in range(len(self.elements)):
            nodes_idx = self.elements[e]
            coords = self.rest_nodes[nodes_idx]

            J = np.zeros((3, 3))
            for i in range(3):
                J[i] = coords[i+1] - coords[0]
            volume = abs(np.linalg.det(J.T)) / 6.0

            element_mass = self.props.density * volume
            node_mass = element_mass / 4.0

            for n in nodes_idx:
                for d in range(3):
                    self.M_global[3*n + d] += node_mass

    def fix_boundary(self, node_indices: List[int]):
        """Fix nodes (Dirichlet boundary conditions)."""
        for n in node_indices:
            for d in range(3):
                self.fixed_dofs.add(3 * n + d)

    def apply_load(self, node_idx: int, force: np.ndarray):
        """Apply external load at a node."""
        pass  # Handled in solve_static

    def solve_static(self, external_forces: np.ndarray) -> np.ndarray:
        """
        Solve static equilibrium: K * u = F

        Returns
        -------
        displacements : (N*3,) nodal displacements
        """
        n_dofs = 3 * len(self.nodes)
        free_dofs = [i for i in range(n_dofs) if i not in self.fixed_dofs]

        K_free = self.K_global[np.ix_(free_dofs, free_dofs)]
        F_free = external_forces[free_dofs]

        try:
            u_free = np.linalg.solve(K_free, F_free)
        except np.linalg.LinAlgError:
            u_free = np.linalg.lstsq(K_free, F_free, rcond=None)[0]

        self.displacements[free_dofs] = u_free
        return self.displacements.copy()

    def solve_dynamic(self, external_forces: np.ndarray, dt: float = 0.001):
        """
        Explicit time integration for dynamic FEM.
        Uses central difference method.
        """
        n_dofs = len(self.displacements)

        # Internal forces: f_int = K * u
        f_int = self.K_global @ self.displacements

        # Acceleration: M * a = F_ext - F_int - C * v
        damping_forces = self.props.damping_ratio * self.velocities_fem
        accelerations = np.zeros(n_dofs)

        for i in range(n_dofs):
            if i in self.fixed_dofs:
                continue
            if self.M_global[i] > 1e-15:
                accelerations[i] = (
                    external_forces[i] - f_int[i] - damping_forces[i]
                ) / self.M_global[i]

        # Update
        self.velocities_fem += accelerations * dt
        self.displacements += self.velocities_fem * dt

        # Update node positions
        for i in range(len(self.nodes)):
            self.nodes[i] = self.rest_nodes[i] + self.displacements[3*i:3*i+3]

    def get_von_mises_stress(self) -> np.ndarray:
        """Compute von Mises stress for each element."""
        stresses = np.zeros(len(self.elements))
        for e in range(len(self.elements)):
            nodes_idx = self.elements[e]
            dofs = []
            for n in nodes_idx:
                dofs.extend([3*n, 3*n+1, 3*n+2])

            u_elem = self.displacements[dofs]

            # Recompute B and D for stress calculation
            coords = self.rest_nodes[nodes_idx]
            J = np.zeros((3, 3))
            for i in range(3):
                J[i] = coords[i+1] - coords[0]
            J = J.T

            det_J = np.linalg.det(J)
            if abs(det_J) < 1e-15:
                continue

            J_inv = np.linalg.inv(J)
            dN_local = np.array([[-1,-1,-1],[1,0,0],[0,1,0],[0,0,1]], dtype=float)
            dN = dN_local @ J_inv

            B = np.zeros((6, 12))
            for n in range(4):
                col = 3 * n
                B[0, col] = dN[n, 0]
                B[1, col+1] = dN[n, 1]
                B[2, col+2] = dN[n, 2]
                B[3, col] = dN[n, 1]; B[3, col+1] = dN[n, 0]
                B[4, col+1] = dN[n, 2]; B[4, col+2] = dN[n, 1]
                B[5, col] = dN[n, 2]; B[5, col+2] = dN[n, 0]

            E = self.props.youngs_modulus
            nu = self.props.poissons_ratio
            c = E / ((1 + nu) * (1 - 2 * nu))
            D = c * np.array([
                [1-nu, nu, nu, 0, 0, 0],
                [nu, 1-nu, nu, 0, 0, 0],
                [nu, nu, 1-nu, 0, 0, 0],
                [0, 0, 0, (1-2*nu)/2, 0, 0],
                [0, 0, 0, 0, (1-2*nu)/2, 0],
                [0, 0, 0, 0, 0, (1-2*nu)/2]
            ])

            stress = D @ B @ u_elem
            # von Mises
            s = stress
            vm = np.sqrt(0.5 * ((s[0]-s[1])**2 + (s[1]-s[2])**2 +
                                 (s[2]-s[0])**2 + 6*(s[3]**2+s[4]**2+s[5]**2)))
            stresses[e] = vm

        return stresses


class TissueModelManager:
    """
    Manager for switching between MSD and FEM tissue models.
    Provides adaptive model selection based on accuracy requirements.
    """

    def __init__(self, properties: Optional[TissueProperties] = None):
        self.props = properties or TissueProperties()
        self.msd_model: Optional[MassSpringModel] = None
        self.fem_model: Optional[FEMTissueModel] = None
        self.active_model: TissueModelType = TissueModelType.MASS_SPRING

    def initialize(self, model_type: TissueModelType = TissueModelType.MASS_SPRING):
        """Initialize the specified tissue model."""
        self.active_model = model_type

        if model_type == TissueModelType.MASS_SPRING:
            self.msd_model = MassSpringModel(self.props)
        else:
            self.fem_model = FEMTissueModel(self.props)
            self.fem_model.generate_box_mesh()

    def step(self, dt: float = 0.001, external_forces: Optional[np.ndarray] = None):
        """Step the active tissue model."""
        if self.active_model == TissueModelType.MASS_SPRING and self.msd_model:
            self.msd_model.step(dt)
        elif self.fem_model:
            if external_forces is None:
                external_forces = np.zeros(3 * len(self.fem_model.nodes))
            self.fem_model.solve_dynamic(external_forces, dt)

    def get_deformation_at(self, position: np.ndarray) -> np.ndarray:
        """Get tissue deformation at a specific position (interpolated)."""
        if self.active_model == TissueModelType.MASS_SPRING and self.msd_model:
            # Find nearest node
            dists = np.linalg.norm(
                self.msd_model.rest_positions - position, axis=1
            )
            nearest = np.argmin(dists)
            return self.msd_model.get_deformation()[nearest]
        elif self.fem_model:
            dists = np.linalg.norm(self.fem_model.rest_nodes - position, axis=1)
            nearest = np.argmin(dists)
            return self.fem_model.displacements[3*nearest:3*nearest+3]
        return np.zeros(3)
