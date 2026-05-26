"""
Module 4: Coverage-dependent lateral interaction models.

Implements:
- First-nearest-neighbor interaction model
- Mean-field lateral interaction correction
- Coverage-dependent activation energy modification
"""
import numpy as np


class LateralInteractionModel:
    """
    Model for adsorbate-adsorbate lateral interactions on catalyst surfaces.
    """
    
    def __init__(self, species, interaction_matrix=None, coordination_number=4):
        """
        Parameters:
            species: list of adsorbate species names
            interaction_matrix: NxN matrix of pairwise interaction energies [eV]
                               positive = repulsive, negative = attractive
            coordination_number: number of nearest neighbors on the surface
        """
        self.species = species
        self.n_species = len(species)
        self.z = coordination_number
        
        if interaction_matrix is not None:
            self.epsilon = np.array(interaction_matrix)
        else:
            self.epsilon = np.zeros((self.n_species, self.n_species))
    
    def mean_field_correction(self, coverages):
        """
        Mean-field correction to binding energies.
        dE_i = sum_j(z * epsilon_ij * theta_j)
        
        Returns correction to binding energy for each species [eV].
        """
        theta = np.array(coverages)
        dE = self.z * self.epsilon @ theta
        return dE
    
    def modified_binding_energy(self, E0, coverages):
        """
        Compute coverage-dependent binding energies.
        E(theta) = E0 + dE(theta)
        """
        dE = self.mean_field_correction(coverages)
        return np.array(E0) + dE
    
    def modified_activation_energy(self, Ea0, coverages, bep_slope=0.5):
        """
        Modify activation energy using BEP relation with coverage correction.
        Ea(theta) = Ea0 + bep_slope * dE(theta)
        
        bep_slope: BEP (Brønsted-Evans-Polanyi) slope (typically 0.3-0.7)
        """
        dE = self.mean_field_correction(coverages)
        Ea_mod = Ea0 + bep_slope * np.sum(np.abs(dE))
        return max(Ea_mod, 0.0)
    
    def differential_heat_of_adsorption(self, E0, coverages, species_idx):
        """
        Compute differential heat of adsorption as function of coverage.
        q_diff = -E0_i - z * sum_j(epsilon_ij * theta_j)
        """
        theta = np.array(coverages)
        q = -E0[species_idx] - self.z * np.dot(self.epsilon[species_idx], theta)
        return q


def create_ft_interaction_matrix():
    """
    Create a representative lateral interaction matrix for
    Fischer-Tropsch key adsorbates on Co(0001).
    
    Species: CO*, H*, C*, O*, CH*, CH2*, CH3*, OH*
    Values in eV (approximate, from DFT literature).
    """
    species = ['CO*', 'H*', 'C*', 'O*', 'CH*', 'CH2*', 'CH3*', 'OH*']
    
    # Interaction energies (eV): positive = repulsive
    epsilon = np.array([
        # CO*    H*     C*     O*     CH*    CH2*   CH3*   OH*
        [0.10,  0.02,  0.05,  0.08,  0.04,  0.03,  0.02,  0.06],  # CO*
        [0.02,  0.01,  0.03,  0.02,  0.02,  0.01,  0.01,  0.02],  # H*
        [0.05,  0.03,  0.12,  0.10,  0.08,  0.06,  0.04,  0.07],  # C*
        [0.08,  0.02,  0.10,  0.15,  0.07,  0.05,  0.03,  0.10],  # O*
        [0.04,  0.02,  0.08,  0.07,  0.06,  0.04,  0.03,  0.05],  # CH*
        [0.03,  0.01,  0.06,  0.05,  0.04,  0.03,  0.02,  0.03],  # CH2*
        [0.02,  0.01,  0.04,  0.03,  0.03,  0.02,  0.02,  0.02],  # CH3*
        [0.06,  0.02,  0.07,  0.10,  0.05,  0.03,  0.02,  0.08],  # OH*
    ])
    
    return species, epsilon
