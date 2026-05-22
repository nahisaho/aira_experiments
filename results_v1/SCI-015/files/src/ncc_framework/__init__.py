"""
NCC Framework — Information-Theoretic Analysis of Neural Correlates of Consciousness
"""
from .iit import PhiCalculator
from .pci import PCISimulator
from .gwt import GlobalWorkspaceAnalyzer
from .clinical import ConsciousnessClassifier
from .utils import generate_anesthesia_data, mutual_information, entropy

__version__ = "1.0.0"
__all__ = [
    "PhiCalculator",
    "PCISimulator",
    "GlobalWorkspaceAnalyzer",
    "ConsciousnessClassifier",
    "generate_anesthesia_data",
    "mutual_information",
    "entropy",
]
