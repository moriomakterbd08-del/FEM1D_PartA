"""1D Finite Element Method package for MSc Thesis."""

from .mesh import Mesh1D
from .quadrature import GaussLegendre1D
from .basis import LinearBasis1D
from .assembly import (
    compute_element_stiffness,
    compute_element_load,
    assemble_poisson_system,
)
from .boundary import (
    apply_homogeneous_dirichlet,
    reconstruct_full_solution,
)
from .solver import solve_poisson_1d
from .errors import FEMErrors, compute_errors
from .manufactured import (
    ManufacturedSolution,
    get_sin_benchmark,
    get_exponential_benchmark,
    get_quadratic_benchmark,
)

__all__ = [
    "Mesh1D",
    "GaussLegendre1D",
    "LinearBasis1D",
    "compute_element_stiffness",
    "compute_element_load",
    "assemble_poisson_system",
    "apply_homogeneous_dirichlet",
    "reconstruct_full_solution",
    "solve_poisson_1d",
    "FEMErrors",
    "compute_errors",
    "ManufacturedSolution",
    "get_sin_benchmark",
    "get_exponential_benchmark",
    "get_quadratic_benchmark",
]
