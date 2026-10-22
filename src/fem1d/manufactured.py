"""Manufactured benchmark solutions for the 1D Poisson model problem."""

from __future__ import annotations
from dataclasses import dataclass
from typing import Callable
import numpy as np


@dataclass
class ManufacturedSolution:
    """Represents a benchmark test problem with known analytical solution."""
    name: str
    description: str
    u: Callable[[np.ndarray], np.ndarray]
    u_prime: Callable[[np.ndarray], np.ndarray]
    u_double_prime: Callable[[np.ndarray], np.ndarray]
    f: Callable[[np.ndarray], np.ndarray]


def get_sin_benchmark() -> ManufacturedSolution:
    """Benchmark 1: u(x) = sin(pi * x)

    -u''(x) = f(x) = pi^2 * sin(pi * x)
    Boundary values: u(0) = 0, u(1) = 0.
    """
    pi = np.pi

    def u(x: np.ndarray) -> np.ndarray:
        return np.sin(pi * x)

    def u_prime(x: np.ndarray) -> np.ndarray:
        return pi * np.cos(pi * x)

    def u_double_prime(x: np.ndarray) -> np.ndarray:
        return -(pi**2) * np.sin(pi * x)

    def f(x: np.ndarray) -> np.ndarray:
        return (pi**2) * np.sin(pi * x)

    return ManufacturedSolution(
        name="sin_pi_x",
        description="u(x) = sin(pi * x)",
        u=u,
        u_prime=u_prime,
        u_double_prime=u_double_prime,
        f=f,
    )


def get_exponential_benchmark() -> ManufacturedSolution:
    """Benchmark 2: u(x) = x * (1 - x) * exp(x)

    u'(x) = (1 - x - x^2) * exp(x)
    u''(x) = (-1 - 3x - x^2) * exp(x)
    -u''(x) = f(x) = (x^2 + 3x + 1) * exp(x)
    Boundary values: u(0) = 0, u(1) = 0.
    """
    def u(x: np.ndarray) -> np.ndarray:
        return x * (1.0 - x) * np.exp(x)

    def u_prime(x: np.ndarray) -> np.ndarray:
        return (1.0 - x - x**2) * np.exp(x)

    def u_double_prime(x: np.ndarray) -> np.ndarray:
        return -x * (x + 3.0) * np.exp(x)

    def f(x: np.ndarray) -> np.ndarray:
        return x * (x + 3.0) * np.exp(x)

    return ManufacturedSolution(
        name="poly_exp",
        description="u(x) = x * (1 - x) * exp(x)",
        u=u,
        u_prime=u_prime,
        u_double_prime=u_double_prime,
        f=f,
    )


def get_quadratic_benchmark() -> ManufacturedSolution:
    """Benchmark 3: u(x) = x * (1 - x)

    u'(x) = 1 - 2x
    u''(x) = -2
    -u''(x) = f(x) = 2 (constant load)
    Boundary values: u(0) = 0, u(1) = 0.
    """
    def u(x: np.ndarray) -> np.ndarray:
        return x * (1.0 - x)

    def u_prime(x: np.ndarray) -> np.ndarray:
        return 1.0 - 2.0 * x

    def u_double_prime(x: np.ndarray) -> np.ndarray:
        return np.full_like(x, -2.0)

    def f(x: np.ndarray) -> np.ndarray:
        return np.full_like(x, 2.0)

    return ManufacturedSolution(
        name="quadratic",
        description="u(x) = x * (1 - x)",
        u=u,
        u_prime=u_prime,
        u_double_prime=u_double_prime,
        f=f,
    )
