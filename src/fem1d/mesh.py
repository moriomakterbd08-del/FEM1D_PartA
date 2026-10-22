"""1D Mesh generation and representation for Finite Element Method."""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class Mesh1D:
    """Represents a 1D finite element mesh on an interval [x_min, x_max].

    Attributes:
        nodes: Coordinates of vertices, shape (num_nodes,).
        elements: Element-to-node connectivity, shape (num_elements, 2).
    """
    nodes: np.ndarray
    elements: np.ndarray

    def __post_init__(self) -> None:
        self.nodes = np.asarray(self.nodes, dtype=float)
        self.elements = np.asarray(self.elements, dtype=int)
        if self.nodes.ndim != 1:
            raise ValueError("nodes must be a 1D array")
        if self.elements.ndim != 2 or self.elements.shape[1] != 2:
            raise ValueError("elements must have shape (num_elements, 2)")
        if not np.all(np.diff(self.nodes) > 0):
            raise ValueError("Node coordinates must be strictly increasing")

    @classmethod
    def uniform(cls, x_min: float = 0.0, x_max: float = 1.0, num_elements: int = 10) -> Mesh1D:
        """Constructs a uniform 1D mesh with given number of elements."""
        if num_elements < 1:
            raise ValueError("num_elements must be at least 1")
        nodes = np.linspace(x_min, x_max, num_elements + 1)
        elements = np.column_stack((np.arange(num_elements), np.arange(1, num_elements + 1)))
        return cls(nodes=nodes, elements=elements)

    @classmethod
    def from_nodes(cls, nodes: np.ndarray) -> Mesh1D:
        """Constructs a 1D mesh from an explicit array of node coordinates."""
        nodes = np.asarray(nodes, dtype=float)
        num_elements = len(nodes) - 1
        elements = np.column_stack((np.arange(num_elements), np.arange(1, num_elements + 1)))
        return cls(nodes=nodes, elements=elements)

    @property
    def num_nodes(self) -> int:
        return len(self.nodes)

    @property
    def num_elements(self) -> int:
        return len(self.elements)

    @property
    def element_lengths(self) -> np.ndarray:
        """Array of length h_e for each element."""
        return self.nodes[self.elements[:, 1]] - self.nodes[self.elements[:, 0]]

    @property
    def h_max(self) -> float:
        """Maximum element diameter (h)."""
        return float(np.max(self.element_lengths))

    @property
    def h_min(self) -> float:
        """Minimum element diameter."""
        return float(np.min(self.element_lengths))

    def get_element_coords(self, elem_idx: int) -> tuple[float, float]:
        """Returns the boundary coordinates (x_left, x_right) of element elem_idx."""
        node_left, node_right = self.elements[elem_idx]
        return float(self.nodes[node_left]), float(self.nodes[node_right])
