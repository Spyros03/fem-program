import numpy as np
from abc import ABC

from femexceptions import *

"""Definition of boundary conditions."""


class BaseSupport(ABC):
    """Base class of a structure's support."""
    def __init__(self, restrictions: list[bool], angles: np.ndarray[np.float64] = None,
                 springs: np.ndarray[np.float64] = None, retreats: np.ndarray[np.float64] = None):
        self.restrictions = restrictions
        self.angles = angles
        self.springs = springs
        self.retreats = retreats

    def get_restrictions(self) -> list[bool]:
        return self.restrictions

    def get_angles(self) -> np.ndarray[np.float64]:
        return self.angles

    def get_springs(self) -> np.ndarray[np.float64]:
        return self.springs

    def get_retreats(self) -> np.ndarray[np.float64]:
        return self.retreats


class PlanarTrussSupport(BaseSupport):
    """A class that resembles a planar truss support."""

    def __init__(self, restrictions: list[bool], angles: np.ndarray[np.float64] = None,
                 springs: np.ndarray[np.float64] = None, retreats: np.ndarray[np.float64] = None):
        super().__init__(restrictions, angles, springs, retreats)


"""Definition of FEM nodes."""


class Node:
    """A class representing a node of a structure."""
    def __init__(self, node_id: int, n_dims: int, coordinates: np.ndarray[np.float64] = None,
                 external_loads: np.ndarray[np.float64] = None, support: BaseSupport = None):
        if not len(coordinates) == n_dims:
            raise FemDimensionError("Planar truss nodes are defined by 2 coordinates.")
        if external_loads is not None:
            if not len(external_loads) == n_dims:
                raise FemDimensionError("Planar truss nodes are loaded in 2 directions.")
        self.node_id = node_id
        self.coordinates = coordinates
        self.p = external_loads  #External loads vector
        self.support = support
        self.n_dims = n_dims
        self.index = self._assign_indexes()
        self.u = None  #Displacement vector
        self.bounded_dofs = self._assign_bounded_dofs()

    def _assign_indexes(self) -> list[int]:
        """Returns the global degrees of freedom that are contained in this node."""
        index = list()
        for idof in range(self.n_dims):
            index.append(self.n_dims * self.node_id + idof)
        return index

    def _assign_bounded_dofs(self) -> list[int]:
        """Returns the global degrees of freedom that are bounded by the support of this node."""
        if self.support is None:
            return list()
        bounded_dofs = list()
        for idof, restriction in enumerate(self.support.get_restrictions()):
            if restriction:
                bounded_dofs.append(self.n_dims * self.node_id + idof)
        return bounded_dofs

    def get_index(self) -> list[int]:
        return self.index

    def get_bounded_dofs(self) -> list[int]:
        return self.bounded_dofs

    def get_springs(self):
        if self.get_support() is None:
            return None
        return self.get_support().get_springs()

    def get_spring_contribution(self, n_nodes: int) -> np.ndarray[np.float64]:
        """Gets the contribution of the elastic support to the global stiffness matrix."""
        km = np.zeros(shape=[self.n_dims * n_nodes, self.n_dims * n_nodes], dtype=np.float64)
        if self.support.get_springs() is None:
            return km
        for idof in range(len(self.support.get_springs())):
            km[self.n_dims * self.node_id + idof, self.n_dims * self.node_id + idof] = self.support.get_springs()[idof]
        return km

    def get_bound_angles(self) -> np.ndarray[np.float64]:
        """Returns the angle of the support of the node."""
        return self.support.get_angles()

    def set_external_load(self, external_loads: np.ndarray[np.float64]):
        self.p = external_loads

    def set_displacements(self, u: np.ndarray[np.float64]):
        """Setter for nodal displacement vector."""
        self.u = u

    def set_support(self, support: BaseSupport):
        """Setter for the support of the node."""
        self.support = support
        self.bounded_dofs = self._assign_bounded_dofs()

    def get_displacements(self):
        return self.u

    def get_coordinates(self) -> np.ndarray[np.float64]:
        return self.coordinates

    def get_external_load(self) -> np.ndarray[np.float64]:
        return self.p

    def get_support(self) -> BaseSupport:
        return self.support
