"""
This module contains classes that represent an element of every type of structure.
"""
import numpy as np
from femnodes import *
from femexceptions import *
from femessentials import *
from femintermediateloads import *
from abc import ABC, abstractmethod


class BaseElement(ABC):
    """A base class for a finite element method element."""
    def __init__(self, element_id: int, n_dims: int, nodes: list[Node], material: Material,
                 properties: BaseProperties, intermediate_loads: list[BaseInterLoad] = None):
        """Constructor that declares and initializes the object's variables"""
        self.element_id = element_id
        if not len(nodes) == 2:
            raise FemDefinitionError("Elements are defined by 2 nodes.")
        self.n_dims = n_dims
        self.nodes = nodes
        self.material = material
        self.properties = properties
        self.intermediate_loads = intermediate_loads
        self.n_element_nodes = len(self.nodes)
        self.index = self._assign_index()
        self.length = self._calculate_length()
        self.transformation_matrix = self._calculate_transformation_matrix()
        self.k = self._calculate_k()
        self.kbar = self._calculate_kbar()

        self.u = None  #temp
        self.ar = self._calculate_ar()
        self.arbar = self._calculate_arbar()
        self.f = None #forces at the end of the element.

    def get_element_id(self) -> int:
        return self.element_id

    def get_element_contribution_to_kg(self, n_nodes: int) -> np.ndarray[np.float64]:
        """Returns the contribution of the element to the global stiffness matrix."""
        kg = np.zeros(shape=[self.n_dims * n_nodes, self.n_dims*n_nodes], dtype=np.float64)
        kg[np.ix_(self.index, self.index)] = self.kbar
        return kg

    def get_element_contribution_to_s(self, n_nodes: int) -> np.ndarray[np.float64]:
        """Returns the contribution of the element to the global consolidation actions vector."""
        s = np.zeros(shape=[self.n_dims * n_nodes, 1], dtype=np.float64)
        s[self.index] = self.arbar
        return s

    def _assign_index(self) -> list[int]:
        """Returns the global degrees of freedom of the element's nodes."""
        return self.nodes[0].get_index() + self.nodes[1].get_index()

    @abstractmethod
    def _calculate_length(self) -> np.ndarray[np.float64]:
        """Abstract method that calculates the length of the element."""
        pass

    @abstractmethod
    def _calculate_transformation_matrix(self) -> np.ndarray[np.float64]:
        """Abstract method that returns the transformation matrix."""
        pass

    @abstractmethod
    def _calculate_k(self) -> np.ndarray[np.float64]:
        """Abstract method that calculates the element's local stiffness matrix."""
        pass

    def _calculate_kbar(self) -> np.ndarray[np.float64]:
        """Returns the element's stiffness matrix transformed to the global coordinates system."""
        return self.transformation_matrix.transpose() @ self.k @ self.transformation_matrix

    def _calculate_ar(self) -> np.ndarray:
        """Gets the consolidation actions of this element from its intermediate loads. """
        if not self.intermediate_loads:
            return np.zeros(shape=[self.n_element_nodes * self.n_dims, 1])
        ar = np.zeros(shape=[self.n_element_nodes * self.n_dims, 1])
        for intermediate_load in self.intermediate_loads:
            ar += intermediate_load.get_consolidation_actions(self)
        return ar

    def _calculate_arbar(self) -> np.ndarray:
        """Returns the consolidation actions of this element transformed to the global coordinates system."""
        return self.transformation_matrix.transpose() @ self.ar

    def set_extreme_displacements(self):
        """Sets the displacements of the element's nodes after the analysis"""
        displacements = np.zeros(shape=[self.n_element_nodes * self.n_dims, 1], dtype=np.float64)
        for inode, node in enumerate(self.nodes):
            for idof in range(self.n_dims):
                displacements[self.n_dims * inode + idof] = node.get_displacements()[idof]
        self.u = displacements

    def calculate_element_forces(self):
        """Calculates forces that are applied at the element's edge."""
        self.set_extreme_displacements()
        self.f = self.ar + self.k @ self.transformation_matrix @ self.u

    def add_intermediate_load(self, intermediate_load: BaseInterLoad):
        """Adds an intermediate load to the element."""
        if self.intermediate_loads is None:
            self.intermediate_loads = []
        self.intermediate_loads.append(intermediate_load)
        self.ar = self._calculate_ar()
        self.arbar = self._calculate_arbar()

    """Getters for the class' variables."""
    def get_properties(self) -> BaseProperties:
        return self.properties

    def get_material(self) -> Material:
        return self.material

    def get_index(self) -> np.ndarray[np.float64]:
        return self.index

    def get_length(self) -> np.float64:
        return self.length

    def get_transformation_matrix(self) -> np.ndarray[np.float64]:
        return self.transformation_matrix

    def get_k(self) -> np.ndarray[np.float64]:
        return self.k

    def get_kbar(self) -> np.ndarray[np.float64]:
        return self.kbar

    def get_ar(self) -> np.ndarray[np.float64]:
        return self.ar

    def get_arbar(self) -> np.ndarray[np.float64]:
        return self.arbar

    def get_nodes(self) -> list[Node]:
        return self.nodes


class PlanarTrussElement(BaseElement):
    """Class of a planar truss element."""
    def __init__(self, element_id: int, n_dims: int, nodes: list[Node], material: Material,
                 properties: PlanarTrussProperties, intermediate_loads: list[BaseInterLoad] = None):
        super().__init__(element_id, n_dims, nodes, material, properties, intermediate_loads)

    def _calculate_length(self) -> np.float64:
        """Implements the abstract method and calculates the length of a planar truss."""
        return np.hypot(self.nodes[1].get_coordinates()[0] - self.nodes[0].get_coordinates()[0],
                        self.nodes[1].get_coordinates()[1] - self.nodes[0].get_coordinates()[1])

    def _calculate_transformation_matrix(self) -> np.ndarray[np.float64]:
        """Calculates the transformation matrix of a planar truss element"""
        dx = self.nodes[1].get_coordinates()[0] - self.nodes[0].get_coordinates()[0]
        dy = self.nodes[1].get_coordinates()[1] - self.nodes[0].get_coordinates()[1]
        sin_phi = dy/self.length
        cos_phi = dx/self.length
        if self.n_dims == 2:
            return np.array([[cos_phi, sin_phi, 0,       0],
                             [-sin_phi, cos_phi, 0,       0],
                             [0,       0,       cos_phi, sin_phi],
                             [0,       0,      -sin_phi, cos_phi]], dtype=np.float64)
        elif self.n_dims == 3:
            return np.array([[cos_phi, sin_phi, 0, 0, 0, 0],
                             [-sin_phi, cos_phi, 0, 0, 0, 0],
                             [0, 0, 1, 0, 0, 0],
                             [0, 0, 0, cos_phi, sin_phi, 0],
                             [0, 0, 0, -sin_phi, cos_phi, 0],
                             [0, 0, 0, 0, 0, 1]], dtype=np.float64)

    def _calculate_k(self) -> np.ndarray[np.float64]:
        """Calculates the stiffness matrix of a planar truss element."""
        E = self.material.get_e_young()
        A = self.properties.get_area()
        L = self.length
        if self.n_dims == 2:
            return E*A/L * np.array([[1, 0, -1, 0],
                                     [0, 0,  0, 0],
                                    [-1, 0,  1, 0],
                                     [0, 0,  0, 0]])
        if self.n_dims == 3:
            return E * A / L * np.array([[1, 0, 0,-1, 0, 0],
                                         [0, 0, 0, 0, 0, 0],
                                         [0, 0, 0, 0, 0, 0],
                                        [-1, 0, 0, 1, 0, 0],
                                         [0, 0, 0, 0, 0, 0],
                                         [0, 0, 0, 0, 0, 0]])


    def calculate_element_forces(self):
        """(Temporary) Returns the axial force of the element."""
        super().calculate_element_forces()
        self.axial_force = -self.f[0]

    def print_internal_forces(self):
        """Temp method to get the outputs."""
        print("Element {:2d} axial load is {:9.4f} kN".format(self.element_id, float(self.axial_force)))


class PlanarBeamElement(BaseElement):

    def __init__(self, element_id: int, n_dims: int, nodes: list[Node], material: Material,
                 properties: PlanarBeamProperties, intermediate_loads: list[BaseInterLoad] = None):
        super().__init__(element_id, n_dims, nodes, material, properties, intermediate_loads)

    def _calculate_length(self) -> np.ndarray[np.float64]:
        """Implements the abstract method and calculates the length of a planar beam."""
        return np.hypot(self.nodes[1].get_coordinates()[0] - self.nodes[0].get_coordinates()[0],
                        self.nodes[1].get_coordinates()[1] - self.nodes[0].get_coordinates()[1])

    def _calculate_transformation_matrix(self) -> np.ndarray[np.float64]:
        """Calculates the transformation matrix of a planar beam element"""
        dx = self.nodes[1].get_coordinates()[0] - self.nodes[0].get_coordinates()[0]
        dy = self.nodes[1].get_coordinates()[1] - self.nodes[0].get_coordinates()[1]
        sin_phi = dy/self.length
        cos_phi = dx/self.length
        return np.array([[cos_phi, sin_phi, 0,        0,      0, 0],
                         [-sin_phi, cos_phi, 0,        0,       0, 0],
                         [0,       0,        1,        0,       0, 0],
                         [0,       0,        0,  cos_phi, sin_phi, 0],
                         [0,       0,        0, -sin_phi, cos_phi, 0],
                         [0,       0,        0,        0,       0, 1]], dtype=np.float64)

    def _calculate_k(self) -> np.ndarray[np.float64]:
        """Calculates the stiffness matrix of a planar truss element."""
        E = self.material.get_e_young()
        A = self.properties.get_area()
        I = self.properties.get_moment_of_inertia()
        L = self.length
        return np.array([[E*A/L,         0,          0,        -E*A/L,         0,          0],
                         [0,       12*E*I/L/L/L,   6*E*I/L/L,      0,-12*E*I/L/L/L,  6*E*I/L/L],
                         [0,          6*E*I/L/L,     4*E*I/L,      0,   -6*E*I/L/L,    2*E*I/L],
                         [-E*A/L,             0,           0,  E*A/L,           0,           0],
                         [0,      -12*E*I/L/L/L,  -6*E*I/L/L,      0, 12*E*I/L/L/L, -6*E*I/L/L],
                         [0,          6*E*I/L/L,     2*E*I/L,      0,   -6*E*I/L/L,    4*E*I/L]], dtype=np.float64)

