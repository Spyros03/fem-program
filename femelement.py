"""
This module contains classes that represent an element of every type of structure.
"""
import numpy as np
from femnodes import *
from femexceptions import *
from femessentials import *
# from femintermediateloads import *
from abc import ABC, abstractmethod


class BaseElement(ABC):
    """A base class for a finite element method element."""
    def __init__(self, element_id: int, n_dims: int, nodes: list[Node], material: Material,
                 properties: BaseProperties, intermediate_loads: list = None):
        """Constructor that declares and initializes the object's variables"""
        self.element_id = element_id
        if not len(nodes) == 2:
            raise FemDefinitionError("Elements are defined by 2 nodes.")
        self.n_dims = n_dims
        self.nodes = nodes
        self.n_nodes = len(nodes)
        self.material = material
        self.properties = properties
        self.intermediate_loads = intermediate_loads
        self.n_element_nodes = len(self.nodes)
        self.index = self._assign_index()
        self.length = self._calculate_length()
        self.angle = None
        self.transformation_matrix = self._calculate_transformation_matrix()
        self.e = self._get_eccentricity()
        self.k = self._calculate_k()
        self.kbar = self._calculate_kbar()
        self.u = None  #temp
        self.solid_node_actions = self._move_solid_node_loads()
        self.ar = self._calculate_ar()
        self.arbar = self._calculate_arbar()
        self.f = None #forces at the end of the element.

    def get_element_id(self) -> int:
        return self.element_id

    def get_element_contribution_to_kg(self, n_dofs: int) -> np.ndarray[np.float64]:
        """Returns the contribution of the element to the global stiffness matrix."""
        kg = np.zeros(shape=[n_dofs, n_dofs], dtype=np.float64)
        kg[np.ix_(self.index, self.index)] = self.kbar
        return kg

    def get_element_contribution_to_s(self, n_dofs: int) -> np.ndarray[np.float64]:
        """Returns the contribution of the element to the global consolidation actions vector."""
        s = np.zeros(shape=[n_dofs, 1], dtype=np.float64)
        s[self.index] = self.arbar + self.solid_node_actions
        return s

    def _assign_index(self):
        """Returns the global degrees of freedom of the element's nodes."""
        return self.nodes[0].get_index() + self.nodes[1].get_index()

    @abstractmethod
    def _calculate_length(self) -> np.ndarray[np.float64]:
        """Abstract method that calculates the length of the element."""
        pass

    @abstractmethod
    def _get_eccentricity(self):
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
        return self.transformation_matrix.transpose() @ self.e.transpose() @ self.k @ self.e @ self.transformation_matrix

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
        self.f = self.ar + self.k @ self.transformation_matrix @ self.e @ self.u

    def _move_solid_node_loads(self) -> np.ndarray:
        pass

    def add_intermediate_load(self, intermediate_load):
        """Adds an intermediate load to the element."""
        if self.intermediate_loads is None:
            self.intermediate_loads = []
        self.intermediate_loads.append(intermediate_load)
        self.ar = self._calculate_ar()
        self.arbar = self._calculate_arbar()
        self.solid_node_actions = self._move_solid_node_loads()

    """Getters for the class' variables."""
    def get_properties(self) -> BaseProperties:
        return self.properties

    def get_material(self) -> Material:
        return self.material

    def get_index(self) -> np.ndarray[np.float64]:
        return self.index

    def get_length(self) -> np.float64:
        return self.length

    def get_angle(self):
        if self.angle is not None:
            return self.angle
        self.angle = 180 / np.pi * np.arctan2(
                                        self.deformable_element_end_coords[1] - self.deformable_element_start_coords[1],
                                        self.deformable_element_end_coords[0] - self.deformable_element_start_coords[0])
        return self.angle

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
                 properties: PlanarTrussProperties, intermediate_loads: list = None):
        super().__init__(element_id, n_dims, nodes, material, properties, intermediate_loads)

    def _calculate_length(self) -> np.float64:
        """Implements the abstract method and calculates the length of a planar truss."""
        return np.hypot(self.nodes[1].get_coordinates()[0] - self.nodes[0].get_coordinates()[0],
                        self.nodes[1].get_coordinates()[1] - self.nodes[0].get_coordinates()[1])

    def _get_eccentricity(self) -> np.ndarray:
        return np.eye(self.n_nodes * self.n_dims, dtype=np.float64)

    def _move_solid_node_loads(self) -> np.ndarray:
        return np.zeros([self.n_element_nodes * self.n_dims, 1], dtype=np.float64)

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
            return E * A / L * np.array([[1, 0, 0, -1, 0, 0],
                                         [0, 0, 0,  0, 0, 0],
                                         [0, 0, 0,  0, 0, 0],
                                        [-1, 0, 0,  1, 0, 0],
                                         [0, 0, 0,  0, 0, 0],
                                         [0, 0, 0,  0, 0, 0]])

    def calculate_element_forces(self):
        """(Temporary) Returns the axial force of the element."""
        super().calculate_element_forces()
        self.axial_force = -self.f[0]

    def print_internal_forces(self):
        """Temp method to get the outputs."""
        print("Element {:2d} axial load is {:9.4f} kN".format(self.element_id, float(self.axial_force)))


class PlanarBeamElement(BaseElement):

    def __init__(self, element_id: int, n_dims: int, nodes: list[Node], material: Material,
                 properties: PlanarBeamProperties, intermediate_loads: list = None,
                 deformable_element_start_coords: np.ndarray[np.float64] = None,
                 deformable_element_end_coords: np.ndarray[np.float64] = None):
        if deformable_element_start_coords is None:
            self.deformable_element_start_coords = nodes[0].get_coordinates()
        else:
            self.deformable_element_start_coords = deformable_element_start_coords
        if deformable_element_end_coords is None:
            self.deformable_element_end_coords = nodes[1].get_coordinates()
        else:
            self.deformable_element_end_coords = deformable_element_end_coords
        super().__init__(element_id, n_dims, nodes, material, properties, intermediate_loads)

    def _calculate_length(self) -> np.ndarray[np.float64]:
        """Implements the abstract method and calculates the length of a planar beam."""
        return np.hypot(self.deformable_element_end_coords[0] - self.deformable_element_start_coords[0],
                        self.deformable_element_end_coords[1] - self.deformable_element_start_coords[1])

    def _get_eccentricity(self):
        transformation_matrix = self.get_transformation_matrix()[np.ix_([0, 1], [0, 1])]
        start_eccentricity_vector = transformation_matrix @ np.array([[self.deformable_element_start_coords[0]
                                                                       - self.nodes[0].get_coordinates()[0]],
                                                                      [self.deformable_element_start_coords[1]
                                                                       - self.nodes[0].get_coordinates()[1]]])
        end_eccentricity_vector = transformation_matrix @ np.array([[self.deformable_element_end_coords[0]
                                                                    - self.nodes[1].get_coordinates()[0]],
                                                                    [self.deformable_element_end_coords[1]
                                                                    - self.nodes[1].get_coordinates()[1]]])
        eccentricity = np.eye(6)
        eccentricity[0, 2] = -start_eccentricity_vector[1]
        eccentricity[1, 2] = start_eccentricity_vector[0]
        eccentricity[3, 5] = -end_eccentricity_vector[1]
        eccentricity[4, 5] = end_eccentricity_vector[0]
        return eccentricity

    def _calculate_transformation_matrix(self) -> np.ndarray[np.float64]:
        """Calculates the transformation matrix of a planar beam element"""
        dx = self.deformable_element_end_coords[0] - self.deformable_element_start_coords[0]
        dy = self.deformable_element_end_coords[1] - self.deformable_element_start_coords[1]
        sin_phi = dy/self.length
        cos_phi = dx/self.length
        return np.array([[cos_phi, sin_phi, 0,        0,      0, 0],
                         [-sin_phi, cos_phi, 0,        0,       0, 0],
                         [0,       0,        1,        0,       0, 0],
                         [0,       0,        0,  cos_phi, sin_phi, 0],
                         [0,       0,        0, -sin_phi, cos_phi, 0],
                         [0,       0,        0,        0,       0, 1]], dtype=np.float64)

    def _move_solid_node_loads(self) -> np.ndarray:
        if not self.intermediate_loads:
            return np.zeros(shape=[self.n_element_nodes * self.n_dims, 1])
        solid_node_actions = np.zeros(shape=[self.n_element_nodes * self.n_dims, 1])
        for intermediate_load in self.intermediate_loads:
            if type(intermediate_load) is PlanarUniformDistributedLoad:
                if intermediate_load.is_acting_on_solid_body:
                    solid_node_actions += intermediate_load.get_solid_node_actions(self)
        return solid_node_actions


    def _calculate_k(self) -> np.ndarray[np.float64]:
        """Calculates the stiffness matrix of a planar truss element."""
        E = self.material.get_e_young()
        A = self.properties.get_area()
        I = self.properties.get_moment_of_inertia()
        L = self.length
        k = np.array([[E*A/L,         0,          0,        -E*A/L,         0,          0],
                         [0,       12*E*I/L/L/L,   6*E*I/L/L,      0, -12*E*I/L/L/L,  6*E*I/L/L],
                         [0,          6*E*I/L/L,     4*E*I/L,      0,    -6*E*I/L/L,    2*E*I/L],
                         [-E*A/L,             0,           0,  E*A/L,            0,           0],
                         [0,      -12*E*I/L/L/L,  -6*E*I/L/L,      0,  12*E*I/L/L/L, -6*E*I/L/L],
                         [0,          6*E*I/L/L,     2*E*I/L,      0,    -6*E*I/L/L,    4*E*I/L]], dtype=np.float64)
        return k


class PlanarBeamElementT1(PlanarBeamElement):
    """Beam o------- hinge at start node."""

    def _calculate_k(self) -> np.ndarray[np.float64]:
        """Calculates the stiffness matrix of a planar truss element."""
        E = self.material.get_e_young()
        A = self.properties.get_area()
        I = self.properties.get_moment_of_inertia()
        L = self.length
        k = np.array([[E*A/L,         0,          0,        -E*A/L,         0,          0],
                         [0,        3*E*I/L/L/L,           0,      0, -3*E*I/L/L/L,   3*E*I/L/L],
                         [0,                  0,           0,      0,            0,           0],
                         [-E*A/L,             0,           0,  E*A/L,            0,           0],
                         [0,       -3*E*I/L/L/L,           0,      0,   3*E*I/L/L/L, -3*E*I/L/L],
                         [0,          3*E*I/L/L,           0,      0,    -3*E*I/L/L,    3*E*I/L]], dtype=np.float64)
        return k


class PlanarBeamElementT2(PlanarBeamElement):
    """Beam --------o hinge at end node."""

    def _calculate_k(self) -> np.ndarray[np.float64]:
        """Calculates the stiffness matrix of a planar truss element."""
        E = self.material.get_e_young()
        A = self.properties.get_area()
        I = self.properties.get_moment_of_inertia()
        L = self.length
        k = np.array([[E*A/L,         0,          0,        -E*A/L,         0,          0],
                         [0,        3*E*I/L/L/L,   3*E*I/L/L,      0,  -3*E*I/L/L/L,          0],
                         [0,          3*E*I/L/L,     3*E*I/L,      0,     -3*E*I/L/L,          0],
                         [-E*A/L,             0,           0,  E*A/L,            0,           0],
                         [0,       -3*E*I/L/L/L,  -3*E*I/L/L,      0,  3*E*I/L/L/L,          0],
                         [0,                  0,           0,       0,             0,          0]], dtype=np.float64)
        return k


class BaseInterLoad(ABC):
    """Base class of intermediate loads."""
    @abstractmethod
    def get_consolidation_actions(self, element):
        pass


class PlanarAxialTemperatureDifference(BaseInterLoad):
    """This class resembles an intermediate load due to a temperature difference."""
    def __init__(self, temp_diff: np.float64):
        self.temp_diff = temp_diff

    def get_consolidation_actions(self, element) -> np.ndarray[np.float64]:
        E = element.get_material().get_e_young()
        a = element.get_material().get_a_thermal()
        A = element.get_properties().get_area()
        DT = self.temp_diff
        if element.n_dims == 2:
            return element.e.transpose() @ np.array([[a * E * A * DT],
                                                     [0],
                                                     [- a * E * A * DT],
                                                     [0]])
        elif element.n_dims == 3:
            return element.e.transpose() @ np.array([[a * E * A * DT],
                                                     [0],
                                                     [0],
                                                     [- a * E * A * DT],
                                                     [0],
                                                     [0]])


class PlanarDefectiveMember(BaseInterLoad):
    """This class resembles an intermediate load due a defective member."""
    def __init__(self, delta: np.float64):
        self.delta = delta

    def get_consolidation_actions(self, element) -> np.ndarray[np.float64]:
        E = element.get_material().get_e_young()
        L = element.get_length()
        A = element.get_properties().get_area()
        d = self.delta
        if element.n_dims == 2:
            return element.e.transpose() @ np.array([[- d * E * A / L],
                                                     [0],
                                                     [d * E * A / L],
                                                     [0]])
        elif element.n_dims == 3:
            return element.e.transpose() @ np.array([[- d * E * A / L],
                                                     [0],
                                                     [0],
                                                     [d * E * A / L],
                                                     [0],
                                                     [0]])


class PlanarUniformDistributedLoad(BaseInterLoad):

    def __init__(self, load_value : np.float64, angle: np.float64 = 0, is_acting_on_solid_body: bool = True):
        """0 angle is a load acting towards the negative of x2 global axis and from there every other angle is
         a clock-wise rotation"""
        self.is_acting_on_solid_body = is_acting_on_solid_body
        self.load_value = load_value
        self.angle = angle

    def get_consolidation_actions(self, element):
        load_angle = element.get_angle() - self.angle
        q1 = self.load_value * np.sin(load_angle * np.pi/180) * np.cos(load_angle * np.pi/180)
        q2 = self.load_value * np.cos(load_angle * np.pi/180)**2
        l = element.get_length()
        if type(element) is PlanarBeamElement:
            deformable_body_actions = element.e.transpose() @ np.array([[q1*l/2],
                                                                        [q2*l/2],
                                                                        [q2*l*l/12],
                                                                        [q1*l/2],
                                                                        [q2*l/2],
                                                                        [-q2*l*l/12]])
        elif type(element) is PlanarBeamElementT1:
            deformable_body_actions = element.e.transpose() @ np.array([[q1 * l / 2],
                                                                        [3 * q2 * l / 8],
                                                                        [0],
                                                                        [q1 * l / 2],
                                                                        [5 * q2 * l / 8],
                                                                        [-q2 * l * l / 8]])
        elif type(element) is PlanarBeamElementT2:
            deformable_body_actions = element.e.transpose() @ np.array([[q1 * l / 2],
                                                                        [5 * q2 * l / 8],
                                                                        [q2 * l * l / 8],
                                                                        [q1 * l / 2],
                                                                        [3 * q2 * l / 8],
                                                                        [0]])
        return deformable_body_actions

    def get_solid_node_actions(self, element):
        start_solid_node = element.get_nodes()[0].get_coordinates()
        start_deformable_node = element.deformable_element_start_coords
        end_solid_node = element.get_nodes()[1].get_coordinates()
        end_deformable_node = element.deformable_element_end_coords
        start_vector = start_deformable_node - start_solid_node
        end_vector = end_deformable_node - end_solid_node
        start_length = np.hypot(start_vector[0], start_vector[1])
        start_angle = np.arctan2(start_vector[1], start_vector[0]) * 180 / np.pi - self.angle
        end_length = np.hypot(end_vector[0], end_vector[1])
        end_angle = np.arctan2(end_vector[1], end_vector[0]) * 180 / np.pi - self.angle
        q1_start = self.load_value * np.sin(start_angle * np.pi/180) * np.cos(start_angle * np.pi/180)
        q2_start = self.load_value * np.cos(start_angle * np.pi/180) ** 2
        q1_end = self.load_value * np.sin(end_angle * np.pi/180) * np.cos(end_angle * np.pi/180)
        q2_end = self.load_value * np.cos(end_angle * np.pi/180) ** 2
        return np.array([[q1_start*start_length],
                         [q2_start*start_length],
                         [q2_start*start_length*start_length/2],
                         [q1_end*end_length],
                         [q2_end*end_length],
                         [q2_end*end_length/2]])


class PlanarLinearTemperatureDifference(BaseInterLoad):

    def __init__(self, temperature_diff: np.float64):
        self.temperature_diff = temperature_diff

    def get_consolidation_actions(self, element):
        DT = self.temperature_diff
        h = element.get_properties().get_height()
        I = element.get_properties().get_moment_of_inertia()
        E = element.get_material().get_e_young()
        a = element.get_material().get_a_thermal()
        if type(element) is PlanarBeamElement:
            return element.e.transpose() @ np.array([[0],
                                                    [0],
                                                    [E*I*a*DT/h],
                                                    [0],
                                                    [0],
                                                    [-E*I*a*DT/h]])
        elif type(element) is PlanarBeamElementT1:
            return element.e.transpose() @ np.array([[0],
                                                     [0],
                                                     [0],
                                                     [0],
                                                     [0],
                                                     [-1.5 * E * I * a * DT / h]])
        elif type(element) is PlanarBeamElementT2:
            return element.e.transpose() @ np.array([[0],
                                                     [0],
                                                     [1.5 * E * I * a * DT / h],
                                                     [0],
                                                     [0],
                                                     [0]])


