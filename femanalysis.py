"""This module provides code that act as an intermediate layer between the user interface and the finite element method
solver."""

import numpy as np
from abc import ABC, abstractmethod

from femnodes import *
from femelement import *
from femessentials import *
from femstructure import *
from femintermediateloads import *
from femdrawstructure import draw_structure


class BaseFemAnalysis(ABC):

    def __init__(self, n_dims: int):
        self.n_dims = n_dims
        self.nodes_count = 0
        self.elements_count = 0
        self.materials_count = 0
        self.properties_count = 0
        self.nodes = list()
        self.elements = list()
        self.materials = dict()
        self.properties = dict()

    @abstractmethod
    def add_node(self):
        pass

    @abstractmethod
    def add_element(self, element_id, start_node_id: int, end_node_id: int, material_name: str, properties_name: str):
        pass

    def add_material(self, material_name: str, e_young: float, a_thermal: float):
        self.materials.update({material_name: Material(material_name, np.float64(e_young), np.float64(a_thermal))})
        self.materials_count += 1

    @abstractmethod
    def add_support(self, node_id: int, restrictions: list[bool], angles: np.ndarray[np.float64],
                    springs: np.ndarray[np.float64], retreats: np.ndarray[np.float64]):
        pass

    def remove_node(self, node_id: int):
        self.nodes.remove(node_id - 1)
        self.nodes_count -= 1

    def remove_element(self, element_id: int):
        self.elements.remove(element_id - 1)
        self.elements_count -= 1

    def remove_materials(self, material_name: str):
        self.materials.pop(material_name)
        self.materials_count -= 1

    def remove_properties(self, material_name: str):
        self.materials.pop(material_name)
        self.properties_count -= 1

    def remove_supports(self, node_id: int):
        self.nodes[node_id - 1].set_support(None)

    @abstractmethod
    def analyze(self):
        pass

    def get_nodes_count(self) -> int:
        return self.nodes_count

    def get_elements_count(self) -> int:
        return self.elements_count

    def get_materials_count(self) -> int:
        return self.materials_count

    def get_properties_count(self) -> int:
        return self.properties_count


class PlanarFemAnalysis(BaseFemAnalysis):

    def __init__(self, n_dims: int):
        super().__init__(n_dims)
        self.structure = None

    def add_node(self, node_id: int, x : float, y : float):
        self.nodes.insert(node_id - 1, Node(node_id - 1, self.n_dims,
                                            np.array([np.float64(x), np.float64(y)], dtype=np.float64)))
        self.nodes_count += 1

    def add_nodal_load(self, node_id: int, px: float, py: float, m=0):
        if self.n_dims == 2:
            self.nodes[node_id - 1].set_external_load(np.array([px, py], dtype=np.float64).transpose())
        elif self.n_dims == 3:
            self.nodes[node_id - 1].set_external_load(np.array([px, py, m], dtype=np.float64).transpose())

    def add_element(self, element_id: int, start_node_id: int, end_node_id: int, material_name: str,
                    properties_name: str, temperature_diff: float = None):
        if isinstance(self.properties[properties_name], PlanarTrussProperties):
            self.elements.insert(element_id, PlanarTrussElement(element_id, self.n_dims, [self.nodes[start_node_id-1],
                                                                self.nodes[end_node_id-1]], self.materials[material_name],
                                                                self.properties[properties_name]))
        elif isinstance(self.properties[properties_name], PlanarBeamProperties):
            self.elements.insert(element_id, PlanarBeamElement(element_id, self.n_dims,
                                                               [self.nodes[start_node_id - 1],
                                                                self.nodes[end_node_id - 1]],
                                                                self.materials[material_name],
                                                                self.properties[properties_name]))

        self.elements_count += 1

    def add_support(self, node_id, restrictions: list[bool], angles: list[float] = None,
                    springs: list[float] = None, retreats: list[float] = None):
        if angles is not None:
            angles = np.array(angles, dtype=np.float64)
        if springs is not None:
            np.array(springs, dtype=np.float64)
        if retreats is not None:
            np.array(retreats, dtype=np.float64)
        self.nodes[node_id - 1].set_support(Support(restrictions, angles, springs, retreats))

    def add_axial_temperature_diff(self, element_id: int, temperature_diff: float):
        self.elements[element_id - 1].add_intermediate_load(
            PlanarAxialTemperatureDifference(np.float64(temperature_diff)))

    def add_linear_temperature_diff(self, element_id: int, temperature_diff: float):
        self.elements[element_id - 1].add_inretmediate_load(
            PlanarLinearTemperatureDifference(np.float64(temperature_diff)))

    def add_intermediate_defective_member(self, element_id: int, delta):
        self.elements[element_id - 1].add_intermediate_load(PlanarDefectiveMember(delta))

    def add_uniform_distributed_load(self, element_id: int, load_value):
        self.elements[element_id - 1].add_intermediate_load(PlanarUniformDistributedLoad(np.float64(load_value)))

    def add_properties(self, name: str, type_of_element: str, area: float, inertia: float = 0):
        if type_of_element == "truss":
            self.properties.update({name: PlanarTrussProperties(name, np.float64(area))})
        elif type_of_element == "beam":
            self.properties.update({name: PlanarBeamProperties(name, np.float64(area), np.float64(inertia))})

    def analyze(self, draw: bool = False):
        self.structure = PlanarStructure(self.nodes, self.elements, self.n_dims)
        self.structure.analyze()
        if draw:
            draw_structure(False, self.elements, self.nodes)
            draw_structure(True, self.elements, self.nodes)

    def print_results(self):
        for element in self.elements:
            element.print_internal_forces()
