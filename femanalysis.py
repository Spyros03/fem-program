"""This module provides code that act as an intermediate layer between the user interface and the finite element method
solver."""

import numpy as np
from abc import ABC, abstractmethod

from femnodes import *
from femelement import *
from femessentials import *
from femstructure import *
from femintermediateloads import *


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
    def add_properties(self, name: str, properties: BaseProperties):
        pass

    # @abstractmethod
    # def add_support(self, restrictions: np.ndarray[np.float64], angles: np.ndarray[np.float64],
    #                 springs: np.ndarray[np.float64], retreats: np.ndarray[np.float64]):
    #     pass

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

    def add_node(self, node_id: int, x:float, y:float):
        self.nodes.insert(node_id - 1, Node(node_id - 1, self.n_dims,
                                            np.array([np.float64(x), np.float64(y)], dtype=np.float64)))
        self.nodes_count += 1

    def add_support(self, node_id: int, support: BaseSupport):
        self.nodes[node_id - 1].set_support(support)

    def add_nodal_load(self, node_id: int, px: float, py: float):
        self.nodes[node_id - 1].set_external_load(np.array([px, py], dtype=np.float64).transpose())

    def add_element(self, element_id: int, start_node_id: int, end_node_id: int, material_name: str,
                    properties_name: str, temperature_diff: float = None):
        self.elements.insert(element_id, PlanarTrussElement(element_id, self.n_dims, [self.nodes[start_node_id-1],
                                                            self.nodes[end_node_id-1]], self.materials[material_name],
                                                            self.properties[properties_name]))
        self.elements_count += 1

    def add_intermediate_temperature_diff(self, element_id: int, temperature_diff):
        self.elements[element_id - 1].add_intermediate_load(PlanarTrussAxialTemperatureDifference(temperature_diff))

    def add_intermediate_defective_member(self, element_id: int, delta):
        self.elements[element_id - 1].add_intermediate_load(PlanarTrussDefectiveMember(delta))

    def add_properties(self, name: str, area: float):
        self.properties.update({name: PlanarTrussProperties(name, np.float64(area))})

    # def add_support(self, restrictions: np.ndarray[np.float64], angles: np.ndarray[np.float64],
    #                 springs: np.ndarray[np.float64], retreats: np.ndarray[np.float64]):
    #     pass

    def analyze(self):
        structure = PlanarStructure(self.nodes, self.elements, self.n_dims)
        structure.analyze()

    def print_results(self):
        for element in self.elements:
            element.print_internal_forces()
