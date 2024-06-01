"""This module contains all kinds of intermediate loads as classes that inherit from the BaseInterLoad class
and implement the get_consolidation_actions method distinct for every type of intermediate load."""

import numpy as np
from abc import ABC, abstractmethod


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
        q1 = self.load_value * np.sin(load_angle) * np.cos(load_angle)
        q2 = self.load_value * np.cos(load_angle)**2
        l = element.get_length()
        deformable_body_actions = element.e.transpose() @ np.array([[q1*l/2],
                                                                    [q2*l/2],
                                                                    [q2*l*l/12],
                                                                    [q1*l/2],
                                                                    [q2*l/2],
                                                                    [q2*l*l/12]])
        if not self.is_acting_on_solid_body:
            return deformable_body_actions
        return deformable_body_actions + self._get_solid_node_actions(element)

    def _get_solid_node_actions(self, element):
        start_solid_node = element.get_nodes()[0].get_coordinates
        start_deformable_node = element.deformalbe_element_start_cords
        end_solid_node = element.get_nodes()[1].get_coordinates
        end_deformable_node = element.deformalbe_element_end_cords
        start_vector = start_deformable_node - start_solid_node
        end_vector = end_deformable_node - end_solid_node
        start_length = np.hypot(start_vector[0], start_vector[1])
        start_angle = np.arctan2(start_vector[1], start_vector[0]) - self.angle
        end_length = np.hypot(end_vector[0], end_vector[1])
        end_angle = np.arctan2(end_vector[1], end_vector[0]) - self.angle
        q1_start = self.load_value * np.sin(start_angle) * np.cos(start_angle)
        q2_start = self.load_value * np.cos(start_angle) ** 2
        q1_end = self.load_value * np.sin(start_angle) * np.cos(start_angle)
        q2_end = self.load_value * np.cos(start_angle) ** 2
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
        return element.e.transpose() @ np.array([[0],
                                                [0],
                                                [E*I*a*DT/h],
                                                [0],
                                                [0],
                                                [-E*I*a*DT/h]])
