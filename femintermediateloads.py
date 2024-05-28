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
            return np.array([[a * E * A * DT],
                            [0],
                            [- a * E * A * DT],
                            [0]])
        elif element.n_dims == 3:
            return np.array([[a * E * A * DT],
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
            return np.array([[- d * E * A / L],
                            [0],
                            [d * E * A / L],
                            [0]])
        elif element.n_dims == 3:
            return np.array([[- d * E * A / L],
                             [0],
                             [0],
                             [d * E * A / L],
                             [0],
                             [0]])


class PlanarUniformDistributedLoad(BaseInterLoad):

    def __init__(self, load_value : np.float64):
        self.load_value = load_value

    def get_consolidation_actions(self, element):
        q = self.load_value
        l = element.get_length()
        return np.array([[0],
                         [q*l/2],
                         [q*l*l/12],
                         [0],
                         [q*l/2],
                         [q*l*l/12]])


class PlanarLinearTemperatureDifference(BaseInterLoad):

    def __init__(self, temperature_diff: np.float64):
        self.temperature_diff = temperature_diff

    def get_consolidation_actions(self, element):
        DT = self.temperature_diff
        h = element.get_properties().get_height()
        I = element.get_properties().get_moment_of_inertia()
        E = element.get_material().get_e_young()
        a = element.get_material().get_a_thermal()
        return np.array([[0],
                         [0],
                         [E*I*a*DT/h],
                         [0],
                         [0],
                         [-E*I*a*DT/h]])
