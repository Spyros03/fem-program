"""This module contains all kinds of intermediate loads as classes that inherit from the BaseInterLoad class
and implement the get_consolidation_actions method distinct for every type of intermediate load."""

import numpy as np
from abc import ABC, abstractmethod


class BaseInterLoad(ABC):
    """Base class of intermediate loads."""
    @abstractmethod
    def get_consolidation_actions(self, element):
        pass


class BasePlanarTrussInterLoad(BaseInterLoad):
    """Base class for intermediate loads for planar truss members."""
    @abstractmethod
    def get_consolidation_actions(self, element) -> np.ndarray[np.float64]:
        pass


class PlanarTrussAxialTemperatureDifference(BasePlanarTrussInterLoad):
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


class PlanarTrussDefectiveMember(BasePlanarTrussInterLoad):
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
