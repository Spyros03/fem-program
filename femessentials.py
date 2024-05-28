"""This module provides classes that contain material properties and geometric properties."""

import numpy as np
from abc import ABC


class Material:
    """This is a class that contains material properties."""
    def __init__(self, name: str, e_young: np.float64, a_thermal: np.float64 = 1e-5):
        self.name = name
        self.e_young = e_young
        self.a_thermal = a_thermal

    def get_e_young(self) -> np.float64:
        return self.e_young

    def get_a_thermal(self) -> np.float64:
        return self.a_thermal


class BaseProperties(ABC):
    """This is a base class that contains geometric properties."""
    def __init__(self, name):
        self.name = name


class PlanarTrussProperties(BaseProperties):
    """Contains the geometric properties(area) of a planar truss structure."""
    def __init__(self, name: str, area: np.float64):
        super().__init__(name)
        self.area = area

    def get_area(self) -> np.float64:
        return self.area


class PlanarBeamProperties(BaseProperties):

    def __init__(self, name: str, area: np.float64, moment_of_inertia: np.float64, height: np.float64 = None):
        super().__init__(name)
        self.area = area
        self.moment_of_inertia = moment_of_inertia
        self.height = height

    def get_area(self) -> np.float64:
        return self.area

    def get_moment_of_inertia(self) -> np.float64:
        return self.moment_of_inertia

    def get_height(self) -> np.float64:
        return self.height
