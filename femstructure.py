"""This module contains classes that resemble structure models."""

from abc import ABC, abstractmethod
import numpy as np

from femnodes import *
from femelement import *


class BaseStructure(ABC):
    """A base class for a structure."""
    def __init__(self, nodes: list[Node], elements: list[BaseElement], n_dims: int):
        """Constructor that declares and initializes all structures variables."""
        self.nodes = nodes
        self.elements = elements
        self.n_dims = n_dims
        self.n_nodes = len(nodes)
        self.n_dofs = self.n_dims * self.n_nodes
        self.kg = self._calculate_kg()
        self.bounded_dofs = self._get_bounded_dofs()
        self.n_bounded_dofs = len(self.bounded_dofs)
        self.n_free = self.n_dofs - self.n_bounded_dofs
        self.V = self._calculate_permutation_vector()
        self.rotations = self._calculate_rotation_matrices()
        self.p_ext = self._get_external_loading()
        self.s = self._calculate_s() #consolidation actions vector.
        self.p = self.p_ext - self.s #nodal loads for the equivalent structure.
        self.p_rot = self.p
        self.u = np.zeros(shape=[self.n_dofs, 1], dtype=np.float64)
        self.u_rot = self.u
        self.km = self._calculate_km()  #modification due to rotated support.
        self.kmm = self._calculate_kmm() #modification due to elastic support.
        self.kmmm = self._calculate_kmmm() #modification to get the free nodes first.

    def _calculate_kg(self) -> np.ndarray[np.float64]:
        """Calculates the global stiffness matrix."""
        kg = np.zeros(shape=[self.n_dofs, self.n_dofs], dtype=np.float64)
        for element in self.elements:
            kg += element.get_element_contribution_to_kg(self.n_nodes)
        return kg

    def _get_bounded_dofs(self) -> list[int]:
        """Gets the degrees of freedom that are bounded."""
        bounded_dofs = list()
        for node in self.nodes:
            bounded_dofs += node.get_bounded_dofs()
        return bounded_dofs

    def _calculate_permutation_vector(self) -> np.ndarray[np.float64]:
        """Calculates the permutation matrix."""
        V0 = np.eye(self.n_dofs)
        all_dofs = np.arange(self.n_dofs).transpose()
        free_dofs = np.setdiff1d(all_dofs, self.bounded_dofs)
        permutation_vector = np.hstack((free_dofs, self.bounded_dofs))
        return V0[permutation_vector, :]

    @abstractmethod
    def _calculate_rotation_matrices(self) -> list[np.ndarray[np.float64]]:
        """Abstract method that returns a list with all the modifications that are needed due a rotated support."""
        pass

    def _get_external_loading(self) -> np.ndarray[np.float64]:
        """Returns the vector of the external loads."""
        p_ext = np.zeros(shape=[self.n_dofs, 1], dtype=np.float64)
        for inode, node in enumerate(self.nodes):
            if node.get_external_load() is None:
                continue
            for idof in range(self.n_dims):
                p_ext[self.n_dims * inode + idof, 0] = node.get_external_load()[idof]
        return p_ext

    def _calculate_s(self) -> np.ndarray[np.float64]:
        """Calculates the consolidation actions."""
        s = np.zeros(shape=[self.n_dofs,1], dtype=np.float64)
        for element in self.elements:
            s += element.get_element_contribution_to_s(self.n_nodes)
        return s

    def _calculate_km(self) -> np.ndarray[np.float64]:
        """Calculates the modified stiffness matrix due to rotated supports."""
        if not self.rotations:
            return self.kg
        km = self.kg
        for rotation in self.rotations:
            km = rotation @ km @ rotation.transpose()
        return km

    def _calculate_kmm(self) -> np.ndarray[np.float64]:
        """Calculates the modified stiffness matrix due to elastic supports."""
        kmm = self.km
        for node in self.nodes:
            if node.get_springs() is None:
                continue
            kmm += node.get_spring_contribution(self.n_nodes)
        return kmm

    def _calculate_kmmm(self) -> np.ndarray[np.float64]:
        """Calculates the permutated stiffness matrix."""
        return self.V @ self.kmm @ self.V.transpose()

    def analyze(self):
        """Solves the structure."""

        #Gets the partial stiffness matrices.
        kff = self.kmmm[0:self.n_free, 0:self.n_free]
        kfs = self.kmmm[0:self.n_free, self.n_free:self.n_dofs]
        ksf = kfs.transpose()
        kss = self.kmmm[self.n_free:self.n_dofs, self.n_free:self.n_dofs]

        #Rotates the external loads due to rotated supports.
        if self.rotations:
            for rotation in self.rotations:
                self.p_rot = rotation @ self.p_rot

        #Permutated displacement and loads vector.
        dm = self.V @ self.u
        pm = self.V @ self.p_rot

        # Supported displacements and free loads (Known variables)
        ds = dm[self.n_free:self.n_dofs, [0]]
        pf = pm[0:self.n_free, [0]]
        #Solves the equations.
        df = np.linalg.inv(kff) @ (pf - kfs @ ds)
        ps = ksf @ df + kss @ ds
        #Displacemnt and external load vectors rotated due to rotated supports.
        self.u_rot = self.V.transpose() @ np.block([[df], [ds]])
        self.p_rot = self.V.transpose() @ np.block([[pf], [ps]])
        #Displacemnt and external load vectors (global coordinate system).
        self.u = self.u_rot
        self.p = self.p_rot
        if self.rotations:
            for rotation in self.rotations[::-1]:
                self.p = rotation.transpose() @ self.p_rot
                self.u = rotation.transpose() @ self.u_rot

        #Exernal loads and support forces of the real structure.
        self.p_ext = self.p + self.s

        self._assign_displacements_to_nodes()
        self._calculate_element_forces()

    def _assign_displacements_to_nodes(self):
        """Assigns the displacements to the nodes."""
        for inode, node in enumerate(self.nodes):
            displacements = np.zeros(shape=[self.n_dims, 1], dtype=np.float64)
            for idof in range(self.n_dims):
                displacements[idof] = self.u[self.n_dims * inode + idof]
            node.set_displacements(displacements)

    def _calculate_element_forces(self):
        """Calculates the element's forces."""
        for element in self.elements:
            element.calculate_element_forces()

    def get_kg(self) -> np.ndarray[np.float64]:
        return self.kg

    def get_external_loading(self) -> np.ndarray[np.float64]:
        return self.p_ext

    def get_s(self) -> np.ndarray[np.float64]:
        return self.s

    def get_km(self) -> np.ndarray[np.float64]:
        return self.km

    def get_kmm(self) -> np.ndarray[np.float64]:
        return self.kmm

    def get_kmmm(self) -> np.ndarray[np.float64]:
        return self.kmmm


class PlanarStructure(BaseStructure):
    """A class for a planar truss structure."""
    def __init__(self, nodes: list[Node], elements: list[PlanarTrussElement], n_dims: int):
        super().__init__(nodes, elements, n_dims)

    def _calculate_rotation_matrices(self) -> list[np.ndarray[np.float64]]:
        """Calculated the rotations matrices due to a rotated support."""
        rotations = list()
        for node in self.nodes:
            if node.get_support() is None:
                continue
            if node.get_bound_angles() is None:
                continue
            r = np.eye(self.n_dofs, dtype=np.float64)
            angle = node.get_support().get_angles()[0]
            cos_angle = np.cos(np.pi/180 * angle)
            sin_angle = np.sin(np.pi/180 * angle)
            r[np.ix_(node.get_index(), node.get_index())] = (
                np.array([[cos_angle, sin_angle], [-sin_angle, cos_angle]]))
            rotations.append(r)
        return rotations