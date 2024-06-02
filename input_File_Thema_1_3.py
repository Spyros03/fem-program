"""Input file erothma 3 gia to thema 1"""

from femanalysis import PlanarFemAnalysis

n_dims = 2

analysis = PlanarFemAnalysis(n_dims)

i = 1
j = 0
k = 0
e_young = 2.1e8
area = 1e-3
Pload = 50
Li = 3 + (i-1) * 0.1
Lj = 4 + (j-1) * 0.1
Lk = 3 + (k-1) * 0.1
a = 1e-5
T = 30
T0 = 10
DT = T - T0


analysis.add_material("1", e_young, a)
analysis.add_properties("1", "truss",  area)

analysis.add_node(1, Li, 0)
analysis.add_node(2, Li, Lk)
analysis.add_node(3, 0, Lk)
analysis.add_node(4, Li, 2*Lk)
analysis.add_node(5, Li + Lj, 2*Lk)
analysis.add_node(6, 2*Li + Lj, Lk)
analysis.add_node(7, Li + Lj, Lk)
analysis.add_node(8, Li + Lj, 0)

analysis.add_element(1, 1, 2, "1", "1")
analysis.add_element(2, 2, 3, "1", "1")
analysis.add_element(3, 3, 4, "1", "1")
analysis.add_element(4, 2, 4, "1", "1")
analysis.add_element(5, 4, 7, "1", "1")
analysis.add_element(6, 4, 5, "1", "1")
analysis.add_element(7, 5, 6, "1", "1")
analysis.add_element(8, 5, 7, "1", "1")
analysis.add_element(9, 2, 5, "1", "1")
analysis.add_element(10, 6, 7, "1", "1")
analysis.add_element(11, 2, 7, "1", "1")
analysis.add_element(12, 2, 8, "1", "1")
analysis.add_element(13, 1, 7, "1", "1")
analysis.add_element(14, 7, 8, "1", "1")
analysis.add_element(15, 1, 8, "1", "1")

analysis.add_support(1, [True, False], springs=[0, 5000])
analysis.add_support(8, [False, True], angles=[30])

analysis.add_nodal_load(4, 0, -2*Pload)
analysis.add_nodal_load(5, 0, -2*Pload)
analysis.add_nodal_load(6, -3*Pload, 0)

analysis.add_axial_temperature_diff(1, DT)
analysis.add_axial_temperature_diff(2, DT)
analysis.add_axial_temperature_diff(3, DT)
analysis.add_axial_temperature_diff(4, DT)
analysis.add_axial_temperature_diff(5, DT)
analysis.add_axial_temperature_diff(6, DT)
analysis.add_axial_temperature_diff(7, DT)
analysis.add_axial_temperature_diff(8, DT)
analysis.add_axial_temperature_diff(9, DT)
analysis.add_axial_temperature_diff(10, DT)
analysis.add_axial_temperature_diff(11, DT)
analysis.add_axial_temperature_diff(12, DT)
analysis.add_axial_temperature_diff(13, DT)
analysis.add_axial_temperature_diff(14, DT)
analysis.add_axial_temperature_diff(15, DT)

analysis.add_intermediate_defective_member(5, delta=0.02)
