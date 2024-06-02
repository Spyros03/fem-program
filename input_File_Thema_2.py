from femanalysis import PlanarFemAnalysis
from femhelperclass import *

n_dims = 3

analysis = PlanarFemAnalysis(n_dims)

i = 1
j = 0
k = 0
e_young1 = 2e7
e_young2 = 2e8
t_internal = 20
t_external = 40
t_initial = 10
Li = 10 + (i - 1) * 0.1

analysis.add_material("beam", e_young1, 1e-5)
analysis.add_material("truss", e_young2, 1e-5)
analysis.add_properties("beam", "beam", get_area(0.20, 0.75), get_moment_of_inertia(0.20, 0.75), 0.75)
analysis.add_properties("column", "beam", get_area(0.20, 0.75), get_moment_of_inertia(0.20, 0.75), 0.75)
analysis.add_properties("truss", "truss", 0.001)

analysis.add_node(1, 0, 0)
analysis.add_node(2, Li, 0)
analysis.add_node(3, Li, - (2 + (j - 1) * 0.1))
analysis.add_node(4, Li, 5 + (j - 1) * 0.1)
analysis.add_node(5, Li/2, 0, free_rotation=True)

analysis.add_element(1, 1, 5, "beam", "beam", start_of_deformed=[1 + (k - 1) * 0.1, 0])
analysis.add_element(2, 2, 5, "beam", "beam")
analysis.add_element(3, 2, 3, "beam", "column")
analysis.add_element(4, 2, 4, "beam", "column")
analysis.add_element(5, 4, 5, "truss", "truss")

analysis.add_uniform_distributed_load(1, 35)
analysis.add_uniform_distributed_load(2, -35)
analysis.add_linear_temperature_diff(1, t_internal - t_external)
analysis.add_linear_temperature_diff(2, t_external - t_internal)
analysis.add_linear_temperature_diff(3, t_internal - t_external)
analysis.add_axial_temperature_diff(1, (t_internal + t_external)/2 - t_initial)
analysis.add_axial_temperature_diff(2, (t_internal + t_external)/2 - t_initial)
analysis.add_axial_temperature_diff(3, (t_internal + t_external)/2 - t_initial)
analysis.add_axial_temperature_diff(4, t_external - t_initial)
analysis.add_axial_temperature_diff(5, t_external - t_initial)


analysis.add_support(1, [False, True, False], angles=[-25], springs=[0, 0, 24000])
analysis.add_support(3, [True, True, True], retreats=[0, -0.002, 0])


