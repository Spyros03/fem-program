from femanalysis import PlanarFemAnalysis

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

analysis.add_material("concrete", 200000000, 1e-5)
analysis.add_properties("concrete", 'beam', 0.004, 0.007)

analysis.add_node(1, 0, 0)
analysis.add_node(2, 0, 2)
analysis.add_node(3, 2, 2)
analysis.add_node(4, 2, 0)

analysis.add_element(1, 1,2, "concrete", "concrete")
analysis.add_element(2, 2, 3, "concrete", "concrete")
analysis.add_element(3, 3, 4, "concrete", "concrete")

analysis.add_support(1, [True, True, True])
analysis.add_support(4, [True, True, True])

analysis.add_nodal_load(3, -200, 0, 0)

