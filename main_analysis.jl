
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import DataFrames
import CSV

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

diamond_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\216_vertices_perfect_diamond.gml"

disorder_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_3\216_vertices_T_0.08_heat_cool_0.1_per_mc_quenched.gml"

disorder_path_2 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_1\216_vertices_T_0.16_heat_cool_0.1_per_mc_quenched.gml"

disorder_path_3 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_1\216_vertices_T_0.24_heat_cool_0.1_per_mc_quenched.gml"


diamond = NG.load_spatial_network_from_gml(diamond_path)
disorder = NG.load_spatial_network_from_gml(disorder_path)
disorder_2 = NG.load_spatial_network_from_gml(disorder_path_2)
disorder_3 = NG.load_spatial_network_from_gml(disorder_path_3)



anisotropy_entropy_diamond = NA.get_anisotropy_entropy_from_bonds(diamond)
anisotropy_entropy_disorder = NA.get_anisotropy_entropy_from_bonds(disorder)
anisotropy_entropy_disorder_2 = NA.get_anisotropy_entropy_from_bonds(disorder_2)
anisotropy_entropy_disorder_3 = NA.get_anisotropy_entropy_from_bonds(disorder_3)

println("anisotropy entropy diamond: ", anisotropy_entropy_diamond)
println("anisotropy entropy disorder: ", anisotropy_entropy_disorder)
println("anisotropy entropy disorder_2: ", anisotropy_entropy_disorder_2)
println("anisotropy entropy disorder_3: ", anisotropy_entropy_disorder_3)

#NG.plot_spatial_network(disorder)