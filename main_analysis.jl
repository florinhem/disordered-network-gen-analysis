
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU


import Plots

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


pore_size_distribution_dict_diamond = NA.get_pore_size_distribution(diamond, sampling_grid_size = 0.4)
pore_size_distribution_dict_disorder = NA.get_pore_size_distribution(disorder, sampling_grid_size = 0.4)
pore_size_distribution_dict_disorder_2 = NA.get_pore_size_distribution(disorder_2, sampling_grid_size = 0.4)
pore_size_distribution_dict_disorder_3 = NA.get_pore_size_distribution(disorder_3, sampling_grid_size = 0.4)

# plot the pore size distribution for the diamond network
Plots.plot(pore_size_distribution_dict_diamond["pore_size_vec"], pore_size_distribution_dict_diamond["pore_size_distribution"], label = "diamond")
Plots.plot!(pore_size_distribution_dict_disorder["pore_size_vec"], pore_size_distribution_dict_disorder["pore_size_distribution"], label = "disorder")
Plots.plot!(pore_size_distribution_dict_disorder_2["pore_size_vec"], pore_size_distribution_dict_disorder_2["pore_size_distribution"], label = "disorder_2")
Plots.plot!(pore_size_distribution_dict_disorder_3["pore_size_vec"], pore_size_distribution_dict_disorder_3["pore_size_distribution"], label = "disorder_3")
Plots.xlabel!("Pore size/d")
Plots.ylabel!("Pore size distribution")