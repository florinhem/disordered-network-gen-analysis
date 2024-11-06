
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

# possible choices of nr_vertices for diamond: 64, 216, 512, 216, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 216 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_1\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\geometrical_models\216_vertices_bond_bending_0.21\run_1\\"

filename = "216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched"

spatial_network = NG.load_spatial_network_from_gml(path * filename * ".gml")

#NG.plot_spatial_network(spatial_network)
NG.plot_spatial_network_2(spatial_network)
