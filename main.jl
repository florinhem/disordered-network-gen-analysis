
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import ProfileView

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007

# julia --threads 20

# load some network
dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"
filename = "216_vertices_T_0.2_heated_for_0.5_steps_quenched"

save_path = raw"..\analysis_data\random_networks\\"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(graph_dict; save_result = true, save_path= save_path*filename)
