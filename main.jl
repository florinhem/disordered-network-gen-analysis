
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007

dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.1]

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*
    "1000_vertices_T_0.1_heated_for_0.5_steps_quenched")

l_max = 12

q_l_total_network_mean_dict = NA.get_q_l_total_network_mean_dict(graph_dict, l_max)

q_l_total_network_mean_vec = NA.convert_q_l_dict_to_vec(q_l_total_network_mean_dict, l_max)
