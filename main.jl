
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

evolution_dict_low_t = GU.load_h5_dict(dict_path*"1000_vertices_T_0.1_heated_for_0.5_steps_quenched_evolution.h5")

evolution_dict_high_t = GU.load_h5_dict(dict_path*"1000_vertices_T_1.0_heated_for_0.5_steps_quenched_evolution.h5")

println(length(evolution_dict_low_t["move_accepted_vec"])/18000 )

println(length(evolution_dict_high_t["move_accepted_vec"])/18000 )