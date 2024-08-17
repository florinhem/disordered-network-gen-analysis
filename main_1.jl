
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import Measurements

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

# julia --threads 23

graph_dict_path = "../structures/random_networks/diamonds/"

structure_dict_path = "../structures/random_networks/binary_structures/diamonds/"

analysis_data_path = "../analysis_data/random_networks/diamonds/"

filename = "216_vertices_perfect_diamond"

NA.get_all_dicts_from_graph_single_file(filename,
    graph_dict_path,
    structure_dict_path,
    analysis_data_path;
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    print_lock = Threads.ReentrantLock())