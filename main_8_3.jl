
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

graph_dicts_path = "../structures/random_networks/512_vertices_bond_bending_0.285/"
structure_dicts_path = "../structures/random_networks/binary_structures/512_vertices_bond_bending_0.285/"
analysis_data_path = "../analysis_data/random_networks/512_vertices_bond_bending_0.285/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
        structure_dicts_path,
        analysis_data_path::String;
        structure_factor_diamond_std_value_ratio = 1,
        spectral_density_diamond_std_value_ratio = 1,
        pore_size_distribution_nr_sampled_voxels = 20000,
        print_progress = true,
        runs_vec = collect(1:5),
        print_lock = Threads.ReentrantLock())