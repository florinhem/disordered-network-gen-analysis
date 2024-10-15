
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs


# possible choices of nr_vertices for diamond: 64, 216, 512, 216, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 216 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

# julia --threads 23

nr_vertices_vec = [216, 512, 1000]


for nr_vertices in nr_vertices_vec

    filename = string(nr_vertices)*"_vertices_T_0.08_heat_cool_0.1_per_mc_quenched"

    for run in 1:5

        spatial_network_path = "../structures/random_networks/"*string(nr_vertices)*"_vertices_bond_bending_0.21/run_"*string(run)*"/"
        structure_dict_path = "../structures/random_networks/binary_structures/"*string(nr_vertices)*"_vertices_bond_bending_0.21/run_"*string(run)*"/"
        analysis_data_path = "../analysis_data/random_networks/"*string(nr_vertices)*"_vertices_bond_bending_0.21/run_"*string(run)*"/"

        println("Nr vertices: ", nr_vertices, ", run: ", run)

        NA.get_all_dicts_from_network_single_file(
            filename,
            spatial_network_path,
            structure_dict_path,
            analysis_data_path;
            bond_radius = 0.35,
            voxel_edge_length = 0.1,
            structure_factor_diamond_std_value_ratio = 1,
            spectral_density_diamond_std_value_ratio = 1,
            pore_size_distribution_nr_sampled_voxels = 20000,
            print_progress = true,
            print_lock = Threads.ReentrantLock())

    end
end

