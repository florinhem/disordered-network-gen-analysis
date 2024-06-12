
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .BinaryDataAnalysis as BDA
import .GeneralUtilities as GU

import Plots
import LinearAlgebra
import Measurements
import MetaGraphsNext
import StatsBase
import Peaks

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

# julia --threads 23


function get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)

    for i in 1:5

        current_graph_dict_path = graph_dict_path*"run_"*string(i)*"\\"

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames

            small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(structure_dict_filename[1:end-13],
            current_graph_dict_path,
            current_analysis_data_path,
            save_result = true)
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\\"

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\\"

get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)
