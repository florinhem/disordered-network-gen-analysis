
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU


import Plots
import Graphs
import MetaGraphsNext
import LinearAlgebra
import Polynomials
import Format
#import Statistics
#import Measurements

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

analysis_data_path = raw"..\analysis_data\neural_network_networks\dia\run_1\\"

order_metric_dict = NA.get_order_metrics_all_files(
    analysis_data_path;
    l_max_steinhardt_q_l = 12,
    save_result = true,
    save_algorithm_parameters_from_filename = true)


analysis_data_path = raw"..\analysis_data\neural_network_networks\ctn\run_1\\"

order_metric_dict = NA.get_order_metrics_all_files(
    analysis_data_path;
    l_max_steinhardt_q_l = 12,
    save_result = true,
    save_algorithm_parameters_from_filename = true)
