
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

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
import Statistics
#import Statistics
#import Measurements

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

neural_network_dataset_path = raw"..\neural_networks\datasets\\"

analysis_data_path = raw"..\analysis_data\neural_network_networks\ctn\\"

order_metrics_dict = GU.load_h5_dict(analysis_data_path*"all_order_metrics.h5")

GU.save_dict_to_h5(order_metrics_dict,
    neural_network_dataset_path*"ctn_all_order_metrics.h5")


analysis_data_path = raw"..\analysis_data\neural_network_networks\dia\\"

order_metrics_dict = GU.load_h5_dict(analysis_data_path*"all_order_metrics.h5")

GU.save_dict_to_h5(order_metrics_dict,
    neural_network_dataset_path*"dia_all_order_metrics.h5")


analysis_data_path = raw"..\analysis_data\neural_network_networks\srs\\"

order_metrics_dict = GU.load_h5_dict(analysis_data_path*"all_order_metrics.h5")

GU.save_dict_to_h5(order_metrics_dict,
    neural_network_dataset_path*"srs_all_order_metrics.h5")
    

analysis_data_path = raw"..\analysis_data\neural_network_networks\lcs\\"

order_metrics_dict = GU.load_h5_dict(analysis_data_path*"all_order_metrics.h5")

GU.save_dict_to_h5(order_metrics_dict,
    neural_network_dataset_path*"lcs_all_order_metrics.h5")
    