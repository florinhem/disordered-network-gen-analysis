
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
#import Statistics
#import Measurements

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

spatial_network_path =  raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\testing\\"
filename = "216_vertices_T_0.11_heat_cool_0.1_per_mc_quenched"
analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\testing\\"

NA.get_all_dicts_from_network_single_file(
    filename,
    spatial_network_path,
    analysis_data_path;
    pore_size_sampling_grid_size = 0.2,
    print_progress = true,
    print_lock = Threads.ReentrantLock())

# load order metrics dict
order_metrics_dict = GU.load_h5_dict(analysis_data_path*filename*"_order_metrics.h5")

for key in keys(order_metrics_dict)
    println(key, ": ", order_metrics_dict[key])
end