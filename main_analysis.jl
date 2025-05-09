
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
#import Statistics
#import Measurements

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

print_lock = Threads.ReentrantLock()

spatial_networks_upper_path = "../structures/random_networks/"
analysis_data_upper_path = "../analysis_data/random_networks/new_order_metrics/"

# get a list of folders in the spatial_networks_path directory that begin with
# "216_vertices"
folder_names = readdir(spatial_networks_upper_path)
folder_names = filter(x -> occursin("216_vertices", x), folder_names)

for folder_name in folder_names
    spatial_networks_path = spatial_networks_upper_path * folder_name * "/"
    analysis_data_path = analysis_data_upper_path * folder_name * "/"

    NA.get_all_dicts_from_networks_multithreading(
    spatial_networks_path,
    analysis_data_path;
    print_progress = true,
    print_lock = print_lock)
end

