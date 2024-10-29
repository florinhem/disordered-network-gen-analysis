
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

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

# calculate diamond correlation function and small scale order metrics

function my_func()
    
    current_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"
    current_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

    # get all files in the current path
    files = readdir(current_path)
    # get vector of all gml files
    gml_files = [file for file in files if endswith(file, ".gml")]
    
    # loop through each file that is a gml file
    for gml_file in gml_files
        filename  = gml_file[1:end-4]
        spatial_network = NG.load_spatial_network_from_gml(current_path*filename*".gml")

        # get correlation functions
        correlation_functions_dict = NA.get_correlation_functions(
            spatial_network;
            distance_histogram_bin_width = 0.02,
            save_result = true,
            save_path = current_analysis_data_path*filename,
            label = nothing)
        # get all order metrics that contain information about small length scales
        small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(
            filename,
            current_path,
            current_analysis_data_path;
            l_max_steinhardt_q_l = 12,
            structure_factor_diamond_std_value_ratio 
                = 1,
            spectral_density_diamond_std_value_ratio 
                = 1,
            save_result = true,
            )
    end
end

my_func()