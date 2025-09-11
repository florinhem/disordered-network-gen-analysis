
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_targeted\dia\1000_vertices\\"

order_metrics_dict = NA.get_order_metrics_all_files(
    analysis_data_path;
    l_max_steinhardt_q_l = 12,
    save_result = true,
    save_algorithm_parameters_from_filename = true)