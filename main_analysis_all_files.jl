
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\local_relaxation\targeted\shell_nr_4\\"

order_metrics_dict = NA.get_order_metrics_all_files(
    analysis_data_path;
    save_result=true,
    save_algorithm_parameters_from_filename=true)
#filename = "all_order_metrics.h5"
#order_metrics_dict = GU.load_h5_dict(load_path*filename)

NA.save_order_metrics_dict_to_csv(order_metrics_dict, analysis_data_path)