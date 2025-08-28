
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Random
import LsqFit
import Plots
import LinearAlgebra
import Measurements

# initial guess
p0_ctn = [0.28, 0.7, 0.25, 0.2, 0.34, 0.05]
p0_dia = [0.28, 0.7, 0.25, 0.2, 0.4, 0.05]
p0_lcs = [0.28, 0.7, 0.25, 0.2, 0.38, 0.05]
p0_srs = [0.28, 0.7, 0.25, 0.2, 0.3, 0.05]


network_type_vec = ["ctn", "dia", "lcs", "srs"]
p0_list = [p0_ctn, p0_dia, p0_lcs, p0_srs]

for (i, network_type) in enumerate(network_type_vec)

    println("")
    println(network_type)

    # load the order metrics dict
    analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_networks\\" * network_type * raw"\\"
    order_metric_dict = GU.load_h5_dict(analysis_data_path * "all_order_metrics.h5")

    # get the r_t dict path
    r_t_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\photonics\tidy3d\simulation_data\neural_network_networks\\" * network_type * raw"\run_1_2_r_t_low_n\\"

    correlations_dict = NA.get_order_peak_height_to_width_correlations(
        order_metric_dict,
        r_t_dict_path,
        analysis_data_path;
        p0=p0_list[i],
        save_results = true)
end