
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import LinearAlgebra
import Plots
import MetaGraphsNext

#load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\local_relaxation\targeted\ctn_pachy\target_6\run\\"
#
#filename = "ctn_beta_9.4000_t_max_15.0000_t_gradient_14.0000_nr_vertices_1792.gml"

#spatial_network = NG.load_spatial_network_from_gml(
#    load_path*filename
#) 
#
#
#correlation_functions_dict = NA.get_correlation_functions(
#    spatial_network;
#    distance_histogram_bin_width = 0.02,
#    periodic_boundary_conditions = true,
#    save_result = false)

#load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\networks\pachy\\"
#
#analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\networks\pachy\\"
#
#filename = "pachy_blue.gml"
#spatial_network = NG.load_spatial_network_from_gml(
#    load_path*filename)
#
#correlation_functions_dict = NA.get_correlation_functions(
#    spatial_network;
#    distance_histogram_bin_width = 0.02,
#    periodic_boundary_conditions = true,
#    save_result = true,
#    save_path = analysis_data_path* 
#        "pachy_blue"
#    )

Plots.plot(
    correlation_functions_dict["vertex_distance_vec"],
    correlation_functions_dict["cumulative_coord_nr_vec"];
    xlabel = "Distance r",
    ylabel = "Cumulative Z",
    bottom_margin = 2Plots.mm,
    xlims = (0, 2),
    ylims = (0, 6)
)
Plots.plot!(
    correlation_functions_dict["vertex_distance_vec"],
    correlation_functions_dict["cumulative_coord_nr_uncoordinated_vec"];
)

# plot a horizontal line at y=1
Plots.hline!([1.0]; linestyle = :dash, color = :black, label = "r = 1.0")