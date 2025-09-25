
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import Format
import StatsBase

network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\local_relaxation\targeted\shell_nr_4\run_1\\"

filename = "lcs_nr_vertices_192_beta_0.2500_t_max_0.2138_t_gradient_0.1140.gml"

spatial_network = NG.load_spatial_network_from_gml(
    network_path*filename)

NG.plot_spatial_network(spatial_network)