
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import Format
import StatsBase

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\crystals\\"

filename = "dia_216_vertices"

evolution_dict = NA.get_evolution_dict(;
            nr_vertices = 216, 
            network_type="dia",
            )

spatial_network = NG.get_periodic_network(evolution_dict)

NG.save_spatial_network_to_gml(
            spatial_network,
            filename;
            save_path = save_path)