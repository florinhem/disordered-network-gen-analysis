
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU


save_path = "../structures/neural_network_targeted/dia/1728_vertices_from_216/"

evolution_dicts_directory_path = "../structures/neural_network_targeted/dia/1728_vertices_from_216/evolution_dicts/"
