
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots

spatial_network_path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\\"
filename=raw"m_BTMC_N=216_T=0.1_Beta=1.0_GradT=0.01_StepsPerT=0.001_Theta_GS=180.0_Trial=1"
structure_dict_path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\structure\\"
analysis_data_path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\analysis_data\\"


NA.get_all_dicts_from_network_single_file(
    filename,
    spatial_network_path,
    structure_dict_path,
    analysis_data_path;
    bond_radius = 0.35,
    voxel_edge_length = 0.25,
    print_progress = true,
    print_lock = Threads.ReentrantLock())



path = raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\pore_size_dist\\"

network_path = structure_dict_path * filename *"_structure.h5"

structure_dict_network = GU.load_h5_dict(network_path)

pore_pixel_radius_array = NA.get_pore_size_distribution(structure_dict_network)

pore_pixel_radius_vec= pore_pixel_radius_array["pore_size_distribution"]

println(pore_pixel_radius_vec)

pore_pixel_radius_filtered_vec = pore_pixel_radius_vec[pore_pixel_radius_vec .> 0.0]
Plots.histogram(pore_pixel_radius_filtered_vec)
Plots.savefig(path*"pore_size_distribution_"* filename *".png")