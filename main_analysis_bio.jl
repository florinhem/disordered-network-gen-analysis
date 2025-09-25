
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU


network_path = "../structures/biological/networks/pachy/"

analysis_data_path = raw"../analysis_data/biological/networks/pachy/"

digital_sphere_mask_path = "../analysis_data/biological/networks/digital_sphere_masks/"

save_filename = "pachy_blue"

spatial_network = NG.load_spatial_network_from_gml(
    network_path*save_filename*".gml")

pore_size_distribution_dict = NA.get_pore_size_distribution(
    spatial_network;
    sampling_grid_size = 0.2,
    max_pore_radius = 3.0,
    periodic_boundary_conditions = false,
    save_result = true,
    save_path = analysis_data_path*save_filename,
    label = nothing,
    digital_sphere_mask_path 
        = digital_sphere_mask_path,
    print_progress = true)


save_filename = "pachy_red"

spatial_network = NG.load_spatial_network_from_gml(
    network_path*save_filename*".gml")

pore_size_distribution_dict = NA.get_pore_size_distribution(
    spatial_network;
    sampling_grid_size = 0.2,
    max_pore_radius = 3.0,
    periodic_boundary_conditions = false,
    save_result = true,
    save_path = analysis_data_path*save_filename,
    label = nothing,
    digital_sphere_mask_path 
        = digital_sphere_mask_path,
    print_progress = true)