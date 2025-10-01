
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

network_path = "../structures/biological/networks/pachy/"
analysis_data_path = "../analysis_data/biological/networks/pachy/"
save_filename = "pachy_blue"

spatial_network = NG.load_spatial_network_from_gml(
    network_path*save_filename*".gml")

consider_bonds = false
maximal_wavevector_int = 5
periodic_boundary_conditions = false
save_result = true
save_path = analysis_data_path*save_filename
print_progress = true
gaussian_filter = true
gaussian_filter_sigma_x = 2*pi/25
gaussian_filter_filtered_data_x_step_length = 2*pi/25

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(
    spatial_network;
    consider_bonds = consider_bonds,
    maximal_wavevector_int = maximal_wavevector_int,
    periodic_boundary_conditions = periodic_boundary_conditions,
    wavevector_array_positive_z = NA.get_wavevector_array_positive_z(spatial_network; 
        maximal_wavevector_int=maximal_wavevector_int,
            periodic_boundary_conditions=periodic_boundary_conditions),
    save_result = save_result,
    save_path = save_path,
    print_progress = print_progress,
    thread_nr = 0,
    print_lock = Threads.ReentrantLock())

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(
        structure_factor_dict;
        consider_bonds = consider_bonds,
        gaussian_filter = gaussian_filter,
        gaussian_filter_sigma_x = gaussian_filter_sigma_x,
        gaussian_filter_filtered_data_x_step_length = gaussian_filter_filtered_data_x_step_length,
        save_result = save_result,
        save_path = save_path)

consider_bonds = true

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(
    spatial_network;
    consider_bonds = consider_bonds,
    maximal_wavevector_int = maximal_wavevector_int,
    periodic_boundary_conditions = periodic_boundary_conditions,
    wavevector_array_positive_z = NA.get_wavevector_array_positive_z(spatial_network; 
        maximal_wavevector_int=maximal_wavevector_int,
            periodic_boundary_conditions=periodic_boundary_conditions),
    save_result = save_result,
    save_path = save_path,
    print_progress = print_progress,
    thread_nr = 0,
    print_lock = Threads.ReentrantLock())

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(
        structure_factor_dict;
        consider_bonds = consider_bonds,
        gaussian_filter = gaussian_filter,
        gaussian_filter_sigma_x = gaussian_filter_sigma_x,
        gaussian_filter_filtered_data_x_step_length = gaussian_filter_filtered_data_x_step_length,
        save_result = save_result,
        save_path = save_path)