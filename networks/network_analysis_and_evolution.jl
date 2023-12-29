"""
these functions are used to evolve and eventually analyze networks
"""


"""
Get a dictionary with all quantities describing the generation
and evolution of a network
"""
function get_network_evolution_dict(;
    nr_vertices::Int64 = 27 , 
    nr_dimensions::Int64 = 3, 
    network_type = "diamond",
    bond_bending_const::Real = 0.285,
    total_energy_fct = get_total_energy_keating,
    nr_max_relaxation_cycles::Int64 = 25,
    reject_during_relaxation_cycle_threshold::Int64 = 10,
    break_at_relative_cluster_energy_change::Float64 = 0.0001,
    shell_nr::Int64 = 5,
    relax_efficiently::Bool = true,
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    print_progress::Bool = false,
    random_evolution_seed = nothing,
    thermal_fluctuations::Bool = false)

    #store all arguments in a dictionary
    network_evolution_dict = Dict(
    "nr_vertices" => nr_vertices, 
    "nr_dimensions" => nr_dimensions, 
    "network_type" => network_type,
    "bond_bending_const" => bond_bending_const,
    "total_energy_fct" => total_energy_fct,
    "nr_max_relaxation_cycles" => nr_max_relaxation_cycles,
    "reject_during_relaxation_cycle_threshold" => reject_during_relaxation_cycle_threshold,
    "break_at_relative_cluster_energy_change" => break_at_relative_cluster_energy_change,
    "shell_nr" => shell_nr,
    "relax_efficiently" => relax_efficiently,
    "relaxation_overshoot_factor_r" => relaxation_overshoot_factor_r,
    "relaxation_optimization_parameter_l" => relaxation_optimization_parameter_l,
    "print_progress" => print_progress,
    "random_evolution_seed" => random_evolution_seed,
    "thermal_fluctuations" => thermal_fluctuations
    )

    return network_evolution_dict

end


"""
Evolve a given network and calculate several order metrics
"""


"""
Generate and evolve multiple networks and calculate the average
of several order metrics over these multiple networks
"""
function get_averaged_order_metrics(;save_result = false,
    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\random_networks\sample_name",
    label = nothing)

    #create dict to save
    order_metrics_dict = Dict("wavevector_array" => wavevector_array,
                        "wavenumber_vec_vec" => wavenumber_vec_vec,
                        "spectral_density_array" => spectral_density_array,
                        "nr_measurements_per_direction" 
                                => autocovariance_fct_dict["nr_measurements_per_direction"],
                        "voxel_edge_length" => autocovariance_fct_dict["voxel_edge_length"],
                        "label" => autocovariance_fct_dict["label"])

    #if desired, adjust voxel edge length and label
    order_metrics_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    #save results if desired
    if save_result
        GU.save_dict_to_h5(copy(order_metrics_dict);
                save_path=save_path*"_order_metrics.h5")

    end

    return order_metrics_dict

end
