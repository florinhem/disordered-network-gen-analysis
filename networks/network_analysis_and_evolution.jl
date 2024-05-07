"""
these functions are used to evolve and eventually analyze networks
"""


"""
Get a dictionary with all quantities describing the generation
and evolution of a network
"""
function get_evolution_dict(;
    nr_vertices::Int64 = 512 , 
    nr_dimensions::Int64 = 3, 
    network_type::String = "diamond",
    bond_bending_const::Float64 = 0.285,
    min_ring_size::Int64 = 3,
    nr_max_relaxation_cycles::Int64 = 25,
    reject_during_relaxation_cycle_threshold::Int64 = 5,
    break_at_relative_cluster_energy_change::Float64 = 0.0001,
    shell_nr::Int64 = 4,
    relax_efficiently::Bool = true,
    relaxation_overshoot_factor_r = 1.5,
    relaxation_optimization_parameter_l = 1,
    inefficient_optimization_method = "newton",
    random_evolution_seed::Int64 = -1,
    thermal_fluctuations::Bool = false,
    temperature_vec::Vector = [2, 1],
    nr_monte_carlo_steps_per_temperature_vec::Vector = [10,10]
    )

    # check if the temperature sequence is given correctly
    if (length(temperature_vec) !== 
        length(nr_monte_carlo_steps_per_temperature_vec))
        @error "temperature_vec and nr_monte_carlo_steps_per_temperature_vec
        must have the same length"
    end

    # store all arguments in a dictionary
    evolution_dict = Dict(
    "nr_vertices" => nr_vertices, 
    "nr_dimensions" => nr_dimensions, 
    "network_type" => network_type,
    "bond_bending_const" => bond_bending_const,
    "min_ring_size" => min_ring_size,
    "nr_max_relaxation_cycles" => nr_max_relaxation_cycles,
    "reject_during_relaxation_cycle_threshold" => reject_during_relaxation_cycle_threshold,
    "break_at_relative_cluster_energy_change" => break_at_relative_cluster_energy_change,
    "shell_nr" => shell_nr,
    "relax_efficiently" => relax_efficiently,
    "relaxation_overshoot_factor_r" => relaxation_overshoot_factor_r,
    "relaxation_optimization_parameter_l" => relaxation_optimization_parameter_l,
    "inefficient_optimization_method" => inefficient_optimization_method,
    "random_evolution_seed" => random_evolution_seed,
    "thermal_fluctuations" => thermal_fluctuations,
    "temperature_vec" => temperature_vec,
    "nr_monte_carlo_steps_per_temperature_vec" => nr_monte_carlo_steps_per_temperature_vec
    )

    return evolution_dict

end


"""
Get temperature sequence for immediate heating and maximal temperature and then
cooling at a constant temperature decrease
"""
function get_temperature_sequence_cooling_gradient(start_temperature = 2;
    end_temperature = 0,
    temperature_decrease_per_monte_carlo_step = 0.5, 
    nr_monte_carlo_steps_per_temperature = 0.01,
    quench::Bool = true )

    # calculate number of monte carlo steps during temperature decrease
    nr_monte_carlo_steps_during_temperature_decrease = ((start_temperature - end_temperature)
    /temperature_decrease_per_monte_carlo_step)

    # create the vector of temperatures
    temperature_vec = (start_temperature 
        .- temperature_decrease_per_monte_carlo_step 
            .* collect(0
                :nr_monte_carlo_steps_per_temperature
                :nr_monte_carlo_steps_during_temperature_decrease))

    # create vector of monte carlo steps per temperature
    nr_monte_carlo_steps_per_temperature_vec = vcat([2], ones(length(temperature_vec)-1) .* nr_monte_carlo_steps_per_temperature )

    # add long quenching time in the end if desired
    if quench
        if end_temperature == 0
            nr_monte_carlo_steps_per_temperature_vec[end] = 50
        else
            push!(temperature_vec, 0)
            push!(nr_monte_carlo_steps_per_temperature_vec, 50)
        end
    end

    return [temperature_vec, nr_monte_carlo_steps_per_temperature_vec]

end


"""
Get temperature sequence for heating and then cooling at
a constant temperature gradient
"""
function get_temperature_sequence_heating_cooling_gradient(maximal_temperature = 2;
    temperature_increase_per_monte_carlo_step = 0.5, 
    nr_monte_carlo_steps_per_temperature = 0.01,
    quench::Bool = true )

    # calculate number of monte carlo steps during temperature increase
    nr_monte_carlo_steps_during_temperature_increase = (maximal_temperature
    /temperature_increase_per_monte_carlo_step)

    # create the vector of temperatures for the temperature increase
    temperature_increase_vec = ( temperature_increase_per_monte_carlo_step 
            .* collect(0
                :nr_monte_carlo_steps_per_temperature
                :nr_monte_carlo_steps_during_temperature_increase))

    # create the vector of temperatures for the temperature decrease
    temperature_decrease_vec = ( maximal_temperature 
        .- temperature_increase_per_monte_carlo_step 
            .* collect(0
                :nr_monte_carlo_steps_per_temperature
                :nr_monte_carlo_steps_during_temperature_increase))[2:end]

    # create total temperatre vector
    temperature_vec = vcat(temperature_increase_vec, temperature_decrease_vec)

    # create vector of monte carlo steps per temperature
    nr_monte_carlo_steps_per_temperature_vec = (
        ones(length(temperature_vec)) .* nr_monte_carlo_steps_per_temperature )

    # add long quenching time in the end if desired
    if quench
        nr_monte_carlo_steps_per_temperature_vec[end] = 50
    end

    return [temperature_vec, nr_monte_carlo_steps_per_temperature_vec]

end


"""
Evolve a given network and calculate several order metrics
"""


"""
Generate and evolve multiple networks and calculate the average
of several order metrics over these multiple networks
"""
function get_averaged_order_metrics(;save_result = false,
    load_path_without_format = raw"..\structures\random_networks\sample_name",
    save_path = raw"..\analysis_data\random_networks\sample_name",
    label = nothing,
    nr_samples = 5)

    # initialize vectors and arrays for order metrics
    bond_length_std_vec = Vector{Fl}

    # loop through samples 
    for i in 1:nr_samples


    end

    # create dict to save
    order_metrics_dict = Dict("wavevector_array" => wavevector_array,
                        "wavenumber_vec_vec" => wavenumber_vec_vec,
                        "spectral_density_array" => spectral_density_array,
                        "nr_measurements_per_direction" 
                                => autocovariance_fct_dict["nr_measurements_per_direction"],
                        "voxel_edge_length" => autocovariance_fct_dict["voxel_edge_length"],
                        "label" => autocovariance_fct_dict["label"])

    # if desired, adjust voxel edge length and label
    order_metrics_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(order_metrics_dict);
                save_path=save_path*"_order_metrics.h5")

    end

    return order_metrics_dict

end
