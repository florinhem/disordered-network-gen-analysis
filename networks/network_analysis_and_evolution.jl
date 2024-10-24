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
    nr_monte_carlo_steps_per_temperature_vec::Vector = [10,10],
    mean_nr_monte_carlo_steps_for_quenching::Float64 = 13.7,
    relax_globally_after_threshold_cycle::Bool = true,
    theta_ground_state::Float64 = 109.5
    )

    # check if the temperature sequence is given correctly
    if (length(temperature_vec) !== 
        length(nr_monte_carlo_steps_per_temperature_vec))
        @error "temperature_vec and nr_monte_carlo_steps_per_temperature_vec
        must have the same length"
    end

    # TODO if mit diamond, cubic, bcc, fcc => coordination nbr (und dann auch in den evol_dict rein)

    # estimate the total number of bond switches
    estimated_nr_bond_switches = 0
    for i in eachindex(temperature_vec)
        if temperature_vec[i] == 0 && nr_monte_carlo_steps_per_temperature_vec[i] == 50
            estimated_nr_bond_switches += 
                mean_nr_monte_carlo_steps_for_quenching * 18 * nr_vertices
                # TODO convert 18 (hardcoded) into 4*3*3/2 with coordination_nr
        else
            estimated_nr_bond_switches += 
                nr_monte_carlo_steps_per_temperature_vec[i] * 18 * nr_vertices
                # TODO convert 18 (hardcoded) into 4*3*3/2 with coordination_nr
        end
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
    "nr_monte_carlo_steps_per_temperature_vec" => nr_monte_carlo_steps_per_temperature_vec,
    "mean_nr_monte_carlo_steps_for_quenching" => mean_nr_monte_carlo_steps_for_quenching,
    "estimated_nr_bond_switches" => estimated_nr_bond_switches,
    "relax_globally_after_threshold_cycle" => relax_globally_after_threshold_cycle,
    "theta_ground_state" => theta_ground_state
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
    nr_monte_carlo_steps_during_temperature_decrease = 
    ((start_temperature - end_temperature)
    /temperature_decrease_per_monte_carlo_step)

    # create the vector of temperatures
    temperature_vec = (start_temperature 
        .- temperature_decrease_per_monte_carlo_step 
            .* collect(0
                :nr_monte_carlo_steps_per_temperature
                :nr_monte_carlo_steps_during_temperature_decrease))

    # create vector of monte carlo steps per temperature
    nr_monte_carlo_steps_per_temperature_vec = vcat([2], 
        ones(length(temperature_vec)-1) .* nr_monte_carlo_steps_per_temperature )

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
a constant temperature gradient. The temperature increase per monte carlo step
is the gradient of heating. The argument 'nr_monte_carlo_steps_per_temperature'
only sets, how finely the temperature is supposed to be sampled.
"""
function get_temperature_sequence_heating_cooling_gradient(
    maximal_temperature = 2;
    temperature_gradient = 0.5, 
    nr_monte_carlo_steps_per_temperature = 0.01,
    quench::Bool = true )

    # calculate number of monte carlo steps during temperature increase
    nr_monte_carlo_steps_during_temperature_increase = (maximal_temperature
    /temperature_gradient)

    # create the vector of temperatures for the temperature increase
    temperature_increase_vec = ( temperature_gradient 
            .* collect(0
                :nr_monte_carlo_steps_per_temperature
                :nr_monte_carlo_steps_during_temperature_increase))

    # create the vector of temperatures for the temperature decrease
    temperature_decrease_vec = ( maximal_temperature 
        .- temperature_gradient 
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
