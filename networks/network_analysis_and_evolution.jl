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

    # estimate the total number of bond switches
    estimated_nr_bond_switches = 0
    for i in eachindex(temperature_vec)
        if temperature_vec[i] == 0 && nr_monte_carlo_steps_per_temperature_vec[i] == 50
            estimated_nr_bond_switches += 
                mean_nr_monte_carlo_steps_for_quenching * 18 * nr_vertices
        else
            estimated_nr_bond_switches += 
                nr_monte_carlo_steps_per_temperature_vec[i] * 18 * nr_vertices
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


"""
Get the melting temperature for a given network type and beta value. The values
are valid for theta_ground_state=180. The melting temperature is defined as the
temperature where the bond switch with the smallest energy change is accepted
with a probability of 0.001. These values were calculated for networks with
similar sizes. The used data is given by
["ctn", 224, 0.0, 180.0, 0.001, 0.0001294330930655222 ],
["ctn", 224, 0.25, 180.0, 0.001, 0.5068377943704439 ],
["ctn", 224, 0.5, 180.0, 0.001, 1.2456941705788351 ],
["ctn", 224, 0.75, 180.0, 0.001, 2.0583613664213662 ],
["ctn", 224, 1.0, 180.0, 0.001, 2.8256795255231606 ],
["dia", 216, 0.0, 180.0, 0.001, 0.00013619791956073077 ],
["dia", 216, 0.25, 180.0, 0.001, 0.35741585656533714 ],
["dia", 216, 0.5, 180.0, 0.001, 0.8893834492799765 ],
["dia", 216, 0.75, 180.0, 0.001, 1.4827028269123552 ],
["dia", 216, 1.0, 180.0, 0.001, 2.3973188103040246 ],
["lcs", 192, 0.0, 180.0, 0.001, 8.054647440611615e-5 ],
["lcs", 192, 0.25, 180.0, 0.001, 0.4029048676264969 ],
["lcs", 192, 0.5, 180.0, 0.001, 0.983114477463025 ],
["lcs", 192, 0.75, 180.0, 0.001, 1.663819575318196 ],
["lcs", 192, 1.0, 180.0, 0.001, 2.82353242478514 ],
["pto", 112, 0.0, 180.0, 0.001, 0.0013268971997985944 ],
["pto", 112, 0.25, 180.0, 0.001, 0.343365472938612 ],
["pto", 112, 0.5, 180.0, 0.001, 0.7280240257823797 ],
["pto", 112, 0.75, 180.0, 0.001, 1.1186689643128152 ],
["pto", 112, 1.0, 180.0, 0.001, 1.5296735456670358 ],
["srd", 270, 0.0, 180.0, 0.001, -0.9593886478273104 ],
["srd", 270, 0.25, 180.0, 0.001, 0.8489788394337111 ],
["srd", 270, 0.5, 180.0, 0.001, 2.266726051364786 ],
["srd", 270, 0.75, 180.0, 0.001, 3.439157715418838 ],
["srd", 270, 1.0, 180.0, 0.001, 4.678955840938107 ],
["srs", 216, 0.0, 180.0, 0.001, 0.00019026075258227854 ],
["srs", 216, 0.25, 180.0, 0.001, 0.45575800978228587 ],
["srs", 216, 0.5, 180.0, 0.001, 1.0842146353235373 ],
["srs", 216, 0.75, 180.0, 0.001, 1.79142389890434 ],
["srs", 216, 1.0, 180.0, 0.001, 2.4612236029485097 ],
"""
function get_melting_temperature(
    network_type::String = "diamond",
    beta::Float64 = 0.0)

    # get the data that was calculated separately 
    beta_vec = collect(0:0.25:1.0)
    t_melt_ctn = [0.0001294330930655222, 0.5068377943704439, 
        1.2456941705788351, 2.0583613664213662, 2.8256795255231606]
    t_melt_dia = [0.00013619791956073077, 0.35741585656533714, 
        0.8893834492799765, 1.4827028269123552, 2.3973188103040246]
    t_melt_lcs = [8.054647440611615e-5, 0.4029048676264969, 
        0.983114477463025, 1.663819575318196, 2.82353242478514]
    t_melt_pto = [0.0013268971997985944, 0.343365472938612, 
        0.7280240257823797, 1.1186689643128152, 1.5296735456670358]
    t_melt_srd = [-0.9593886478273104, 0.8489788394337111, 
        2.266726051364786, 3.439157715418838, 4.678955840938107]
    t_melt_srs = [0.00019026075258227854, 0.45575800978228587, 
        1.0842146353235373, 1.79142389890434, 2.4612236029485097]

    # interpolate linearly the melting temperature for the given beta value
    if network_type == "diamond" || network_type == "dia"
        t_melt = t_melt_dia
    elseif network_type == "ctn"
        t_melt = t_melt_ctn
    elseif network_type == "lcs"
        t_melt = t_melt_lcs
    elseif network_type == "pto"
        t_melt = t_melt_pto
    elseif network_type == "srd"
        t_melt = t_melt_srd
    elseif network_type == "srs"
        t_melt = t_melt_srs
    else
        @error "network type not recognized"
    end 

    # get the melting temperature for the given beta value
    if beta < 0.0 || beta > 1.0
        @error "beta must be between 0 and 1"
    end
    if beta == 0.0
        t_melt_beta = t_melt[1]
    elseif beta == 1.0
        t_melt_beta = t_melt[end]
    else
        # get the index of the melting temperature for the given beta value
        index = findfirst(x -> x >= beta, beta_vec)
        if index == nothing
            @error "beta value not found"
        end
        # interpolate linearly the melting temperature for the given beta value
        t_melt_beta = (t_melt[index-1] 
            + (t_melt[index] - t_melt[index-1]) 
            * (beta - beta_vec[index-1]) 
            / (beta_vec[index] - beta_vec[index-1]))
    end
    return t_melt_beta
end


"""
Get the melting temperature for a given network type and beta value. The values
are valid for theta_ground_state=180. The melting temperature is defined as the
temperature where the bond switch with the smallest energy change is accepted
with a probability of 0.06. The melting temperature data was calculated for
mixed network sizes which is why the use of the function 
get_melting_temperature() is recommended. The data used in this function is
data = [
["ctn", 224, 0.0, 180.0, 0.06, 0.00031779651885594445 ],
["ctn", 224, 0.25, 180.0, 0.06, 1.1511070269351324 ],
["ctn", 224, 0.5, 180.0, 0.06, 2.7958110652423622 ],
["ctn", 224, 0.75, 180.0, 0.06, 4.585660340210036 ],
["ctn", 224, 1.0, 180.0, 0.06, 6.415852631785852 ],
["dia", 64, 0.0, 180.0, 0.06, 0.0009204175214628252 ],
["dia", 64, 0.25, 180.0, 0.06, 0.3061062214612626 ],
["dia", 64, 0.5, 180.0, 0.06, 0.6897997003855234 ],
["dia", 64, 0.75, 180.0, 0.06, 1.1801646529256444 ],
["dia", 64, 1.0, 180.0, 0.06, 1.8723883922952942 ],
["lcs", 192, 0.0, 180.0, 0.06, 0.00019776541351310378 ],
["lcs", 192, 0.25, 180.0, 0.06, 1.071663395930038 ],
["lcs", 192, 0.5, 180.0, 0.06, 2.5884136096353902 ],
["lcs", 192, 0.75, 180.0, 0.06, 4.313285672070859 ],
["lcs", 192, 1.0, 180.0, 0.06, 7.128334985976234 ],
["pto", 112, 0.0, 180.0, 0.06, 0.0032579250096590583 ],
["pto", 112, 0.25, 180.0, 0.06, 0.7813861085162169 ],
["pto", 112, 0.5, 180.0, 0.06, 1.5409201350892336 ],
["pto", 112, 0.75, 180.0, 0.06, 2.2759438445283595 ],
["pto", 112, 1.0, 180.0, 0.06, 3.3613117404334374 ],
["srd", 80, 0.0, 180.0, 0.06, -0.7407408623269777 ],
["srd", 80, 0.25, 180.0, 0.06, 0.3300682873059035 ],
["srd", 80, 0.5, 180.0, 0.06, 1.1013913412535 ],
["srd", 80, 0.75, 180.0, 0.06, 1.830137437807958 ],
["srd", 80, 1.0, 180.0, 0.06, 2.9236032169832975 ],
["srs", 64, 0.0, 180.0, 0.06, 0.0013561265418390608 ],
["srs", 64, 0.25, 180.0, 0.06, 0.34187449273557846 ],
["srs", 64, 0.5, 180.0, 0.06, 0.7734640806305035 ],
["srs", 64, 0.75, 180.0, 0.06, 1.2002240438551928 ],
["srs", 64, 1.0, 180.0, 0.06, 1.707019829012384 ],
]
"""
function get_melting_temperature_mixed_network_sizes(
    network_type::String = "diamond",
    beta::Float64 = 0.0)

    # get the data that was calculated separately 
    beta_vec = collect(0:0.25:1.0)
    t_melt_dia = [0.0009204175214628252, 0.3061062214612626, 
        0.6897997003855234, 1.1801646529256444, 1.8723883922952942]
    t_melt_ctn = [0.00031779651885594445, 1.1511070269351324, 
        2.7958110652423622, 4.585660340210036, 6.415852631785852]
    t_melt_lcs = [0.00019776541351310378, 1.071663395930038, 
        2.5884136096353902, 4.313285672070859, 7.128334985976234]
    t_melt_pto = [0.0032579250096590583, 0.7813861085162169, 
        1.5409201350892336, 2.2759438445283595, 3.3613117404334374]
    t_melt_srd = [-0.7407408623269777, 0.3300682873059035, 
        1.1013913412535, 1.830137437807958, 2.9236032169832975]
    t_melt_srs = [0.0013561265418390608, 0.34187449273557846, 
        0.7734640806305035, 1.2002240438551928, 1.707019829012384]

    # interpolate linearly the melting temperature for the given beta value
    if network_type == "diamond" || network_type == "dia"
        t_melt = t_melt_dia
    elseif network_type == "ctn"
        t_melt = t_melt_ctn
    elseif network_type == "lcs"
        t_melt = t_melt_lcs
    elseif network_type == "pto"
        t_melt = t_melt_pto
    elseif network_type == "srd"
        t_melt = t_melt_srd
    elseif network_type == "srs"
        t_melt = t_melt_srs
    else
        @error "network type not recognized"
    end 

    # get the melting temperature for the given beta value
    if beta < 0.0 || beta > 1.0
        @error "beta must be between 0 and 1"
    end
    if beta == 0.0
        t_melt_beta = t_melt[1]
    elseif beta == 1.0
        t_melt_beta = t_melt[end]
    else
        # get the index of the melting temperature for the given beta value
        index = findfirst(x -> x >= beta, beta_vec)
        if index == nothing
            @error "beta value not found"
        end
        # interpolate linearly the melting temperature for the given beta value
        t_melt_beta = (t_melt[index-1] 
            + (t_melt[index] - t_melt[index-1]) 
            * (beta - beta_vec[index-1]) 
            / (beta_vec[index] - beta_vec[index-1]))
    end
    return t_melt_beta
end