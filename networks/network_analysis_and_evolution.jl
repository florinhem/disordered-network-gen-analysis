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
similar sizes. The following data is for networks that were relaxed globally
coming with the problem that the melting temperatures change drastically for
larger networks:
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
["srs", 216, 0.0, 180.0, 0.001, 0.00019026075258227854 ],
["srs", 216, 0.25, 180.0, 0.001, 0.45575800978228587 ],
["srs", 216, 0.5, 180.0, 0.001, 1.0842146353235373 ],
["srs", 216, 0.75, 180.0, 0.001, 1.79142389890434 ],
["srs", 216, 1.0, 180.0, 0.001, 2.4612236029485097 ],
["bcu_cn_5_6_7_8", 432, 0.0, 180.0, 0.001, true, 4, 0.44726145342442325 ],
["bcu_cn_5_6_7_8", 432, 0.25, 180.0, 0.001, true, 4, 100.65980650498858 ],
["bcu_cn_5_6_7_8", 432, 0.5, 180.0, 0.001, true, 4, 262.6291684535291 ],
["bcu_cn_5_6_7_8", 432, 0.75, 180.0, 0.001, true, 4, 431.6794774973545 ],
["bcu_cn_5_6_7_8", 432, 1.0, 180.0, 0.001, true, 4, 577.6926948220212 ],
["pcu_cn_4_5_6", 216, 0.0, 180.0, 0.001, true, 4, 0.0003545596706363815 ],
["pcu_cn_4_5_6", 216, 0.25, 180.0, 0.001, true, 4, 1.5104535301442776 ],
["pcu_cn_4_5_6", 216, 0.5, 180.0, 0.001, true, 4, 20.535890429251697 ],
["pcu_cn_4_5_6", 216, 0.75, 180.0, 0.001, true, 4, 41.871921656928805 ],
["pcu_cn_4_5_6", 216, 1.0, 180.0, 0.001, true, 4, 67.45750150259808 ],

The following data is for networks that were relaxed up to the 4th neighbor 
shell:
["ctn", 224, 0.0, 180.0, 0.001, false, 4, 0.00041499170346413945 ],
["ctn", 224, 0.25, 180.0, 0.001, false, 4, 0.1252534628370132 ],
["ctn", 224, 0.5, 180.0, 0.001, false, 4, 0.30497680344164985 ],
["ctn", 224, 0.75, 180.0, 0.001, false, 4, 0.4352128146806752 ],
["ctn", 224, 1.0, 180.0, 0.001, false, 4, 0.6287240577847529 ],
["dia", 216, 0.0, 180.0, 0.001, false, 4, 0.0003252933056548639 ],
["dia", 216, 0.25, 180.0, 0.001, false, 4, 0.11916487389418548 ],
["dia", 216, 0.5, 180.0, 0.001, false, 4, 0.28976530223293046 ],
["dia", 216, 0.75, 180.0, 0.001, false, 4, 0.4545309351335605 ],
["dia", 216, 1.0, 180.0, 0.001, false, 4, 0.7733782038178849 ],
["lcs", 192, 0.0, 180.0, 0.001, false, 4, 0.00017379958338485042 ],
["lcs", 192, 0.25, 180.0, 0.001, false, 4, 0.14250906823844808 ],
["lcs", 192, 0.5, 180.0, 0.001, false, 4, 0.34873519151914834 ],
["lcs", 192, 0.75, 180.0, 0.001, false, 4, 0.5638818094183685 ],
["lcs", 192, 1.0, 180.0, 0.001, false, 4, 0.9224890576183211 ],
["srs", 216, 0.0, 180.0, 0.001, false, 4, 0.0006744918502233107 ],
["srs", 216, 0.25, 180.0, 0.001, false, 4, 0.04589951551091327 ],
["srs", 216, 0.5, 180.0, 0.001, false, 4, 0.1020805728671242 ],
["srs", 216, 0.75, 180.0, 0.001, false, 4, 0.16119291348132522 ],
["srs", 216, 1.0, 180.0, 0.001, false, 4, 0.24773319499055632 ],
["bcu_cn_5_6_7_8", 432, 0.0, 180.0, 0.001, false, 4, 0.2108305865744253 ],
["bcu_cn_5_6_7_8", 432, 0.25, 180.0, 0.001, false, 4, 27.370479812098893 ],
["bcu_cn_5_6_7_8", 432, 0.5, 180.0, 0.001, false, 4, 75.58737112513226 ],
["bcu_cn_5_6_7_8", 432, 0.75, 180.0, 0.001, false, 4, 117.53668306629757 ],
["bcu_cn_5_6_7_8", 432, 1.0, 180.0, 0.001, false, 4, 166.69221669496523 ],
["pcu_cn_4_5_6", 216, 0.0, 180.0, 0.001, false, 4, 0.0003996287259513749 ],
["pcu_cn_4_5_6", 216, 0.25, 180.0, 0.001, false, 4, 0.44142359212153215 ],
["pcu_cn_4_5_6", 216, 0.5, 180.0, 0.001, false, 4, 8.537230402733845 ],
["pcu_cn_4_5_6", 216, 0.75, 180.0, 0.001, false, 4, 17.27276876274271 ],
["pcu_cn_4_5_6", 216, 1.0, 180.0, 0.001, false, 4, 28.542212916214314 ],

The following data is for networks that were relaxed up to the 3rd neighbor 
shell where one expects the smallest dependency on system size:
["ctn", 224, 0.0, 180.0, 0.001, false, 3, 0.0006011552485005689 ],
["ctn", 224, 0.25, 180.0, 0.001, false, 3, 0.056329651723466774 ],
["ctn", 224, 0.5, 180.0, 0.001, false, 3, 0.12025892877483091 ],
["ctn", 224, 0.75, 180.0, 0.001, false, 3, 0.1817910450163424 ],
["ctn", 224, 1.0, 180.0, 0.001, false, 3, 0.2837122928706454 ],
["dia", 216, 0.0, 180.0, 0.001, false, 3, 0.0004500125750497054 ],
["dia", 216, 0.25, 180.0, 0.001, false, 3, 0.07214908414976968 ],
["dia", 216, 0.5, 180.0, 0.001, false, 3, 0.1662209382225799 ],
["dia", 216, 0.75, 180.0, 0.001, false, 3, 0.2640566152053155 ],
["dia", 216, 1.0, 180.0, 0.001, false, 3, 0.4195846024594645 ],
["lcs", 192, 0.0, 180.0, 0.001, false, 3, 0.00025383038096432386 ],
["lcs", 192, 0.25, 180.0, 0.001, false, 3, 0.07345060880056713 ],
["lcs", 192, 0.5, 180.0, 0.001, false, 3, 0.16290517832054138 ],
["lcs", 192, 0.75, 180.0, 0.001, false, 3, 0.28595208776627684 ],
["lcs", 192, 1.0, 180.0, 0.001, false, 3, 0.4296392047793182 ],
["srs", 216, 0.0, 180.0, 0.001, false, 3, 0.0006835069104931692 ],
["srs", 216, 0.25, 180.0, 0.001, false, 3, 0.024276394585506206 ],
["srs", 216, 0.5, 180.0, 0.001, false, 3, 0.06254822987957036 ],
["srs", 216, 0.75, 180.0, 0.001, false, 3, 0.08630480253631828 ],
["srs", 216, 1.0, 180.0, 0.001, false, 3, 0.11034305336366997 ],
["bcu_cn_5_6_7_8", 432, 0.0, 180.0, 0.001, false, 3, 0.0949361861248105 ],
["bcu_cn_5_6_7_8", 432, 0.25, 180.0, 0.001, false, 3, 12.070760027826928 ],
["bcu_cn_5_6_7_8", 432, 0.5, 180.0, 0.001, false, 3, 29.429401371933043 ],
["bcu_cn_5_6_7_8", 432, 0.75, 180.0, 0.001, false, 3, 47.62558728446081 ],
["bcu_cn_5_6_7_8", 432, 1.0, 180.0, 0.001, false, 3, 69.75965974214769 ],
["pcu_cn_4_5_6", 216, 0.0, 180.0, 0.001, false, 3, 0.0006445035631994268 ],
["pcu_cn_4_5_6", 216, 0.25, 180.0, 0.001, false, 3, 0.3610517424915666 ],
["pcu_cn_4_5_6", 216, 0.5, 180.0, 0.001, false, 3, 3.7569276899170836 ],
["pcu_cn_4_5_6", 216, 0.75, 180.0, 0.001, false, 3, 8.182869022170127 ],
["pcu_cn_4_5_6", 216, 1.0, 180.0, 0.001, false, 3, 13.242477453061841 ],
"""
function get_melting_temperature(
    network_type::String = "diamond",
    beta::Float64 = 0.0;
    relax_globally_after_threshold_cycle::Bool = false,
    shell_nr::Int64 = 4,)

    # get the data that was calculated separately 
    beta_vec = collect(0:0.25:1.0)

    if relax_globally_after_threshold_cycle
        t_melt_ctn = [0.0001294330930655222, 0.5068377943704439, 
        1.2456941705788351, 2.0583613664213662, 2.8256795255231606]
        t_melt_dia = [0.00013619791956073077, 0.35741585656533714, 
            0.8893834492799765, 1.4827028269123552, 2.3973188103040246]
        t_melt_lcs = [8.054647440611615e-5, 0.4029048676264969, 
            0.983114477463025, 1.663819575318196, 2.82353242478514]
        t_melt_srs = [0.00019026075258227854, 0.45575800978228587, 
            1.0842146353235373, 1.79142389890434, 2.4612236029485097]
        t_melt_bcu_cn_5_6_7_8 = [0.44726145342442325, 100.65980650498858, 
            262.6291684535291, 431.6794774973545, 577.6926948220212]
        t_melt_pcu_cn_4_5_6 = [0.0003545596706363815, 1.5104535301442776, 
            20.535890429251697, 41.871921656928805, 67.45750150259808]
    
    else
        if shell_nr != 3 && shell_nr != 4
            @error "shell_nr must be 3 or 4"
        elseif shell_nr == 4
            t_melt_ctn = [0.00041499170346413945, 0.1252534628370132, 
                0.30497680344164985, 0.4352128146806752, 0.6287240577847529]
            t_melt_dia = [0.0003252933056548639, 0.11916487389418548, 
                0.28976530223293046, 0.4545309351335605, 0.7733782038178849]
            t_melt_lcs = [0.00017379958338485042, 0.14250906823844808, 
                0.34873519151914834, 0.5638818094183685, 0.9224890576183211]
            t_melt_srs = [0.0006744918502233107, 0.04589951551091327, 
                0.1020805728671242, 0.16119291348132522, 0.24773319499055632]
            t_melt_bcu_cn_5_6_7_8 = [0.2108305865744253, 27.370479812098893, 
                75.58737112513226, 117.53668306629757, 166.69221669496523]
            t_melt_pcu_cn_4_5_6 = [0.0003996287259513749, 0.44142359212153215, 
                8.537230402733845, 17.27276876274271, 28.542212916214314]

        else 
            t_melt_ctn = [0.0006011552485005689, 0.056329651723466774, 
                0.12025892877483091, 0.1817910450163424, 0.2837122928706454]
            t_melt_dia = [0.0004500125750497054, 0.07214908414976968, 
                0.1662209382225799, 0.2640566152053155, 0.4195846024594645]
            t_melt_lcs = [0.00025383038096432386, 0.07345060880056713, 
                0.16290517832054138, 0.28595208776627684, 0.4296392047793182]
            t_melt_srs = [0.0006835069104931692, 0.024276394585506206, 
                0.06254822987957036, 0.08630480253631828, 0.11034305336366997]
            t_melt_bcu_cn_5_6_7_8 = [0.0949361861248105, 12.070760027826928, 
                29.429401371933043, 47.62558728446081, 69.75965974214769]
            t_melt_pcu_cn_4_5_6 = [0.0006445035631994268, 0.3610517424915666, 
                3.7569276899170836, 8.182869022170127, 13.242477453061841]
        end
    end
    

    # interpolate linearly the melting temperature for the given beta value
    if network_type == "diamond" || network_type == "dia"
        t_melt = t_melt_dia
    elseif network_type == "ctn"
        t_melt = t_melt_ctn
    elseif network_type == "lcs"
        t_melt = t_melt_lcs
    elseif network_type == "srs"
        t_melt = t_melt_srs
    elseif network_type == "bcu_cn_5_6_7_8"
        t_melt = t_melt_bcu_cn_5_6_7_8
    elseif network_type == "pcu_cn_4_5_6"
        t_melt = t_melt_pcu_cn_4_5_6
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


