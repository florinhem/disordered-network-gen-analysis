
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Format

nr_samples = 500

network_type_vec = ["bcu_cn_5_6_7_8"]
nr_vertices_vec = [ 432]


theta_ground_state = 180.0
shell_nr = 4
relax_globally_after_threshold_cycle = false

for (nr_vertices, network_type) in zip(nr_vertices_vec, network_type_vec)

    println("Generating $nr_vertices vertices of type $network_type")

    # choose random beta values for the samples between 5 and 10
    beta_vec = 5 .+ 5 .* rand(nr_samples)

    # get the melting temperature for the beta values
    t_melt_vec = [NA.get_melting_temperature(network_type, beta; relax_globally_after_threshold_cycle=relax_globally_after_threshold_cycle, shell_nr=shell_nr) for beta in beta_vec]

    # get random values of t_max between t_melt/2 and 2*t_melt 
    t_max_vec = t_melt_vec .* (1/2 .+ 3/2 .* rand(nr_samples))

    # get random values of the heating/cooling gradient between t_melt/4 and 
    # 2*t_melt
    t_gradient_vec = t_melt_vec .* (1/4 .+ 7/4 .* rand(nr_samples))

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\local_relaxation\random\\"*network_type*raw"\evolution_dicts_3\\"

    for i in 1:nr_samples
        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(t_max_vec[i],
            temperature_gradient = t_gradient_vec[i], 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,
            temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = beta_vec[i], network_type = network_type,
            theta_ground_state = theta_ground_state,
            relax_globally_after_threshold_cycle = relax_globally_after_threshold_cycle,
            shell_nr = shell_nr,)

        filename = Format.format(network_type*"_beta_{1:.4f}_t_max_{2:.4f}_t_gradient_{3:.4f}", beta_vec[i], t_max_vec[i], t_gradient_vec[i])
        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")
    end
end
