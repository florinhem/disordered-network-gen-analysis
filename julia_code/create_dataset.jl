# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis
# of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Format

# generate evolution dictionaries for all network types

nr_samples = 100

network_type_vec = ["ctn", "dia", "lcs", "srs", "bcu_cn_5_6_7_8", 
    "pcu_cn_4_5_6"]
nr_vertices_vec = [ 224, 216, 192, 216, 432, 216]

theta_ground_state = 180.0
shell_nr = 4
relax_globally_after_threshold_cycle = false

evolution_dicts_path = "../data/structures/evolution_dicts/"

for (nr_vertices, network_type) in zip(nr_vertices_vec, network_type_vec)

    println("Generating evolution dicts of type $network_type")

    # choose random beta values for the samples between 1 and 5
    beta_vec =  10.0 .* rand(nr_samples)

    # get the melting temperature for the beta values
    t_melt_vec = [NA.get_melting_temperature(network_type, beta; 
        relax_globally_after_threshold_cycle
            =relax_globally_after_threshold_cycle, 
        shell_nr=shell_nr) 
        for beta in beta_vec]

    # get random values of t_max between 1.5*t_melt and 2*t_melt 
    t_max_vec = t_melt_vec .* (1/2 .+ 3/2 .* rand(nr_samples))

    # get random values of the heating/cooling gradient between t_melt/4 and 
    # 2*t_melt
    t_gradient_vec = t_melt_vec .* (1/4 .+ 7/4 .* rand(nr_samples))

    for i in 1:nr_samples
        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = (
            NA.get_temperature_sequence_heating_cooling_gradient(t_max_vec[i],
                temperature_gradient = t_gradient_vec[i], 
                nr_monte_carlo_steps_per_temperature = 0.01,
                quench = true ))

        evolution_dict = NA.get_evolution_dict(;
            nr_vertices = nr_vertices,
            temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec 
                = nr_monte_carlo_steps_per_temperature_vec, 
            min_ring_size = 3,
            bond_bending_const = beta_vec[i], 
            network_type = network_type,
            theta_ground_state = theta_ground_state,
            relax_globally_after_threshold_cycle 
                = relax_globally_after_threshold_cycle,
            shell_nr = shell_nr,)

        filename = Format.format(network_type
            *"_beta_{1:.4f}_t_max_{2:.4f}_t_gradient_{3:.4f}", 
            beta_vec[i], t_max_vec[i], t_gradient_vec[i])
        GU.save_dict_to_h5(evolution_dict, 
            evolution_dicts_path*filename*"_evolution.h5")
    end
end

# generate networks

spatial_network_path = "../data/structures/" 

print_every_nr_attempted_bond_switches = 200
print_progress = true
save_network_after_each_temperature = false
further_evolve_previous_networks = false
runs_vec = [1] 
random_evolution_seed = -1
print_lock = Threads.ReentrantLock()

NG.generate_spatial_networks_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_path,
    spatial_network_path;
    print_every_nr_attempted_bond_switches
        =print_every_nr_attempted_bond_switches,
    print_progress=print_progress,
    random_evolution_seed=random_evolution_seed,
    save_network_after_each_temperature=save_network_after_each_temperature,
    further_evolve_previous_networks=further_evolve_previous_networks,
    runs_vec=runs_vec,
    print_lock=print_lock)


# analyze networks

analysis_data_path = "../data/analysis_data/" 

NA.get_all_dicts_from_networks_multithreading(
    spatial_network_path,
    analysis_data_path;
    print_progress = print_progress,
    runs_vec = runs_vec,
    print_lock = print_lock)


# for each network type, create one dictionary with the order metrics of all 
# networks

for network_type in network_type_vec
    analysis_data_path_network_type = analysis_data_path *
        network_type * "/"

    order_metrics_dict = NA.get_order_metrics_all_files(
        analysis_data_path_network_type;
        save_result=true,
        save_algorithm_parameters_from_filename=true)
    
    # create corresponding CSV file
    NA.save_order_metrics_dict_to_csv(order_metrics_dict, 
        analysis_data_path_network_type)
end


# create arrays of melting temperatures for ctn networks with different bond
# bending constants. This will be used to access the melting temperature
# for order metric prediction or data plotting.

bond_bending_const_vec = collect(0.0:0.2:10.0)
t_melt_vec = [NA.get_melting_temperature("ctn", beta; 
        relax_globally_after_threshold_cycle=false, shell_nr=4) 
    for beta in bond_bending_const_vec]

t_gradient_vec = 0.5 .* t_melt_vec
t_max_over_t_melt_vec = collect(0.5:0.02:2.0)
t_max_arr = t_max_over_t_melt_vec' .* t_melt_vec

# save the arrays to a h5 file
save_path = "../data/analysis_data/ctn/" 

GU.save_dict_to_h5(Dict(
    "bond_bending_const_vec" => bond_bending_const_vec,
    "t_melt_vec" => t_melt_vec,
    "t_gradient_vec" => t_gradient_vec,
    "t_max_over_t_melt_vec" => t_max_over_t_melt_vec,
    "t_max_arr" => t_max_arr,
), save_path*"ctn_t_melt_vec.h5")