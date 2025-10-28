
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Format



theta_ground_state = 180.0
shell_nr = 4
relax_globally_after_threshold_cycle = false

nr_vertices = 224 #1792
network_type = "ctn"


# choose random beta values for the samples between 0 and 1
beta = 9.4
t_max = 15.0
t_gradient = 14.0

# get random values of the heating/cooling gradient between t_melt/4 and 
# 2*t_melt
save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\local_relaxation\targeted\ctn_pachy\target_6\evolution_dicts\\"

temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(t_max,
        temperature_gradient = t_gradient, 
        nr_monte_carlo_steps_per_temperature = 0.01,
        quench = true )

evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,
    temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
    bond_bending_const = beta, network_type = network_type,
    theta_ground_state = theta_ground_state,
    relax_globally_after_threshold_cycle = relax_globally_after_threshold_cycle,
    shell_nr = shell_nr,)
filename = Format.format(network_type*"_beta_{1:.4f}_t_max_{2:.4f}_t_gradient_{3:.4f}_nr_vertices_{4}", beta, t_max, t_gradient, nr_vertices)

GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")
