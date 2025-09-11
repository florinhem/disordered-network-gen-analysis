
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Format

beta = 0.25
nr_vertices = 1000
network_type = "dia"
theta_ground_state = 180.0

t_max_216 = 0.4592
t_gradient_216 = 0.2449

t_melt_216 = 0.8912301132026756
t_melt_1000 = 5.520317556297713

t_max = t_max_216 * (t_melt_1000 / t_melt_216)
t_gradient = t_gradient_216 * (t_melt_1000 / t_melt_216)

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\neural_network_targeted\dia\1000_vertices_old_generation_fct\evolution_dicts_3\\"

temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(t_max,
    temperature_gradient = t_gradient, 
    nr_monte_carlo_steps_per_temperature = 0.01,
    quench = true )
evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,     
    temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
    bond_bending_const = beta, network_type = network_type,
    theta_ground_state = theta_ground_state,)
filename = Format.format("dia_beta_{1:.4f}_t_max_{2:.4f}_t_gradient_{3:.4f}", beta, t_max, t_gradient)
GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")