
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Format
import Random

nr_samples = 24

theta_ground_state = 180.0
shell_nr = 4
relax_globally_after_threshold_cycle = false

nr_vertices = 1792 #1792
network_type = "ctn"

# choose random beta values for the samples between 0 and 1
beta = 6.3359
t_max = 8.6070
t_gradient = 11.7170



# Function to generate uniformly random values in [0.9*value, 1.1*value]
function uniform_random_vec(value::Float64, n::Int=100)
    lo, hi = 0.9 * value, 1.1 * value
    return Random.rand(lo:0.0001:hi, n)
end

# Your values
beta        = 6.3359
t_max       = 8.6070
t_gradient  = 11.7170

# Generate lists
beta_vec       = uniform_random_vec(beta, nr_samples)
t_max_vec      = uniform_random_vec(t_max, nr_samples)
t_gradient_vec = uniform_random_vec(t_gradient, nr_samples)

# get random values of the heating/cooling gradient between t_melt/4 and 
# 2*t_melt
save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\local_relaxation\targeted\ctn_pachy\target_8\evolution_dicts\\"


#GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")
for i in 1:nr_samples
    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(t_max_vec[i],
        temperature_gradient = t_gradient_vec[i], 
        nr_monte_carlo_steps_per_temperature = 0.01,
        quench = true )

    # remove the last entry of the nr_monte_carlo_steps_per_temperature_vec
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec[1:end-1]
    # attach 100 times the number 0.5 to the nr_monte_carlo_steps_per_temperature_vec
    nr_monte_carlo_steps_per_temperature_vec = vcat(nr_monte_carlo_steps_per_temperature_vec, fill(0.5, 100))

    # attach 99 times the last temperature to the temperature_vec
    temperature_vec = vcat(temperature_vec, fill(0.0, 99))

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