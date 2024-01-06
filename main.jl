
#include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

#import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import Plots
import Measurements

evolution_dict = NA.get_evolution_dict(;nr_vertices = 250 ,temperature_vec = [3, 0],
nr_monte_carlo_steps_per_temperature_vec = [1,1])

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence(graph_dict,
        evolution_dict;
        print_progress = true,
        random_evolution_seed = 2)

structure_factor_dict = NA.get_structure_factor_isotrope_by_wavenumber_vec(
        graph_dict)

hyperuniformity_parameter = NA.get_effective_hyperuniformity_parameter(structure_factor_dict)
println("hyperuniformity parameter: "*string(hyperuniformity_parameter))

local_nr_variance_dict = NA.get_local_nr_variance_by_window_radius_vec(graph_dict;
structure_factor_dict = structure_factor_dict)

Plots.plot(local_nr_variance_dict["window_radius_vec"], 
local_nr_variance_dict["local_nr_variance_vec"] ./ local_nr_variance_dict["window_radius_vec"].^3)

NG.plot_network(graph_dict)