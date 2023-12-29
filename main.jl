
#include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

#import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import Plots
import Measurements

evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [5, 0.1],
nr_monte_carlo_steps_per_temperature_vec = [5,5])

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence(graph_dict,
        evolution_dict;
        print_progress = true,
        random_evolution_seed = 2)

NG.plot_network(graph_dict)