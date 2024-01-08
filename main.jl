
#include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

#import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import Plots
import Measurements

evolution_dict = NA.get_evolution_dict(;nr_vertices = 512 ,temperature_vec = [0.2,0],
nr_monte_carlo_steps_per_temperature_vec = [2,20])

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 100)

Plots.plot(collect(1:length(total_energy_vec)), total_energy_vec)