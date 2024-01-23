
#include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

#import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU
import Plots

#possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr n

nr_samples = 5

evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 )

evolution_dict["temperature_vec"] = zeros(nr_samples) .+ 0.125
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = ones(Int64, nr_samples)
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][1] = 3 

graph_dict = NG.get_periodic_network(evolution_dict)


graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200,
    save_network_after_each_step = true,
    filename = "1000_vertices_T_0.125",)