
#include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

#import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU
import Plots

#possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr n


#evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [0.05],
#nr_monte_carlo_steps_per_temperature_vec = [1])
#
#graph_dict = NG.get_periodic_network(evolution_dict)

evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = 3

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
        move_accepted_vec=move_accepted_vec,
        total_energy_vec=total_energy_vec,
    print_progress = true,
    print_every_nr_attempted_bond_switches = 20)

Plots.plot(collect(1:length(total_energy_vec)), total_energy_vec)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "64_vertices_T_0.05"

save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)
        