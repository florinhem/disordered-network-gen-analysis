
#include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

#import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 ,temperature_vec = [1,0],
nr_monte_carlo_steps_per_temperature_vec = [1,30])

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200)

Plots.plot(collect(1:length(total_energy_vec)), total_energy_vec)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

filename = "1000_vertices_T_1_quenched"

save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\\"

NG.save_mesh_from_network(graph_dict_to_save, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict_to_save,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)