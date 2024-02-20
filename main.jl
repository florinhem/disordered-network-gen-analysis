
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU
import Plots

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural


evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [1, 0],
nr_monte_carlo_steps_per_temperature_vec = [1, 50], min_ring_size = 4)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "216_vertices_T_1_quenched"

save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\with_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 512 ,temperature_vec = [1, 0],
nr_monte_carlo_steps_per_temperature_vec = [1, 50], min_ring_size = 4)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "512_vertices_T_1_quenched"

save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\with_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)

