
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU
import Plots
import LaTeXStrings as Latex # to display latex symbols in plot labels
import NaNStatistics

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007

# 


temperatures = [0.3, 0.4, 0.5]

for temperature in temperatures

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_cooling_gradient(temperature;
    temperature_decrease_per_monte_carlo_step = 0.1,
    quench = true )

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_cool_0.1_per_mc_quenched"

    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
    temperature_increase_per_monte_carlo_step = 0.1,
    quench = true )

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_cool_0.1_per_mc_quenched"

    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end