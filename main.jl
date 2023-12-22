
#include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

#import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import Plots

graph_dict = NG.get_periodic_network( ; nr_vertices = 500 , 
                            nr_dimensions = 3, 
                            network_type = "diamond")

temperature = 10
shell_nr = 4

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network(graph_dict,
    20, 
    temperature; 
    nr_max_relaxation_cycles = 25,
        break_at_relative_cluster_energy_change = 0.0001,
        reject_during_relaxation_cycle_threshold = 10,
        relax_efficiently = true,
        shell_nr = shell_nr,
    print_progress = true,
    random_evolution_seed = 3,
    thermal_fluctuations = false)

wavenumber_vec, structure_factor_vec = NA.get_structure_factor_isotrope_by_wavenumber_vec(graph_dict)

Plots.plot(wavenumber_vec, structure_factor_vec)