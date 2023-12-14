"""
In this code analyses whether a 3d binary structure is hyperuniform.
This is done by calculating the local volume fraction variance as a function
measuring window size.

The calculations are taken from
10.1016/j.physrep.2018.03.001, pages 13, 18, 30
10.1103/PhysRevE.96.032909

"""


#include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

#import my module that contains all functions for the analysis of binary structure data
import .NetworkGeneration as NG

temperature = 1

graph_dict = NG.get_periodic_network( ; nr_atoms = 128 , 
        nr_dimensions = 3, 
        network_type = "diamond",
        temperature = temperature,
        bond_bending_const = 0.285,
        thermal_fluctuations= true)

shell_nr = 4

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network(graph_dict,
    256, 
    temperature; 
    nr_max_relaxation_cycles = 25,
        break_at_relative_cluster_energy_change = 0.001,
        reject_during_relaxation_cycle_threshold = 10,
        relax_efficiently = true,
        shell_nr = shell_nr,
    print_progress = true,
    random_evolution_seed = 3,
    thermal_fluctuations = true)


NG.plot_network(graph_dict)

#plot_dict = NG.cut_bonds_out_of_supercell!(graph_dict)
#NG.save_mesh_from_network(plot_dict, "diamond_small_t20")
