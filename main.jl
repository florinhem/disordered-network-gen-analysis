
#include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

#import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import Plots
import Measurements

nr_graphs = 3

y_arr = Array{Real}(undef, 28, nr_graphs)

for i in 1:nr_graphs

    graph_dict = NG.get_periodic_network( ; nr_vertices = 1000 , 
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
        print_progress = false,
        random_evolution_seed = i,
        thermal_fluctuations = false)

    wavenumber_vec, structure_factor_vec = NA.get_structure_factor_isotrope_by_wavenumber_vec(graph_dict)

    y_arr[:,i] = structure_factor_vec

    println(string(i))

end

structure_factor_average_vec = NG.get_multiple_graph_average_vec(
    y_arr)

Plots.plot(wavenumber_vec, Measurements.value.(structure_factor_average_vec), ribbon = Measurements.uncertainty.(structure_factor_average_vec) )