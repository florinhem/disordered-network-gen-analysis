
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots

function plot_keating_energy()

    # plot the theoretical energy functions. E_harmonic is the Taylor of E_str
    r_norm=collect(0:0.05:1.5)
    E_str=3/16 * (r_norm.^2 .-1).^2
    E_harmonic=3/16 * 4 *(r_norm .-1).^2

    Plots.plot(r_norm,E_str,label="E_str")
    Plots.plot!(r_norm,E_harmonic,label="E_harmonic")

    
    # create the network
    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216, network_type="diamond", bond_bending_const=0.285, min_ring_size=3)
    spatial_network = NG.get_periodic_network(evolution_dict)


    # heat and cool down
    maximal_temperature=0.1

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = 
        NA.get_temperature_sequence_heating_cooling_gradient(
            maximal_temperature;
            temperature_gradient = 10.0, 
            nr_monte_carlo_steps_per_temperature = 2/(18*216),
            quench = false)

    evolution_dict["temperature_vec"] = temperature_vec
    evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = nr_monte_carlo_steps_per_temperature_vec

    spatial_network, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(
        spatial_network,
        evolution_dict;
        print_progress = true,
        print_every_nr_attempted_bond_switches = 10)

    println("Heat cool down finished, now plotting")


    # prepare for scatter plotting
    r=[]
    E=[]

    for bond in MetaGraphsNext.edge_labels(spatial_network)
        append!(r,sqrt(spatial_network[bond...]["distance_squared"]))
        append!(E,NG.local_bond_stretching_energy_keating(spatial_network, bond))
    end

    Plots.scatter!(r,E,label="Calculated points")

    println(r)
    println(E)

    Plots.savefig(raw"./my_networks/E_str/E_str_20_maxT="
        *"$(maximal_temperature)"
        *".png")

end

plot_keating_energy()
