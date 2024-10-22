
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
import Base.Iterator

function get_bond_length_AND_keating_energy_per_vertex_after_evolve_network(
    ;maximal_temperature=1)

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216, network_type="diamond", bond_bending_const=0.285, min_ring_size=3)
    spatial_network = NG.get_periodic_network(evolution_dict)
    NG.plot_spatial_network(spatial_network)

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = 
        NA.get_temperature_sequence_heating_cooling_gradient(
            maximal_temperature;
            temperature_gradient = 0.5, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = false)

    evolution_dict["temperature_vec"] = temperature_vec
    evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = nr_monte_carlo_steps_per_temperature_vec

    spatial_network, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(
        spatial_network,
        evolution_dict;
        print_progress = true,
        print_every_nr_attempted_bond_switches = 10)

    #NG.plot_spatial_network(spatial_network)

    #=
    save_path = raw".\my_networks\\"
    filename = "Test3"

    NG.save_spatial_network_to_gml(
        spatial_network,
        filename;
        evolution_dict = evolution_dict,
        save_path = save_path)
    =#

    #Now do here a standard deviation of bond
    bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)

    #Get the Keating energy
    total_energy_keating=NG.get_total_energy_keating(spatial_network)
    energy_keating_per_vertex=total_energy_keating/evolution_dict["nr_vertices"]

    return bond_length_std,energy_keating_per_vertex
end


x=[]
y=[]

nr_trials_per_temperature=4
temperature_range=0:1:2
temperature_array=collect(Base.Iterator.take(Base.Iterator.cycle(temperature_range),nr_trials_per_temperature*length(temperature_range)))


for i in eachindex(temperature_array)

    maximal_temperature=temperature_array[i]
    bond_length_std,energy_keating_per_vertex=get_bond_length_AND_keating_energy_per_vertex_after_evolve_network(
        ;maximal_temperature=maximal_temperature)
    push!(x,energy_keating_per_vertex)
    push!(y,bond_length_std)
end



#Create a point in the graph
#x=[0,energy_keating_per_vertex]
#y=[0,bond_length_std]

#Plot where we are in the graph

P=Plots.scatter(x, y,
    xlabel="Keating energy per vertex", 
    ylabel="Bond length std",  
    label="",
    marker_z=temperature_array,
    color = :jet,
    colorbar_title="Maximal temperature",
    flip_axis=false
    )

#savefig(P,raw".\my_networks\KE_VS_BLSTD\KE_VS_BLSTD_3_N=216_beta=0_285_T=0-1-2_trials=4.png")

#Do this all in a function => iterate over temperature (color?) and over beta (what we want to know) and stop at different Keating energies(?)