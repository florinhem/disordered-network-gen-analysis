# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")    #*#

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
Plots.plotlyjs()
import .Threads
import Statistics
import LinearAlgebra

function save_multiple_N_T_trials_beta_gml(
    ;
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature_array,
    network_type,
    save_path,
    filename_start
    )
    
    println(Threads.nthreads())


    nr_vertices=nr_vertices_array[1]
    maximal_temperature=maximal_temperature_array[1]
    bond_bending_const=bond_bending_const_array[1]
    temperature_gradient=temperature_gradient_array[1]
    nr_monte_carlo_steps_per_temperature=nr_monte_carlo_steps_per_temperature_array[1]
    theta_ground_state=theta_ground_state_array[1]
    trial=nr_trials_per_temperature_array[1]

    println("$nr_vertices"*", "*
    "$maximal_temperature"*", "*
    "$bond_bending_const"*", "*
    "$temperature_gradient"*", "*
    "$nr_monte_carlo_steps_per_temperature"*", "*
    "$theta_ground_state"*", "*
    "$trial" )

    evolution_dict = NA.get_evolution_dict(;
        nr_vertices = nr_vertices, 
        network_type=network_type, 
        bond_bending_const=bond_bending_const, 
        min_ring_size=3,
        theta_ground_state=theta_ground_state
    )

    # ==== Ordered 

    spatial_network = NG.get_periodic_network(evolution_dict)

    filename = (filename_start
        *"_NW="*"$network_type"
        *"_N="*"$nr_vertices"
        *"_T="*"$maximal_temperature"
        *"_Beta="*"$bond_bending_const"
        *"_GradT="*"$temperature_gradient"
        *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
        *"_Theta_GS="*"$theta_ground_state"
        *"_Trial="*"$trial"
        )

    NG.save_spatial_network_to_gml(
        spatial_network,
        filename;
        evolution_dict = evolution_dict,
        save_path = save_path)

    NG.get_spatial_network_for_simulation!(
        spatial_network,
        vector_out_of_supercell_length = 1,
        duplicate_bonds_close_to_supercell_edge = true,
        save_result = true,
        filename = filename*"_o",
        save_path = save_path)

    #=
    plot1=NG.plot_spatial_network_2(spatial_network)
    Plots.xlabel!("x")
    Plots.ylabel!("y")
    Plots.zlabel!("z")
    display(plot1)

    println("sigma_L, $((NA.get_bond_length_std(spatial_network))[1])")
    println("sigma_A, $((NA.get_bond_angle_std(spatial_network))[1])")
    =#

    # ==== Disordered

    # force one bond switch and then save the network

    #=
    switched_chain = NG.get_random_chain(
        spatial_network; 
        declined_chains = [],
        remaining_chains = [],
        min_ring_size = evolution_dict["min_ring_size"])

    println("switched_chain", switched_chain)
    =#

    spatial_network = NG.switch_chain!(spatial_network, (3, 11, 27, 7))
    
    #spatial_network = NG.switch_chain!(spatial_network, switched_chain)
    #=
    spatial_network, move_accepted = NG.monte_carlo_move!(
        spatial_network, evolution_dict, maximal_temperature; switched_chain)

    move_accepted_vec=[move_accepted]
        =#

    total_energy_vec=NG.get_total_energy_keating(spatial_network)
    
    #=
    plot1=NG.plot_spatial_network_2(spatial_network)
    Plots.xlabel!("x")
    Plots.ylabel!("y")
    Plots.zlabel!("z")
    display(plot1)

    println("sigma_L, $((NA.get_bond_length_std(spatial_network))[1])")
    println("sigma_A, $((NA.get_bond_angle_std(spatial_network))[1])")
        =#

    evolution_dict["total_energy_vec"] = total_energy_vec
    #evolution_dict["move_accepted_vec"] = move_accepted_vec

    
	
    
    NG.save_spatial_network_to_gml(
        spatial_network,
        filename;
        evolution_dict = evolution_dict,
        save_path = save_path)

    NG.get_spatial_network_for_simulation!(
        spatial_network,
        vector_out_of_supercell_length = 1,
        duplicate_bonds_close_to_supercell_edge = true,
        save_result = true,
        filename = filename*"_d",
        save_path = save_path)
    
    
    return
end


save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array=[28*1^3],
    maximal_temperature_array=[10.0],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    network_type="ctn",
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/disordered_pto_ctn/",
    filename_start="MC=1_R=No_Q=No_1"
)