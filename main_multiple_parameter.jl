
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
import .Threads

function save_multiple_N_T_trials_beta_gml(
    ;
    nr_vertices_array,
    maximal_temperature_array,
    nr_trials_per_temperature,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    save_path,
    filename_start
    )
    
    println(Threads.nthreads())

    Threads.@threads for (nr_vertices,maximal_temperature,i,bond_bending_const,temperature_gradient,nr_monte_carlo_steps_per_temperature,theta_ground_state) in collect(
        Iterators.product(nr_vertices_array,maximal_temperature_array,1:nr_trials_per_temperature,bond_bending_const_array,temperature_gradient_array,nr_monte_carlo_steps_per_temperature_array,theta_ground_state_array))
                
        println("$nr_vertices"*", "*
		"$maximal_temperature"*", "*
		"$i"*", "*
		"$bond_bending_const"*", "*
		"$temperature_gradient"*", "*
                "$nr_monte_carlo_steps_per_temperature"*", "*
                "$theta_ground_state" )
    
        evolution_dict = NA.get_evolution_dict(;
            nr_vertices = nr_vertices, 
            network_type="diamond", 
            bond_bending_const=bond_bending_const, 
            min_ring_size=3,
            theta_ground_state=theta_ground_state
            )
        spatial_network = NG.get_periodic_network(evolution_dict)
    
        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = 
            NA.get_temperature_sequence_heating_cooling_gradient(
                maximal_temperature;
                temperature_gradient = temperature_gradient, 
                nr_monte_carlo_steps_per_temperature = nr_monte_carlo_steps_per_temperature,
                quench = false)

        evolution_dict["temperature_vec"] = temperature_vec
        evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = nr_monte_carlo_steps_per_temperature_vec

        total_energy_vec::Vector{Float64}=[]
        move_accepted_vec::Vector{Bool}=[]

        spatial_network, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(
            spatial_network,
            evolution_dict;
            total_energy_vec = total_energy_vec,
            move_accepted_vec= move_accepted_vec,
            print_progress = true,
            print_every_nr_attempted_bond_switches = 1000)

        evolution_dict["total_energy_vec"] = total_energy_vec
        evolution_dict["move_accepted_vec"] = move_accepted_vec

        filename = (filename_start
            *"_N="*"$nr_vertices"
            *"_T="*"$maximal_temperature"
            *"_Trial="*"$i"
            *"_Beta="*"$bond_bending_const"
            *"_Theta_GS="*"$theta_ground_state"
            *"_GradT="*"$temperature_gradient"
            *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
            )
	
        NG.save_spatial_network_to_gml(
            spatial_network,
            filename;
            evolution_dict = evolution_dict,
            save_path = save_path)
                    
    end
end

try
    save_multiple_N_T_trials_beta_gml(;
        nr_vertices_array=[216],
        maximal_temperature_array=[0.1],
        nr_trials_per_temperature=1,
        bond_bending_const_array=[0.21,0.285,0.36],
        temperature_gradient_array=[0.025,0.05,0.1,0.2,0.4],
        nr_monte_carlo_steps_per_temperature_array=[0.0025,0.005,0.01,0.02,0.04],
        theta_ground_state_array=[110.0,180.0],
        save_path ="/home/glauserv/Documents/GitLinux/GitF/code_photonic_structures/simulations/multiple_parameters/",
        filename_start="multiple_BTMC_"
    )
catch e
    error_msg = sprint(showerror, e)
    st = sprint((io,v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
    @warn "Trouble doing things:\n$(error_msg)\n$(st)"
    println("Trouble doing things:\n$(error_msg)\n$(st)")
end
