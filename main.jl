
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

function save_multiple_N_T_trials_beta_gml(
    ;
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature_array,
    edge_length_unit_cell_array,
    network_type,
    quench,
    save_path,
    filename_start
    )
    
    println(Threads.nthreads())

    Iter=collect(Iterators.product(
        nr_vertices_array,
        maximal_temperature_array,
        bond_bending_const_array,
        temperature_gradient_array,
        nr_monte_carlo_steps_per_temperature_array,
        theta_ground_state_array,
        nr_trials_per_temperature_array,
        edge_length_unit_cell_array))

    Threads.@threads for (
        nr_vertices,
        maximal_temperature,
        bond_bending_const,
        temperature_gradient,
        nr_monte_carlo_steps_per_temperature,
        theta_ground_state,
        trial,
        edge_length_unit_cell) in Iter
                
        println("$nr_vertices"*", "*
		"$maximal_temperature"*", "*
		"$bond_bending_const"*", "*
		"$temperature_gradient"*", "*
        "$nr_monte_carlo_steps_per_temperature"*", "*
        "$theta_ground_state"*", "*
        "$trial"*", "*
        "$edge_length_unit_cell" )

        #edges=NG.get_edges_ctn(edge_length_unit_cell)

        #break
    
        evolution_dict = NA.get_evolution_dict(;
            nr_vertices = nr_vertices, 
            network_type=network_type, 
            bond_bending_const=bond_bending_const, 
            min_ring_size=3,
            theta_ground_state=theta_ground_state
            )
        spatial_network = NG.get_periodic_network(
            evolution_dict,
            edge_length_unit_cell)

        println("sigma_L, $((NA.get_bond_length_std(spatial_network))[1])")
        println("sigma_A, $((NA.get_bond_angle_std(spatial_network))[1])")

        nr_accepted_moves=0
        nr_attempted_moves=100

        spatial_network_copy=deepcopy(spatial_network)
        evolution_dict_copy=deepcopy(evolution_dict)

        E=NG.get_total_energy_keating(spatial_network_copy)
        println("E, $E")
        
        for i in 1:nr_attempted_moves
            spatial_network=deepcopy(spatial_network_copy)
            evolution_dict=deepcopy(evolution_dict_copy)
            
            move=NG.monte_carlo_move!(
                spatial_network, 
                evolution_dict,
                maximal_temperature)

            #println("move[2], $(move[2])")

            if move[2]
                nr_accepted_moves+=1
            end
        end
        println(nr_accepted_moves)
        println("l=$(edge_length_unit_cell), Ratio: $(nr_accepted_moves/nr_attempted_moves)")
                    
    end
end

#network_type="dia"     #diamond
#network_type="srs"      #gyroid
network_type="ctn"
network_type="srd"

save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array=[216],
    maximal_temperature_array=[0.1],
    bond_bending_const_array=[0.1125],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    edge_length_unit_cell_array=[0],
    network_type=network_type,
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/networks_dia_srd_ctn/",     
    filename_start="m_$(network_type)_3"
)