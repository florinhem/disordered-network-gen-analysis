
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")    #*#

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import .Threads
import Plots
Plots.pyplot()

function save_multiple_N_T_trials_beta_gml(
    ;
    network_type_array,
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature_array,
    quench,
    save_path,
    filename_start
    )
    
    println(Threads.nthreads())

    Iter=collect(Iterators.product(
        network_type_array,
        nr_vertices_array,
        maximal_temperature_array,
        bond_bending_const_array,
        temperature_gradient_array,
        nr_monte_carlo_steps_per_temperature_array,
        theta_ground_state_array,
        nr_trials_per_temperature_array
        ))

    Threads.@threads for (
        network_type,
        nr_vertices,
        maximal_temperature,
        bond_bending_const,
        temperature_gradient,
        nr_monte_carlo_steps_per_temperature,
        theta_ground_state,
        trial) in Iter
                
        println(
        "$network_type"*", "*
        "$nr_vertices"*", "*
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

        spatial_network = NG.get_periodic_network(evolution_dict)

        #=
        plot1=NG.plot_spatial_network_2(spatial_network)
        Plots.xlabel!("x")
        Plots.ylabel!("y")
        Plots.zlabel!("z")
        display(plot1)
        =#

        println("sigma_L, $((NA.get_bond_length_std(spatial_network))[1])")
        println("sigma_A, $((NA.get_bond_angle_std(spatial_network))[1])")
    
        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = 
            NA.get_temperature_sequence_heating_cooling_gradient(
                maximal_temperature;
                temperature_gradient = temperature_gradient, 
                nr_monte_carlo_steps_per_temperature = nr_monte_carlo_steps_per_temperature,
                quench = quench) 

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

        #=
        plot1=NG.plot_spatial_network_2(spatial_network)
        Plots.xlabel!("x")
        Plots.ylabel!("y")
        Plots.zlabel!("z")
        display(plot1)
        =#

        #break

        println("nbr acc moves, $(length(move_accepted_vec)), $(sum(move_accepted_vec))")
        println("sigma_L, $((NA.get_bond_length_std(spatial_network))[1])")
        println("sigma_A, $((NA.get_bond_angle_std(spatial_network))[1])")
        println("evolve_network_temperature_sequence end")

        evolution_dict["total_energy_vec"] = total_energy_vec
        evolution_dict["move_accepted_vec"] = move_accepted_vec

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
                    
    end
end

save_multiple_N_T_trials_beta_gml(;
    network_type_array=["dia"],
    nr_vertices_array=[8*3^3],
    maximal_temperature_array=[0.5,1.0,1.5],
    bond_bending_const_array=[0.0,1.0],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_5"
)

#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["pto"],
    nr_vertices_array=[14*2^3],
    maximal_temperature_array=[0.35],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_4"
)
    =#
#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["dia"],
    nr_vertices_array=[8*3^3],
    maximal_temperature_array=[0.44],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_4"
)
    =#
#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["srd"],
    nr_vertices_array=[10*3^3],
    maximal_temperature_array=[0.29],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_4"
)=#

#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["srd"],
    nr_vertices_array=[10*3^3],
    maximal_temperature_array=[0.25],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_4"
)=#

#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["srd"],
    nr_vertices_array=[10*3^3],
    maximal_temperature_array=[0.1],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_3"
)=#

#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["srs"],
    nr_vertices_array=[8*3^3],
    maximal_temperature_array=[0.19],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_3"
)
    =#

#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["srs"],
    nr_vertices_array=[8*3^3],
    maximal_temperature_array=[0.21],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_3"
)
    =#

#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["dia"],
    nr_vertices_array=[8*3^3],
    maximal_temperature_array=[0.5],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_3"
)
    =#

#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["srd"],
    nr_vertices_array=[10*3^3],
    maximal_temperature_array=[0.5],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_3"
)
    =#



#=

save_multiple_N_T_trials_beta_gml(;
    network_type_array=["ctn"],
    nr_vertices_array=[28*3^3],
    maximal_temperature_array=[1.0],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_3"
)

=#
#=

save_multiple_N_T_trials_beta_gml(;
    network_type_array=["pto"],
    nr_vertices_array=[14*2^3],
    maximal_temperature_array=[0.07, 0.14,0.28,0.56],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_2"
)
    =#


#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["srs"],
    nr_vertices_array=[8*3^3],
    maximal_temperature_array=[0.17],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",     
    filename_start="mts_2"
)
    =#







#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["pto"],
    nr_vertices_array=[14*3^3],
    maximal_temperature_array=[0.45],
    bond_bending_const_array=[0.1],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/ft_1/",     
    filename_start="test_1"
)
=#
#=
save_multiple_N_T_trials_beta_gml(;
    network_type_array=["srd"],#["pto", "lcs"],
    nr_vertices_array=[18],#[216],
    maximal_temperature_array=[0.0001],
    bond_bending_const_array=[0.1],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/ft_1/",     
    filename_start="test_1"
)

save_multiple_N_T_trials_beta_gml(;
    network_type_array=["dia", "srs", "srd", "ctn"],
    nr_vertices_array=[216],
    maximal_temperature_array=[0.1,0.15,0.2],
    bond_bending_const_array=[0.1,0.05,0.025],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    quench=true,
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/ft_1/",     
    filename_start="ft_1"
)
=#