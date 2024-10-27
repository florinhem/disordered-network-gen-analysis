
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
import Colors

function scatter_plot_for_mulitple_gml(;
    nr_vertices_array,
    maximal_temperature_array,
    nr_trials_per_temperature,
    bond_bending_const_array,
    temperature_gradient,
    nr_monte_carlo_steps_per_temperature,
    theta_ground_state_array,
    save_path,
    filename_start,
    plot_save_path,
    plot_filename_start,
    shape_array)

    # test before we begin
    @assert length(nr_vertices_array)>=1
    @assert length(maximal_temperature_array)>=1
    @assert nr_trials_per_temperature>=1
    @assert length(maximal_temperature_array)>=1
    @assert length(theta_ground_state_array)===length(shape_array)

    # generate color for the different temperatures
    if length(maximal_temperature_array)===1
        color_array=[:blue]
    else 
        color_array=range(Colors.colorant"blue", stop=Colors.colorant"red", length=length(maximal_temperature_array))
        println(color_array)
    end
    
    @assert length(maximal_temperature_array)===length(color_array)

    P=Plots.scatter(
        title="Bond length std vs Keating energy per vertex"
    )

    for k in eachindex(nr_vertices_array)

        nr_vertices=nr_vertices_array[k]
        #shape=shape_array[k]

        for j in eachindex(maximal_temperature_array)

            maximal_temperature=maximal_temperature_array[j]
            color=color_array[j]

            for i in 1:nr_trials_per_temperature

                for m in eachindex(bond_bending_const_array)

                    bond_bending_const=bond_bending_const_array[m]
                       
                    for n in eachindex(theta_ground_state_array)

                        theta_ground_state=theta_ground_state_array[n]
			            shape=shape_array[n]

                        println("$nr_vertices"*", "*"$maximal_temperature"*", "*"$i"*", "*"$bond_bending_const"*", "*"$theta_ground_state" )

                
                        filename = (filename_start
                            *"_N="*"$nr_vertices"
                            *"_T="*"$maximal_temperature"
                            *"_Trial="*"$i"
                            *"_Beta="*"$bond_bending_const"
                            *"_Theta_GS="*"$theta_ground_state"
                            *"_GradT="*"$temperature_gradient"
                            *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
                            *".gml"
                            )

                        total_path=save_path*filename

                        spatial_network=NG.load_spatial_network_from_gml(total_path)

                        #bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)
                        bond_angle_std, bond_angle_vec = NA.get_bond_angle_std(spatial_network)


                        total_energy_keating=NG.get_total_energy_keating(spatial_network)
                        energy_keating_per_vertex=total_energy_keating/spatial_network[]["nr_vertices"]

                        P=Plots.scatter!(
                            P,
                            [energy_keating_per_vertex],
                            #[bond_length_std],
                            [bond_angle_std],
                            xlabel="Keating energy per vertex",
                            #ylabel="Bond length std",
                            ylabel="Bond angle std",
                            markercolor=color,
                            markershape=shape,
                            label=false,
                            cbar=true,
                            show=true)
                    end
                end
            end
        end
    end

    minimum_nr_vertices=minimum(nr_vertices_array)
    maximum_nr_vertices=maximum(nr_vertices_array)

    minimum_temperature=minimum(maximal_temperature_array)
    maximum_temperature=maximum(maximal_temperature_array)
    
    minimum_bond_bending_const=minimum(bond_bending_const_array)
    maximum_bond_bending_const=maximum(bond_bending_const_array)

    minimum_theta=minimum(theta_ground_state_array)
    maximum_theta=maximum(theta_ground_state_array)


    plot_filename = (plot_filename_start
        *"_N="*"$minimum_nr_vertices" * "-" * "$maximum_nr_vertices"
        *"_T="*"$minimum_temperature" * "-" * "$maximum_temperature"
        *"_Trials="*"$nr_trials_per_temperature"
        *"_Beta="*"$minimum_bond_bending_const" * "-" * "$maximum_bond_bending_const"
        *"_Theta_GS="*"$minimum_theta" * "-" * "$maximum_theta"
        *"_GradT="*"$temperature_gradient"
        *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
        *".png")
    plot_total_path=plot_save_path*plot_filename

    Plots.savefig(P,plot_total_path)
end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216,512],
    maximal_temperature_array=[0.1,0.135,0.17,0.205,0.24],
    nr_trials_per_temperature=1,
    bond_bending_const_array=[0.285],
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state_array=[110.0,180.0],
    #save_path = "/home/glauserv/Documents/GitLinux/GitF/code_photonic_structures/simulations/multiple_parameters/",      
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "multiple_p_quench_false_theta_array",    #"multiple_p",
    #plot_save_path = "/home/glauserv/Documents/GitLinux/GitF/code_photonic_structures/simulations/analysis_plot/",      
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "multiple_p_angle",        #"multiple_p",
    shape_array=[:circle,:rect]
)
