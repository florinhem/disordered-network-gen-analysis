
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
import Glob

function scatter_plot_for_mulitple_gml(;
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature,
    save_path,
    filename_start,
    plot_save_path,
    plot_filename_start,
    markershape_array,
    markersize_array,
    markerstrokewidth_array)

    # test before we begin
    @assert length(nr_vertices_array)>=1
    @assert length(maximal_temperature_array)>=1
    @assert length(bond_bending_const_array)>=1
    @assert length(temperature_gradient_array)>=1
    @assert length(nr_monte_carlo_steps_per_temperature_array)>=1
    @assert length(theta_ground_state_array)>=1
    @assert nr_trials_per_temperature>=1
    @assert length(markershape_array)>=1

    @assert length(theta_ground_state_array)===length(markershape_array)


    # generate color for the different temperatures
    if length(temperature_gradient_array)===1
        markercolor_array=[:blue]
    else 
        markercolor_array=range(Colors.colorant"blue", stop=Colors.colorant"red", length=length(temperature_gradient_array))
    end
    
    @assert length(temperature_gradient_array)===length(markercolor_array)

    # store array with all the paths in the directory to check if we have
    # this .h5 and .gml to be able to plot it.
    path_array=Glob.glob(filename_start*"*.gml",save_path)

    P=Plots.scatter(
        title="Bond length std vs Bond angle std",
        xlabel="Bond length std",
        ylabel="Bond angle std"
    )

    for k in eachindex(nr_vertices_array)

        nr_vertices=nr_vertices_array[k]

        for j in eachindex(maximal_temperature_array)

            maximal_temperature=maximal_temperature_array[j]

            for m in eachindex(bond_bending_const_array)

                bond_bending_const=bond_bending_const_array[m]
                markerstrokewidth=markerstrokewidth_array[m]

                for n in eachindex(temperature_gradient_array)

                    temperature_gradient=temperature_gradient_array[n]
                    markercolor=markercolor_array[n]
                    
                    for o in eachindex(nr_monte_carlo_steps_per_temperature_array)

                        nr_monte_carlo_steps_per_temperature=nr_monte_carlo_steps_per_temperature_array[o]
                        markersize=markersize_array[o]
                    
                        for p in eachindex(theta_ground_state_array)
                            
                            theta_ground_state=theta_ground_state_array[p]
                            markershape=markershape_array[p]

                            for i in 1:nr_trials_per_temperature
                                
                                #=
                                println("$nr_vertices"*", "*
                                    "$maximal_temperature"*", "*
                                    "$i"*", "*
                                    "$bond_bending_const"*", "*
                                    "$temperature_gradient"*", "*
                                    "$nr_monte_carlo_steps_per_temperature"*", "*
                                    "$theta_ground_state" )
                                =#
                        
                                filename = (filename_start
                                    *"_N="*"$nr_vertices"
                                    *"_T="*"$maximal_temperature"
                                    *"_Trial="*"$i"
                                    # TODO Take trial after theta_GS
                                    *"_Beta="*"$bond_bending_const"
                                    *"_GradT="*"$temperature_gradient"
                                    *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
                                    *"_Theta_GS="*"$theta_ground_state"
                                    *".gml"
                                    )

                                total_path=save_path*filename

                                if(total_path in path_array)
                                    #println("scatter done")

                                    spatial_network=NG.load_spatial_network_from_gml(total_path)

                                    bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)
                                    bond_angle_std, bond_angle_vec = NA.get_bond_angle_std(spatial_network)

                                    println("["
                                        #*"$nr_vertices"*","
                                        #*"$maximal_temperature"*","
                                        *"$bond_bending_const"*","
                                        *"$temperature_gradient"*","
                                        *"$nr_monte_carlo_steps_per_temperature"*","
                                        *"$theta_ground_state"*","
                                        #*"$i"*","
                                        *"$bond_length_std"*","
                                        *"$bond_angle_std"
                                        *"],")

                                    Plots.scatter!(
                                        P,
                                        [bond_length_std],
                                        [bond_angle_std],
                                        markercolor=markercolor,
                                        markershape=markershape,
                                        markersize=markersize,
                                        markerstrokewidth=markerstrokewidth,
                                        legend=false,
                                        cbar=true,
                                        show=true)
                                else
                                    #println("file not in directory")
                                end
                            end
                        end
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

    minimum_temperature_gradient=minimum(temperature_gradient_array)
    maximum_temperature_gradient=maximum(temperature_gradient_array)

    minimum_nr_monte_carlo_steps_per_temperature=minimum(nr_monte_carlo_steps_per_temperature_array)
    maximum_nr_monte_carlo_steps_per_temperature=maximum(nr_monte_carlo_steps_per_temperature_array)

    minimum_theta=minimum(theta_ground_state_array)
    maximum_theta=maximum(theta_ground_state_array)


    plot_filename = (plot_filename_start
        *"_N="*"$minimum_nr_vertices" * "-" * "$maximum_nr_vertices"
        *"_T="*"$minimum_temperature" * "-" * "$maximum_temperature"
        *"_Beta="*"$minimum_bond_bending_const" * "-" * "$maximum_bond_bending_const"
        *"_GradT="*"$minimum_temperature_gradient" * "-" * "$maximum_temperature_gradient"
        *"_StepsPerT="*"$minimum_nr_monte_carlo_steps_per_temperature" * "-" * "$maximum_nr_monte_carlo_steps_per_temperature"
        *"_Theta_GS="*"$minimum_theta" * "-" * "$maximum_theta"
        *"_Trials="*"$nr_trials_per_temperature"
        *".png")
    plot_total_path=plot_save_path*plot_filename

    Plots.savefig(P,plot_total_path)
end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216],
    maximal_temperature_array=[0.1],
    nr_trials_per_temperature=1,
    bond_bending_const_array=[0.135,0.21,0.285,0.36,0.435],                #[0.135,0.21,0.285,0.36,0.435],
    temperature_gradient_array=[0.001,0.01,0.1,1,10],                      #[0.001,0.01,0.1,1,10],
    nr_monte_carlo_steps_per_temperature_array=[0.0001,0.001,0.01,0.1,1],  #[0.0001,0.001,0.01,0.1,1],
    theta_ground_state_array=[110.0,180.0],   
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "multiple_BTMC",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "multiple_BTMC_g_1",
    markershape_array=[:circle,:rect],                              #[:circle,:rect],       #[:circle,:rect#=,:star5,:cross,:+=#],
    markersize_array=1.5 .*[1,2,3,4,5],                             #1.5 .*[1,2,3,4,5],
    markerstrokewidth_array=[1,2,3,4,5]                             #[1,2,3,4,5]
)
