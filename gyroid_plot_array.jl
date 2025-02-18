
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
import DataFrames
import LaTeXStrings
import IterTools
import Combinatorics
import Statistics
using StatsPlots

function scatter_plot_for_mulitple_gml(;
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature_array,
    save_path,
    filename_start,
    plot_save_path,
    plot_filename_start)

    # test before we begin
    @assert length(nr_vertices_array)>=1
    @assert length(maximal_temperature_array)>=1
    @assert length(bond_bending_const_array)>=1
    @assert length(temperature_gradient_array)>=1
    @assert length(nr_monte_carlo_steps_per_temperature_array)>=1
    @assert length(theta_ground_state_array)>=1
    @assert length(nr_trials_per_temperature_array)>=1

    path_array=Glob.glob(filename_start*"*",save_path)

    for k in eachindex(nr_vertices_array)

        nr_vertices=nr_vertices_array[k]

        for j in eachindex(maximal_temperature_array)

            maximal_temperature=maximal_temperature_array[j]

            for m in eachindex(bond_bending_const_array)

                bond_bending_const=bond_bending_const_array[m]

                for n in eachindex(temperature_gradient_array)

                    temperature_gradient=temperature_gradient_array[n]
                    
                    for o in eachindex(nr_monte_carlo_steps_per_temperature_array)

                        nr_monte_carlo_steps_per_temperature=nr_monte_carlo_steps_per_temperature_array[o]
                    
                        for p in eachindex(theta_ground_state_array)
                            
                            theta_ground_state=theta_ground_state_array[p]

                            for i in eachindex(nr_trials_per_temperature_array)

                                trial=nr_trials_per_temperature_array[i]
                                
                                filename = (filename_start
                                    *"_N="*"$nr_vertices"
                                    *"_T="*"$maximal_temperature"
                                    *"_Beta="*"$bond_bending_const"
                                    *"_GradT="*"$temperature_gradient"
                                    *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
                                    *"_Theta_GS="*"$theta_ground_state"
                                    *"_Trial="*"$trial"
                                    )

                                total_path=save_path*filename

                                if(total_path*".gml" in path_array)
                                    println(filename)

                                    spatial_network=NG.load_spatial_network_from_gml(total_path*".gml")
                                    println(spatial_network)

                                    plot1=NG.plot_spatial_network_2(spatial_network)
                                    display(plot1)

                                else
                                    println("file not in directory")
                                end
                            end
                        end
                    end
                end
            end
        end
    end

   
end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[8*8],          
    maximal_temperature_array=[0.1],
    bond_bending_const_array=[0.15],
    temperature_gradient_array=[1],     
    nr_monte_carlo_steps_per_temperature_array=[0.1],    
    theta_ground_state_array=[120.0],
    nr_trials_per_temperature_array=[1],
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/multiple_parameters/",     
    filename_start="m_a_g_1",  
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_a_g_1"
)
