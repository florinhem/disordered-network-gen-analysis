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
using StatsPlots

function scatter_plot_for_mulitple_gml(;
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature_array,
    simulation_path,
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

    # store array with all the paths in the directory to check if we have
    # this .h5 and .gml to be able to plot it.
    spatial_network_path=simulation_path * raw"melting_temp_search\\"
    digital_sphere_masks_dict_path=simulation_path * raw"digital_sphere_masks\\"
    pore_size_dist_data_path=simulation_path * raw"pore_size_dist\\"

    spatial_network_path_array=Glob.glob(filename_start*"*",spatial_network_path)
    pore_size_dist_data_path_array=Glob.glob(filename_start*"*",pore_size_dist_data_path)

    # for data storage 
    data::Vector{Vector{Float64}}=[]

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

                                nr_trials_per_temperature=nr_trials_per_temperature_array[i]
                                
                                filename = (filename_start
                                    *"_N="*"$nr_vertices"
                                    *"_T="*"$maximal_temperature"
                                    *"_Beta="*"$bond_bending_const"
                                    *"_GradT="*"$temperature_gradient"
                                    *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
                                    *"_Theta_GS="*"$theta_ground_state"
                                    *"_Trial="*"$nr_trials_per_temperature"
                                    )

                                total_path=spatial_network_path*filename

                                if(total_path*".gml" in spatial_network_path_array)
                                  
                                    println(filename)

                                    spatial_network=NG.load_spatial_network_from_gml(total_path*".gml")

                                    bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)
                                    bond_angle_std, bond_angle_vec = NA.get_bond_angle_std(spatial_network)

                                    pore_size_dist_data_file = pore_size_dist_data_path * filename *"_pore_size_distribution.h5"
                                    if (!(pore_size_dist_data_file in pore_size_dist_data_path_array))

                                        pore_pixel_radius_array = 
                                            NA.get_pore_size_distribution(
                                                spatial_network;
                                                sampling_grid_size = 0.2,
                                                save_result = true,
                                                save_path = pore_size_dist_data_path*filename,
                                                label = nothing,
                                                digital_sphere_mask_path 
                                                    = raw"simulations\digital_sphere_masks\\",
                                                print_progress = true,
                                                thread_nr = 0,
                                                print_lock = Threads.ReentrantLock()
                                            )

                                    else
                                        println("Already calculated structure")
                                    end
                                    #save_path = simulation_path * raw"pore_size_dist\\"

                                    println("A")
                                    pore_size_dict = GU.load_h5_dict(pore_size_dist_data_file)
                                    println("B")
                                    println("pore_size_dict, $pore_size_dict")
                                    
                                    println("C")
                                    #pore_pixel_radius_vec= pore_pixel_radius_array["pore_size_distribution"]
                                    #println(pore_pixel_radius_vec)
                                    pore_size_distribution_second_moment=NA.get_pore_size_distribution_second_moment(pore_size_dict)
                                    println(pore_size_distribution_second_moment)

                                    

                                    println("bond_length_std, $bond_length_std")
                                    println("bond_angle_std, $bond_angle_std")
                                    println("pore_size_distribution_second_moment, $pore_size_distribution_second_moment")

                                    evolution_dict = GU.load_h5_dict(total_path*"_evolution.h5")
                                    accepted_moves=sum(evolution_dict["move_accepted_vec"])

                                    data=push!(data,
                                        [
                                            #nr_vertices,
                                            maximal_temperature,
                                            bond_bending_const,
                                            temperature_gradient,
                                            nr_monte_carlo_steps_per_temperature,
                                            theta_ground_state,
                                            nr_trials_per_temperature,
                                            bond_length_std,
                                            bond_angle_std,
                                            accepted_moves,
                                            pore_size_distribution_second_moment
                                        ])
                            
                                else
                                    println("file not in directory")
                                    println("total_path, $total_path")
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    
    #println(data)
    #println(typeof(data))

    # data cleaning

    matrix=mapreduce(permutedims, vcat, data)

    df=DataFrames.DataFrame(matrix,:auto)

    DataFrames.rename!(df, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10] .=>  [:MaxT, :Beta, :GradT, :MCsteps, :Theta, :Trial, :BondLenghtStd, :BondAngleStd, :AcceptedMoves, :PoreSizeDistSecMoment])



    

    df.AcceptedMovesLog10 = log10.(df.AcceptedMoves .+ 1)
   
    #println(df)

    fontsize=11

    # plot

    df.Marker = map(
        beta -> beta == 0.0 ? :utriangle :
        beta == 0.25 ? :rect :
        beta == 0.5 ?  :star5 :
        beta == 0.75 ? :star6 :
        beta == 1.0 ? :circle :
        :auto, df.Beta)

    println(df)

    A=@df df Plots.scatter(
        #left_margin=5Plots.PlotMeasures.mm,
        layout = (1,2),
        #size = (1200, 800),
        xtickfont = font(fontsize),  # Set x-axis tick font size
        ytickfont = font(fontsize),   # Set y-axis tick font size
        #=
        xlim=[-0.005,0.18+0.005],
        xticks = (
            0:0.025:0.175, 
            ["0","0.025","0.05","0.075","0.1", "0.125", "0.15", "0.175"]
        ),
        =#
        xguidefont = font(fontsize),
        yguidefont = font(fontsize),
        )
        

    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        color_palette=[:blue,:turquoise1,:green,:yellow,:gold,:orange,:red,:maroon],
        legendtitle =LaTeXStrings.L"T_\mathrm{max}",
        group=:MaxT,
        markershape=:Marker,
        xlabel=LaTeXStrings.L"\sigma_\mathrm{length} / d",
        ylabel=LaTeXStrings.L"\sigma_\mathrm{angle} / \mathrm{rad}",
        subplot=1)

    @df df StatsPlots.scatter!(
        :PoreSizeDistSecMoment, 
        :BondAngleStd, 
        color_palette=[:blue,:turquoise1,:green,:yellow,:gold,:orange,:red,:maroon],
        legendtitle =LaTeXStrings.L"T_\mathrm{max}",
        group=:MaxT,
        markershape=:Marker,
        xlabel=LaTeXStrings.L"r_\mathrm{crit-pore-size}", 
        ylabel=LaTeXStrings.L"\sigma_\mathrm{angle} / \mathrm{rad}",
        subplot=2)    


        
   

    # plot path and name

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

    minimum_nr_trials_per_temperature=minimum(nr_trials_per_temperature_array)
    maximum_nr_trials_per_temperature=maximum(nr_trials_per_temperature_array)

    plot_filename = (plot_filename_start
        *"_N="*"$minimum_nr_vertices" * "-" * "$maximum_nr_vertices"
        *"_T="*"$minimum_temperature" * "-" * "$maximum_temperature"
        *"_Beta="*"$minimum_bond_bending_const" * "-" * "$maximum_bond_bending_const"
        *"_GradT="*"$minimum_temperature_gradient" * "-" * "$maximum_temperature_gradient"
        *"_StepsPerT="*"$minimum_nr_monte_carlo_steps_per_temperature" * "-" * "$maximum_nr_monte_carlo_steps_per_temperature"
        *"_Theta_GS="*"$minimum_theta" * "-" * "$maximum_theta"
        *"_Trials="*"$minimum_nr_trials_per_temperature" * "-" * "$maximum_nr_trials_per_temperature"
        *".png")

    plot_total_path=plot_save_path*plot_filename
    Plots.savefig(A,plot_total_path)

end

network_type="dia"

scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216],
    maximal_temperature_array=[0.5,1.0,1.5],
    bond_bending_const_array=[0.0,0.25,0.5,0.75,1.0],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1], 
    simulation_path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\\",
    filename_start = "mts_5_NW=$(network_type)",
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "thesis_pore_$(network_type)_2"
)