
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

    # store array with all the paths in the directory to check if we have
    # this .h5 and .gml to be able to plot it.
    path_array=Glob.glob(filename_start*"*",save_path)

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

                            for i in eachindex(nr_trials_per_temperature)

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

                                    stretching_energy=NG.get_stretching_energy_keating(spatial_network)
                                    bending_energy=NG.get_bending_energy_keating(spatial_network)
                                    total_energy=NG.get_total_energy_keating(spatial_network)
                                    ratio_energy=stretching_energy/bending_energy

                                    #=display(bending_energy)
                                    display(stretching_energy)
                                    display(total_energy)
                                    display(ratio_energy)=#

                                    @assert stretching_energy+bending_energy===total_energy



                                    
                                    bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)
                                    bond_angle_std, bond_angle_vec = NA.get_bond_angle_std(spatial_network)

                                    evolution_dict = GU.load_h5_dict(total_path*"_evolution.h5")
                                    accepted_moves=sum(evolution_dict["move_accepted_vec"])
                                    

                                    #=
                                    println("["
                                        #*"$nr_vertices"*","
                                        *"$maximal_temperature"*","
                                        *"$bond_bending_const"*","
                                        *"$temperature_gradient"*","
                                        #*"$nr_monte_carlo_steps_per_temperature"*","
                                        *"$theta_ground_state"*","
                                        #*"$trial"*","
                                        *"$bond_length_std"*","
                                        *"$bond_angle_std"*","
                                        *"$accepted_moves"
                                        *"],")
                                    =#

                                    data=push!(data,
                                        [
                                            nr_vertices,
                                            maximal_temperature,
                                            bond_bending_const,
                                            temperature_gradient,
                                            nr_monte_carlo_steps_per_temperature,
                                            theta_ground_state,
                                            #trial,
                                            bond_length_std,
                                            bond_angle_std,
                                            accepted_moves,
                                            stretching_energy,
                                            bending_energy,
                                            total_energy,
                                            ratio_energy
                                        ])
                                    
                                    

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

    
    #println(data)
    #println(typeof(data))

    # data cleaning

    matrix=mapreduce(permutedims, vcat, data)

    df=DataFrames.DataFrame(matrix,:auto)

    DataFrames.rename!(df, 
        [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13] 
        .=>  [:NrVertices, :MaxT, :Beta, :GradT, :MCsteps, :Theta, :BondLenghtStd, 
            :BondAngleStd, :AccMoves, :StretchEnergy, :BendEnergy ,
            :TotalEnergy, :RatioEnergy])

    df.AccMovesLg10 = log10.(df.AccMoves .+ 1)
    #println(df)

    #Cleaning
    df=filter(row -> row.AccMoves > 250, df)

    #Calculate first definition "a" with sqrt
    df_calculation=deepcopy(df)
    df_calculation=filter(row -> row.Theta === 110.0, df_calculation)
    NrAngle=6*df_calculation.NrVertices
    df_calculation.a1=sqrt.(8/3*df_calculation.BendEnergy ./ NrAngle)
    deltaCos=cosd(180)-cosd(110)
    df_calculation.betaRatioPred1=(df_calculation.a1 ./ (df_calculation.a1 .- deltaCos)).^2


    #Calculate second definition "a" with 4 terms
    
    x1=(1 .+ df_calculation.BondLenghtStd).^2
    x2=(1 .- df_calculation.BondLenghtStd).^2
    y1=cosd.(110 .+ df_calculation.BondAngleStd/pi*180) .- cosd(110)
    y2=cosd.(110 .- df_calculation.BondAngleStd/pi*180) .- cosd(110)
    display(x1)
    display(x2)
    display(y1)
    display(y2)

    df_calculation.a2=1/4 .* (
        (x1 .- y1) .^2
        .+ (x2 .- y1) .^2
        .+ (x1 .- y2) .^2
        .+ (x2 .- y2) .^2
    )

    df_calculation.betaRatioPred2=(df_calculation.a2 ./ (df_calculation.a2 .- deltaCos)).^2
        
    display(df_calculation)

    mean_df_calculation = DataFrames.DataFrame([DataFrames.mean(df_calculation[!, col]) for col in DataFrames.names(df_calculation)]', DataFrames.names(df_calculation))
    display("Means of all columns:")
    display(round.(mean_df_calculation, digits=3))


end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216],
    maximal_temperature_array=[0.1,0.125,0.15,0.175,0.2],
    bond_bending_const_array=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[110.0,180.0], #[100.0,110.0,130.0,150.0,180.0],
    nr_trials_per_temperature=1, 
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_rad_",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_rad_br_1_"
)
