
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
    nr_trials_per_temperature,
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
    @assert nr_trials_per_temperature>=1

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

                            for i in 1:nr_trials_per_temperature
                                
                                filename = (filename_start
                                    *"_N="*"$nr_vertices"
                                    *"_T="*"$maximal_temperature"
                                    *"_Beta="*"$bond_bending_const"
                                    *"_GradT="*"$temperature_gradient"
                                    *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
                                    *"_Theta_GS="*"$theta_ground_state"
                                    *"_Trial="*"$i"
                                    )

                                total_path=save_path*filename

                                if(total_path*".gml" in path_array && 
                                    #filename!="m_BTMC_N=216_T=0.17_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.1_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1"
                                    )
                                    #println("scatter done")
                                    println(filename)

                                    spatial_network=NG.load_spatial_network_from_gml(total_path*".gml")

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
                                        #*"$i"*","
                                        *"$bond_length_std"*","
                                        *"$bond_angle_std"*","
                                        *"$accepted_moves"
                                        *"],")
                                    =#

                                    data=push!(data,
                                        [
                                            #nr_vertices,
                                            maximal_temperature,
                                            bond_bending_const,
                                            temperature_gradient,
                                            nr_monte_carlo_steps_per_temperature,
                                            theta_ground_state,
                                            #i,
                                            bond_length_std,
                                            bond_angle_std,
                                            accepted_moves
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

    DataFrames.rename!(df, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8] .=>  [:MaxT, :Beta, :GradT, :MCsteps, :Theta, :BondLenghtStd, :BondAngleStd, :AcceptedMoves])



    println(df)


    df.BondLenghtStdDiff.=0.0
    df.BondAngleStdDiff.=0.0

    for (i, row) in enumerate( eachrow( df ) ) 
        println(i)
        
        if i<size(df,1)
            if ((df[i,"MaxT"]===df[i+1,"MaxT"]) &&
                (df[i,"Beta"]===df[i+1,"Beta"]) && 
                (df[i,"GradT"]===df[i+1,"GradT"]) && 
                (df[i,"MCsteps"]===df[i+1,"MCsteps"]) && 
                (df[i,"Theta"]!=df[i+1,"Theta"]))

                if df[i,"Theta"]===180.0 && df[i+1,"Theta"]===110.0

                    df[i,"BondLenghtStdDiff"]=df[i,"BondLenghtStd"]-df[i+1,"BondLenghtStd"]
                    df[i,"BondAngleStdDiff"]=df[i,"BondAngleStd"]-df[i+1,"BondAngleStd"]

                elseif df[i+1,"Theta"]===180.0 && df[i,"Theta"]===110.0

                    df[i,"BondLenghtStdDiff"]=df[i+1,"BondLenghtStd"]-df[i,"BondLenghtStd"]
                    df[i,"BondAngleStdDiff"]=df[i+1,"BondAngleStd"]-df[i,"BondAngleStd"]
                
                else
                    println("i+1: Theta not 180 and 110")
                end
            elseif (i>1 && (
                (df[i,"MaxT"]===df[i-1,"MaxT"]) && 
                (df[i,"Beta"]===df[i-1,"Beta"]) && 
                (df[i,"GradT"]===df[i-1,"GradT"]) && 
                (df[i,"MCsteps"]===df[i-1,"MCsteps"]) && 
                (df[i,"Theta"]!=df[i-1,"Theta"])))

                if df[i,"Theta"]===180.0 && df[i-1,"Theta"]===110.0

                    df[i,"BondLenghtStdDiff"]=df[i,"BondLenghtStd"]-df[i-1,"BondLenghtStd"]
                    df[i,"BondAngleStdDiff"]=df[i,"BondAngleStd"]-df[i-1,"BondAngleStd"]

                elseif df[i-1,"Theta"]===180.0 && df[i,"Theta"]===110.0

                    df[i,"BondLenghtStdDiff"]=df[i-1,"BondLenghtStd"]-df[i,"BondLenghtStd"]
                    df[i,"BondAngleStdDiff"]=df[i-1,"BondAngleStd"]-df[i,"BondAngleStd"]
                
                else
                    println("i-1: Theta not 180 and 110")
                end
            else
                println("Beta, GradT, MCsteps not the same OR Theta the same")
            end
        else
            if ((df[i,"MaxT"]===df[i-1,"MaxT"]) &&
                (df[i,"Beta"]===df[i-1,"Beta"]) && 
                (df[i,"GradT"]===df[i-1,"GradT"]) && 
                (df[i,"MCsteps"]===df[i-1,"MCsteps"]) && 
                (df[i,"Theta"]!=df[i-1,"Theta"]))

                if df[i,"Theta"]===180.0 && df[i-1,"Theta"]===110.0

                    df[i,"BondLenghtStdDiff"]=df[i,"BondLenghtStd"]-df[i-1,"BondLenghtStd"]
                    df[i,"BondAngleStdDiff"]=df[i,"BondAngleStd"]-df[i-1,"BondAngleStd"]

                elseif df[i-1,"Theta"]===180.0 && df[i,"Theta"]===110.0

                    df[i,"BondLenghtStdDiff"]=df[i-1,"BondLenghtStd"]-df[i,"BondLenghtStd"]
                    df[i,"BondAngleStdDiff"]=df[i-1,"BondAngleStd"]-df[i,"BondAngleStd"]
                
                else
                    println("i-1: Theta not 180 and 110")
                end
            end
        end
    end


    df.BondLenghtStdDiffLog10 = log10.(df.BondLenghtStdDiff .+ 1)
    df.BondAngleStdDiffLog10 = log10.(df.BondAngleStdDiff .+ 1)



    df.Shape = @. ifelse(df.Theta == 110.0, :diamond, :rect); df

    df.GradTLog10 = log10.(df.GradT)
    df.AcceptedMovesLog10 = log10.(df.AcceptedMoves .+ 1)
   
    #println(df)


    # plot

    A=@df df Plots.scatter(
        #plot_title=plot_title,
        plot_title=LaTeXStrings.LaTeXString("Bond length std vs Bond angle std\n\$^\\textrm{"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "Tmax="*"$(maximal_temperature_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])"*
            "}\$"),
        top_margin=5Plots.PlotMeasures.mm,
        titlefont = font(6),
        layout = (2,2),
        size = (900, 700),
        colorbar=:top,
        legend=false,
        xlabel="Bond length std / d", 
        ylabel="Bond angle std / rad"
        )

    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        colorbar_title="MaxT",
        zcolor=:MaxT,
        color=cgrad(:roma, rev=true),
        shape=:Shape,
        subplot=1)

    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        colorbar_title="Beta",
        zcolor=:Beta,
        color=cgrad(:roma, rev=true),
        shape=:Shape,
        subplot=2)

    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        colorbar_title="Theta",
        zcolor=:Theta,
        color=cgrad(:roma, rev=true),
        shape=:Shape,
        subplot=3)
    
    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        colorbar_title="AcceptedMovesLog10",
        zcolor=:AcceptedMovesLog10,
        color=cgrad(:roma, rev=true),
        shape=:Shape,
        subplot=4)







    # plot

    B=@df df Plots.scatter(
        #plot_title=plot_title,
        plot_title=LaTeXStrings.LaTeXString("Bond length std diff log10 vs Bond angle std diff log10 \n\$^\\textrm{"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "Tmax="*"$(maximal_temperature_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])"*
            "}\$"),
        top_margin=5Plots.PlotMeasures.mm,
        titlefont = font(6),
        layout = (2,2),
        size = (900, 700),
        colorbar=:top,
        legend=false,
        xlabel="Bond length std diff log10 / d", 
        ylabel="Bond angle std diff log10 / rad"
        )

    @df df StatsPlots.scatter!(
        :BondLenghtStdDiffLog10, 
        :BondAngleStdDiffLog10, 
        colorbar_title="MaxT",
        zcolor=:MaxT,
        color=cgrad(:roma, rev=true),
        shape=:Shape,
        subplot=1)

    @df df StatsPlots.scatter!(
        :BondLenghtStdDiffLog10, 
        :BondAngleStdDiffLog10,
        colorbar_title="Beta",
        zcolor=:Beta,
        color=cgrad(:roma, rev=true),
        shape=:Shape,
        subplot=2)

    @df df StatsPlots.scatter!(
        :BondLenghtStdDiffLog10, 
        :BondAngleStdDiffLog10,
        colorbar_title="Theta",
        zcolor=:Theta,
        color=cgrad(:roma, rev=true),
        shape=:Shape,
        subplot=3)
    
    @df df StatsPlots.scatter!(
        :BondLenghtStdDiffLog10, 
        :BondAngleStdDiffLog10,
        colorbar_title="AcceptedMovesLog10",
        zcolor=:AcceptedMovesLog10,
        color=cgrad(:roma, rev=true),
        shape=:Shape,
        subplot=4)    





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

    plot_filename = (plot_filename_start
        *"_N="*"$minimum_nr_vertices" * "-" * "$maximum_nr_vertices"
        *"_T="*"$minimum_temperature" * "-" * "$maximum_temperature"
        *"_Beta="*"$minimum_bond_bending_const" * "-" * "$maximum_bond_bending_const"
        *"_GradT="*"$minimum_temperature_gradient" * "-" * "$maximum_temperature_gradient"
        *"_StepsPerT="*"$minimum_nr_monte_carlo_steps_per_temperature" * "-" * "$maximum_nr_monte_carlo_steps_per_temperature"
        *"_Theta_GS="*"$minimum_theta" * "-" * "$maximum_theta"
        *"_Trials="*"$nr_trials_per_temperature"
        *".png")

    plot_filename_diff = (plot_filename_start
        *"d_"
        *"_N="*"$minimum_nr_vertices" * "-" * "$maximum_nr_vertices"
        *"_T="*"$minimum_temperature" * "-" * "$maximum_temperature"
        *"_Beta="*"$minimum_bond_bending_const" * "-" * "$maximum_bond_bending_const"
        *"_GradT="*"$minimum_temperature_gradient" * "-" * "$maximum_temperature_gradient"
        *"_StepsPerT="*"$minimum_nr_monte_carlo_steps_per_temperature" * "-" * "$maximum_nr_monte_carlo_steps_per_temperature"
        *"_Theta_GS="*"$minimum_theta" * "-" * "$maximum_theta"
        *"_Trials="*"$nr_trials_per_temperature"
        *".png")

    plot_total_path=plot_save_path*plot_filename
    plot_total_path_diff=plot_save_path*plot_filename_diff

    Plots.savefig(A,plot_total_path)

    Plots.savefig(B,plot_total_path_diff)

end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216],
    maximal_temperature_array=[0.1,0.135,0.17,0.205,0.24],
    bond_bending_const_array=[0,0.21,0.285,0.36,0.5],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[110.0,180.0],
    nr_trials_per_temperature=1, 
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_BTMC_q_t_",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_BTMC_matrix_q_t_5_"
)
