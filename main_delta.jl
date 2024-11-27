
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
import CairoMakie
import FileIO
using StatsPlots
import PlotlyJS


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
                                    filename!="m_BTMC_q_t__N=216_T=0.1_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=180.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.1_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" #&&
                                    #filename!="m_BTMC_q_t__N=216_T=0.17_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=0.0_Trial=1" &&
                                    #filename!="m_BTMC_q_t__N=216_T=0.205_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=0.0_Trial=1"
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
                                    #=
                                    data=push!(data,
                                        [
                                            #nr_vertices,
                                            maximal_temperature,
                                            bond_bending_const,
                                            temperature_gradient,
                                            nr_monte_carlo_steps_per_temperature,
                                            theta_ground_state,
                                            #i,
                                            0,
                                            0,
                                            0
                                        ])
                                    =#
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    
    #println(data)
    println(typeof(data))

    # data cleaning

    matrix=mapreduce(permutedims, vcat, data)
    println(typeof(matrix))

    df=DataFrames.DataFrame(matrix,:auto)

    println(typeof(df))

    DataFrames.rename!(df, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8] .=>  [:MaxT, :Beta, :GradT, :MCsteps, :Theta, :BondLenghtStd, :BondAngleStd, :AcceptedMoves])



    println(df)


    df.BondLenghtStdDiff.=0.0
    df.BondAngleStdDiff.=0.0
    df.BondLenghtStdRatio.=0.0
    df.BondAngleStdRatio.=0.0

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
                    df[i,"BondLenghtStdRatio"]=df[i,"BondLenghtStd"]/df[i+1,"BondLenghtStd"]
                    df[i,"BondAngleStdRatio"]=df[i,"BondAngleStd"]/df[i+1,"BondAngleStd"]

                elseif df[i+1,"Theta"]===180.0 && df[i,"Theta"]===110.0

                    df[i,"BondLenghtStdDiff"]=df[i+1,"BondLenghtStd"]-df[i,"BondLenghtStd"]
                    df[i,"BondAngleStdDiff"]=df[i+1,"BondAngleStd"]-df[i,"BondAngleStd"]
                    df[i,"BondLenghtStdRatio"]=df[i+1,"BondLenghtStd"]/df[i,"BondLenghtStd"]
                    df[i,"BondAngleStdRatio"]=df[i+1,"BondAngleStd"]/df[i,"BondAngleStd"]
                
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
                    df[i,"BondLenghtStdRatio"]=df[i,"BondLenghtStd"]/df[i-1,"BondLenghtStd"]
                    df[i,"BondAngleStdRatio"]=df[i,"BondAngleStd"]/df[i-1,"BondAngleStd"]

                elseif df[i-1,"Theta"]===180.0 && df[i,"Theta"]===110.0

                    df[i,"BondLenghtStdDiff"]=df[i-1,"BondLenghtStd"]-df[i,"BondLenghtStd"]
                    df[i,"BondAngleStdDiff"]=df[i-1,"BondAngleStd"]-df[i,"BondAngleStd"]
                    df[i,"BondLenghtStdRatio"]=df[i-1,"BondLenghtStd"]/df[i,"BondLenghtStd"]
                    df[i,"BondAngleStdRatio"]=df[i-1,"BondAngleStd"]/df[i,"BondAngleStd"]
                
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
                    df[i,"BondLenghtStdRatio"]=df[i,"BondLenghtStd"]/df[i-1,"BondLenghtStd"]
                    df[i,"BondAngleStdRatio"]=df[i,"BondAngleStd"]/df[i-1,"BondAngleStd"]

                elseif df[i-1,"Theta"]===180.0 && df[i,"Theta"]===110.0

                    df[i,"BondLenghtStdDiff"]=df[i-1,"BondLenghtStd"]-df[i,"BondLenghtStd"]
                    df[i,"BondAngleStdDiff"]=df[i-1,"BondAngleStd"]-df[i,"BondAngleStd"]
                    df[i,"BondLenghtStdRatio"]=df[i-1,"BondLenghtStd"]/df[i,"BondLenghtStd"]
                    df[i,"BondAngleStdRatio"]=df[i-1,"BondAngleStd"]/df[i,"BondAngleStd"]
                
                else
                    println("i-1: Theta not 180 and 110")
                end
            end
        end
    end

    d=1
    df.Delta1=sqrt.(df.BondLenghtStdDiff.^2 ./ d.^2 + df.BondAngleStdDiff.^2)
    df.Delta2=sqrt.( (1 .- df.BondLenghtStdRatio).^2 + (1 .- df.BondAngleStdRatio).^2)

    df.Shape = @. ifelse(df.Theta == 180.0, :rect, 
                    ifelse(df.Theta == 110.0, :diamond, 
                    ifelse(df.Theta == 90.0, :circle, :x))); df

    df.GradTLog10 = log10.(df.GradT)
    df.AcceptedMovesLog10 = log10.(df.AcceptedMoves .+ 1)
   
    #println(df)


    # plot

    A=@df df Plots.histogram(
        :Delta1,
        plot_title=LaTeXStrings.LaTeXString("Number of networks with a Delta1\n\$^\\textrm{"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "GradT="*"$(temperature_gradient_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])"*
            "}\$"),
    
        top_margin=10Plots.PlotMeasures.mm,
        titlefont = font(6),
        legend=false,
        xlabel="Delta1 / rad", 
        ylabel="Number of networks"
        )

    B=@df df Plots.histogram(
        :Delta2,
        plot_title=LaTeXStrings.LaTeXString("Number of networks with a Delta2\n\$^\\textrm{"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "GradT="*"$(temperature_gradient_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])"*
            "}\$"),
        
        top_margin=10Plots.PlotMeasures.mm,
        titlefont = font(6),
        legend=false,
        xlabel="Delta2 / rad", 
        ylabel="Number of networks"
        )

    
    # D1
    # take only Theta=180 and then take the delta dataframe array, convert to matrix and shuffle matrix
    df180=df[(df.Theta .=== 180.0), :]
    Delta1_array=Array(df180[:, 13])
    Delta1_matrix=reshape(Delta1_array,length(bond_bending_const_array),length(maximal_temperature_array))
    Delta1_permuted=permutedims(Delta1_matrix,[2,1])

    # H1

    Delta1_trace = PlotlyJS.heatmap(
        x=bond_bending_const_array,
        y=maximal_temperature_array,
        z=Delta1_permuted
    )

    Delta1_layout = PlotlyJS.Layout(
        title="Delta1 for different MaxT and Beta", 
        xaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=bond_bending_const_array,
            tickvals=bond_bending_const_array
        ),
        yaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=maximal_temperature_array,
            tickvals=maximal_temperature_array
        ),
        xaxis_title="Beta",
        yaxis_title="MaxT",
    )

    C=PlotlyJS.plot(Delta1_trace,Delta1_layout)

    # D2

    Delta2_array=Array(df180[:, 14])
    Delta2_matrix=reshape(Delta2_array,length(bond_bending_const_array),length(maximal_temperature_array))
    Delta2_permuted=permutedims(Delta2_matrix,[2,1])

    # H2

    Delta2_trace = PlotlyJS.heatmap(
        x=bond_bending_const_array,
        y=maximal_temperature_array,
        z=Delta2_permuted
    )

    Delta2_layout = PlotlyJS.Layout(
        title="Delta1 for different MaxT and Beta", 
        xaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=bond_bending_const_array,
            tickvals=bond_bending_const_array
        ),
        yaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=maximal_temperature_array,
            tickvals=maximal_temperature_array
        ),
        xaxis_title="Beta",
        yaxis_title="MaxT",
    )

    D=PlotlyJS.plot(Delta2_trace,Delta2_layout)




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
 
    plot_filename_end=(
        "_N="*"$minimum_nr_vertices" * "-" * "$maximum_nr_vertices"
        *"_T="*"$minimum_temperature" * "-" * "$maximum_temperature"
        *"_Beta="*"$minimum_bond_bending_const" * "-" * "$maximum_bond_bending_const"
        *"_GradT="*"$minimum_temperature_gradient" * "-" * "$maximum_temperature_gradient"
        *"_StepsPerT="*"$minimum_nr_monte_carlo_steps_per_temperature" * "-" * "$maximum_nr_monte_carlo_steps_per_temperature"
        *"_Theta_GS="*"$minimum_theta" * "-" * "$maximum_theta"
        *"_Trials="*"$nr_trials_per_temperature"
        *".png")

    plot_filename_A = (
        plot_filename_start
        *"_Hist_D1"
        *plot_filename_end)


    plot_filename_B = (
        plot_filename_start
        *"_Hist_D2"
        *plot_filename_end)

    plot_filename_C = (
        plot_filename_start
        *"_Heat_D1"
        *plot_filename_end)

    plot_filename_D = (
        plot_filename_start
        *"_Heat_D2"
        *plot_filename_end)

    plot_total_path_A=plot_save_path*plot_filename_A
    plot_total_path_B=plot_save_path*plot_filename_B
    plot_total_path_C=plot_save_path*plot_filename_C
    plot_total_path_D=plot_save_path*plot_filename_D

    Plots.savefig(A,plot_total_path_A)
    Plots.savefig(B,plot_total_path_B)
    PlotlyJS.savefig(C,plot_total_path_C)
    PlotlyJS.savefig(D,plot_total_path_D)


end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216],
    maximal_temperature_array=[0.1,0.125,0.15,0.175,0.2],
    bond_bending_const_array=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[110.0,180.0],
    nr_trials_per_temperature=1,
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_BTMC_q_t_",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_delta_9_"
)
