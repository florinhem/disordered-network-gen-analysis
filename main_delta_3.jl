
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
import Statistics



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
                                    filename!="m_BTMC_q_t__N=216_T=0.125_Beta=0.25_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.125_Beta=0.3_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.125_Beta=0.3_GradT=0.1_StepsPerT=0.01_Theta_GS=180.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.15_Beta=0.3_GradT=0.1_StepsPerT=0.01_Theta_GS=180.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.15_Beta=0.25_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.15_Beta=0.25_GradT=0.1_StepsPerT=0.01_Theta_GS=180.0_Trial=1"
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


    #-----------------------
    # Delta2
    #-----------------------

    Delta2Matrix=zeros(size(df,1),size(df,1),3)
    #println(Delta2Matrix)
    for i in 1:size(df,1)
        for j in 1:size(df,1)
            delta2=sqrt((1-df[i,6]/df[j,6])^2+(1-df[i,7]/df[j,7])^2)
            if df[i,5]===df[j,5] || i===j #i>=j ||
                Delta2Matrix[i,j,3]=Inf64
            elseif df[i,5]===110.0 && df[j,5]===180.0
                Delta2Matrix[i,j,1]=df[i,2]-df[j,2]
                Delta2Matrix[i,j,2]=df[i,1]-df[j,1]
                Delta2Matrix[i,j,3]=delta2
            elseif df[i,5]===180.0 && df[j,5]===110.0
                Delta2Matrix[i,j,1]=df[j,2]-df[i,2]
                Delta2Matrix[i,j,2]=df[j,1]-df[i,1]
                Delta2Matrix[i,j,3]=delta2
            end
        end
    end
    #display(Delta2Matrix)

    # SCATTER
    Delta2MatrixRound=zeros(size(df,1),size(df,1),3)
    display(Delta2Matrix)
    for i in 1:size(Delta2Matrix,1)
        for j in 1:size(Delta2Matrix,2)
            for k in 1:size(Delta2Matrix,3)
                Delta2MatrixRound[i,j,k]=round(Delta2Matrix[i,j,k],digits=3)
            end
        end
    end
    display(Delta2MatrixRound)

    


    x=Delta2MatrixRound[:,:,1] |> unique |> sort
    x=filter(!iszero,x)
    y=Delta2MatrixRound[:,:,2] |> unique |> sort
    y=filter(!iszero,y)
    z=Delta2MatrixRound[:,:,3] |> vec |> sort
    z=filter(!isinf,z)

    display(x)
    display(y)
    display(z)

    #ratio_small=0.025
    ratio_small=0.75
    #ratio_small=1
    number_of_smallest=ceil(Int,ratio_small*size(z,1))
    max=maximum(first(z, number_of_smallest), dims = 1)[1]
    display("max: $max")
    


    Counts=zeros(size(x,1),size(y,1))
    PossibleCounts=zeros(size(x,1),size(y,1))
    for i in 1:size(Delta2MatrixRound,1)
        for j in 1:size(Delta2MatrixRound,2)
            T=Delta2MatrixRound[i,j,1]
            B=Delta2MatrixRound[i,j,2]
            D=Delta2MatrixRound[i,j,3]
            #display(T)
            #display(B)
            #display(D)
           
            for k in 1:size(x,1)
                for m in 1:size(y,1)
                    if T===x[k] && B===y[m]
                        if D<=max
                            Counts[k,m]+=1
                        end
                        if D!=Inf64
                            #display("$T,$B,$D")
                            PossibleCounts[k,m]+=1
                        end
                    end
                end
            end
        end
    end
    display(Counts)
    display(PossibleCounts)

    NormalizedCounts=zeros(size(x,1),size(y,1))
    for i in 1:size(x,1)
        for j in 1:size(y,1)
            NormalizedCounts[i,j]=Counts[i,j]/PossibleCounts[i,j]
        end
    end

    display(NormalizedCounts)


    




    

    
    
   

    # HISTOGRAM
    Delta2Sort=Delta2MatrixRound[Delta2MatrixRound[:,:,3].<=max,:]

    # plot
    Plot_H2=Plots.histogram(
        Delta2Sort[:,3],
        title=LaTeXStrings.LaTeXString(
            "Distribution of smallest $(ratio_small*100)% of Delta2 \n\$^\\textrm{"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "GradT="*"$(temperature_gradient_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])"*
            "}\$"),
        xlabel=LaTeXStrings.LaTeXString("\$ \\Delta_2 \$"),
        ylabel="Number of networks",
        legend=false
    )

    #display(Plot_H2)

    




    # HEATMAP
  
    NormalizedCounts_permuted=permutedims(NormalizedCounts,[2,1])
    display(NormalizedCounts_permuted)
    
    # plot
    Delta2_trace=PlotlyJS.heatmap(
        x=x,
        y=y,
        z=NormalizedCounts_permuted,
        vline=0
    )

    #PlotlyJS.line(0,0,-1,1)
    #vline=PlotlyJS.add_vline!(x=0)
    

    Delta2_layout = PlotlyJS.Layout(
        title=
            "Number of instances of smallest $(ratio_small*100)% of Delta2, <br>"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "GradT="*"$(temperature_gradient_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])",
        xaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=x,
            tickvals=x
        ),
        yaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=y,
            tickvals=y
        ),
        xaxis_title="Delta beta",
        yaxis_title="Delta MaxT",
        autosize=false
        
    )

    Plot_Heat2=PlotlyJS.plot(Delta2_trace,Delta2_layout)


    PlotlyJS.add_shape!(Plot_Heat2, PlotlyJS.line(
        x0=0, 
        x1=0, 
        y0=minimum(y)-0.0125,
        y1=maximum(y)+0.0125,
        line=PlotlyJS.attr(color=:lightgreen, width=3),
    ))

    PlotlyJS.add_shape!(Plot_Heat2, PlotlyJS.line(
        x0=minimum(x)-0.025, 
        x1=maximum(x)+0.025, 
        y0=0,
        y1=0,
        line=PlotlyJS.attr(color=:lightgreen, width=3),
    ))

    # statistics
    x_mean=Statistics.mean(Delta2Sort[:,1])
    y_mean=Statistics.mean(Delta2Sort[:,2])
    x_std=Statistics.std(Delta2Sort[:,1])
    y_std=Statistics.std(Delta2Sort[:,2])

    if isnan(x_std)
        x_std=0
    end
    if isnan(y_std)
        y_std=0
    end

    println("$x_mean, $y_mean, $x_std, $y_std")

    PlotlyJS.add_shape!(Plot_Heat2, PlotlyJS.line(
        x0=x_mean-x_std, 
        x1=x_mean+x_std,
        y0=y_mean,
        y1=y_mean,
        line=PlotlyJS.attr(color=:red, width=3),
    ))

    PlotlyJS.add_shape!(Plot_Heat2, PlotlyJS.line(
        x0=x_mean, 
        x1=x_mean,
        y0=y_mean-y_std,
        y1=y_mean+y_std,
        line=PlotlyJS.attr(color=:red, width=3),
    ))

    PlotlyJS.add_shape!(Plot_Heat2, PlotlyJS.circle(
        x0=x_mean*0.95, 
        x1=x_mean*1.05,
        y0=y_mean*0.95,
        y1=y_mean*1.05,
        line=PlotlyJS.attr(color=:red, width=3),
    ))
    #display(Plot_Heat2)

    








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

   

    plot_filename_B = (
        plot_filename_start
        *"_Histogram"
        *plot_filename_end)

    plot_filename_C = (
        plot_filename_start
        *"_Heatmap"
        *plot_filename_end)

    
    plot_total_path_B=plot_save_path*plot_filename_B
    plot_total_path_C=plot_save_path*plot_filename_C

    
    Plots.savefig(Plot_H2,plot_total_path_B)
    PlotlyJS.savefig(Plot_Heat2,plot_total_path_C)


end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216],
    maximal_temperature_array=[0.13,0.14,0.15,0.16,0.17, 0.125,0.135,0.145,0.155,0.165],                   #[0.1,0.125,0.15,0.175,0.2],
    bond_bending_const_array=[0.24,0.26,0.28,0.30,0.32, 0.23,0.25,0.27,0.29,0.31],     #[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[110.0,180.0],
    nr_trials_per_temperature=1,
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_BTMC_q_t_",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_d3_f9_"
)
