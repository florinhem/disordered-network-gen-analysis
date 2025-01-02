
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
import LsqFit


function linear_f(x, p)
    a1, b1 = p
    display("x= $x")
    display("a1= $a1")
    display("b1= $b1")
    display("a1 * x .+ b1= $(a1 * x .+ b1)")
    return a1 * x .+ b1
end


function scatter_plot_for_mulitple_gml(;
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature,
    theta_compare,
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

                                if(total_path*".gml" in path_array #=&& 
                                    filename!="m_BTMC_q_t__N=216_T=0.125_Beta=0.25_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1"=#
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
    #println(typeof(data))

    # data cleaning

    matrix=mapreduce(permutedims, vcat, data)
    #println(typeof(matrix))

    df=DataFrames.DataFrame(matrix,:auto)

    #println(typeof(df))

    DataFrames.rename!(df, [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8] .=>  [:MaxT, :Beta, :GradT, :MCsteps, :Theta, :BondLenghtStd, :BondAngleStd, :AcceptedMoves])


    #-----------------------
    # Delta2
    #-----------------------

    Delta2Matrix=zeros(size(df,1),size(df,1),7)
    #println(Delta2Matrix)
    for i in 1:size(df,1)
        for j in 1:size(df,1)
            delta2=sqrt((1-df[i,6]/df[j,6])^2+(1-df[i,7]/df[j,7])^2)
            Delta2Matrix[i,j,1]=df[i,2]
            Delta2Matrix[i,j,2]=df[j,2]
            Delta2Matrix[i,j,3]=df[i,1]
            Delta2Matrix[i,j,4]=df[j,1]
            Delta2Matrix[i,j,5]=delta2
            if  i===j || df[i,5]===df[j,5]
                Delta2Matrix[i,j,5]=Inf64
            end
            Delta2Matrix[i,j,6]=df[i,5]
            Delta2Matrix[i,j,7]=df[j,5]
        end
    end
    #display(Delta2Matrix)

    Delta2MatrixRound=zeros(size(df,1),size(df,1),7)
    #display(Delta2Matrix)
    for i in 1:size(Delta2Matrix,1)
        for j in 1:size(Delta2Matrix,2)
            for k in 1:size(Delta2Matrix,3)
                Delta2MatrixRound[i,j,k]=round(Delta2Matrix[i,j,k],digits=3)
            end
        end
    end
    #display(Delta2MatrixRound)




    b1=Delta2MatrixRound[:,:,1] |> unique |> sort
    b1=filter(!iszero,b1)
    b2=Delta2MatrixRound[:,:,2] |> unique |> sort
    b2=filter(!iszero,b2)
    t1=Delta2MatrixRound[:,:,3] |> unique |> sort
    t1=filter(!iszero,t1)
    t2=Delta2MatrixRound[:,:,4] |> unique |> sort
    t2=filter(!iszero,t2)
    d=Delta2MatrixRound[:,:,5] |> vec |> sort
    d=filter(!isinf,d)
    theta1=Delta2MatrixRound[:,:,6] |> unique |> sort
    theta1=filter(!iszero,theta1)
    theta2=Delta2MatrixRound[:,:,7] |> unique |> sort
    theta2=filter(!iszero,theta2)

    #display(b1)
    #display(b2)
    #display(t1)
    #display(t2)
    #display(d)
    #display(theta1)
    #display(theta2)

    ratio_small=0.025
    #ratio_small=0.2
    #ratio_small=0.5
    number_of_smallest=ceil(Int,ratio_small*size(d,1))
    max=maximum(first(d, number_of_smallest), dims = 1)[1]
    #display("max: $max")
    


    Countsb1b2=zeros(size(b1,1),size(b2,1))
    PossibleCountsb1b2=zeros(size(b1,1),size(b2,1))
    Countst1t2=zeros(size(t1,1),size(t2,1))
    PossibleCountst1t2=zeros(size(t1,1),size(t2,1))
    for i in 1:size(Delta2MatrixRound,1)
        for j in 1:size(Delta2MatrixRound,2)
            B1=Delta2MatrixRound[i,j,1]
            B2=Delta2MatrixRound[i,j,2]
            T1=Delta2MatrixRound[i,j,3]
            T2=Delta2MatrixRound[i,j,4]
            D=Delta2MatrixRound[i,j,5]
            Theta1=Delta2MatrixRound[i,j,6]
            Theta2=Delta2MatrixRound[i,j,7]
            #display("$B1, $B2, $T1, $T2, $D")
            #display(B)
            #display(D)
           
            for k in 1:size(b1,1)
                for m in 1:size(b2,1)
                    if B1===b1[k] && B2===b2[m] && Theta1===180.0 && Theta2===theta_compare
                        if D<=max
                            #display("$B1, $B2")
                            Countsb1b2[k,m]+=1
                        end
                        if D!=Inf64
                            #display("$T,$B,$D")
                            PossibleCountsb1b2[k,m]+=1
                        end
                    end
                end
            end

            for k in 1:size(t1,1)
                for m in 1:size(t2,1)
                    if T1===t1[k] && T2===t2[m] && Theta1===180.0 && Theta2===theta_compare && B2-B1===0.1
                        if D<=max
                            display("$B1, $B2")
                            Countst1t2[k,m]+=1
                        end
                        if D!=Inf64
                            #display("$T,$B,$D")
                            PossibleCountst1t2[k,m]+=1
                        end
                    end
                end
            end
        end
    end
    #display(Countsb1b2)
    #display(PossibleCountsb1b2)
    #display(Countst1t2)
    #display(PossibleCountst1t2)

    NormalizedCountsb1b2=zeros(size(b1,1),size(b2,1))
    for i in 1:size(b1,1)
        for j in 1:size(b2,1)
            NormalizedCountsb1b2[i,j]=Countsb1b2[i,j]/PossibleCountsb1b2[i,j]
        end
    end

    #display(NormalizedCountsb1b2)

    NormalizedCountst1t2=zeros(size(t1,1),size(t2,1))
    for i in 1:size(t1,1)
        for j in 1:size(t2,1)
            NormalizedCountst1t2[i,j]=Countst1t2[i,j]/PossibleCountst1t2[i,j]
        end
    end

    #display(NormalizedCountst1t2)



    



    
    # HISTOGRAM
    Delta2Sort=Delta2MatrixRound[Delta2MatrixRound[:,:,5].<=max,:]
    #=
    display("Delta2Sort:")
    println("[")
    for i in 1:size(Delta2Sort,1)
        print("[")
        for j in 1:size(Delta2Sort,2)-1
            print("$(Delta2Sort[i,j]), ")
        end
        print("$(Delta2Sort[i,size(Delta2Sort,2)])")
        println("],")
    end
    println("]")
    =#

    # plot
    Plot_H2=Plots.histogram(
        Delta2Sort[:,5],
        title=LaTeXStrings.LaTeXString(
            "Distribution of smallest $(ratio_small*100)% of Delta \n\$^\\textrm{"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "GradT="*"$(temperature_gradient_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])"*
            "}\$"),
        xlabel=LaTeXStrings.LaTeXString("\$ \\Delta \$"),
        ylabel="Number of networks",
        legend=false
    )

    #display(Plot_H2)

    




    # HEATMAP
  
    #NormalizedCounts_permutedb1b2=permutedims(NormalizedCountsb1b2,[2,1])
    #display(NormalizedCounts_permutedb1b2)
    
    
    # Plot_Heatb1b2
    Delta2_traceb1b2=PlotlyJS.heatmap(
        x=b1,
        y=b2,
        z=NormalizedCountsb1b2,
        vline=0
    )
    

    Delta2_layoutb1b2 = PlotlyJS.Layout(
        title=
            "Number of instances of smallest $(ratio_small*100)% of Delta, <br>"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "GradT="*"$(temperature_gradient_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])",
        xaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=b1,
            tickvals=b1
        ),
        yaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=b2,
            tickvals=b2
        ),
        xaxis_title="Beta 2, Theta $(theta_compare)°",
        yaxis_title="Beta 1, Theta 180°",
        autosize=false
        
    )

    

   

    #Linear fit
    x=[]
    y=[]
    w=[]
    for i in 1:size(NormalizedCountsb1b2,1)
        for j in 1:size(NormalizedCountsb1b2,2)
            #display(NormalizedCounts_permutedb1b2[i,j])
            #display(b1[i])
            #display(b2[j])
            if NormalizedCountsb1b2[i,j]>0.0
                append!(y,b1[i])
                append!(x,b2[j])
                #print(NormalizedCounts_permutedb1b2[i,j])
                append!(w,1 ./ NormalizedCountsb1b2[i,j] .^2)
            end
        end
    end
    display("x:")
    display(x)
    display("y")
    display(y)
    display("w:")
    display(w)
    
    p0=[1,0.1]

    fit=LsqFit.curve_fit(linear_f,x,y,w,p0)

    cf=LsqFit.coef(fit)

    display("cf:")
    display(cf)

    #display("cf[1] .+ cf[2]*x= $(cf[1]*x .+ cf[2])")

    #Prepare plot
    x=[minimum(b1),maximum(b1)]
    y=cf[1]*x .+ cf[2]

    df_plot = DataFrames.DataFrame(
        x=x,
        y=y,
    )

    Delta2_traceb1b2_line=PlotlyJS.scatter(
        df_plot,
        x=:x, 
        y=:y, 
        mode="lines",
        line = PlotlyJS.attr(color = "green")
        ) 

    display("finished")

    Plot_Heatb1b2=PlotlyJS.plot([Delta2_traceb1b2,Delta2_traceb1b2_line],Delta2_layoutb1b2)
# ------------------------------------------------


    #NormalizedCounts_permutedt1t2=permutedims(NormalizedCountst1t2,[2,1])
    #display(NormalizedCounts_permutedt1t2)

    # Plot_Heatt1t2
    Delta2_tracet1t2=PlotlyJS.heatmap(
        x=t1,
        y=t2,
        z=NormalizedCountst1t2,
        vline=0
    )
    

    Delta2_layoutt1t2 = PlotlyJS.Layout(
        title=
            "Number of instances of smallest $(ratio_small*100)% of Delta, <br>"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "GradT="*"$(temperature_gradient_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])" *", "*
            "B2-B1=0.1",
        xaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=t1,
            tickvals=t1
        ),
        yaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=t2,
            tickvals=t2
        ),
        xaxis_title="MaxT 2, Theta $(theta_compare)°",
        yaxis_title="MaxT 1, Theta 180°",
        autosize=false
        
    )

    Plot_Heatt1t2=PlotlyJS.plot(Delta2_tracet1t2,Delta2_layoutt1t2)


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
        *"_HeatmapB"
        *plot_filename_end)

    plot_filename_D = (
        plot_filename_start
        *"_HeatmapT"
        *plot_filename_end)

    
    plot_total_path_B=plot_save_path*plot_filename_B
    plot_total_path_C=plot_save_path*plot_filename_C
    plot_total_path_D=plot_save_path*plot_filename_D

    
    Plots.savefig(Plot_H2,plot_total_path_B)
    PlotlyJS.savefig(Plot_Heatb1b2,plot_total_path_C)
    PlotlyJS.savefig(Plot_Heatt1t2,plot_total_path_D)


end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216],
    maximal_temperature_array=[0.1,0.125,0.15,0.175,0.2], 
    bond_bending_const_array=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[110.0,180.0],                 #change 110 to 100
    nr_trials_per_temperature=1,
    theta_compare=110.0,                                    #change 110 to 100
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_rad_",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_r_1_"
)
