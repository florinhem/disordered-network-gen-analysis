


###

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
    shape_array)

    @assert length(nr_vertices_array)===length(shape_array)

    color_array=range(Colors.colorant"blue", stop=Colors.colorant"red", length=length(maximal_temperature_array))
    println(color_array)

    @assert length(maximal_temperature_array)===length(color_array)

    P=Plots.scatter(
        title="Bond length std vs Keating energy per vertex"
    )

    for k in eachindex(nr_vertices_array)

        nr_vertices=nr_vertices_array[k]
        shape=shape_array[k]

        for j in eachindex(maximal_temperature_array)

            maximal_temperature=maximal_temperature_array[j]
            color=color_array[j]

            for i in 1:nr_trials_per_temperature

                for m in eachindex(bond_bending_const_array)

                    bond_bending_const=bond_bending_const_array[m]

                    println("$nr_vertices"*", "*"$maximal_temperature"*", "*"$i"*", "*"$bond_bending_const" )

                    save_path = raw".\my_networks\multiple_gml\\"
                    filename = ("multiple_gml_10"
                        *"_N="*"$nr_vertices"
                        *"_T="*"$maximal_temperature"
                        *"_trial="*"$i"
                        *"_beta="*"$bond_bending_const"
                        *".gml")
                    total_path=save_path*filename

                    spatial_network=NG.load_spatial_network_from_gml(total_path)

                    #bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)

                    total_energy_keating=NG.get_total_energy_keating(spatial_network)
                    energy_keating_per_vertex=total_energy_keating/spatial_network[]["nr_vertices"]


                    nr_angles_per_vertex=spatial_network[]["coordination_nr"]*(spatial_network[]["coordination_nr"]-1)/2
                    println("---")
                    println("E:"*"$energy_keating_per_vertex")
                    println("beta:"*"$bond_bending_const")
                    #println("angles:"*"$nr_angles_per_vertex")
                    println(energy_keating_per_vertex/(2+bond_bending_const*nr_angles_per_vertex))
                    println()
                    #=
                    P=Plots.scatter!(
                        P,
                        [energy_keating_per_vertex],
                        [bond_length_std],
                        xlabel="Keating energy per vertex",
                        ylabel="Bond length std",
                        markercolor=color,
                        markershape=shape,
                        label=false,
                        cbar=true,
                        show=true)
                    =#
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
    #=
    plot_save_path = raw".\my_networks\multiple_gml\\"
    plot_filename = ("multiple_gml_23"
        *"_N="*"$minimum_nr_vertices" * "-" * "$maximum_nr_vertices"
        *"_T="*"$minimum_temperature" * "-" * "$maximum_temperature"
        *"_trials="*"$nr_trials_per_temperature"
        *"_beta="*"$minimum_bond_bending_const" * "-" * "$maximum_bond_bending_const"
        *".png")
    plot_total_path=plot_save_path*plot_filename

    
    Plots.savefig(P,plot_total_path)
    =#
end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216,512],
    maximal_temperature_array=[0.1,0.2],
    nr_trials_per_temperature=1,
    bond_bending_const_array=[0.21,0.285,0.36],
    shape_array=[:circle,:rect]
)




###

using Plots

# Sample data
x = 1:10
y = rand(10)
z = rand(10)

# Create a scatter plot with a colorbar at the top
scatter = Plots.scatter(x, y, z, zcolor=z, size=(400, 300), label="", colorbar=true)

# Adjust the layout to position the colorbar at the top
Plots.plot!(scatter, layout=(2, 1), size=(600, 400), colorbar=:left)

###

include("structure_analysis_modules.jl")
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

#path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\data_old\test\216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_evolution.h5"
path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\\"

filename=raw"multiple_p_quench_false_theta_array_N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=180.0_GradT=0.1_StepsPerT=0.01"

#spatial_network = NG.load_spatial_network_from_gml(path*filename)

evolution_dict = GU.load_h5_dict(path*filename*"_evolution.h5")

println(sum(evolution_dict["move_accepted_vec"]))


###

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

                                if(total_path*".gml" in path_array #&& 
                                    #filename!="m_BTMC_N=216_T=0.17_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" &&
                                    #filename!="m_BTMC_q_t__N=216_T=0.1_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=180.0_Trial=1" &&
                                    #filename!="m_BTMC_q_t__N=216_T=0.1_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" #&&
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
                Delta2Matrix[i,j,1]=df[i,2]-df[j,2]     #Delta beta of network 110° minus 180°
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

    ratio_small=0.025
    number_of_smallest=ceil(Int,ratio_small*size(df,1)^2)
    max=maximum(first(Delta2Matrix[:,:,3] |> vec |> sort, number_of_smallest), dims = 1)
    #display("max: $max")
    Delta2Sort=Delta2Matrix[Delta2Matrix[:,:,3].<max,:]
    display(Delta2Sort)

    x_mean=Statistics.mean(Delta2Sort[:,1])
    y_mean=Statistics.mean(Delta2Sort[:,2])
    x_std=Statistics.std(Delta2Sort[:,1])
    y_std=Statistics.std(Delta2Sort[:,2])
   
    # HISTOGRAM
    
    # plot
    Plot_H2=Plots.histogram(
        Delta2Sort[:,3],
        title=LaTeXStrings.LaTeXString(
            "Distribution of smallest $(ratio_small*100)% of delta \n\$^\\textrm{"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "GradT="*"$(temperature_gradient_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])"*
            "}\$"),
        xlabel=LaTeXStrings.LaTeXString("\$ \\delta \$"),
        ylabel="Number of networks",
        legend=false
    )

    #display(Plot_H2)

    # HEATMAP

    fontsize=30

    # calculate counts
    Delta2Round=zeros(size(Delta2Sort))

    for i in 1:size(Delta2Sort,1)
        for j in 1:size(Delta2Sort,2)
            Delta2Round[i,j]=round(Delta2Sort[i,j],digits=3)
        end
    end

    x=Delta2Round[:,1] |> unique |> sort
    y=Delta2Round[:,2] |> unique |> sort
    
    Counts=zeros(size(x,1),size(y,1))
    for k in 1:size(Delta2Round,1)
        T=Delta2Round[k,1]
        B=Delta2Round[k,2]
        for i in 1:size(x,1)
            for j in 1:size(y,1)
                if T===x[i] && B===y[j]
                    Counts[i,j]+=1
                end
            end
        end
    end

    display(Counts)
  
    Counts_permuted=permutedims(Counts,[2,1])
    
    # plot
    Delta2_trace=PlotlyJS.heatmap(
        x=x,
        y=y,
        z=Counts_permuted,
        vline=0,
        colorbar = PlotlyJS.attr(
            tickfont = PlotlyJS.attr(size = fontsize)
        )
    )

    Delta2_layout = PlotlyJS.Layout(
        xaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=x,
            tickvals=x,
            tickfont = PlotlyJS.attr(size = fontsize),
            title = PlotlyJS.attr(
                #text = "Delta beta",
                font = PlotlyJS.attr(size = fontsize)
            ),
            domain = [0.52, 1]
        ),
        yaxis = PlotlyJS.attr(
            tickmode = "array",
            ticktext=y,
            tickvals=y,
            tickfont = PlotlyJS.attr(size = fontsize),
            title = PlotlyJS.attr(
                #text = "Delta MaxT",
                font = PlotlyJS.attr(size = fontsize)
            ),
            scaleanchor = "x",  # Link the scale of the x-axis to the y-axis
            scaleratio = 2
            
        )
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

    PlotlyJS.add_shape!(Plot_Heat2, PlotlyJS.line(
        x0=x_mean-x_std, 
        x1=x_mean+x_std,
        y0=y_mean,
        y1=y_mean,
        line=PlotlyJS.attr(color=:black, width=3),
    ))

    PlotlyJS.add_shape!(Plot_Heat2, PlotlyJS.line(
        x0=x_mean, 
        x1=x_mean,
        y0=y_mean-y_std,
        y1=y_mean+y_std,
        line=PlotlyJS.attr(color=:black, width=3),
    ))

    PlotlyJS.add_shape!(Plot_Heat2, PlotlyJS.circle(
        x0=x_mean-0.025*x_std, 
        x1=x_mean+0.025*x_std,
        y0=y_mean-0.025*y_std,
        y1=y_mean+0.025*y_std,
        line=PlotlyJS.attr(color=:black, width=3),
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
        *"_B="*"$minimum_bond_bending_const" * "-" * "$maximum_bond_bending_const"
        *"_GradT="*"$minimum_temperature_gradient" * "-" * "$maximum_temperature_gradient"
        *"_MC="*"$minimum_nr_monte_carlo_steps_per_temperature" * "-" * "$maximum_nr_monte_carlo_steps_per_temperature"
        *"_Theta="*"$minimum_theta" * "-" * "$maximum_theta"
        *"_i="*"$nr_trials_per_temperature"
        *"_frac="*"$ratio_small"
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
    maximal_temperature_array=[0.1,0.125,0.15,0.175,0.2],
    bond_bending_const_array=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[110.0,180.0],
    nr_trials_per_temperature=1,
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_rad_",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_rad_d5_png_3_"
)



###

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

                                if(total_path*".gml" in path_array #=&& 
                                    filename!="m_BTMC_q_t__N=216_T=0.125_Beta=0.25_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.125_Beta=0.3_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.125_Beta=0.3_GradT=0.1_StepsPerT=0.01_Theta_GS=180.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.15_Beta=0.3_GradT=0.1_StepsPerT=0.01_Theta_GS=180.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.15_Beta=0.25_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.15_Beta=0.25_GradT=0.1_StepsPerT=0.01_Theta_GS=180.0_Trial=1"=#
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

    Delta2Matrix=zeros(size(df,1),size(df,1),5)
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
        end
    end
    display(Delta2Matrix)

    # SCATTER
    Delta2MatrixRound=zeros(size(df,1),size(df,1),5)
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

    #display(b1)
    #display(b2)
    #display(t1)
    #display(t2)
    #display(d)
    

    #ratio_small=0.025
    ratio_small=0.025
    #ratio_small=1
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
            #display("$B1, $B2, $T1, $T2, $D")
            #display(B)
            #display(D)
           
            for k in 1:size(b1,1)
                for m in 1:size(b2,1)
                    if B1===b1[k] && B2===b2[m]
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
                    if T1===t1[k] && T2===t2[m]
                        if D<=max
                            #display("$B1, $B2")
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

    # plot
    Plot_H2=Plots.histogram(
        Delta2Sort[:,5],
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
  
    NormalizedCounts_permutedb1b2=permutedims(NormalizedCountsb1b2,[2,1])
    #display(NormalizedCounts_permutedb1b2)
    
    
    # Plot_Heatb1b2
    Delta2_traceb1b2=PlotlyJS.heatmap(
        x=b1,
        y=b2,
        z=NormalizedCounts_permutedb1b2,
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
        xaxis_title="Beta 1",
        yaxis_title="Beta 2",
        autosize=false
        
    )

    Plot_Heatb1b2=PlotlyJS.plot(Delta2_traceb1b2,Delta2_layoutb1b2)


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

    #println("$x_mean, $y_mean, $x_std, $y_std")

    PlotlyJS.add_shape!(Plot_Heatb1b2, PlotlyJS.line(
        x0=x_mean-x_std, 
        x1=x_mean+x_std,
        y0=y_mean,
        y1=y_mean,
        line=PlotlyJS.attr(color=:red, width=3),
    ))

    PlotlyJS.add_shape!(Plot_Heatb1b2, PlotlyJS.line(
        x0=x_mean, 
        x1=x_mean,
        y0=y_mean-y_std,
        y1=y_mean+y_std,
        line=PlotlyJS.attr(color=:red, width=3),
    ))

    PlotlyJS.add_shape!(Plot_Heatb1b2, PlotlyJS.circle(
        x0=x_mean-0.025*x_std, 
        x1=x_mean+0.025*x_std,
        y0=y_mean-0.025*y_std,
        y1=y_mean+0.025*y_std,
        line=PlotlyJS.attr(color=:red, width=3),
    ))
    



# ------------------------------------------------


    NormalizedCounts_permutedt1t2=permutedims(NormalizedCountst1t2,[2,1])
    #display(NormalizedCounts_permutedt1t2)

    # Plot_Heatt1t2
    Delta2_tracet1t2=PlotlyJS.heatmap(
        x=t1,
        y=t2,
        z=NormalizedCounts_permutedt1t2,
        vline=0
    )
    

    Delta2_layoutt1t2 = PlotlyJS.Layout(
        title=
            "Number of instances of smallest $(ratio_small*100)% of Delta, <br>"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "GradT="*"$(temperature_gradient_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])",
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
        xaxis_title="MaxT 1",
        yaxis_title="MaxT 2",
        autosize=false
        
    )

    Plot_Heatt1t2=PlotlyJS.plot(Delta2_tracet1t2,Delta2_layoutt1t2)


    # statistics
    x_mean=Statistics.mean(Delta2Sort[:,3])
    y_mean=Statistics.mean(Delta2Sort[:,4])
    x_std=Statistics.std(Delta2Sort[:,3])
    y_std=Statistics.std(Delta2Sort[:,4])

    if isnan(x_std)
        x_std=0
    end
    if isnan(y_std)
        y_std=0
    end

    #println("$x_mean, $y_mean, $x_std, $y_std")

    PlotlyJS.add_shape!(Plot_Heatt1t2, PlotlyJS.line(
        x0=x_mean-x_std, 
        x1=x_mean+x_std,
        y0=y_mean,
        y1=y_mean,
        line=PlotlyJS.attr(color=:red, width=3),
    ))

    PlotlyJS.add_shape!(Plot_Heatt1t2, PlotlyJS.line(
        x0=x_mean, 
        x1=x_mean,
        y0=y_mean-y_std,
        y1=y_mean+y_std,
        line=PlotlyJS.attr(color=:red, width=3),
    ))

    PlotlyJS.add_shape!(Plot_Heatt1t2, PlotlyJS.circle(
        x0=x_mean-0.025*x_std, 
        x1=x_mean+0.025*x_std,
        y0=y_mean-0.025*y_std,
        y1=y_mean+0.025*y_std,
        line=PlotlyJS.attr(color=:red, width=3),
    ))


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
    maximal_temperature_array=[0.1,0.125,0.15,0.175,0.2],                       #[0.13,0.14,0.15,0.16,0.17, 0.125,0.135,0.145,0.155,0.165],  
    bond_bending_const_array=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5],         #[0.24,0.26,0.28,0.30,0.32, 0.23,0.25,0.27,0.29,0.31],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[110.0,180.0],
    nr_trials_per_temperature=1,
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_rad_",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_rad_d4_5_"
)


###

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



###


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

    #ratio_small=0.24
    ratio_small=1
    number_of_smallest=ceil(Int,ratio_small*size(df,1)^2)
    max=maximum(first(Delta2Matrix[:,:,3] |> vec |> sort, number_of_smallest), dims = 1)
    #display("max: $max")
    Delta2Sort=Delta2Matrix[Delta2Matrix[:,:,3].<max,:]
    display(Delta2Sort)

    # plot
    Plot_S2=Plots.scatter(
        Delta2Sort[:,1],
        Delta2Sort[:,2],
        title=LaTeXStrings.LaTeXString(
            "Distribution of smallest $(ratio_small*100)% of Delta2 \n\$^\\textrm{"*
            "N="*"$(nr_vertices_array[1])" *", "*
            "GradT="*"$(temperature_gradient_array[1])" *", "*
            "MCSteps="*"$(nr_monte_carlo_steps_per_temperature_array[1])"*
            "}\$"),
        xlabel=LaTeXStrings.LaTeXString("\$ \\Delta \\beta \$"),
        ylabel=LaTeXStrings.LaTeXString("\$ \\Delta T_{max} \$"),
        legend=false
    )

    # statistics
    #display(Delta2Sort[:,1])
    #println(Delta2Sort[:,2])

    x_mean=Statistics.mean(Delta2Sort[:,1])
    y_mean=Statistics.mean(Delta2Sort[:,2])
    x_std=Statistics.std(Delta2Sort[:,1])
    y_std=Statistics.std(Delta2Sort[:,2])
    # plot one red point
    Plots.scatter!(
        [x_mean],
        [y_mean],
        xerr=x_std,
        yerr=y_std,
        color=:red
        )

    #display(Plot_S2)


    # HISTOGRAM
    
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

    # calculate counts
    Delta2Round=zeros(size(Delta2Sort))

    for i in 1:size(Delta2Sort,1)
        for j in 1:size(Delta2Sort,2)
            Delta2Round[i,j]=round(Delta2Sort[i,j],digits=3)
        end
    end

    x=Delta2Round[:,1] |> unique |> sort
    y=Delta2Round[:,2] |> unique |> sort
    
    Counts=zeros(size(x,1),size(y,1))
    for k in 1:size(Delta2Round,1)
        T=Delta2Round[k,1]
        B=Delta2Round[k,2]
        for i in 1:size(x,1)
            for j in 1:size(y,1)
                if T===x[i] && B===y[j]
                    Counts[i,j]+=1
                end
            end
        end
    end

    display(Counts)
  
    Counts_permuted=permutedims(Counts,[2,1])
    
    # plot
    Delta2_trace=PlotlyJS.heatmap(
        x=x,
        y=y,
        z=Counts_permuted,
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

    plot_filename_A = (
        plot_filename_start
        *"_Scatter"
        *plot_filename_end)

    plot_filename_B = (
        plot_filename_start
        *"_Histogram"
        *plot_filename_end)

    plot_filename_C = (
        plot_filename_start
        *"_Heatmap"
        *plot_filename_end)

    plot_total_path_A=plot_save_path*plot_filename_A
    plot_total_path_B=plot_save_path*plot_filename_B
    plot_total_path_C=plot_save_path*plot_filename_C

    Plots.savefig(Plot_S2,plot_total_path_A)
    Plots.savefig(Plot_H2,plot_total_path_B)
    PlotlyJS.savefig(Plot_Heat2,plot_total_path_C)


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
    plot_filename_start = "m_delta2_6_"
)



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
    plot_filename_start = "m_delta_10_"
)


###

# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
import LinearAlgebra
import Statistics

function plot_streching_energy(;
        nr_vertices,
        maximal_temperature,
        bond_bending_const)

    # plot the theoretical and taylor function around the equilibrium
    r_theoretical=collect(0:0.05:1.5)
    r_equilibrium=1
    E_str=3/16 * (r_theoretical.^2 .-r_equilibrium).^2
    E_taylor=3/16 * 4 *(r_theoretical .-r_equilibrium).^2

    Plots.plot(r_theoretical,E_str,label="E_str",legend=:topleft)
    Plots.plot!(r_theoretical,E_taylor,label="E_taylor")
    Plots.plot!([r_equilibrium], seriestype="vline", label="Equilibrium length", color=:blue)

    
    # create the network
    evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices, network_type="diamond", bond_bending_const=bond_bending_const, min_ring_size=3)
    spatial_network = NG.get_periodic_network(evolution_dict)


    # heat and cool down
    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = 
        NA.get_temperature_sequence_heating_cooling_gradient(
            maximal_temperature;
            temperature_gradient = 0.5, 
            nr_monte_carlo_steps_per_temperature = 0.1,
            quench = false)

    evolution_dict["temperature_vec"] = temperature_vec
    evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = nr_monte_carlo_steps_per_temperature_vec

    spatial_network, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(
        spatial_network,
        evolution_dict;
        print_progress = true,
        print_every_nr_attempted_bond_switches = 10)


    # prepare for scatter plotting
    r=[]
    E=[]

    for bond in MetaGraphsNext.edge_labels(spatial_network)
        append!(r,sqrt(spatial_network[bond...]["distance_squared"]))
        append!(E,NG.local_bond_stretching_energy_keating(spatial_network, bond))
    end

    println(r)
    println(E)


    # prepare for std around equilibrium
    r_std=Statistics.std(r)

    Plots.plot!([r_equilibrium-r_std], seriestype="vline", label=false, color=:red)
    Plots.plot!([r_equilibrium+r_std], seriestype="vline", label=false, color=:red)
    
    
    Plots.scatter!(r,E,xlabel="bond length / d",ylabel="streching energy / (α d^2)",label="Measured")

    # save picture
    save_path = raw".\my_networks\E_str\\"
    filename = ("E_str_27_"
        *"_N="*"$nr_vertices"
        *"_T="*"$maximal_temperature"
        *"_beta="*"$bond_bending_const"
        *".png")

    total_path=save_path*filename

    Plots.savefig(total_path)

end






function plot_bending_energy(;
        nr_vertices,
        maximal_temperature,
        bond_bending_const)
    
    # SIMULATED
    # create the network
    evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices, network_type="diamond", bond_bending_const=bond_bending_const, min_ring_size=3)
    spatial_network = NG.get_periodic_network(evolution_dict)


    # heat and cool down
    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = 
        NA.get_temperature_sequence_heating_cooling_gradient(
            maximal_temperature;
            temperature_gradient = 0.5, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = false)

    evolution_dict["temperature_vec"] = temperature_vec
    evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = nr_monte_carlo_steps_per_temperature_vec

    spatial_network, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(
        spatial_network,
        evolution_dict;
        print_progress = false,
        print_every_nr_attempted_bond_switches = 10)

    #NG.plot_spatial_network(spatial_network)

    # prepare for scatter plotting
    θ=[]
    E=[]

    for vertex_label in MetaGraphsNext.labels(spatial_network)

        neighbor_label_vec::Vector{Int64} = collect(MetaGraphsNext.neighbor_labels(
        spatial_network, 
        vertex_label))

        for j in 1:spatial_network[]["coordination_nr"]
            a=sign(neighbor_label_vec[j] - vertex_label).* 
                spatial_network[vertex_label, neighbor_label_vec[j]]["vector"]

            for k in j+1:spatial_network[]["coordination_nr"]
                b=sign(neighbor_label_vec[k] - vertex_label).* 
                    spatial_network[vertex_label, neighbor_label_vec[k]]["vector"]
                
                append!(θ,acos(LinearAlgebra.dot(a,b)/(LinearAlgebra.norm(a)*LinearAlgebra.norm(b))))
                append!(E,3/8*spatial_network[]["bond_bending_const"]*(LinearAlgebra.dot(a,b) + 1/3)^2)
                #=
                if(acos(LinearAlgebra.dot(a,b)/(LinearAlgebra.norm(a)*LinearAlgebra.norm(b)))>2.5)
                    println(acos(LinearAlgebra.dot(a,b)/(LinearAlgebra.norm(a)*LinearAlgebra.norm(b))))
                    println(3/8*spatial_network[]["bond_bending_const"]*(LinearAlgebra.dot(a,b) + 1/3)^2)
                    println(a)
                    println(b)
                    println(LinearAlgebra.dot(a,b))
                end
                =#
            end
        end
    end

    #println(θ)
    #println(E)


    # PLOT
    E_min=minimum(E)
    E_max=maximum(E)
    E_range=E_max-E_min
    E_start=E_min-E_range*0.05
    E_end=E_max+E_range*0.05

    Plots.scatter(
        θ,
        E,
        xlabel="bond angle in rad",
        ylabel="bending energy / (α d^2)",
        ylimits=(E_start,E_end),
        label="Measured", 
        legend=:topleft,
        markersize=3, 
        markerstrokewidth=1, 
        color=:violet)


    # statistics of theta
    θ_mean=Statistics.mean(θ)
    θ_std=Statistics.std(θ)
    
    Plots.plot!([θ_mean], seriestype="vline", label="θ_mean", color=:green)
    Plots.plot!([θ_mean-θ_std], seriestype="vline", label="θ_mean±θ_std", color=:blue)
    Plots.plot!([θ_mean+θ_std], seriestype="vline", label=false, color=:blue)

    
    # THEORY
    # plot the theoretical and taylor function around the equilibrium
    nr_steps=100
    theta_min=minimum(θ)
    theta_max=maximum(θ)
    theta_range=theta_max-theta_min
    theta_step=theta_range/nr_steps
    theta_start=theta_min-theta_range*0.05
    theta_end=theta_max+theta_range*0.05
    theta_theoretical=collect(theta_start:theta_step:theta_end)
    theta_equilibrium=acos(-1/3)
    r_norm=1
    E_bend=3/8 * bond_bending_const * (r_norm.^2*cos.(theta_theoretical) .+1/3).^2
    second_taylor_constant=3/8 * bond_bending_const * 2 * (sin(theta_equilibrium)^2-cos(theta_equilibrium)^2-1/3*cos(theta_equilibrium)) 
    E_taylor=1/2 * second_taylor_constant *(theta_theoretical .- theta_equilibrium).^2

    Plots.plot!(theta_theoretical,E_bend,label="E_bend", color=:red)
    Plots.plot!(theta_theoretical,E_taylor,label="E_taylor", color=:orange)
    Plots.plot!([theta_equilibrium], seriestype="vline", label="θ_eq", color=:yellow)



    # SAVE
    # save picture
    save_path = raw".\my_networks\E_bend\\"
    filename = ("E_bend_31"
        *"_N="*"$nr_vertices"
        *"_T="*"$maximal_temperature"
        *"_beta="*"$bond_bending_const"
        *".png")

    total_path=save_path*filename

    Plots.savefig(total_path)

end


#call functions
#=
plot_streching_energy(
    nr_vertices=216,
    maximal_temperature=0.2,
    bond_bending_const=0.285
)
=#

plot_bending_energy(
    nr_vertices=216,
    maximal_temperature=0.1,
    bond_bending_const=0.85)


###


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
import Statistics

function plot_stretch_vs_bend(;
    nr_vertices_array,
    maximal_temperature_array,
    nr_trials_per_temperature,
    bond_bending_const,
    shape_array)

    @assert length(nr_vertices_array)===length(shape_array)

    color_array=range(Colors.colorant"blue", stop=Colors.colorant"red", length=length(maximal_temperature_array))
    println(color_array)

    @assert length(maximal_temperature_array)===length(color_array)

    P=Plots.scatter(
        title="Stretching energy vs bending energy"
    )

    for k in eachindex(nr_vertices_array)

        nr_vertices=nr_vertices_array[k]
        shape=shape_array[k]

        for j in eachindex(maximal_temperature_array)

            maximal_temperature=maximal_temperature_array[j]
            color=color_array[j]

            for i in 1:nr_trials_per_temperature

                println("$nr_vertices"*", "*"$maximal_temperature"*", "*"$i")

                save_path = raw".\my_networks\multiple_gml\\"
                filename = ("multiple_gml_13"
                    *"_N="*"$nr_vertices"
                    *"_T="*"$maximal_temperature"
                    *"_trial="*"$i"
                    *"_beta="*"$bond_bending_const"
                    *".gml")
                total_path=save_path*filename

                spatial_network=NG.load_spatial_network_from_gml(total_path)


                # E_stretch
                E_stretch=[]
                for bond in MetaGraphsNext.edge_labels(spatial_network)
                    append!(E_stretch,NG.local_bond_stretching_energy_keating(
                        spatial_network, bond))
                end

                # E_bend
                E_bend=[]
                for vertex in MetaGraphsNext.labels(spatial_network)
                    append!(E_bend,NG.local_bond_bending_energy_keating(
                        spatial_network, vertex))
                end

                # Prepare plot
                x=Statistics.mean(E_stretch)
                x_err=Statistics.std(E_stretch)
                y=Statistics.mean(E_bend)
                y_err=Statistics.std(E_bend)

                println(x)
                println(x_err)
                println(y)
                println(y_err)
                println("Ratio: "*"$(y_err/x_err)")

                P=Plots.scatter!(
                    P,
                    [x],
                    [y],
                    xerr=x_err,
                    yerr=y_err,
                    xlabel="Stretching energy with 1 std errorbar",
                    ylabel="Bending energy with 1 std errorbar",
                    markercolor=color,
                    markershape=shape,
                    label=false,
                    cbar=true,
                    show=true,
                    aspect_ratio = :equal
                    )

            end
        end
    end

    minimum_temperature=minimum(maximal_temperature_array)
    maximum_temperature=maximum(maximal_temperature_array)
    
    minimum_nr_vertices=minimum(nr_vertices_array)
    maximum_nr_vertices=maximum(nr_vertices_array)

    plot_save_path = raw".\my_networks\stretch_vs_bend\\"
    plot_filename = ("stretch_vs_bend_6"
        *"_N="*"$minimum_nr_vertices" * "-" * "$maximum_nr_vertices"
        *"_T="*"$minimum_temperature" * "-" * "$maximum_temperature"
        *"_trials="*"$nr_trials_per_temperature"
        *"_beta="*"$bond_bending_const"
        *".png")
    plot_total_path=plot_save_path*plot_filename

    
    Plots.savefig(P,plot_total_path)
end


plot_stretch_vs_bend(
    nr_vertices_array=[216,512],
    maximal_temperature_array=[0.1,0.2,0.3,0.4],
    nr_trials_per_temperature=1,
    bond_bending_const=0.285,
    shape_array=[:circle,:rect]
)

###


# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
import LinearAlgebra
import Statistics

function plot_streching_energy(;
    filename,
    characteristics,
    r_equilibrium)

    # load spatial network
    path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\\"
    type=raw".gml"

    spatial_network = NG.load_spatial_network_from_gml(path*filename*characteristics*type)


    # prepare for scatter plotting
    r=[]
    E=[]

    for bond in MetaGraphsNext.edge_labels(spatial_network)
        append!(r,sqrt(spatial_network[bond...]["distance_squared"]))
        append!(E,NG.local_bond_stretching_energy_keating(spatial_network, bond))
    end


    # THEORY
    # plot the theoretical and taylor function around the equilibrium
    nr_steps=100
    length_min=minimum(r)
    length_max=maximum(r)
    length_range=length_max-length_min
    length_step=length_range/nr_steps
    length_start=length_min-length_range*0.05
    length_end=length_max+length_range*0.05
    length_theoretical=collect(length_start:length_step:length_end)


    # plot histogram
    Plots.histogram!(
        r, 
        label="Measured",
        xlabel="bond length / d",
        ylabel="streching energy / (α d^2)", 
        normalize=:pdf,
        color=:violet,
        alpha = 0.3
        )

    # prepare for std around equilibrium
    r_mean=Statistics.mean(r)
    r_std=Statistics.std(r)

    # plot gaussian over histogram
    fit_r_E=1/(sqrt(2*pi)*r_std) .* exp.(-1/2 .* ((length_theoretical .- r_mean) ./ r_std) .^2)

    Plots.plot!(
        length_theoretical,
        fit_r_E,
        label="Gaussian",
        color=:red
        )

    Plots.plot!([r_equilibrium], seriestype="vline", label="Equilibrium length", color=:blue)
    Plots.plot!([r_mean], seriestype="vline", label=false, color=:green)
    Plots.plot!([r_mean-r_std], seriestype="vline", label=false, color=:blue)
    Plots.plot!([r_mean+r_std], seriestype="vline", label=false, color=:blue)
    

    

    


    # save picture
    save_path = raw".\simulations\metric_E_str\\"
    save_filename = ("metric_together_histo_E_str_4_"
        *characteristics
        *".png")

    save_total_path=save_path*save_filename

    Plots.savefig(save_total_path)

end






function plot_bending_energy(;
    filename,
    characteristics,
    theta_equilibrium)
    
    # load spatial network
    path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\\"
    type=raw".gml"

    spatial_network = NG.load_spatial_network_from_gml(path*filename*characteristics*type)
    bond_bending_const=spatial_network[]["bond_bending_const"]

    # prepare for scatter plotting
    θ=[]
    E=[]

    for vertex_label in MetaGraphsNext.labels(spatial_network)

        neighbor_label_vec::Vector{Int64} = collect(MetaGraphsNext.neighbor_labels(
        spatial_network, 
        vertex_label))

        for j in 1:spatial_network[]["coordination_nr"]
            a=sign(neighbor_label_vec[j] - vertex_label).* 
                spatial_network[vertex_label, neighbor_label_vec[j]]["vector"]

            for k in j+1:spatial_network[]["coordination_nr"]
                b=sign(neighbor_label_vec[k] - vertex_label).* 
                    spatial_network[vertex_label, neighbor_label_vec[k]]["vector"]
                
                append!(θ,acos(LinearAlgebra.dot(a,b)/(LinearAlgebra.norm(a)*LinearAlgebra.norm(b))))
                append!(E,3/8*bond_bending_const*(LinearAlgebra.dot(a,b) + 1/3)^2)
            end
        end
    end

 
    # THEORY
    # plot the theoretical and taylor function around the equilibrium
    nr_steps=100
    theta_min=minimum(θ)
    theta_max=maximum(θ)
    theta_range=theta_max-theta_min
    theta_step=theta_range/nr_steps
    theta_start=theta_min-theta_range*0.05
    theta_end=theta_max+theta_range*0.05
    theta_theoretical=collect(theta_start:theta_step:theta_end)

    
    # plot histogram
    Plots.histogram!(
        θ, 
        label="Measured",
        xlabel="bond angle in rad",
        ylabel="bending energy / (α d^2)", 
        normalize=:pdf,
        color=:violet,
        alpha = 0.3
        )


    # statistics of theta
    θ_mean=Statistics.mean(θ)
    θ_std=Statistics.std(θ)

    # plot gaussian over histogram
    fit_θ_E=1/(sqrt(2*pi)*θ_std) .* exp.(-1/2 .* ((theta_theoretical .- θ_mean) ./ θ_std) .^2)

    Plots.plot!(
        theta_theoretical,
        fit_θ_E,
        label="Gaussian",
        color=:red
        )

    Plots.plot!([theta_equilibrium], seriestype="vline", label="Equilibrium length", color=:blue)
    Plots.plot!([θ_mean], seriestype="vline", label="θ_mean", color=:green)
    Plots.plot!([θ_mean-θ_std], seriestype="vline", label="θ_mean±θ_std", color=:blue)
    Plots.plot!([θ_mean+θ_std], seriestype="vline", label=false, color=:blue)


    # SAVE
    # save picture
    save_path = raw".\simulations\metric_E_bend\\"
    save_filename = ("metric_together_histo_E_bend_4_"
        *characteristics
        *".png")

    save_total_path=save_path*save_filename

    Plots.savefig(save_total_path)

end



function together_hist_str()
    plot_streching_energy(;
        filename=raw"multiple_p_quench_false_theta_array_",
        characteristics=raw"N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=110.0_GradT=0.1_StepsPerT=0.01",
        r_equilibrium=1
    )

    plot_streching_energy(;
        filename=raw"multiple_p_quench_false_theta_array_",
        characteristics=raw"N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=180.0_GradT=0.1_StepsPerT=0.01",
        r_equilibrium=1
    )
end

function together_hist_bend()
    plot_bending_energy(;
        filename=raw"multiple_p_quench_false_theta_array_",
        characteristics=raw"N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=110.0_GradT=0.1_StepsPerT=0.01",
        theta_equilibrium=110.0/360.0*2*pi
    )

    plot_bending_energy(;
        filename=raw"multiple_p_quench_false_theta_array_",
        characteristics=raw"N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=180.0_GradT=0.1_StepsPerT=0.01",
        theta_equilibrium=180.0/360.0*2*pi
    )
end


#call functions
#=
plot_streching_energy(;
    filename=raw"multiple_p_quench_false_theta_array_",
    characteristics=raw"N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=110.0_GradT=0.1_StepsPerT=0.01",
    r_equilibrium=1
)

#Plots.closeall()


plot_bending_energy(;
    filename=raw"multiple_p_quench_false_theta_array_",
    characteristics=raw"N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=110.0_GradT=0.1_StepsPerT=0.01",
    theta_equilibrium=110.0/360.0*2*pi
)
=#

Plots.plot()
#together_hist_str()
together_hist_bend()



###


# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
import .Threads

function save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array,
    maximal_temperature_array,
    nr_trials_per_temperature,
    bond_bending_const_array)
    
    println(Threads.nthreads())

    Threads.@threads for (nr_vertices,maximal_temperature,i,bond_bending_const) in collect(
        Iterators.product(nr_vertices_array,maximal_temperature_array,1:nr_trials_per_temperature,bond_bending_const_array))
                
        println("$nr_vertices"*", "*"$maximal_temperature"*", "*"$i"*", "*"$bond_bending_const" )
        

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices, network_type="diamond", bond_bending_const=bond_bending_const, min_ring_size=3)
        spatial_network = NG.get_periodic_network(evolution_dict)
    
        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = 
            NA.get_temperature_sequence_heating_cooling_gradient(
                maximal_temperature;
                temperature_gradient = 10.0, 
                nr_monte_carlo_steps_per_temperature = 2/(18*216),
                quench = false)

        evolution_dict["temperature_vec"] = temperature_vec
        evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = nr_monte_carlo_steps_per_temperature_vec

        spatial_network, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(
            spatial_network,
            evolution_dict;
            print_progress = true,
            print_every_nr_attempted_bond_switches = 10000)

        save_path = raw".\my_networks\multiple_gml\\"
        filename = ("multiple_gml_13"
            *"_N="*"$nr_vertices"
            *"_T="*"$maximal_temperature"
            *"_trial="*"$i"
            *"_beta="*"$bond_bending_const")

    

        NG.save_spatial_network_to_gml(
            spatial_network,
            filename;
            evolution_dict = evolution_dict,
            save_path = save_path)
                    
             
    end
end


save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array=[216,512],
    maximal_temperature_array=[0.1,0.2,0.3,0.4],
    nr_trials_per_temperature=1,
    bond_bending_const_array=[0.285]
)
#=
nr_vertices_array=[216,512],
    maximal_temperature_array=[0.1,0.2,0.4,0.8,1.6],
    nr_trials_per_temperature=1,
    bond_bending_const_array=[0.21,0.285,0.36]
=#
#=
nr_vertices_array=[216],
    maximal_temperature_array=[1.6],
    nr_trials_per_temperature=1,
    bond_bending_const_array=[0.285]
=#
#=
nr_vertices_array=[216,512],
    maximal_temperature_array=[0.1,0.2,0.3,0.4],
    nr_trials_per_temperature=1,
    bond_bending_const_array=[0.285]
=#


###

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
                                    filename!="m_BTMC_q_t__N=216_T=0.1_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.17_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=0.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.205_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=0.0_Trial=1"
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



    df.Shape = @. ifelse(df.Theta == 180.0, :rect, 
                    ifelse(df.Theta == 110.0, :diamond, 
                    ifelse(df.Theta == 90.0, :circle, :x))); df

    df.GradTLog10 = log10.(df.GradT)
    df.AcceptedMovesLog10 = log10.(df.AcceptedMoves .+ 1)
   
    #println(df)

    fontsize=20
    legend_fontsize=fontsize
    legendtitle_fontsize=legend_fontsize
    colorbar_title_fontsize=legend_fontsize
    guidefontsize=fontsize+5
    
    

    # plot

    A=@df df Plots.scatter(
        left_margin=5Plots.PlotMeasures.mm,
        layout = (2,2),
        size = (900, 700),
        xtickfont = font(fontsize),  # Set x-axis tick font size
        ytickfont = font(fontsize),   # Set y-axis tick font size
        xlim=[-0.005,0.18+0.005],
        xticks = (
            0:0.05:0.15, 
            ["0","0.05","0.1", "0.15"]
        ),
        xlabel=LaTeXStrings.L"\sigma_\mathrm{length} / d",
        ylabel=LaTeXStrings.L"\sigma_\mathrm{angle} / \mathrm{rad}",
        #xguidefont = font(fontsize*1.5),
        #yguidefont = font(fontsize*1.5),
        guidefont=font(guidefontsize),
        legendfontsize = legend_fontsize,
        legendtitlefontsize = legendtitle_fontsize,
        legend=:bottomright,
        #legendborderwidth = 0,
        legendborderalpha = 0,
        legendspacing = 0,  
        legendmargin = 0,
        legendmarkersize =5
        )

    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        color_palette=[:blue,:turquoise1,:yellow,:orange,:red],
        legendtitle =LaTeXStrings.L"T_\mathrm{max}",
        
        group=:MaxT,
        subplot=1)

    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        color_palette=[:purple3,:royalblue,:turquoise1, :green,:greenyellow,:yellow,:gold,:coral,:firebrick,:maroon],
        legendtitle =LaTeXStrings.L"\beta",
        group=:Beta,
        subplot=2)

    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        color_palette=[:grey80,:gray20],
        legendtitle =LaTeXStrings.L"\theta_\mathrm{eq}",
        group=:Theta,
        subplot=3)
    
    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        colorbar_title=LaTeXStrings.L"\mathrm{Log_{10}} \left( N_\mathrm{Acc} \right)",
        #colorbar_titlefont = font(colorbar_title_fontsize),
        colorbar_titlefontsize = colorbar_title_fontsize,
        zcolor=:AcceptedMovesLog10,
        color=cgrad([:blue,:red]),
        colorbar=:top,
        legend=false,
        xlim=[-0.005,0.1375],
        xticks = (
            0:0.025:0.125, 
            ["0","0.025","0.05","0.075","0.1", "0.125"]
        ),
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

    plot_total_path=plot_save_path*plot_filename
    Plots.savefig(A,plot_total_path)

end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216],
    maximal_temperature_array=[0.1,0.125,0.15,0.175,0.2],
    bond_bending_const_array=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[110.0,180.0],
    nr_trials_per_temperature=1, 
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_rad_",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_rad_ma_png_2_"
)


###

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
                                    filename!="m_BTMC_q_t__N=216_T=0.1_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.17_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=0.0_Trial=1" &&
                                    filename!="m_BTMC_q_t__N=216_T=0.205_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=0.0_Trial=1"
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



    df.Shape = @. ifelse(df.Theta == 180.0, :rect, 
                    ifelse(df.Theta == 110.0, :diamond, 
                    ifelse(df.Theta == 90.0, :circle, :x))); df

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
    theta_ground_state_array=[0.0,90.0,110.0,180.0],
    nr_trials_per_temperature=1, 
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_BTMC_q_t_",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_BTMC_matrix_q_t_8_"
)



###

# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Colors
import Glob
import DataFrames
import LaTeXStrings
using StatsPlots


function scatter_plot_for_multiple_gml(;
    data,
    plot_save_path,
    plot_filename_start
    )

    matrix=mapreduce(permutedims, vcat, data)

    df=DataFrames.DataFrame(matrix,:auto)

    DataFrames.rename!(df, [:x1, :x2, :x3, :x4, :x5, :x6, :x7] .=>  [:Beta, :GradT, :MCsteps, :Theta, :BondLenghtStd, :BondAngleStd, :AcceptedMoves])

    df.BondLenghtStdDiff.=0.0
    df.BondAngleStdDiff.=0.0

    for (i, row) in enumerate( eachrow( df ) ) 
        println(i)
        
        if i<size(df,1)
            if ((df[i,"Beta"]===df[i+1,"Beta"]) && 
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
            elseif ((df[i,"Beta"]===df[i-1,"Beta"]) && 
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
            else
                println("Beta, GradT, MCsteps not the same OR Theta the same")
            end
        else
            if ((df[i,"Beta"]===df[i-1,"Beta"]) && 
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


    df.Shape = @. ifelse(df.Theta == 110.0, :diamond, :rect); df
    df.GradTLog10 = log10.(df.GradT)
    df.MCstepsLog10 = log10.(df.MCsteps)
    df.AcceptedMovesLog10 = log10.(df.AcceptedMoves .+ 1)
    df.BondLenghtStdDiffLog10 = log10.(df.BondLenghtStdDiff .+ 1)
    df.BondAngleStdDiffLog10 = log10.(df.BondAngleStdDiff .+ 1)
    df.SurfaceLog10 = log10.(df.MCsteps ./ df.GradT)
   

    println(df)

    plot_title="Bond length std vs Bond angle std"

    A=@df df Plots.scatter(
        #plot_title=plot_title,
        plot_title=LaTeXStrings.LaTeXString("Bond length std vs Bond angle std\n\$^\\textrm{"*
            "N=216, "*
            "MaxT=0.1"*
            "}\$"),
        top_margin=5Plots.PlotMeasures.mm,
        titlefont = font(6),
        layout = (2,2),
        size = (900, 700),
        colorbar=:top,
        legend=false,
        xlabel="Bond length std / d", 
        ylabel="Bond angle std / rad",
        right_margin=5Plots.PlotMeasures.mm
        )

    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        colorbar_title="Beta",
        zcolor=:Beta,
        color=cgrad(:roma, rev=true),
        shape=:Shape,
        subplot=1)

    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        colorbar_title="GradTLog10",
        zcolor=:GradTLog10,
        color=cgrad(:roma, rev=true),
        shape=:Shape,
        subplot=2)

    @df df StatsPlots.scatter!(
        :BondLenghtStd, 
        :BondAngleStd, 
        colorbar_title=" \nMCstepsLog10",
        zcolor=:MCstepsLog10,
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

    plot_total_path=(plot_save_path
        *plot_filename_start
        *plot_title
        *".png")

    Plots.savefig(A,plot_total_path)


    
    # DIFF



    plot_title="Bond length std diff vs Bond angle std diff"

    B=@df df Plots.plot(
        plot_title=plot_title,
        titlefont = font(6),
        layout = (2,2),
        size = (1000, 700))

    @df df StatsPlots.scatter!(
        :BondLenghtStdDiffLog10, 
        :BondAngleStdDiffLog10, 
        colorbar_title="Beta",
        colorbar_tickfontsize = 5,
        zcolor=:Beta,
        color=cgrad(:roma, rev=true),
        #shape=:Shape,
        legend=false,
        subplot=1)
    
    @df df StatsPlots.scatter!(
        :BondLenghtStdDiffLog10, 
        :BondAngleStdDiffLog10, 
        colorbar_title="GradTLog10",
        colorbar_tickfontsize = 5,
        zcolor=:GradTLog10,
        color=cgrad(:roma, rev=true),
        #shape=:Shape,
        legend=false,
        subplot=2)

    @df df StatsPlots.scatter!(
        :BondLenghtStdDiffLog10, 
        :BondAngleStdDiffLog10, 
        colorbar_title="MCstepsLog10",
        colorbar_tickfontsize = 5,
        zcolor=:MCstepsLog10,
        color=cgrad(:roma, rev=true),
        #shape=:Shape,
        legend=false,
        subplot=3)
    #=
    @df df StatsPlots.scatter!(
        :BondLenghtStdDiffLog10, 
        :BondAngleStdDiffLog10, 
        colorbar_title="SurfaceLog10",
        colorbar_tickfontsize = 5,
        zcolor=:SurfaceLog10,
        color=cgrad(:roma, rev=true),
        #shape=:Shape,
        legend=false,
        subplot=9)
    =#
    @df df StatsPlots.scatter!(
        :BondLenghtStdDiffLog10, 
        :BondAngleStdDiffLog10, 
        colorbar_title="AcceptedMovesLog10",
        colorbar_tickfontsize = 5,
        zcolor=:AcceptedMovesLog10,
        color=cgrad(:roma, rev=true),
        #shape=:Shape,
        legend=false,
        subplot=4)

    plot_total_path=(plot_save_path
        *plot_filename_start
        *plot_title
        *".png")

    Plots.savefig(B,plot_total_path)

    println("Finished")
end

data=
[
    [0.0,100.0,0.0005,110.0,0.004626424076958542,0.1842324686232803,5],
    [0.0,100.0,0.0005,180.0,0.006219872731802025,0.19531572045827497,5],
    [0.0,100.0,0.001,110.0,0.002651781866545889,0.13868983172684338,3],
    [0.0,100.0,0.001,180.0,0.006496935874876855,0.17819070885422514,3],
    [0.0,100.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.0,100.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.0,100.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.0,100.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.0,10.0,0.0005,110.0,0.016564431365386764,0.3357886773940933,30],
    [0.0,10.0,0.0005,180.0,0.015208028545476635,0.315916934263737,34],
    [0.0,10.0,0.001,110.0,0.016513947408933517,0.3294658524639375,31],
    [0.0,10.0,0.001,180.0,0.01763606051969618,0.31311810020202796,25],
    [0.0,10.0,0.01,110.0,0.017806186490104136,0.2906090079899177,20],
    [0.0,10.0,0.01,180.0,0.019105526533857782,0.3308778653030486,28],
    [0.0,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.0,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.0,1.0,0.0005,110.0,0.04404973843573433,0.5837653164561262,189],
    [0.0,1.0,0.0005,180.0,0.047848462243526464,0.6061327576706804,225],
    [0.0,1.0,0.001,110.0,0.04327525236852168,0.5954030088450593,243],
    [0.0,1.0,0.001,180.0,0.049531646511830674,0.5808467277553878,224],
    [0.0,1.0,0.01,110.0,0.041129378467737454,0.5767390495573411,228],
    [0.0,1.0,0.01,180.0,0.04620543948343673,0.5826473542522266,227],
    [0.0,1.0,0.1,110.0,0.033031100173689994,0.5663175252411269,230],
    [0.0,1.0,0.1,180.0,0.03593422206253258,0.5766154288283338,221],
    [0.0,0.1,0.0005,110.0,0.0027106344363123816,0.7326490257883934,2387],
    [0.0,0.1,0.0005,180.0,0.002412277377183146,0.7205060812564397,2379],
    [0.0,0.1,0.001,110.0,0.001385464020536888,0.709216818897079,2597],
    [0.0,0.1,0.001,180.0,0.002333011017511657,0.6877226942875097,2430],
    [0.0,0.1,0.01,110.0,0.005432903679512262,0.7434525230586142,2477],
    [0.0,0.1,0.01,180.0,0.00048290346645791355,0.7152497656534547,2615],
    [0.0,0.1,0.1,110.0,0.004120183094563972,0.7246143182870131,2428],
    [0.0,0.1,0.1,180.0,0.0032885250965832146,0.7541431799906408,2366],
    [0.0,0.01,0.0005,110.0,0.0008786429158044023,0.7642145863573229,44703],
    [0.0,0.01,0.0005,180.0,0.00025395185425300295,0.7882211001589146,43130],
    [0.0,0.01,0.001,110.0,0.00012004108399042846,0.7245973304443243,44837],
    [0.0,0.01,0.001,180.0,9.678393546421083e-5,0.7469049361889326,43835],
    [0.0,0.01,0.01,110.0,0.00014668194361288202,0.7455120000347742,42024],
    [0.0,0.01,0.01,180.0,0.00012960139826509332,0.7624184285216029,41318],
    [0.0,0.01,0.1,110.0,0.00021467022961704264,0.7588369105610075,43453],
    [0.0,0.01,0.1,180.0,0.0001328901954700292,0.7702589950817739,42427],
    [0.1425,100.0,0.0005,110.0,0.02047604129415642,0.07420478236372059,1],
    [0.1425,100.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.1425,100.0,0.001,110.0,0.007308019219609829,0.045345311307515734,1],
    [0.1425,100.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.1425,100.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.1425,100.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.1425,100.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.1425,100.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.1425,10.0,0.0005,110.0,0.007308461549498563,0.045344608577493464,1],
    [0.1425,10.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.1425,10.0,0.001,110.0,0.02327881820210157,0.11448042302234333,5],
    [0.1425,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.1425,10.0,0.01,110.0,0.023288351057625705,0.09794405750616053,3],
    [0.1425,10.0,0.01,180.0,0.014844207764056361,0.08989230410776343,4],
    [0.1425,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.1425,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.1425,1.0,0.0005,110.0,0.03158638391624253,0.18560174769811133,27],
    [0.1425,1.0,0.0005,180.0,0.049360159966287814,0.27245624962245907,50],
    [0.1425,1.0,0.001,110.0,0.042803170513221,0.22697111578234563,32],
    [0.1425,1.0,0.001,180.0,0.04118171108596398,0.2326692073405377,35],
    [0.1425,1.0,0.01,110.0,0.04771485629620374,0.27677498504785647,49],
    [0.1425,1.0,0.01,180.0,0.033297470019568086,0.1926950730940322,27],
    [0.1425,1.0,0.1,110.0,0.04886058064982134,0.27969327845033365,50],
    [0.1425,1.0,0.1,180.0,0.05325295227228932,0.30343392626829935,69],
    [0.1425,0.1,0.0005,110.0,0.06761717148588525,0.42085984726168807,646],
    [0.1425,0.1,0.0005,180.0,0.0610931869722474,0.4018094076569822,738],
    [0.1425,0.1,0.001,110.0,0.06446351951832385,0.41473581238825025,604],
    [0.1425,0.1,0.001,180.0,0.05510569878043807,0.4069766233201634,739],
    [0.1425,0.1,0.01,110.0,0.06852200106088832,0.40160211778659144,645],
    [0.1425,0.1,0.01,180.0,0.05587910405162278,0.39469319230015,721],
    [0.1425,0.1,0.1,110.0,0.06558753298543431,0.41231790189680695,667],
    [0.1425,0.1,0.1,180.0,0.06159389243388003,0.3962338834271377,742],
    [0.1425,0.01,0.0005,110.0,0.05602784023350045,0.3278374552564045,7368],
    [0.1425,0.01,0.0005,180.0,0.04769321771481373,0.4170531846831608,14183],
    [0.1425,0.01,0.001,110.0,0.055498054903738395,0.39635662586895776,8538],
    [0.1425,0.01,0.001,180.0,0.04512902809229777,0.4206179805355421,13137],
    [0.1425,0.01,0.01,110.0,0.05492614796519643,0.3356235044590366,7093],
    [0.1425,0.01,0.01,180.0,0.046166462050800085,0.4235244702684507,14603],
    [0.1425,0.01,0.1,110.0,0.059933583751347795,0.37881749365965306,7412],
    [0.1425,0.01,0.1,180.0,0.04671570058920258,0.42397394940540123,14395],
    [0.285,100.0,0.0005,110.0,0.011705892412006737,0.04056363651435685,1],
    [0.285,100.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.285,100.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.285,100.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.285,100.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.285,100.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.285,100.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.285,100.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.285,10.0,0.0005,110.0,0.011702070955184442,0.040568775052264396,1],
    [0.285,10.0,0.0005,180.0,0.029437469776391473,0.08266487648536916,2],
    [0.285,10.0,0.001,110.0,0.016590546662329616,0.05698868019862161,2],
    [0.285,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.285,10.0,0.01,110.0,0.011703348516466941,0.04056805604746311,1],
    [0.285,10.0,0.01,180.0,0.028403634179098743,0.08378561156734662,2],
    [0.285,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.285,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.285,1.0,0.0005,110.0,0.040217732779356995,0.1238921893637239,10],
    [0.285,1.0,0.0005,180.0,0.03539363341728674,0.13711782158795421,15],
    [0.285,1.0,0.001,110.0,0.03660418010063556,0.1070937753428414,11],
    [0.285,1.0,0.001,180.0,0.03120278663622889,0.10755904932950623,9],
    [0.285,1.0,0.01,110.0,0.04449793850429103,0.14689269835879576,15],
    [0.285,1.0,0.01,180.0,0.01163509233289244,0.040303175985411145,5],
    [0.285,1.0,0.1,110.0,0.061071970579252156,0.1862434548816739,25],
    [0.285,1.0,0.1,180.0,0.04843567795192965,0.1387945541846145,14],
    [0.285,0.1,0.0005,110.0,0.08486020409282612,0.3102926716375029,282],
    [0.285,0.1,0.0005,180.0,0.06574057857495988,0.2857489561595977,244],
    [0.285,0.1,0.001,110.0,0.06494759314446726,0.2231226796642656,196],
    [0.285,0.1,0.001,180.0,0.06849132153902895,0.24964922632212364,194],
    [0.285,0.1,0.01,110.0,0.07664576225411741,0.2518578310339166,209],
    [0.285,0.1,0.01,180.0,0.0716564866839261,0.2757548619684896,256],
    [0.285,0.1,0.1,110.0,0.08253087155250682,0.29526005939548233,202],
    [0.285,0.1,0.1,180.0,0.05330993831986317,0.19905924244958623,165],
    [0.285,0.01,0.0005,110.0,0.08070515049430328,0.2892999978079795,2969],
    [0.285,0.01,0.0005,180.0,0.07004497988867374,0.2591668517864453,4214],
    [0.285,0.01,0.001,110.0,0.07634384860633442,0.2960962912230227,2906],
    [0.285,0.01,0.001,180.0,0.06005795730958101,0.27357726727921255,4051],
    [0.285,0.01,0.01,110.0,0.07992957646801474,0.28205207410319255,2745],
    [0.285,0.01,0.01,180.0,0.06555689367189359,0.2634426233136765,3788],
    [0.285,0.01,0.1,110.0,0.08064548912882502,0.2896146752489243,2715],
    [0.285,0.01,0.1,180.0,0.07181036277054521,0.26433459669573456,3876],
    [0.5,100.0,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,100.0,0.0005,180.0,0.016473105058986626,0.03674549212509805,1],
    [0.5,100.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,100.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,100.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,100.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,100.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,100.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,10.0,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,10.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,10.0,0.01,180.0,0.016471659456487667,0.036746296304094116,1],
    [0.5,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,1.0,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,1.0,0.0005,180.0,0.016601687033206268,0.03656641598611299,2],
    [0.5,1.0,0.001,110.0,0.03290584048045969,0.07386914070326657,4],
    [0.5,1.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [0.5,1.0,0.01,110.0,0.016512100404591133,0.036790961724148384,1],
    [0.5,1.0,0.01,180.0,0.02386772254084038,0.052641887785486466,2],
    [0.5,1.0,0.1,110.0,0.04675068877360048,0.1079381368676503,3],
    [0.5,1.0,0.1,180.0,0.016474815205035343,0.03674340389264217,1],
    [0.5,0.1,0.0005,110.0,0.031457467827875056,0.06965687145693322,34],
    [0.5,0.1,0.0005,180.0,0.03955289908213162,0.09281099210237856,25],
    [0.5,0.1,0.001,110.0,0.023484224581852717,0.05005219695124578,17],
    [0.5,0.1,0.001,180.0,0.027944775663527197,0.06382038839684076,26],
    [0.5,0.1,0.01,110.0,0.04760498677840435,0.09833911856200878,12],
    [0.5,0.1,0.01,180.0,0.0362169077479976,0.08404560705881536,20],
    [0.5,0.1,0.1,110.0,0.016725555083932746,0.036492002722743735,11],
    [0.5,0.1,0.1,180.0,0.03848838264324312,0.08473368274397011,34],
    [0.5,0.01,0.0005,110.0,0.11445416982648741,0.24591656987149318,806],
    [0.5,0.01,0.0005,180.0,0.0018162665691991663,0.0018307382009918403,441],
    [0.5,0.01,0.001,110.0,0.12020870097181274,0.27280831099412545,987],
    [0.5,0.01,0.001,180.0,0.06997580461803161,0.1992933877280043,1283],
    [0.5,0.01,0.01,110.0,0.12255980843825075,0.2710680959060793,877],
    [0.5,0.01,0.01,180.0,0.06243526777597644,0.1613134229901907,798],
    [0.5,0.01,0.1,110.0,0.14136327103924257,0.2635646087035402,896],
    [0.5,0.01,0.1,180.0,0.08167883887500126,0.21668508878520087,1018],
    [1.0,100.0,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,100.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,100.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,100.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,100.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,100.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,100.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,100.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,10.0,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,10.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,10.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,1.0,0.0005,110.0,0.023412421950557625,0.03324851870787081,1],
    [1.0,1.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,1.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,1.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,1.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,1.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,1.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,1.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.1,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.1,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.1,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.1,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.1,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.1,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.1,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.1,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.01,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.01,0.0005,180.0,0.002131903769235298,0.0010152835088053076,9],
    [1.0,0.01,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.01,0.001,180.0,0.0024833483923236313,0.0011590762698844561,6],
    [1.0,0.01,0.01,110.0,0.0035760558000548495,0.0018621618683401734,9],
    [1.0,0.01,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
    [1.0,0.01,0.1,110.0,0.0032856543734298025,0.0016443314804323857,7],
    [1.0,0.01,0.1,180.0,0.0023281943151224004,0.0010259545568706058,17]
]

scatter_plot_for_multiple_gml(;
    data,    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_BTMC_sub_2_"
)


###

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
    path_array=Glob.glob(filename_start*"*",save_path)

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
                                    *"_Beta="*"$bond_bending_const"
                                    *"_GradT="*"$temperature_gradient"
                                    *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
                                    *"_Theta_GS="*"$theta_ground_state"
                                    *"_Trial="*"$i"
                                    )

                                total_path=save_path*filename

                                if(total_path*".gml" in path_array)
                                    #println("scatter done")

                                    spatial_network=NG.load_spatial_network_from_gml(total_path*".gml")

                                    bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)
                                    bond_angle_std, bond_angle_vec = NA.get_bond_angle_std(spatial_network)

                                    evolution_dict = GU.load_h5_dict(total_path*"_evolution.h5")
                                    accepted_moves=sum(evolution_dict["move_accepted_vec"])

                                    println("["
                                        #*"$nr_vertices"*","
                                        #*"$maximal_temperature"*","
                                        *"$bond_bending_const"*","
                                        *"$temperature_gradient"*","
                                        *"$nr_monte_carlo_steps_per_temperature"*","
                                        *"$theta_ground_state"*","
                                        #*"$i"*","
                                        *"$bond_length_std"*","
                                        *"$bond_angle_std"*","
                                        *"$accepted_moves"
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
    bond_bending_const_array=[0,0.1425,0.285,0.5,1],
    temperature_gradient_array=[0.01,0.1,1,10,100], #[100,10,1,0.1,0.01],
    nr_monte_carlo_steps_per_temperature_array=[0.0005,0.001,0.01,0.1],
    theta_ground_state_array=[110.0,180.0],
    nr_trials_per_temperature=1,   
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_BTMC",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_BTMC_T_1__",
    markershape_array=[:circle,:rect],                              #[:circle,:rect],       #[:circle,:rect#=,:star5,:cross,:+=#],
    markersize_array=1.5 .*[1,2,3,4],                             #1.5 .*[1,2,3,4,5],
    markerstrokewidth_array=[1,2,3,4,5]                             #[1,2,3,4,5]
)


###


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
    shape_array,
    plot_metric)

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

    if plot_metric==="stretching"
        title="Bond length std vs Keating energy per vertex"
        ylabel="Bond length std"
    elseif plot_metric==="bending"
        title="Bond angle std vs Keating energy per vertex"
        ylabel="Bond angle std"
    elseif plot_metric==="dihedral"
        title="Dihedral angle std vs Keating energy per vertex"
        ylabel="Dihedral angle std"
    end

    P=Plots.scatter(
        title=title,
        xlabel="Keating energy per vertex",
        ylabel=ylabel
    )

    for k in eachindex(nr_vertices_array)

        nr_vertices=nr_vertices_array[k]

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


                        if plot_metric==="stretching"
                            y, bond_length_vec = NA.get_bond_length_std(spatial_network)
                            
                        elseif plot_metric==="bending"
                            y, bond_angle_vec = NA.get_bond_angle_std(spatial_network)
                        elseif plot_metric==="dihedral"
                            y, dihedral_angle_vec = NA.get_dihedral_angle_std(spatial_network)
                        end
                        #bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)
                        #bond_angle_std, bond_angle_vec = NA.get_bond_angle_std(spatial_network)
                        #dihedral_angle_std, dihedral_angle_vec = get_dihedral_angle_std(spatial_network)


                        total_energy_keating=NG.get_total_energy_keating(spatial_network)
                        energy_keating_per_vertex=total_energy_keating/spatial_network[]["nr_vertices"]

                        P=Plots.scatter!(
                            P,
                            [energy_keating_per_vertex],
                            [y],
                            markercolor=color,
                            markershape=shape,
                            legend=false,
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
        *"_"*plot_metric
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
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "multiple_p_quench_false_theta_array",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "multiple_SBD_1",
    shape_array=[:circle,:rect],
    #plot_metric="stretching" 
    #plot_metric="bending"
    plot_metric="dihedral"
)




###
#Pairplot

#%%

import PyPlot 
using PyCall
import Pandas
import DataFrames

@pyimport seaborn as sns
@pyimport pandas as pd 
@pyimport matplotlib.pyplot as plt 

sns.set_theme(style="ticks")

# data
x=[



[0.0,100.0,0.0005,110.0,0.004626424076958542,0.1842324686232803,5],
[0.0,100.0,0.0005,180.0,0.006219872731802025,0.19531572045827497,5],
[0.0,100.0,0.001,110.0,0.002651781866545889,0.13868983172684338,3],
[0.0,100.0,0.001,180.0,0.006496935874876855,0.17819070885422514,3],
[0.0,100.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.0,100.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.0,100.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.0,100.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.0,10.0,0.0005,110.0,0.016564431365386764,0.3357886773940933,30],
[0.0,10.0,0.0005,180.0,0.015208028545476635,0.315916934263737,34],
[0.0,10.0,0.001,110.0,0.016513947408933517,0.3294658524639375,31],
[0.0,10.0,0.001,180.0,0.01763606051969618,0.31311810020202796,25],
[0.0,10.0,0.01,110.0,0.017806186490104136,0.2906090079899177,20],
[0.0,10.0,0.01,180.0,0.019105526533857782,0.3308778653030486,28],
[0.0,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.0,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.0,1.0,0.0005,110.0,0.04404973843573433,0.5837653164561262,189],
[0.0,1.0,0.0005,180.0,0.047848462243526464,0.6061327576706804,225],
[0.0,1.0,0.001,110.0,0.04327525236852168,0.5954030088450593,243],
[0.0,1.0,0.001,180.0,0.049531646511830674,0.5808467277553878,224],
[0.0,1.0,0.01,110.0,0.041129378467737454,0.5767390495573411,228],
[0.0,1.0,0.01,180.0,0.04620543948343673,0.5826473542522266,227],
[0.0,1.0,0.1,110.0,0.033031100173689994,0.5663175252411269,230],
[0.0,1.0,0.1,180.0,0.03593422206253258,0.5766154288283338,221],
[0.0,0.1,0.0005,110.0,0.0027106344363123816,0.7326490257883934,2387],
[0.0,0.1,0.0005,180.0,0.002412277377183146,0.7205060812564397,2379],
[0.0,0.1,0.001,110.0,0.001385464020536888,0.709216818897079,2597],
[0.0,0.1,0.001,180.0,0.002333011017511657,0.6877226942875097,2430],
[0.0,0.1,0.01,110.0,0.005432903679512262,0.7434525230586142,2477],
[0.0,0.1,0.01,180.0,0.00048290346645791355,0.7152497656534547,2615],
[0.0,0.1,0.1,110.0,0.004120183094563972,0.7246143182870131,2428],
[0.0,0.1,0.1,180.0,0.0032885250965832146,0.7541431799906408,2366],
[0.0,0.01,0.0005,110.0,0.0008786429158044023,0.7642145863573229,44703],
[0.0,0.01,0.0005,180.0,0.00025395185425300295,0.7882211001589146,43130],
[0.0,0.01,0.001,110.0,0.00012004108399042846,0.7245973304443243,44837],
[0.0,0.01,0.001,180.0,9.678393546421083e-5,0.7469049361889326,43835],
[0.0,0.01,0.01,110.0,0.00014668194361288202,0.7455120000347742,42024],
[0.0,0.01,0.01,180.0,0.00012960139826509332,0.7624184285216029,41318],
[0.0,0.01,0.1,110.0,0.00021467022961704264,0.7588369105610075,43453],
[0.0,0.01,0.1,180.0,0.0001328901954700292,0.7702589950817739,42427],
[0.1425,100.0,0.0005,110.0,0.02047604129415642,0.07420478236372059,1],
[0.1425,100.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.1425,100.0,0.001,110.0,0.007308019219609829,0.045345311307515734,1],
[0.1425,100.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.1425,100.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.1425,100.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.1425,100.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.1425,100.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.1425,10.0,0.0005,110.0,0.007308461549498563,0.045344608577493464,1],
[0.1425,10.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.1425,10.0,0.001,110.0,0.02327881820210157,0.11448042302234333,5],
[0.1425,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.1425,10.0,0.01,110.0,0.023288351057625705,0.09794405750616053,3],
[0.1425,10.0,0.01,180.0,0.014844207764056361,0.08989230410776343,4],
[0.1425,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.1425,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.1425,1.0,0.0005,110.0,0.03158638391624253,0.18560174769811133,27],
[0.1425,1.0,0.0005,180.0,0.049360159966287814,0.27245624962245907,50],
[0.1425,1.0,0.001,110.0,0.042803170513221,0.22697111578234563,32],
[0.1425,1.0,0.001,180.0,0.04118171108596398,0.2326692073405377,35],
[0.1425,1.0,0.01,110.0,0.04771485629620374,0.27677498504785647,49],
[0.1425,1.0,0.01,180.0,0.033297470019568086,0.1926950730940322,27],
[0.1425,1.0,0.1,110.0,0.04886058064982134,0.27969327845033365,50],
[0.1425,1.0,0.1,180.0,0.05325295227228932,0.30343392626829935,69],
[0.1425,0.1,0.0005,110.0,0.06761717148588525,0.42085984726168807,646],
[0.1425,0.1,0.0005,180.0,0.0610931869722474,0.4018094076569822,738],
[0.1425,0.1,0.001,110.0,0.06446351951832385,0.41473581238825025,604],
[0.1425,0.1,0.001,180.0,0.05510569878043807,0.4069766233201634,739],
[0.1425,0.1,0.01,110.0,0.06852200106088832,0.40160211778659144,645],
[0.1425,0.1,0.01,180.0,0.05587910405162278,0.39469319230015,721],
[0.1425,0.1,0.1,110.0,0.06558753298543431,0.41231790189680695,667],
[0.1425,0.1,0.1,180.0,0.06159389243388003,0.3962338834271377,742],
[0.1425,0.01,0.0005,110.0,0.05602784023350045,0.3278374552564045,7368],
[0.1425,0.01,0.0005,180.0,0.04769321771481373,0.4170531846831608,14183],
[0.1425,0.01,0.001,110.0,0.055498054903738395,0.39635662586895776,8538],
[0.1425,0.01,0.001,180.0,0.04512902809229777,0.4206179805355421,13137],
[0.1425,0.01,0.01,110.0,0.05492614796519643,0.3356235044590366,7093],
[0.1425,0.01,0.01,180.0,0.046166462050800085,0.4235244702684507,14603],
[0.1425,0.01,0.1,110.0,0.059933583751347795,0.37881749365965306,7412],
[0.1425,0.01,0.1,180.0,0.04671570058920258,0.42397394940540123,14395],
[0.285,100.0,0.0005,110.0,0.011705892412006737,0.04056363651435685,1],
[0.285,100.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.285,100.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.285,100.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.285,100.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.285,100.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.285,100.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.285,100.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.285,10.0,0.0005,110.0,0.011702070955184442,0.040568775052264396,1],
[0.285,10.0,0.0005,180.0,0.029437469776391473,0.08266487648536916,2],
[0.285,10.0,0.001,110.0,0.016590546662329616,0.05698868019862161,2],
[0.285,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.285,10.0,0.01,110.0,0.011703348516466941,0.04056805604746311,1],
[0.285,10.0,0.01,180.0,0.028403634179098743,0.08378561156734662,2],
[0.285,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.285,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.285,1.0,0.0005,110.0,0.040217732779356995,0.1238921893637239,10],
[0.285,1.0,0.0005,180.0,0.03539363341728674,0.13711782158795421,15],
[0.285,1.0,0.001,110.0,0.03660418010063556,0.1070937753428414,11],
[0.285,1.0,0.001,180.0,0.03120278663622889,0.10755904932950623,9],
[0.285,1.0,0.01,110.0,0.04449793850429103,0.14689269835879576,15],
[0.285,1.0,0.01,180.0,0.01163509233289244,0.040303175985411145,5],
[0.285,1.0,0.1,110.0,0.061071970579252156,0.1862434548816739,25],
[0.285,1.0,0.1,180.0,0.04843567795192965,0.1387945541846145,14],
[0.285,0.1,0.0005,110.0,0.08486020409282612,0.3102926716375029,282],
[0.285,0.1,0.0005,180.0,0.06574057857495988,0.2857489561595977,244],
[0.285,0.1,0.001,110.0,0.06494759314446726,0.2231226796642656,196],
[0.285,0.1,0.001,180.0,0.06849132153902895,0.24964922632212364,194],
[0.285,0.1,0.01,110.0,0.07664576225411741,0.2518578310339166,209],
[0.285,0.1,0.01,180.0,0.0716564866839261,0.2757548619684896,256],
[0.285,0.1,0.1,110.0,0.08253087155250682,0.29526005939548233,202],
[0.285,0.1,0.1,180.0,0.05330993831986317,0.19905924244958623,165],
[0.285,0.01,0.0005,110.0,0.08070515049430328,0.2892999978079795,2969],
[0.285,0.01,0.0005,180.0,0.07004497988867374,0.2591668517864453,4214],
[0.285,0.01,0.001,110.0,0.07634384860633442,0.2960962912230227,2906],
[0.285,0.01,0.001,180.0,0.06005795730958101,0.27357726727921255,4051],
[0.285,0.01,0.01,110.0,0.07992957646801474,0.28205207410319255,2745],
[0.285,0.01,0.01,180.0,0.06555689367189359,0.2634426233136765,3788],
[0.285,0.01,0.1,110.0,0.08064548912882502,0.2896146752489243,2715],
[0.285,0.01,0.1,180.0,0.07181036277054521,0.26433459669573456,3876],
[0.5,100.0,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,100.0,0.0005,180.0,0.016473105058986626,0.03674549212509805,1],
[0.5,100.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,100.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,100.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,100.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,100.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,100.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,10.0,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,10.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,10.0,0.01,180.0,0.016471659456487667,0.036746296304094116,1],
[0.5,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,1.0,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,1.0,0.0005,180.0,0.016601687033206268,0.03656641598611299,2],
[0.5,1.0,0.001,110.0,0.03290584048045969,0.07386914070326657,4],
[0.5,1.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[0.5,1.0,0.01,110.0,0.016512100404591133,0.036790961724148384,1],
[0.5,1.0,0.01,180.0,0.02386772254084038,0.052641887785486466,2],
[0.5,1.0,0.1,110.0,0.04675068877360048,0.1079381368676503,3],
[0.5,1.0,0.1,180.0,0.016474815205035343,0.03674340389264217,1],
[0.5,0.1,0.0005,110.0,0.031457467827875056,0.06965687145693322,34],
[0.5,0.1,0.0005,180.0,0.03955289908213162,0.09281099210237856,25],
[0.5,0.1,0.001,110.0,0.023484224581852717,0.05005219695124578,17],
[0.5,0.1,0.001,180.0,0.027944775663527197,0.06382038839684076,26],
[0.5,0.1,0.01,110.0,0.04760498677840435,0.09833911856200878,12],
[0.5,0.1,0.01,180.0,0.0362169077479976,0.08404560705881536,20],
[0.5,0.1,0.1,110.0,0.016725555083932746,0.036492002722743735,11],
[0.5,0.1,0.1,180.0,0.03848838264324312,0.08473368274397011,34],
[0.5,0.01,0.0005,110.0,0.11445416982648741,0.24591656987149318,806],
[0.5,0.01,0.0005,180.0,0.0018162665691991663,0.0018307382009918403,441],
[0.5,0.01,0.001,110.0,0.12020870097181274,0.27280831099412545,987],
[0.5,0.01,0.001,180.0,0.06997580461803161,0.1992933877280043,1283],
[0.5,0.01,0.01,110.0,0.12255980843825075,0.2710680959060793,877],
[0.5,0.01,0.01,180.0,0.06243526777597644,0.1613134229901907,798],
[0.5,0.01,0.1,110.0,0.14136327103924257,0.2635646087035402,896],
[0.5,0.01,0.1,180.0,0.08167883887500126,0.21668508878520087,1018],
[1.0,100.0,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,100.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,100.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,100.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,100.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,100.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,100.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,100.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,10.0,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,10.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,10.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,1.0,0.0005,110.0,0.023412421950557625,0.03324851870787081,1],
[1.0,1.0,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,1.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,1.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,1.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,1.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,1.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,1.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.1,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.1,0.0005,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.1,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.1,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.1,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.1,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.1,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.1,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.01,0.0005,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.01,0.0005,180.0,0.002131903769235298,0.0010152835088053076,9],
[1.0,0.01,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.01,0.001,180.0,0.0024833483923236313,0.0011590762698844561,6],
[1.0,0.01,0.01,110.0,0.0035760558000548495,0.0018621618683401734,9],
[1.0,0.01,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15,0],
[1.0,0.01,0.1,110.0,0.0032856543734298025,0.0016443314804323857,7],
[1.0,0.01,0.1,180.0,0.0023281943151224004,0.0010259545568706058,17]
]

# data cleaning

matrix=mapreduce(permutedims, vcat, x)
jdf=DataFrames.DataFrame(matrix,:auto)

DataFrames.rename!(jdf, [:x1, :x2, :x3, :x4, :x5, :x6, :x7] .=>  [:Beta, :GradT, :MCsteps, :Theta, :BondLenghtStd, :BondAngleStd, :AcceptedMoves])

jdf.BetaLog2 = log10.(jdf.Beta)
jdf.GradTLog10 = log10.(jdf.GradT)
jdf.MCstepsLog10 = log10.(jdf.MCsteps)
jdf.AcceptedMovesLog10 = log10.(jdf.AcceptedMoves)

pdf=Pandas.DataFrame(jdf)

# plot

vars=[:BetaLog2, :GradTLog10, :MCstepsLog10, :AcceptedMovesLog10, :BondLenghtStd, :BondAngleStd]
g = sns.PairGrid(pdf,hue="Theta",palette=sns.color_palette("tab10"),vars = vars)
g.add_legend() #title="Theta",fontsize=14, bbox_to_anchor=(1.5,1)
g.map_upper(sns.scatterplot)
g.map_lower(sns.regplot)
g.map_diag(sns.histplot)

# save

path=raw".\simulations\analysis_plot\\"
filename="pairgrid_19"
total_path=path*filename

plt.savefig(total_path)











#%%


using DataFrames, StatsPlots

df = DataFrames.DataFrame(randn(100, 3), [:Tenure, :Balance, :CreditScore])

N = DataFrames.ncol(df)

@df df StatsPlots.histogram(cols(1:N); layout=N)

# Try to make a bond length angle plot, but with always the color tells us what changes => 8 subplots?

#%%

using DataFrames, StatsPlots

df = DataFrames.DataFrame(randn(100, 3), [:Tenure, :Balance, :CreditScore])

N = DataFrames.ncol(df)

@df df StatsPlots.histogram(cols(1:N); layout=N)











#%%

import PyPlot 
using PyCall
import Pandas
import DataFrames

@pyimport seaborn as sns
@pyimport pandas as pd 
@pyimport matplotlib.pyplot as plt 

sns.set_theme(style="ticks")

# data
x=[



    [0.135,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.01,110.0,0.06401338557732181,0.4060359551967697],
    [0.135,0.1,0.01,180.0,0.05601616236847719,0.41441432448271637],
    [0.135,0.1,0.1,110.0,0.061021038271114415,0.4017501099339115],
    [0.135,0.1,0.1,180.0,0.05790225273741632,0.4082923625999096],
    [0.135,0.1,1.0,110.0,0.05897650192383255,0.4104274400477412],
    [0.135,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,0.001,110.0,0.042772092411191245,0.2491108770309085],
    [0.135,1.0,0.001,180.0,0.03516704002128376,0.231927891539188],
    [0.135,1.0,0.01,110.0,0.04534523097905104,0.24346082832206614],
    [0.135,1.0,0.01,180.0,0.039450907125217774,0.23283584896706738],
    [0.135,1.0,0.1,110.0,0.0513763408470471,0.2777466848434963],
    [0.135,1.0,0.1,180.0,0.04617737510212192,0.24376432020623962],
    [0.135,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.001,110.0,0.029486503383955923,0.11636377503522224],
    [0.135,10.0,0.001,180.0,0.023680147947620137,0.11761914640609153],
    [0.135,10.0,0.01,110.0,0.014582277814301785,0.09048097276452079],
    [0.135,10.0,0.01,180.0,0.009589334343224454,0.06305994133486469],
    [0.135,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.1,180.0,0.05963159674363879,0.3277816597244015],
    [0.21,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.1,0.001,110.0,0.07091891548046385,0.3404057617828189],
    [0.21,0.1,0.01,110.0,0.06813133697371909,0.332960711051228],
    [0.21,0.1,0.01,180.0,0.06285080077101268,0.34238991375194816],
    [0.21,0.1,0.1,110.0,0.07515351783491563,0.34872603045153594],
    [0.21,0.1,0.1,180.0,0.062093685718115334,0.3473716182448396],
    [0.21,0.1,1.0,110.0,0.08030152475293591,0.35514811068133506],
    [0.21,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,0.001,110.0,0.03946255851224546,0.16048730987901297],
    [0.21,1.0,0.001,180.0,0.03778667986180233,0.1367606380556472],
    [0.21,1.0,0.01,110.0,0.04524079219157094,0.18863155314768074],
    [0.21,1.0,0.01,180.0,0.024659510937369695,0.11332417786930574],
    [0.21,1.0,0.1,110.0,0.0526742357804406,0.22416645793245832],
    [0.21,1.0,0.1,180.0,0.04430625683546078,0.1806697362975857],
    [0.21,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.001,110.0,0.00955822655835486,0.042743082489513244],
    [0.21,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.01,110.0,0.009556540659477994,0.042745390095866115],
    [0.21,10.0,0.01,180.0,0.009437557881115406,0.04244814254695355],
    [0.21,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.1,180.0,0.06468369957102785,0.2760068630870718],
    [0.285,0.01,1.0,110.0,0.08249820918238371,0.2834942741953725],
    [0.285,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.1,0.001,110.0,0.08001153603871601,0.26792380065206794],
    [0.285,0.1,0.001,180.0,0.06427847485375804,0.2377611647289092],
    [0.285,0.1,0.01,110.0,0.08700429074430302,0.30641276227284037],
    [0.285,0.1,0.01,180.0,0.06400031417499383,0.2579974253553776],
    [0.285,0.1,0.1,110.0,0.06981030610951264,0.2507019069839894],
    [0.285,0.1,0.1,180.0,0.05274836646038725,0.21840089946952432],
    [0.285,0.1,1.0,110.0,0.08500959026300058,0.3121860789865292],
    [0.285,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,0.001,110.0,0.016482797871082706,0.055094444432573136],
    [0.285,1.0,0.001,180.0,0.023607762935759048,0.08424307494241305],
    [0.285,1.0,0.01,110.0,0.016420057462119685,0.056514604215375255],
    [0.285,1.0,0.01,180.0,0.025751312624126475,0.08904247069917819],
    [0.285,1.0,0.1,110.0,0.06062750788622491,0.2006928547974821],
    [0.285,1.0,0.1,180.0,0.03256236861693596,0.10678637292109114],
    [0.285,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.01,110.0,0.020195763557142046,0.06849297493772595],
    [0.285,10.0,0.01,180.0,0.011609533140198082,0.04044463858766386],
    [0.285,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.1,180.0,0.08360229710952409,0.2627644110643584],
    [0.36,0.01,1.0,110.0,0.09575810585223143,0.2764219440731812],
    [0.36,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.1,0.001,110.0,0.0731359938273293,0.20293821935155515],
    [0.36,0.1,0.001,180.0,0.04864516244293532,0.15224415553899875],
    [0.36,0.1,0.01,110.0,0.08784682233469933,0.22944925414224968],
    [0.36,0.1,0.01,180.0,0.05985313415684717,0.2066815944504365],
    [0.36,0.1,0.1,110.0,0.06036738605467378,0.17233075664898892],
    [0.36,0.1,0.1,180.0,0.039467762407835665,0.11549374892544934],
    [0.36,0.1,1.0,110.0,0.09624535372177485,0.2799101501149262],
    [0.36,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.001,110.0,0.03731551435349158,0.10456722494280188],
    [0.36,1.0,0.001,180.0,0.013492498381847583,0.03887475717575668],
    [0.36,1.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.01,180.0,0.013574192524228739,0.038738277328439936],
    [0.36,1.0,0.1,110.0,0.03895254476528791,0.10640954516459276],
    [0.36,1.0,0.1,180.0,0.028280826832015982,0.08144372491319456],
    [0.36,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.01,180.0,0.07879594971601717,0.23396920638186364],
    [0.435,0.01,0.1,110.0,0.11177413494982522,0.26781213984153257],
    [0.435,0.01,0.1,180.0,0.07982332216222571,0.23482174230900055],
    [0.435,0.01,1.0,110.0,0.10404336410294569,0.2653633751069886],
    [0.435,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.1,0.001,110.0,0.0514064157848007,0.11829627953551232],
    [0.435,0.1,0.001,180.0,0.024863487764014657,0.06315739336480754],
    [0.435,0.1,0.01,110.0,0.025841337025231764,0.06486046199439431],
    [0.435,0.1,0.01,180.0,0.03492304739777853,0.08697522836176026],
    [0.435,0.1,0.1,110.0,0.03738125408829439,0.09143436097440404],
    [0.435,0.1,0.1,180.0,0.03297336040443855,0.08218000709150271],
    [0.435,0.1,1.0,110.0,0.08652066536889622,0.21970131477371324],
    [0.435,0.1,1.0,180.0,0.05947106152287292,0.16182636931731134],
    [0.435,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.001,110.0,0.036602559829868035,0.08971152601751382],
    [0.435,1.0,0.001,180.0,0.0013058232073252154,0.0018719903367372267],
    [0.435,1.0,0.01,110.0,0.015227611796538492,0.03750268386910032],
    [0.435,1.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.1,110.0,0.01996487793673578,0.04930644716403539],
    [0.435,1.0,0.1,180.0,0.029502838607312336,0.07606843198185774],
    [0.435,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.001,110.0,0.0152253136763512,0.03767997016309791],
    [0.435,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15]
]

# data cleaning

matrix=mapreduce(permutedims, vcat, x)
jdf=DataFrames.DataFrame(matrix,:auto)

DataFrames.rename!(jdf, [:x1, :x2, :x3, :x4, :x5, :x6] .=>  [:Beta, :GradT, :MCsteps, :Theta, :BondLenghtStd, :BondAngleStd])

jdf.GradTLog10 = log10.(jdf.GradT)
jdf.MCstepsLog10 = log10.(jdf.MCsteps)

pdf=Pandas.DataFrame(jdf)

# plot

vars=[:Beta, :GradTLog10, :MCstepsLog10, :BondLenghtStd, :BondAngleStd]
g = sns.PairGrid(pdf,hue="Theta",palette=sns.color_palette("tab10"),vars = vars)
g.add_legend(title="Theta",fontsize=14, bbox_to_anchor=(1.5,1))
g.map_upper(sns.histplot)
g.map_lower(sns.regplot)
g.map_diag(sns.histplot)

# save

path=raw".\simulations\analysis_plot\\"
filename="pairgrid_12"
total_path=path*filename

plt.savefig(total_path)























#%%


import PyPlot 
using PyCall
import Pandas

@pyimport seaborn as sns
@pyimport pandas as pd 

sns.set_theme(style="ticks")

# data
x=[

    [0.135,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.01,110.0,0.06401338557732181,0.4060359551967697],
    [0.135,0.1,0.01,180.0,0.05601616236847719,0.41441432448271637],
    [0.135,0.1,0.1,110.0,0.061021038271114415,0.4017501099339115],
    [0.135,0.1,0.1,180.0,0.05790225273741632,0.4082923625999096],
    [0.135,0.1,1.0,110.0,0.05897650192383255,0.4104274400477412],
    [0.135,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,0.001,110.0,0.042772092411191245,0.2491108770309085],
    [0.135,1.0,0.001,180.0,0.03516704002128376,0.231927891539188],
    [0.135,1.0,0.01,110.0,0.04534523097905104,0.24346082832206614],
    [0.135,1.0,0.01,180.0,0.039450907125217774,0.23283584896706738],
    [0.135,1.0,0.1,110.0,0.0513763408470471,0.2777466848434963],
    [0.135,1.0,0.1,180.0,0.04617737510212192,0.24376432020623962],
    [0.135,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.001,110.0,0.029486503383955923,0.11636377503522224],
    [0.135,10.0,0.001,180.0,0.023680147947620137,0.11761914640609153],
    [0.135,10.0,0.01,110.0,0.014582277814301785,0.09048097276452079],
    [0.135,10.0,0.01,180.0,0.009589334343224454,0.06305994133486469],
    [0.135,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.1,180.0,0.05963159674363879,0.3277816597244015],
    [0.21,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.1,0.001,110.0,0.07091891548046385,0.3404057617828189],
    [0.21,0.1,0.01,110.0,0.06813133697371909,0.332960711051228],
    [0.21,0.1,0.01,180.0,0.06285080077101268,0.34238991375194816],
    [0.21,0.1,0.1,110.0,0.07515351783491563,0.34872603045153594],
    [0.21,0.1,0.1,180.0,0.062093685718115334,0.3473716182448396],
    [0.21,0.1,1.0,110.0,0.08030152475293591,0.35514811068133506],
    [0.21,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,0.001,110.0,0.03946255851224546,0.16048730987901297],
    [0.21,1.0,0.001,180.0,0.03778667986180233,0.1367606380556472],
    [0.21,1.0,0.01,110.0,0.04524079219157094,0.18863155314768074],
    [0.21,1.0,0.01,180.0,0.024659510937369695,0.11332417786930574],
    [0.21,1.0,0.1,110.0,0.0526742357804406,0.22416645793245832],
    [0.21,1.0,0.1,180.0,0.04430625683546078,0.1806697362975857],
    [0.21,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.001,110.0,0.00955822655835486,0.042743082489513244],
    [0.21,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.01,110.0,0.009556540659477994,0.042745390095866115],
    [0.21,10.0,0.01,180.0,0.009437557881115406,0.04244814254695355],
    [0.21,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.1,180.0,0.06468369957102785,0.2760068630870718],
    [0.285,0.01,1.0,110.0,0.08249820918238371,0.2834942741953725],
    [0.285,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.1,0.001,110.0,0.08001153603871601,0.26792380065206794],
    [0.285,0.1,0.001,180.0,0.06427847485375804,0.2377611647289092],
    [0.285,0.1,0.01,110.0,0.08700429074430302,0.30641276227284037],
    [0.285,0.1,0.01,180.0,0.06400031417499383,0.2579974253553776],
    [0.285,0.1,0.1,110.0,0.06981030610951264,0.2507019069839894],
    [0.285,0.1,0.1,180.0,0.05274836646038725,0.21840089946952432],
    [0.285,0.1,1.0,110.0,0.08500959026300058,0.3121860789865292],
    [0.285,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,0.001,110.0,0.016482797871082706,0.055094444432573136],
    [0.285,1.0,0.001,180.0,0.023607762935759048,0.08424307494241305],
    [0.285,1.0,0.01,110.0,0.016420057462119685,0.056514604215375255],
    [0.285,1.0,0.01,180.0,0.025751312624126475,0.08904247069917819],
    [0.285,1.0,0.1,110.0,0.06062750788622491,0.2006928547974821],
    [0.285,1.0,0.1,180.0,0.03256236861693596,0.10678637292109114],
    [0.285,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.01,110.0,0.020195763557142046,0.06849297493772595],
    [0.285,10.0,0.01,180.0,0.011609533140198082,0.04044463858766386],
    [0.285,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.1,180.0,0.08360229710952409,0.2627644110643584],
    [0.36,0.01,1.0,110.0,0.09575810585223143,0.2764219440731812],
    [0.36,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.1,0.001,110.0,0.0731359938273293,0.20293821935155515],
    [0.36,0.1,0.001,180.0,0.04864516244293532,0.15224415553899875],
    [0.36,0.1,0.01,110.0,0.08784682233469933,0.22944925414224968],
    [0.36,0.1,0.01,180.0,0.05985313415684717,0.2066815944504365],
    [0.36,0.1,0.1,110.0,0.06036738605467378,0.17233075664898892],
    [0.36,0.1,0.1,180.0,0.039467762407835665,0.11549374892544934],
    [0.36,0.1,1.0,110.0,0.09624535372177485,0.2799101501149262],
    [0.36,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.001,110.0,0.03731551435349158,0.10456722494280188],
    [0.36,1.0,0.001,180.0,0.013492498381847583,0.03887475717575668],
    [0.36,1.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.01,180.0,0.013574192524228739,0.038738277328439936],
    [0.36,1.0,0.1,110.0,0.03895254476528791,0.10640954516459276],
    [0.36,1.0,0.1,180.0,0.028280826832015982,0.08144372491319456],
    [0.36,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.01,180.0,0.07879594971601717,0.23396920638186364],
    [0.435,0.01,0.1,110.0,0.11177413494982522,0.26781213984153257],
    [0.435,0.01,0.1,180.0,0.07982332216222571,0.23482174230900055],
    [0.435,0.01,1.0,110.0,0.10404336410294569,0.2653633751069886],
    [0.435,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.1,0.001,110.0,0.0514064157848007,0.11829627953551232],
    [0.435,0.1,0.001,180.0,0.024863487764014657,0.06315739336480754],
    [0.435,0.1,0.01,110.0,0.025841337025231764,0.06486046199439431],
    [0.435,0.1,0.01,180.0,0.03492304739777853,0.08697522836176026],
    [0.435,0.1,0.1,110.0,0.03738125408829439,0.09143436097440404],
    [0.435,0.1,0.1,180.0,0.03297336040443855,0.08218000709150271],
    [0.435,0.1,1.0,110.0,0.08652066536889622,0.21970131477371324],
    [0.435,0.1,1.0,180.0,0.05947106152287292,0.16182636931731134],
    [0.435,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.001,110.0,0.036602559829868035,0.08971152601751382],
    [0.435,1.0,0.001,180.0,0.0013058232073252154,0.0018719903367372267],
    [0.435,1.0,0.01,110.0,0.015227611796538492,0.03750268386910032],
    [0.435,1.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.1,110.0,0.01996487793673578,0.04930644716403539],
    [0.435,1.0,0.1,180.0,0.029502838607312336,0.07606843198185774],
    [0.435,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.001,110.0,0.0152253136763512,0.03767997016309791],
    [0.435,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15]
]

matrix=mapreduce(permutedims, vcat, x)
jdf=DataFrames.DataFrame(matrix,:auto)

rename!(jdf, [:x1, :x2, :x3, :x4, :x5, :x6] .=>  [:Beta, :GradT, :MCsteps, :Theta, :BondLenghtStd, :BondAngleStd])

pdf=Pandas.DataFrame(jdf)

# plot
sns_plot=sns.pairplot(pdf, hue="Theta", kind="reg", diag_kind="hist",palette="hls")

# save
path=raw".\simulations\analysis_plot\\"
filename="pairplot_8"
total_path=path*filename
plt.savefig(total_path)
















#%%

using DataFrames, CairoMakie

function pairplot(df; stride=1, colormap=:thermal, resolution=(1200,1200))

    dim = size(df,2) # how many colums there are in the dataframe
    idxs = 1:stride:size(df,1)
    colorant = range(0, 1, length=length(idxs))

    pp_theme = Attributes(
        Axis = (
            aspect = 1,
            topspinevisible = false,
            rightspinevisible = false,
        ),
        Scatter = (
            colormap = colormap, # try :thermal, :darkrainbow
            markersize = 6
        )
    )

    f = with_theme(pp_theme) do
        f = CairoMakie.Figure(size=resolution)

        for i in 1:dim, j in 1:dim

            ax = CairoMakie.Axis(f[i, j])
            CairoMakie.scatter!(df[idxs,j], df[idxs,i], color = colorant)

            if i==dim
                ax.xticklabelsvisible = true
                ax.xlabel = DataFrames.names(df)[j]
            end
            if j==1
                ax.yticklabelsvisible = true
                ax.ylabel = DataFrames.names(df)[i]
            end
        end
        f
    end
end

x=[
    [1,2,3],
    [0.3,0.1,0.5],
    [3,4,1]
]

matrix=mapreduce(permutedims, vcat, x)

#df=convert(DataFrames.DataFrame,matrix)
df=DataFrames.DataFrame(matrix,:auto)


pairplot(df)




#%%




import Seaborn
import Plots
import Images

#dfs = rand(50, 5) * 100;

sbplot = Seaborn.heatmap([1 1.5; 3 2.8])
#sbplot = Seaborn.pairplot(dfs)
Seaborn.savefig("SBplot.png")
sbplot_Plots = Plots.plot(Images.load("SBplot.png")) # Load in and make it a Plots plot object

hist = Plots.bar([1,1.5, 2.,7])

heat_hist_plot = Plots.plot(sbplot_Plots, hist, layout = @layout([A{0.8h}; B]))

#%%

import PyPlot 
using PyCall
import Pandas

@pyimport seaborn as sns
@pyimport pandas as pd 

sns.set_theme(style="ticks")

# data
x=[
    [1,2,3],
    [0.3,0.1,0.5],
    [3,4,1]
]

matrix=mapreduce(permutedims, vcat, x)
jdf=DataFrames.DataFrame(matrix,:auto)
pdf=Pandas.DataFrame(jdf)

# plot
sns_plot=sns.pairplot(pdf)

# save
path=raw".\simulations\analysis_plot\\"
filename="pairplot_1"
total_path=path*filename
plt.savefig(total_path)

#%%


import PyPlot 
import PyCall

PyCall.@pyimport seaborn as sns

x=randn(1000)
plot=sns.distplot(x)
plot

#%%

import CairoMakie
import PairPlots
import DataFrames

#=
N = 100_000
α = [2randn(N÷2) .+ 6; randn(N÷2)]
β = [3randn(N÷2); 2randn(N÷2)]
γ = randn(N)
δ = β .+ 0.6randn(N)

df = DataFrames.DataFrame(;α, β, γ, δ)
=#
x=[

[0.135,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.01,110.0,0.06401338557732181,0.4060359551967697],
    [0.135,0.1,0.01,180.0,0.05601616236847719,0.41441432448271637],
    [0.135,0.1,0.1,110.0,0.061021038271114415,0.4017501099339115],
    [0.135,0.1,0.1,180.0,0.05790225273741632,0.4082923625999096],
    [0.135,0.1,1.0,110.0,0.05897650192383255,0.4104274400477412],
    [0.135,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,0.001,110.0,0.042772092411191245,0.2491108770309085],
    [0.135,1.0,0.001,180.0,0.03516704002128376,0.231927891539188],
    [0.135,1.0,0.01,110.0,0.04534523097905104,0.24346082832206614],
    [0.135,1.0,0.01,180.0,0.039450907125217774,0.23283584896706738],
    [0.135,1.0,0.1,110.0,0.0513763408470471,0.2777466848434963],
    [0.135,1.0,0.1,180.0,0.04617737510212192,0.24376432020623962],
    [0.135,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.001,110.0,0.029486503383955923,0.11636377503522224],
    [0.135,10.0,0.001,180.0,0.023680147947620137,0.11761914640609153],
    [0.135,10.0,0.01,110.0,0.014582277814301785,0.09048097276452079],
    [0.135,10.0,0.01,180.0,0.009589334343224454,0.06305994133486469],
    [0.135,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.1,180.0,0.05963159674363879,0.3277816597244015],
    [0.21,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.1,0.001,110.0,0.07091891548046385,0.3404057617828189],
    [0.21,0.1,0.01,110.0,0.06813133697371909,0.332960711051228],
    [0.21,0.1,0.01,180.0,0.06285080077101268,0.34238991375194816],
    [0.21,0.1,0.1,110.0,0.07515351783491563,0.34872603045153594],
    [0.21,0.1,0.1,180.0,0.062093685718115334,0.3473716182448396],
    [0.21,0.1,1.0,110.0,0.08030152475293591,0.35514811068133506],
    [0.21,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,0.001,110.0,0.03946255851224546,0.16048730987901297],
    [0.21,1.0,0.001,180.0,0.03778667986180233,0.1367606380556472],
    [0.21,1.0,0.01,110.0,0.04524079219157094,0.18863155314768074],
    [0.21,1.0,0.01,180.0,0.024659510937369695,0.11332417786930574],
    [0.21,1.0,0.1,110.0,0.0526742357804406,0.22416645793245832],
    [0.21,1.0,0.1,180.0,0.04430625683546078,0.1806697362975857],
    [0.21,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.001,110.0,0.00955822655835486,0.042743082489513244],
    [0.21,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.01,110.0,0.009556540659477994,0.042745390095866115],
    [0.21,10.0,0.01,180.0,0.009437557881115406,0.04244814254695355],
    [0.21,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.1,180.0,0.06468369957102785,0.2760068630870718],
    [0.285,0.01,1.0,110.0,0.08249820918238371,0.2834942741953725],
    [0.285,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.1,0.001,110.0,0.08001153603871601,0.26792380065206794],
    [0.285,0.1,0.001,180.0,0.06427847485375804,0.2377611647289092],
    [0.285,0.1,0.01,110.0,0.08700429074430302,0.30641276227284037],
    [0.285,0.1,0.01,180.0,0.06400031417499383,0.2579974253553776],
    [0.285,0.1,0.1,110.0,0.06981030610951264,0.2507019069839894],
    [0.285,0.1,0.1,180.0,0.05274836646038725,0.21840089946952432],
    [0.285,0.1,1.0,110.0,0.08500959026300058,0.3121860789865292],
    [0.285,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,0.001,110.0,0.016482797871082706,0.055094444432573136],
    [0.285,1.0,0.001,180.0,0.023607762935759048,0.08424307494241305],
    [0.285,1.0,0.01,110.0,0.016420057462119685,0.056514604215375255],
    [0.285,1.0,0.01,180.0,0.025751312624126475,0.08904247069917819],
    [0.285,1.0,0.1,110.0,0.06062750788622491,0.2006928547974821],
    [0.285,1.0,0.1,180.0,0.03256236861693596,0.10678637292109114],
    [0.285,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.01,110.0,0.020195763557142046,0.06849297493772595],
    [0.285,10.0,0.01,180.0,0.011609533140198082,0.04044463858766386],
    [0.285,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.1,180.0,0.08360229710952409,0.2627644110643584],
    [0.36,0.01,1.0,110.0,0.09575810585223143,0.2764219440731812],
    [0.36,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.1,0.001,110.0,0.0731359938273293,0.20293821935155515],
    [0.36,0.1,0.001,180.0,0.04864516244293532,0.15224415553899875],
    [0.36,0.1,0.01,110.0,0.08784682233469933,0.22944925414224968],
    [0.36,0.1,0.01,180.0,0.05985313415684717,0.2066815944504365],
    [0.36,0.1,0.1,110.0,0.06036738605467378,0.17233075664898892],
    [0.36,0.1,0.1,180.0,0.039467762407835665,0.11549374892544934],
    [0.36,0.1,1.0,110.0,0.09624535372177485,0.2799101501149262],
    [0.36,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.001,110.0,0.03731551435349158,0.10456722494280188],
    [0.36,1.0,0.001,180.0,0.013492498381847583,0.03887475717575668],
    [0.36,1.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.01,180.0,0.013574192524228739,0.038738277328439936],
    [0.36,1.0,0.1,110.0,0.03895254476528791,0.10640954516459276],
    [0.36,1.0,0.1,180.0,0.028280826832015982,0.08144372491319456],
    [0.36,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.01,180.0,0.07879594971601717,0.23396920638186364],
    [0.435,0.01,0.1,110.0,0.11177413494982522,0.26781213984153257],
    [0.435,0.01,0.1,180.0,0.07982332216222571,0.23482174230900055],
    [0.435,0.01,1.0,110.0,0.10404336410294569,0.2653633751069886],
    [0.435,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.1,0.001,110.0,0.0514064157848007,0.11829627953551232],
    [0.435,0.1,0.001,180.0,0.024863487764014657,0.06315739336480754],
    [0.435,0.1,0.01,110.0,0.025841337025231764,0.06486046199439431],
    [0.435,0.1,0.01,180.0,0.03492304739777853,0.08697522836176026],
    [0.435,0.1,0.1,110.0,0.03738125408829439,0.09143436097440404],
    [0.435,0.1,0.1,180.0,0.03297336040443855,0.08218000709150271],
    [0.435,0.1,1.0,110.0,0.08652066536889622,0.21970131477371324],
    [0.435,0.1,1.0,180.0,0.05947106152287292,0.16182636931731134],
    [0.435,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.001,110.0,0.036602559829868035,0.08971152601751382],
    [0.435,1.0,0.001,180.0,0.0013058232073252154,0.0018719903367372267],
    [0.435,1.0,0.01,110.0,0.015227611796538492,0.03750268386910032],
    [0.435,1.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.1,110.0,0.01996487793673578,0.04930644716403539],
    [0.435,1.0,0.1,180.0,0.029502838607312336,0.07606843198185774],
    [0.435,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.001,110.0,0.0152253136763512,0.03767997016309791],
    [0.435,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15]
]

#matrix=reduce(hcat,x)'
matrix=mapreduce(permutedims, vcat, x)

#df=convert(DataFrames.DataFrame,matrix)
df=DataFrames.DataFrame(matrix,:auto)

#PairPlots.pairplot(df)
PairPlots.pairplot(
    df => (
        PairPlots.Scatter(),
        PairPlots.MarginHist()#,
        #PairPlots.MarginConfidenceLimits()
    )
)






#%%

import Plots
import StatsPlots

function pairplot(df, title="Pairplot")
    rows, cols = size(df)
    plots = []
	p_mat = Matrix(undef,cols,cols);
    for i = 1:cols, j = 1:cols
		p_mat[i,j] = i>=j ? (label = :°, blank = false) : (label = :_, blank = true)
        if i>j
			subplot = Plots.scatter(df[:,i], df[:,j], legend = false,markersize=14/cols, 
			alpha=0.4, markerstrokecolor=nothing, grid=nothing, 
			tickfontsize=round(32/cols), tick_direction=:in, 
			xrot=45,showaxis=(i==cols ? (j==1 ? true : :x) : (j==1 ? :y : false)),
			xticks=(i==cols ? :auto : nothing), yticks=(j==1 ? :auto : nothing))
			push!(plots,subplot)
		elseif i==j
			subplot = Plots.histogram(df[:,i], normalize=true, legend=false, alpha=0.4, ticks=nothing, showaxis=(i==cols ? (j==1 ? true : :x) : (j==1 ? :y : false)),xticks=(i==cols ? :auto : nothing), 
				tick_direction=:in, tickfontcolor="white", xrot=45
			)
			StatsPlots.density!(df[:,i], grid=nothing, tickfontsize=round(32/cols), trim=true,
			xlimits=(minimum(df[:,i]),maximum(df[:,i])))
			push!(plots,subplot)
		end
    end
    return Plots.plot(plots..., layout=p_mat, margin=0Plots.PlotMeasures.mm, plot_title=title, plot_titlevspan=0.05)
end

dfs = rand(50, 5) * 100;
pairplot(dfs)

#%%

x=[
    [0.135,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,0.1,0.01,110.0,0.06401338557732181,0.4060359551967697],
    [0.135,0.1,0.01,180.0,0.05601616236847719,0.41441432448271637],
    [0.135,0.1,0.1,110.0,0.061021038271114415,0.4017501099339115],
    [0.135,0.1,0.1,180.0,0.05790225273741632,0.4082923625999096],
    [0.135,0.1,1.0,110.0,0.05897650192383255,0.4104274400477412],
    [0.135,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,0.001,110.0,0.042772092411191245,0.2491108770309085],
    [0.135,1.0,0.001,180.0,0.03516704002128376,0.231927891539188],
    [0.135,1.0,0.01,110.0,0.04534523097905104,0.24346082832206614],
    [0.135,1.0,0.01,180.0,0.039450907125217774,0.23283584896706738],
    [0.135,1.0,0.1,110.0,0.0513763408470471,0.2777466848434963],
    [0.135,1.0,0.1,180.0,0.04617737510212192,0.24376432020623962],
    [0.135,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.001,110.0,0.029486503383955923,0.11636377503522224],
    [0.135,10.0,0.001,180.0,0.023680147947620137,0.11761914640609153],
    [0.135,10.0,0.01,110.0,0.014582277814301785,0.09048097276452079],
    [0.135,10.0,0.01,180.0,0.009589334343224454,0.06305994133486469],
    [0.135,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.135,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.01,0.1,180.0,0.05963159674363879,0.3277816597244015],
    [0.21,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,0.1,0.001,110.0,0.07091891548046385,0.3404057617828189],
    [0.21,0.1,0.01,110.0,0.06813133697371909,0.332960711051228],
    [0.21,0.1,0.01,180.0,0.06285080077101268,0.34238991375194816],
    [0.21,0.1,0.1,110.0,0.07515351783491563,0.34872603045153594],
    [0.21,0.1,0.1,180.0,0.062093685718115334,0.3473716182448396],
    [0.21,0.1,1.0,110.0,0.08030152475293591,0.35514811068133506],
    [0.21,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,0.001,110.0,0.03946255851224546,0.16048730987901297],
    [0.21,1.0,0.001,180.0,0.03778667986180233,0.1367606380556472],
    [0.21,1.0,0.01,110.0,0.04524079219157094,0.18863155314768074],
    [0.21,1.0,0.01,180.0,0.024659510937369695,0.11332417786930574],
    [0.21,1.0,0.1,110.0,0.0526742357804406,0.22416645793245832],
    [0.21,1.0,0.1,180.0,0.04430625683546078,0.1806697362975857],
    [0.21,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.001,110.0,0.00955822655835486,0.042743082489513244],
    [0.21,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.01,110.0,0.009556540659477994,0.042745390095866115],
    [0.21,10.0,0.01,180.0,0.009437557881115406,0.04244814254695355],
    [0.21,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.21,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.01,0.1,180.0,0.06468369957102785,0.2760068630870718],
    [0.285,0.01,1.0,110.0,0.08249820918238371,0.2834942741953725],
    [0.285,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,0.1,0.001,110.0,0.08001153603871601,0.26792380065206794],
    [0.285,0.1,0.001,180.0,0.06427847485375804,0.2377611647289092],
    [0.285,0.1,0.01,110.0,0.08700429074430302,0.30641276227284037],
    [0.285,0.1,0.01,180.0,0.06400031417499383,0.2579974253553776],
    [0.285,0.1,0.1,110.0,0.06981030610951264,0.2507019069839894],
    [0.285,0.1,0.1,180.0,0.05274836646038725,0.21840089946952432],
    [0.285,0.1,1.0,110.0,0.08500959026300058,0.3121860789865292],
    [0.285,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,0.001,110.0,0.016482797871082706,0.055094444432573136],
    [0.285,1.0,0.001,180.0,0.023607762935759048,0.08424307494241305],
    [0.285,1.0,0.01,110.0,0.016420057462119685,0.056514604215375255],
    [0.285,1.0,0.01,180.0,0.025751312624126475,0.08904247069917819],
    [0.285,1.0,0.1,110.0,0.06062750788622491,0.2006928547974821],
    [0.285,1.0,0.1,180.0,0.03256236861693596,0.10678637292109114],
    [0.285,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.01,110.0,0.020195763557142046,0.06849297493772595],
    [0.285,10.0,0.01,180.0,0.011609533140198082,0.04044463858766386],
    [0.285,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.285,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.01,0.1,180.0,0.08360229710952409,0.2627644110643584],
    [0.36,0.01,1.0,110.0,0.09575810585223143,0.2764219440731812],
    [0.36,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,0.1,0.001,110.0,0.0731359938273293,0.20293821935155515],
    [0.36,0.1,0.001,180.0,0.04864516244293532,0.15224415553899875],
    [0.36,0.1,0.01,110.0,0.08784682233469933,0.22944925414224968],
    [0.36,0.1,0.01,180.0,0.05985313415684717,0.2066815944504365],
    [0.36,0.1,0.1,110.0,0.06036738605467378,0.17233075664898892],
    [0.36,0.1,0.1,180.0,0.039467762407835665,0.11549374892544934],
    [0.36,0.1,1.0,110.0,0.09624535372177485,0.2799101501149262],
    [0.36,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.001,110.0,0.03731551435349158,0.10456722494280188],
    [0.36,1.0,0.001,180.0,0.013492498381847583,0.03887475717575668],
    [0.36,1.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,0.01,180.0,0.013574192524228739,0.038738277328439936],
    [0.36,1.0,0.1,110.0,0.03895254476528791,0.10640954516459276],
    [0.36,1.0,0.1,180.0,0.028280826832015982,0.08144372491319456],
    [0.36,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.36,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.001,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.001,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.01,0.01,180.0,0.07879594971601717,0.23396920638186364],
    [0.435,0.01,0.1,110.0,0.11177413494982522,0.26781213984153257],
    [0.435,0.01,0.1,180.0,0.07982332216222571,0.23482174230900055],
    [0.435,0.01,1.0,110.0,0.10404336410294569,0.2653633751069886],
    [0.435,0.1,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.1,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,0.1,0.001,110.0,0.0514064157848007,0.11829627953551232],
    [0.435,0.1,0.001,180.0,0.024863487764014657,0.06315739336480754],
    [0.435,0.1,0.01,110.0,0.025841337025231764,0.06486046199439431],
    [0.435,0.1,0.01,180.0,0.03492304739777853,0.08697522836176026],
    [0.435,0.1,0.1,110.0,0.03738125408829439,0.09143436097440404],
    [0.435,0.1,0.1,180.0,0.03297336040443855,0.08218000709150271],
    [0.435,0.1,1.0,110.0,0.08652066536889622,0.21970131477371324],
    [0.435,0.1,1.0,180.0,0.05947106152287292,0.16182636931731134],
    [0.435,1.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.001,110.0,0.036602559829868035,0.08971152601751382],
    [0.435,1.0,0.001,180.0,0.0013058232073252154,0.0018719903367372267],
    [0.435,1.0,0.01,110.0,0.015227611796538492,0.03750268386910032],
    [0.435,1.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,0.1,110.0,0.01996487793673578,0.04930644716403539],
    [0.435,1.0,0.1,180.0,0.029502838607312336,0.07606843198185774],
    [0.435,1.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,1.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.0001,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.0001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.001,110.0,0.0152253136763512,0.03767997016309791],
    [0.435,10.0,0.001,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.01,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.01,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.1,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,0.1,180.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,1.0,110.0,3.9130083816065697e-16,1.2468865458354986e-15],
    [0.435,10.0,1.0,180.0,3.9130083816065697e-16,1.2468865458354986e-15]
]

matrix=reduce(hcat,x)'
pairplot(matrix)





###
#bond length bond angle standard deviation plots energy

# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
import LinearAlgebra
import Statistics

function plot_streching_energy(;
        filename,
        characteristics)

    # plot the theoretical and taylor function around the equilibrium
    r_theoretical=collect(0:0.05:1.5)
    r_equilibrium=1
    E_str=3/16 * (r_theoretical.^2 .-r_equilibrium).^2
    E_taylor=3/16 * 4 *(r_theoretical .-r_equilibrium).^2

    Plots.plot(r_theoretical,E_str,label="E_str",legend=:topleft)
    Plots.plot!(r_theoretical,E_taylor,label="E_taylor")
    Plots.plot!([r_equilibrium], seriestype="vline", label="Equilibrium length", color=:blue)

    # load spatial network
    path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\\"
    #filename=raw"multiple_p_quench_false_theta_array_"
    #characteristics=raw"N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=180.0_GradT=0.1_StepsPerT=0.01"
    type=raw".gml"

    spatial_network = NG.load_spatial_network_from_gml(path*filename*characteristics*type)


    # prepare for scatter plotting
    r=[]
    E=[]

    for bond in MetaGraphsNext.edge_labels(spatial_network)
        append!(r,sqrt(spatial_network[bond...]["distance_squared"]))
        append!(E,NG.local_bond_stretching_energy_keating(spatial_network, bond))
    end

    println(r)
    println(E)


    # prepare for std around equilibrium
    r_std=Statistics.std(r)

    Plots.plot!([r_equilibrium-r_std], seriestype="vline", label=false, color=:red)
    Plots.plot!([r_equilibrium+r_std], seriestype="vline", label=false, color=:red)
    
    
    Plots.scatter!(r,E,xlabel="bond length / d",ylabel="streching energy / (α d^2)",label="Measured")

    # save picture
    save_path = raw".\simulations\metric_E_str\\"
    save_filename = ("metric_E_str_2_"
        *characteristics
        *".png")

    save_total_path=save_path*save_filename

    Plots.savefig(save_total_path)

end






function plot_bending_energy(;
        filename,
        characteristics,
        theta_equilibrium=180.0)
    
    # load spatial network
    path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\\"
    type=raw".gml"

    spatial_network = NG.load_spatial_network_from_gml(path*filename*characteristics*type)
    bond_bending_const=spatial_network[]["bond_bending_const"]

    # prepare for scatter plotting
    θ=[]
    E=[]

    for vertex_label in MetaGraphsNext.labels(spatial_network)

        neighbor_label_vec::Vector{Int64} = collect(MetaGraphsNext.neighbor_labels(
        spatial_network, 
        vertex_label))

        for j in 1:spatial_network[]["coordination_nr"]
            a=sign(neighbor_label_vec[j] - vertex_label).* 
                spatial_network[vertex_label, neighbor_label_vec[j]]["vector"]

            for k in j+1:spatial_network[]["coordination_nr"]
                b=sign(neighbor_label_vec[k] - vertex_label).* 
                    spatial_network[vertex_label, neighbor_label_vec[k]]["vector"]
                #println(a)
                #println(LinearAlgebra.norm(a))

                #println(LinearAlgebra.dot(a,b))
                

                
                append!(θ,acos(LinearAlgebra.dot(a,b)/(LinearAlgebra.norm(a)*LinearAlgebra.norm(b))))
                append!(E,3/8*bond_bending_const*(LinearAlgebra.dot(a,b) + 1/3)^2)
                #println(LinearAlgebra.dot(a,b))
                #=
                if(acos(LinearAlgebra.dot(a,b)/(LinearAlgebra.norm(a)*LinearAlgebra.norm(b)))>2.5)
                    println(acos(LinearAlgebra.dot(a,b)/(LinearAlgebra.norm(a)*LinearAlgebra.norm(b))))
                    println(3/8*spatial_network[]["bond_bending_const"]*(LinearAlgebra.dot(a,b) + 1/3)^2)
                    println(a)
                    println(b)
                    println(LinearAlgebra.dot(a,b))
                end
                =#
            end
        end
    end

    #println(θ)
    #println(E)


    # PLOT
    E_min=minimum(E)
    E_max=maximum(E)
    E_range=E_max-E_min
    E_start=E_min-E_range*0.05
    E_end=E_max+E_range*0.05

    Plots.scatter(
        θ,
        E,
        xlabel="bond angle in rad",
        ylabel="bending energy / (α d^2)",
        ylimits=(E_start,E_end),
        label="Measured", 
        legend=:topleft,
        markersize=3, 
        markerstrokewidth=1, 
        color=:violet)


    # statistics of theta
    θ_mean=Statistics.mean(θ)
    θ_std=Statistics.std(θ)
    
    Plots.plot!([θ_mean], seriestype="vline", label="θ_mean", color=:green)
    Plots.plot!([θ_mean-θ_std], seriestype="vline", label="θ_mean±θ_std", color=:blue)
    Plots.plot!([θ_mean+θ_std], seriestype="vline", label=false, color=:blue)

    
    # THEORY
    # plot the theoretical and taylor function around the equilibrium
    nr_steps=100
    theta_min=minimum(θ)
    theta_max=maximum(θ)
    theta_range=theta_max-theta_min
    theta_step=theta_range/nr_steps
    theta_start=theta_min-theta_range*0.05
    theta_end=theta_max+theta_range*0.05
    theta_theoretical=collect(theta_start:theta_step:theta_end)

    r_norm=1
    E_bend=3/8 * bond_bending_const * (r_norm.^2*cos.(theta_theoretical) .- cos(theta_equilibrium)).^2
    #second_taylor_constant=3/8 * bond_bending_const * 2 * (sin(theta_equilibrium)^2-cos(theta_equilibrium)^2-1/3*cos(theta_equilibrium)) 
    #second_taylor_constant=3/8 * bond_bending_const * 2 * (sin(theta_equilibrium)^2-cos(theta_equilibrium)^2+cos(theta_equilibrium)*cos(theta_equilibrium)) 
    second_taylor_constant=3/8 * bond_bending_const * 2 * sin(theta_equilibrium)^2
    E_taylor=1/2 * second_taylor_constant *(theta_theoretical .- theta_equilibrium).^2

    Plots.plot!(theta_theoretical,E_bend,label="E_bend", color=:red)
    Plots.plot!(theta_theoretical,E_taylor,label="E_taylor", color=:orange)
    Plots.plot!([theta_equilibrium], seriestype="vline", label="θ_eq", color=:yellow)



    # SAVE
    # save picture
    save_path = raw".\simulations\metric_E_bend\\"
    save_filename = ("metric_E_bend_2_"
        *characteristics
        *".png")

    save_total_path=save_path*save_filename

    Plots.savefig(save_total_path)

end


#call functions

filename=raw"multiple_p_quench_false_theta_array_"
characteristics=raw"N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=110.0_GradT=0.1_StepsPerT=0.01"

plot_streching_energy(;
    filename=filename,
    characteristics=characteristics
)


plot_bending_energy(;
    filename=filename,
    characteristics=characteristics,
    theta_equilibrium=110.0/360.0*2*pi
    #theta_equilibrium=179.9/360.0*2*pi
)


###

include("structure_analysis_modules.jl")
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
using MetaGraphsNext


path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\\"
filename=raw"m_rad__N=216_T=0.2_Beta=0.3_GradT=0.1_StepsPerT=0.01_Theta_GS=110.0_Trial=1.gml"
spatial_network = NG.load_spatial_network_from_gml(path*filename)

#=
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216, network_type="diamond", bond_bending_const=0.285, min_ring_size=3)
spatial_network = NG.get_periodic_network(evolution_dict)
=#
function plot_SN(;
    filename::String,
    spatial_network::MetaGraphsNext.MetaGraph
    )

    P=Plots.plot()

    for vertex in MetaGraphsNext.labels(spatial_network)
        x=spatial_network[vertex]["position"][1]
        y=spatial_network[vertex]["position"][2]
        z=spatial_network[vertex]["position"][3]
        println(x)
        Plots.plot!([x],[y],[z],seriestype=:scatter,legend=false,camera = (50, 40))
    end

    plot_save_path = raw".\simulations\analysis_plot\\"
    plot_filename_start = "m_plotNW_a2_"

    plot_total_path=(plot_save_path
        *plot_filename_start
        *filename
        *".png")

    Plots.savefig(P,plot_total_path)

end

plot_SN(
    filename=filename,
    spatial_network=spatial_network
)


###

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
    shape_array)

    @assert length(nr_vertices_array)===length(shape_array)

    color_array=range(Colors.colorant"blue", stop=Colors.colorant"red", length=length(maximal_temperature_array))
    println(color_array)

    @assert length(maximal_temperature_array)===length(color_array)

    P=Plots.scatter(
        title="Bond length std vs Keating energy per vertex"
    )

    for k in eachindex(nr_vertices_array)

        nr_vertices=nr_vertices_array[k]
        shape=shape_array[k]

        for j in eachindex(maximal_temperature_array)

            maximal_temperature=maximal_temperature_array[j]
            color=color_array[j]

            for i in 1:nr_trials_per_temperature

                for m in eachindex(bond_bending_const_array)

                    bond_bending_const=bond_bending_const_array[m]

                    println("$nr_vertices"*", "*"$maximal_temperature"*", "*"$i"*", "*"$bond_bending_const" )

                    save_path = raw".\my_networks\multiple_gml\\"
                    filename = ("multiple_gml_10"
                        *"_N="*"$nr_vertices"
                        *"_T="*"$maximal_temperature"
                        *"_trial="*"$i"
                        *"_beta="*"$bond_bending_const"
                        *".gml")
                    total_path=save_path*filename

                    spatial_network=NG.load_spatial_network_from_gml(total_path)

                    bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)

                    total_energy_keating=NG.get_total_energy_keating(spatial_network)
                    energy_keating_per_vertex=total_energy_keating/spatial_network[]["nr_vertices"]

                    P=Plots.scatter!(
                        P,
                        [energy_keating_per_vertex],
                        [bond_length_std],
                        xlabel="Keating energy per vertex",
                        ylabel="Bond length std",
                        markercolor=color,
                        markershape=shape,
                        label=false,
                        cbar=true,
                        show=true)
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

    plot_save_path = raw".\my_networks\multiple_gml\\"
    plot_filename = ("multiple_gml_23"
        *"_N="*"$minimum_nr_vertices" * "-" * "$maximum_nr_vertices"
        *"_T="*"$minimum_temperature" * "-" * "$maximum_temperature"
        *"_trials="*"$nr_trials_per_temperature"
        *"_beta="*"$minimum_bond_bending_const" * "-" * "$maximum_bond_bending_const"
        *".png")
    plot_total_path=plot_save_path*plot_filename

    
    Plots.savefig(P,plot_total_path)
end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[216,512],
    maximal_temperature_array=[0.1,0.2,0.4,0.8,1.6],
    nr_trials_per_temperature=1,
    bond_bending_const_array=[0.21,0.285,0.36],
    shape_array=[:circle,:rect]
)


###

# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots

spatial_network_path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\\"
filename=raw"m_BTMC_N=216_T=0.1_Beta=1.0_GradT=0.01_StepsPerT=0.001_Theta_GS=180.0_Trial=1"
structure_dict_path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\structure\\"
analysis_data_path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\analysis_data\\"


NA.get_all_dicts_from_network_single_file(
    filename,
    spatial_network_path,
    structure_dict_path,
    analysis_data_path;
    bond_radius = 0.35,
    voxel_edge_length = 0.25,
    print_progress = true,
    print_lock = Threads.ReentrantLock())



path = raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\pore_size_dist\\"

network_path = structure_dict_path * filename *"_structure.h5"

structure_dict_network = GU.load_h5_dict(network_path)

pore_pixel_radius_array = NA.get_pore_size_distribution(structure_dict_network)

pore_pixel_radius_vec= pore_pixel_radius_array["pore_size_distribution"]

println(pore_pixel_radius_vec)

pore_pixel_radius_filtered_vec = pore_pixel_radius_vec[pore_pixel_radius_vec .> 0.0]
Plots.histogram(pore_pixel_radius_filtered_vec)
Plots.savefig(path*"pore_size_distribution_"* filename *".png")

###

import Plots
Plots.scatter([0,1,2], [6, 10, 2],  xerr = [0.1, 0.33, 0.2],yerr = [1, 2, 3], msc = 1)

###

using StatsPlots # no need for `using Plots` as that is reexported here
gr(size=(400,300))

#%%

using DataFrames
#using IndexedTables
df = DataFrame(a = 1:10, b = 10 .* rand(10), c = 10 .* rand(10))
@df df plot(:a, [:b :c], colour = [:red :blue])
@df df scatter(:a, :b, markersize = 4 .* log.(:c .+ 0.1))
t = table(1:10, rand(10), names = [:a, :b]) # IndexedTable
@df t scatter(2 .* :b)

#%%

@df df plot(:a, cols(2:3), colour = [:red :blue])

#%%

s = :b
@df df plot(:a, cols(s))

#%%

df[:red] = rand(10)
@df df plot(:a, [:b :c], colour = ^([:red :blue]))

#%%


###


"""
Copy this code and paste it into the terminal with julia
"""


# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

# import the benchmark module for using @btime
using BenchmarkTools

# import the plot module
using Plots

# prepare short tests
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216, network_type="diamond", bond_bending_const=0.285, min_ring_size=3)
spatial_network = NG.get_periodic_network(evolution_dict)

@benchmark NG.get_total_energy_keating2(spatial_network)

@benchmark NG.get_total_energy_keating(spatial_network)

#@benchmark NG.get_total_energy_keating2(spatial_network)

println("end")


###

include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and 
# analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import Random
import LaTeXStrings as Latex
import Measurements
import Polynomials
import FFTW

fontsize=18

Plots.gr()
Plots.default(grid=false, 
legend = true, 
dpi=250,
xtickfontsize=fontsize,
ytickfontsize=fontsize,
xguidefontsize=fontsize,
yguidefontsize=fontsize,
legendfontsize=fontsize,
bottom_margin = 3Plots.mm,
linewidth=3, 
thickness_scaling = 1,
framestyle = :box)

# functions to have pi ticks
function pitick(start, stop, denom; mode=:text)
    a = Int(cld(start, 2*π/denom))
    b = Int(fld(stop, 2*π/denom))
    tick = range(a*2*π/denom, b*2*π/denom; step=2*π/denom)
    ticklabel = piticklabel.( 2 .* (a:b) .// denom, Val(mode))
    tick, ticklabel
end

function piticklabel(x::Rational, ::Val{:text})
    iszero(x) && return "0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return S * N * "π"
    S * N * "π/" * repr(d)
end

function piticklabel(x::Rational, ::Val{:latex})
    iszero(x) && return Latex.L"0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return Latex.L"%$S%$N\pi"
    Latex.L"%$S\frac{%$N\pi}{%$d}"
end


"""
    get_tickslogscale(lims; skiplog=false)
Return a tuple (ticks, ticklabels) for the axis limit `lims`
where multiples of 10 are major ticks with label and minor ticks have no label
skiplog argument should be set to true if `lims` is already in log scale.
"""
function get_tickslogscale(lims::Tuple{T, T}; skiplog::Bool=false) where {T<:AbstractFloat}
    mags = if skiplog
        # if the limits are already in log scale
        floor.(lims)
    else
        floor.(log10.(lims))
    end
    rlims = if skiplog; 10 .^(lims) else lims end

    total_tickvalues = []
    total_ticknames = []

    rgs = range(mags..., step=1)
    for (i, m) in enumerate(rgs)
        if m >= 0
            tickvalues = range(Int(10^m), Int(10^(m+1)); step=Int(10^m))
            ticknames  = vcat([string(round(Int, 10^(m)))],
                              ["" for i in 2:9],
                              [string(round(Int, 10^(m+1)))])
        else
            tickvalues = range(10^m, 10^(m+1); step=10^m)
            ticknames  = vcat([string(10^(m))], ["" for i in 2:9], [string(10^(m+1))])
        end

        if i==1
            # lower bound
            indexlb = findlast(x->x<rlims[1], tickvalues)
            if isnothing(indexlb); indexlb=1 end
        else
            indexlb = 1
        end
        if i==length(rgs)
            # higher bound
            indexhb = findfirst(x->x>rlims[2], tickvalues)
            if isnothing(indexhb); indexhb=10 end
        else
            # do not take the last index if not the last magnitude
            indexhb = 9
        end

        total_tickvalues = vcat(total_tickvalues, tickvalues[indexlb:indexhb])
        total_ticknames = vcat(total_ticknames, ticknames[indexlb:indexhb])
    end
    return (total_tickvalues, total_ticknames)
end

"""
    fancylogscale!(p; forcex=false, forcey=false)
Transform the ticks to log scale for the axis with scale=:log10.
forcex and forcey can be set to true to force the transformation
if the variable is already expressed in log10 units.
"""
function fancylogscale!(p::Plots.Subplot; forcex::Bool=false, forcey::Bool=false)
    kwargs = Dict()
    for (ax, force, lims) in zip((:x, :y), (forcex, forcey), (Plots.xlims, Plots.ylims))
        axis = Symbol("$(ax)axis")
        ticks = Symbol("$(ax)ticks")

        if force || p.attr[axis][:scale] == :log10
            # Get limits of the plot and convert to Float
            ls = float.(lims(p))
            ts = if force
                (vals, labs) = get_tickslogscale(ls; skiplog=true)
                (log10.(vals), labs)
            else
                get_tickslogscale(ls)
            end
            kwargs[ticks] = ts
        end
    end

    if length(kwargs) > 0
        Plots.plot!(p; kwargs...)
    end
    p
end
fancylogscale!(p::Plots.Plot; kwargs...) = (fancylogscale!(p.subplots[1]; kwargs...); return p)
fancylogscale!(; kwargs...) = fancylogscale!(Plots.plot!(); kwargs...)

path = raw"..\..\presentations\material\\"

# Set the size of the checkerboard
n_rows, n_cols = 8, 8  # 8x8 checkerboard

# Create an empty grid to hold the colors
checkerboard = [Random.rand() for i in 1:n_rows, j in 1:n_cols]

colormap = Plots.cgrad(:roma)

checkerboard_colors = [colormap[checkerboard[i,j]] for i in 1:n_rows, j in 1:n_cols]

# Plot the checkerboard
Plots.plot(Plots.heatmap(1:n_rows, 1:n_cols, checkerboard_colors, aspect_ratio=:equal))
Plots.savefig(path*"checkerboard.png")

colormap = Plots.cgrad(:roma, scale=:lin)
checkerboard_fft = FFTW.fft(checkerboard)
checkerboard_fft_normalized = abs.(checkerboard_fft) ./ 4
checkerboard_fft_colors = [colormap[abs(checkerboard_fft_normalized[i,j])] for i in 1:n_rows, j in 1:n_cols]
Plots.plot(Plots.heatmap(1:n_rows, 1:n_cols, checkerboard_fft_colors, aspect_ratio=:equal))
Plots.savefig(path*"checkerboard_fft.png")


###



"""
This module contains all functions with tests of functions or modules
"""
module MyTests

# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

# import the text module for using @test
using Test

#import the module MetaGraphsNext
using MetaGraphsNext



# prepare short tests
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216, network_type="diamond", bond_bending_const=0.285, min_ring_size=3)
spatial_network = NG.get_periodic_network(evolution_dict)

# testing with @test
@test true
@test spatial_network[]["supercell_edge_length"]===6.9282032302755105






try

    #Energy test
    for vertex in MetaGraphsNext.labels(spatial_network)
        println(vertex)
        println(NG.local_bond_bending_energy_keating(spatial_network, vertex))
        println(NG.local_bond_bending_energy_keating2(spatial_network, vertex))
        @test NG.local_bond_bending_energy_keating(spatial_network, vertex)===
              NG.local_bond_bending_energy_keating2(spatial_network, vertex)
        println("yes test has passed") 
    end
catch e
           error_msg = sprint(showerror, e)
           st = sprint((io,v) -> show(io, "text/plain", v), stacktrace(catch_backtrace()))
           @warn "Trouble doing things:\n$(error_msg)\n$(st)"
           println("Trouble doing things:\n$(error_msg)\n$(st)")
end





# testing multiple tests with table (passed, failed, total, time)
@testset "All tests" begin

    @testset "First" begin
        @test false | true
        @test false | true
    end

    @testset "Second" begin
        # is run, ok
        @test true 
    end
end




end




###

import plotly.express as px
from sklearn.decomposition import PCA
from sklearn import datasets
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import plotly.graph_objects as go

list0_25=[
[0.05, 0.05, 0.1, 0.1, 0.264],
[0.1, 0.05, 0.1, 0.1, 0.215],
[0.05, 0.05, 0.125, 0.1, 0.23],
[0.1, 0.05, 0.125, 0.1, 0.251],
[0.05, 0.05, 0.15, 0.1, 0.235],
[0.1, 0.05, 0.15, 0.1, 0.276],
[0.05, 0.05, 0.175, 0.1, 0.2],
[0.1, 0.05, 0.175, 0.1, 0.232],
[0.05, 0.05, 0.2, 0.1, 0.305],
[0.05, 0.05, 0.15, 0.1, 0.246],
[0.05, 0.05, 0.175, 0.1, 0.252],
[0.05, 0.05, 0.2, 0.1, 0.322],
[0.1, 0.1, 0.1, 0.1, 0.307],
[0.15, 0.1, 0.1, 0.1, 0.199],
[0.2, 0.1, 0.1, 0.1, 0.227],
[0.1, 0.1, 0.125, 0.1, 0.261],
[0.15, 0.1, 0.125, 0.1, 0.102],
[0.2, 0.1, 0.125, 0.1, 0.194],
[0.1, 0.1, 0.15, 0.1, 0.233],
[0.15, 0.1, 0.15, 0.1, 0.141],
[0.2, 0.1, 0.15, 0.1, 0.312],
[0.1, 0.1, 0.175, 0.1, 0.243],
[0.15, 0.1, 0.175, 0.1, 0.162],
[0.2, 0.1, 0.175, 0.1, 0.249],
[0.1, 0.1, 0.2, 0.1, 0.333],
[0.15, 0.1, 0.2, 0.1, 0.259],
[0.05, 0.1, 0.1, 0.1, 0.267],
[0.05, 0.1, 0.125, 0.1, 0.325],
[0.1, 0.1, 0.125, 0.1, 0.281],
[0.05, 0.1, 0.15, 0.1, 0.296],
[0.1, 0.1, 0.175, 0.1, 0.285],
[0.15, 0.15, 0.1, 0.1, 0.282],
[0.2, 0.15, 0.1, 0.1, 0.279],
[0.25, 0.15, 0.1, 0.1, 0.302],
[0.3, 0.15, 0.1, 0.1, 0.344],
[0.15, 0.15, 0.125, 0.1, 0.169],
[0.2, 0.15, 0.125, 0.1, 0.238],
[0.25, 0.15, 0.125, 0.1, 0.227],
[0.1, 0.15, 0.15, 0.1, 0.344],
[0.15, 0.15, 0.15, 0.1, 0.082],
[0.2, 0.15, 0.15, 0.1, 0.123],
[0.25, 0.15, 0.15, 0.1, 0.234],
[0.15, 0.15, 0.175, 0.1, 0.141],
[0.2, 0.15, 0.175, 0.1, 0.074],
[0.15, 0.15, 0.2, 0.1, 0.152],
[0.2, 0.15, 0.2, 0.1, 0.216],
[0.25, 0.15, 0.2, 0.1, 0.223],
[0.1, 0.15, 0.1, 0.1, 0.236],
[0.1, 0.15, 0.125, 0.1, 0.268],
[0.1, 0.15, 0.15, 0.1, 0.267],
[0.1, 0.15, 0.175, 0.1, 0.237],
[0.25, 0.2, 0.1, 0.1, 0.257],
[0.3, 0.2, 0.1, 0.1, 0.274],
[0.15, 0.2, 0.125, 0.1, 0.318],
[0.2, 0.2, 0.125, 0.1, 0.32],
[0.25, 0.2, 0.125, 0.1, 0.065],
[0.3, 0.2, 0.125, 0.1, 0.066],
[0.35, 0.2, 0.125, 0.1, 0.189],
[0.4, 0.2, 0.125, 0.1, 0.212],
[0.15, 0.2, 0.15, 0.1, 0.343],
[0.2, 0.2, 0.15, 0.1, 0.151],
[0.25, 0.2, 0.15, 0.1, 0.043],
[0.3, 0.2, 0.15, 0.1, 0.16],
[0.35, 0.2, 0.15, 0.1, 0.26],
[0.4, 0.2, 0.15, 0.1, 0.25],
[0.2, 0.2, 0.175, 0.1, 0.188],
[0.25, 0.2, 0.175, 0.1, 0.193],
[0.3, 0.2, 0.175, 0.1, 0.253],
[0.35, 0.2, 0.175, 0.1, 0.298],
[0.2, 0.2, 0.2, 0.1, 0.164],
[0.25, 0.2, 0.2, 0.1, 0.098],
[0.1, 0.2, 0.1, 0.1, 0.291],
[0.1, 0.2, 0.15, 0.1, 0.337],
[0.1, 0.2, 0.175, 0.1, 0.325],
[0.25, 0.25, 0.125, 0.1, 0.273],
[0.3, 0.25, 0.125, 0.1, 0.17],
[0.35, 0.25, 0.125, 0.1, 0.325],
[0.4, 0.25, 0.125, 0.1, 0.327],
[0.2, 0.25, 0.15, 0.1, 0.277],
[0.25, 0.25, 0.15, 0.1, 0.236],
[0.3, 0.25, 0.15, 0.1, 0.102],
[0.35, 0.25, 0.15, 0.1, 0.096],
[0.4, 0.25, 0.15, 0.1, 0.23],
[0.2, 0.25, 0.175, 0.1, 0.312],
[0.25, 0.25, 0.175, 0.1, 0.085],
[0.3, 0.25, 0.175, 0.1, 0.06],
[0.35, 0.25, 0.175, 0.1, 0.168],
[0.4, 0.25, 0.175, 0.1, 0.173],
[0.5, 0.25, 0.175, 0.1, 0.331],
[0.2, 0.25, 0.2, 0.1, 0.208],
[0.25, 0.25, 0.2, 0.1, 0.302],
[0.3, 0.25, 0.2, 0.1, 0.076],
[0.35, 0.25, 0.2, 0.1, 0.13],
[0.4, 0.25, 0.2, 0.1, 0.209],
[0.5, 0.25, 0.2, 0.1, 0.317],
[0.2, 0.25, 0.1, 0.1, 0.314],
[0.3, 0.25, 0.1, 0.1, 0.187],
[0.35, 0.25, 0.1, 0.1, 0.27],
[0.25, 0.3, 0.1, 0.1, 0.21],
[0.3, 0.3, 0.1, 0.1, 0.136],
[0.35, 0.3, 0.125, 0.1, 0.235],
[0.4, 0.3, 0.125, 0.1, 0.241],
[0.5, 0.3, 0.15, 0.1, 0.275],
[0.3, 0.3, 0.1, 0.1, 0.125],
[0.35, 0.3, 0.1, 0.1, 0.209],
[0.25, 0.35, 0.1, 0.1, 0.269],
[0.3, 0.35, 0.1, 0.1, 0.2],
[0.35, 0.35, 0.125, 0.1, 0.198],
[0.4, 0.35, 0.125, 0.1, 0.183],
[0.5, 0.35, 0.15, 0.1, 0.176],
[0.4, 0.5, 0.1, 0.1, 0.325],
[0.5, 0.5, 0.125, 0.1, 0.327],
[0.05, 0.05, 0.1, 0.125, 0.296],
[0.1, 0.05, 0.1, 0.125, 0.245],
[0.05, 0.05, 0.125, 0.125, 0.248],
[0.1, 0.05, 0.125, 0.125, 0.265],
[0.05, 0.05, 0.15, 0.125, 0.275],
[0.1, 0.05, 0.15, 0.125, 0.239],
[0.05, 0.05, 0.175, 0.125, 0.22],
[0.1, 0.05, 0.175, 0.125, 0.21],
[0.05, 0.05, 0.2, 0.125, 0.324],
[0.1, 0.05, 0.2, 0.125, 0.304],
[0.05, 0.05, 0.1, 0.125, 0.276],
[0.05, 0.05, 0.125, 0.125, 0.323],
[0.05, 0.05, 0.15, 0.125, 0.181],
[0.05, 0.05, 0.175, 0.125, 0.163],
[0.05, 0.05, 0.2, 0.125, 0.218],
[0.1, 0.1, 0.1, 0.125, 0.222],
[0.15, 0.1, 0.1, 0.125, 0.212],
[0.2, 0.1, 0.1, 0.125, 0.273],
[0.1, 0.1, 0.125, 0.125, 0.174],
[0.15, 0.1, 0.125, 0.125, 0.198],
[0.2, 0.1, 0.125, 0.125, 0.262],
[0.1, 0.1, 0.15, 0.125, 0.124],
[0.15, 0.1, 0.15, 0.125, 0.278],
[0.1, 0.1, 0.175, 0.125, 0.128],
[0.15, 0.1, 0.175, 0.125, 0.27],
[0.1, 0.1, 0.2, 0.125, 0.258],
[0.05, 0.1, 0.1, 0.125, 0.298],
[0.1, 0.1, 0.125, 0.125, 0.205],
[0.05, 0.1, 0.15, 0.125, 0.336],
[0.1, 0.1, 0.15, 0.125, 0.316],
[0.1, 0.1, 0.175, 0.125, 0.206],
[0.15, 0.15, 0.1, 0.125, 0.327],
[0.2, 0.15, 0.1, 0.125, 0.32],
[0.25, 0.15, 0.1, 0.125, 0.312],
[0.15, 0.15, 0.125, 0.125, 0.22],
[0.2, 0.15, 0.125, 0.125, 0.281],
[0.25, 0.15, 0.125, 0.125, 0.19],
[0.3, 0.15, 0.125, 0.125, 0.277],
[0.15, 0.15, 0.15, 0.125, 0.117],
[0.2, 0.15, 0.15, 0.125, 0.077],
[0.25, 0.15, 0.15, 0.125, 0.181],
[0.15, 0.15, 0.175, 0.125, 0.164],
[0.2, 0.15, 0.175, 0.125, 0.061],
[0.15, 0.15, 0.2, 0.125, 0.138],
[0.2, 0.15, 0.2, 0.125, 0.141],
[0.25, 0.15, 0.2, 0.125, 0.198],
[0.1, 0.15, 0.1, 0.125, 0.113],
[0.15, 0.15, 0.1, 0.125, 0.195],
[0.1, 0.15, 0.125, 0.125, 0.205],
[0.15, 0.15, 0.125, 0.125, 0.269],
[0.1, 0.15, 0.15, 0.125, 0.159],
[0.15, 0.15, 0.15, 0.125, 0.306],
[0.1, 0.15, 0.175, 0.125, 0.18],
[0.15, 0.15, 0.175, 0.125, 0.237],
[0.1, 0.15, 0.2, 0.125, 0.297],
[0.15, 0.15, 0.2, 0.125, 0.295],
[0.25, 0.2, 0.1, 0.125, 0.323],
[0.3, 0.2, 0.1, 0.125, 0.344],
[0.15, 0.2, 0.125, 0.125, 0.317],
[0.25, 0.2, 0.125, 0.125, 0.133],
[0.3, 0.2, 0.125, 0.125, 0.104],
[0.35, 0.2, 0.125, 0.125, 0.27],
[0.4, 0.2, 0.125, 0.125, 0.289],
[0.15, 0.2, 0.15, 0.125, 0.265],
[0.2, 0.2, 0.15, 0.125, 0.112],
[0.25, 0.2, 0.15, 0.125, 0.087],
[0.3, 0.2, 0.15, 0.125, 0.162],
[0.35, 0.2, 0.15, 0.125, 0.264],
[0.4, 0.2, 0.15, 0.125, 0.3],
[0.15, 0.2, 0.175, 0.125, 0.309],
[0.2, 0.2, 0.175, 0.125, 0.155],
[0.25, 0.2, 0.175, 0.125, 0.192],
[0.3, 0.2, 0.175, 0.125, 0.242],
[0.35, 0.2, 0.175, 0.125, 0.322],
[0.15, 0.2, 0.2, 0.125, 0.268],
[0.2, 0.2, 0.2, 0.125, 0.052],
[0.25, 0.2, 0.2, 0.125, 0.167],
[0.3, 0.2, 0.2, 0.125, 0.307],
[0.1, 0.2, 0.1, 0.125, 0.24],
[0.15, 0.2, 0.1, 0.125, 0.287],
[0.1, 0.2, 0.125, 0.125, 0.318],
[0.1, 0.2, 0.15, 0.125, 0.288],
[0.1, 0.2, 0.175, 0.125, 0.287],
[0.15, 0.2, 0.175, 0.125, 0.316],
[0.25, 0.25, 0.125, 0.125, 0.24],
[0.3, 0.25, 0.125, 0.125, 0.131],
[0.35, 0.25, 0.125, 0.125, 0.283],
[0.4, 0.25, 0.125, 0.125, 0.285],
[0.2, 0.25, 0.15, 0.125, 0.264],
[0.25, 0.25, 0.15, 0.125, 0.206],
[0.3, 0.25, 0.15, 0.125, 0.054],
[0.35, 0.25, 0.15, 0.125, 0.06],
[0.4, 0.25, 0.15, 0.125, 0.188],
[0.2, 0.25, 0.175, 0.125, 0.301],
[0.25, 0.25, 0.175, 0.125, 0.032],
[0.3, 0.25, 0.175, 0.125, 0.03],
[0.35, 0.25, 0.175, 0.125, 0.135],
[0.4, 0.25, 0.175, 0.125, 0.151],
[0.5, 0.25, 0.175, 0.125, 0.325],
[0.2, 0.25, 0.2, 0.125, 0.209],
[0.25, 0.25, 0.2, 0.125, 0.268],
[0.3, 0.25, 0.2, 0.125, 0.132],
[0.35, 0.25, 0.2, 0.125, 0.139],
[0.4, 0.25, 0.2, 0.125, 0.194],
[0.5, 0.25, 0.2, 0.125, 0.308],
[0.1, 0.25, 0.1, 0.125, 0.34],
[0.15, 0.25, 0.1, 0.125, 0.236],
[0.2, 0.25, 0.1, 0.125, 0.069],
[0.3, 0.25, 0.1, 0.125, 0.305],
[0.35, 0.25, 0.1, 0.125, 0.314],
[0.15, 0.25, 0.125, 0.125, 0.217],
[0.2, 0.25, 0.125, 0.125, 0.148],
[0.25, 0.25, 0.125, 0.125, 0.305],
[0.15, 0.25, 0.15, 0.125, 0.271],
[0.2, 0.25, 0.15, 0.125, 0.228],
[0.15, 0.25, 0.175, 0.125, 0.207],
[0.2, 0.25, 0.175, 0.125, 0.234],
[0.15, 0.25, 0.2, 0.125, 0.224],
[0.2, 0.25, 0.2, 0.125, 0.226],
[0.25, 0.25, 0.2, 0.125, 0.297],
[0.3, 0.3, 0.125, 0.125, 0.29],
[0.3, 0.3, 0.15, 0.125, 0.226],
[0.35, 0.3, 0.15, 0.125, 0.164],
[0.4, 0.3, 0.15, 0.125, 0.274],
[0.25, 0.3, 0.175, 0.125, 0.203],
[0.3, 0.3, 0.175, 0.125, 0.162],
[0.35, 0.3, 0.175, 0.125, 0.185],
[0.4, 0.3, 0.175, 0.125, 0.152],
[0.5, 0.3, 0.175, 0.125, 0.203],
[0.2, 0.3, 0.2, 0.125, 0.344],
[0.3, 0.3, 0.2, 0.125, 0.143],
[0.35, 0.3, 0.2, 0.125, 0.076],
[0.4, 0.3, 0.2, 0.125, 0.147],
[0.5, 0.3, 0.2, 0.125, 0.203],
[0.15, 0.3, 0.1, 0.125, 0.288],
[0.2, 0.3, 0.1, 0.125, 0.063],
[0.25, 0.3, 0.1, 0.125, 0.194],
[0.35, 0.3, 0.1, 0.125, 0.337],
[0.15, 0.3, 0.125, 0.125, 0.254],
[0.2, 0.3, 0.125, 0.125, 0.109],
[0.25, 0.3, 0.125, 0.125, 0.147],
[0.15, 0.3, 0.15, 0.125, 0.29],
[0.2, 0.3, 0.15, 0.125, 0.145],
[0.25, 0.3, 0.15, 0.125, 0.195],
[0.15, 0.3, 0.175, 0.125, 0.255],
[0.2, 0.3, 0.175, 0.125, 0.131],
[0.25, 0.3, 0.175, 0.125, 0.287],
[0.15, 0.3, 0.2, 0.125, 0.252],
[0.2, 0.3, 0.2, 0.125, 0.168],
[0.25, 0.3, 0.2, 0.125, 0.139],
[0.3, 0.35, 0.125, 0.125, 0.329],
[0.3, 0.35, 0.15, 0.125, 0.27],
[0.35, 0.35, 0.15, 0.125, 0.206],
[0.4, 0.35, 0.15, 0.125, 0.289],
[0.25, 0.35, 0.175, 0.125, 0.249],
[0.3, 0.35, 0.175, 0.125, 0.214],
[0.35, 0.35, 0.175, 0.125, 0.207],
[0.4, 0.35, 0.175, 0.125, 0.171],
[0.5, 0.35, 0.175, 0.125, 0.167],
[0.3, 0.35, 0.2, 0.125, 0.216],
[0.35, 0.35, 0.2, 0.125, 0.13],
[0.4, 0.35, 0.2, 0.125, 0.152],
[0.5, 0.35, 0.2, 0.125, 0.173],
[0.2, 0.35, 0.1, 0.125, 0.232],
[0.3, 0.35, 0.1, 0.125, 0.201],
[0.35, 0.35, 0.1, 0.125, 0.167],
[0.5, 0.35, 0.125, 0.125, 0.267],
[0.25, 0.35, 0.2, 0.125, 0.321],
[0.3, 0.4, 0.15, 0.125, 0.342],
[0.35, 0.4, 0.15, 0.125, 0.279],
[0.4, 0.4, 0.15, 0.125, 0.294],
[0.25, 0.4, 0.175, 0.125, 0.325],
[0.3, 0.4, 0.175, 0.125, 0.303],
[0.35, 0.4, 0.175, 0.125, 0.241],
[0.4, 0.4, 0.175, 0.125, 0.212],
[0.5, 0.4, 0.175, 0.125, 0.114],
[0.35, 0.4, 0.2, 0.125, 0.237],
[0.4, 0.4, 0.2, 0.125, 0.179],
[0.5, 0.4, 0.2, 0.125, 0.129],
[0.2, 0.4, 0.1, 0.125, 0.262],
[0.3, 0.4, 0.1, 0.125, 0.206],
[0.35, 0.4, 0.1, 0.125, 0.159],
[0.5, 0.4, 0.125, 0.125, 0.247],
[0.25, 0.4, 0.2, 0.125, 0.328],
[0.4, 0.5, 0.125, 0.125, 0.326],
[0.5, 0.5, 0.15, 0.125, 0.024],
[0.05, 0.05, 0.1, 0.15, 0.198],
[0.1, 0.05, 0.1, 0.15, 0.258],
[0.05, 0.05, 0.125, 0.15, 0.168],
[0.1, 0.05, 0.125, 0.15, 0.316],
[0.05, 0.05, 0.15, 0.15, 0.167],
[0.05, 0.05, 0.175, 0.15, 0.136],
[0.1, 0.05, 0.175, 0.15, 0.333],
[0.05, 0.05, 0.2, 0.15, 0.244],
[0.05, 0.05, 0.1, 0.15, 0.308],
[0.05, 0.05, 0.15, 0.15, 0.201],
[0.05, 0.05, 0.175, 0.15, 0.219],
[0.05, 0.05, 0.2, 0.15, 0.291],
[0.1, 0.1, 0.1, 0.15, 0.289],
[0.15, 0.1, 0.1, 0.15, 0.215],
[0.2, 0.1, 0.1, 0.15, 0.252],
[0.1, 0.1, 0.125, 0.15, 0.243],
[0.15, 0.1, 0.125, 0.15, 0.139],
[0.2, 0.1, 0.125, 0.15, 0.224],
[0.1, 0.1, 0.15, 0.15, 0.186],
[0.15, 0.1, 0.15, 0.15, 0.165],
[0.1, 0.1, 0.175, 0.15, 0.202],
[0.15, 0.1, 0.175, 0.15, 0.166],
[0.2, 0.1, 0.175, 0.15, 0.296],
[0.1, 0.1, 0.2, 0.15, 0.279],
[0.15, 0.1, 0.2, 0.15, 0.277],
[0.05, 0.1, 0.1, 0.15, 0.225],
[0.1, 0.1, 0.1, 0.15, 0.247],
[0.05, 0.1, 0.125, 0.15, 0.218],
[0.1, 0.1, 0.125, 0.15, 0.117],
[0.05, 0.1, 0.15, 0.15, 0.287],
[0.1, 0.1, 0.15, 0.15, 0.199],
[0.05, 0.1, 0.175, 0.15, 0.314],
[0.1, 0.1, 0.175, 0.15, 0.139],
[0.05, 0.1, 0.2, 0.15, 0.315],
[0.1, 0.1, 0.2, 0.15, 0.262],
[0.15, 0.15, 0.125, 0.15, 0.249],
[0.2, 0.15, 0.125, 0.15, 0.313],
[0.25, 0.15, 0.125, 0.15, 0.221],
[0.3, 0.15, 0.125, 0.15, 0.29],
[0.15, 0.15, 0.15, 0.15, 0.11],
[0.2, 0.15, 0.15, 0.15, 0.114],
[0.25, 0.15, 0.15, 0.15, 0.207],
[0.15, 0.15, 0.175, 0.15, 0.138],
[0.2, 0.15, 0.175, 0.15, 0.107],
[0.15, 0.15, 0.2, 0.15, 0.089],
[0.2, 0.15, 0.2, 0.15, 0.144],
[0.25, 0.15, 0.2, 0.15, 0.233],
[0.1, 0.15, 0.1, 0.15, 0.126],
[0.15, 0.15, 0.1, 0.15, 0.078],
[0.1, 0.15, 0.125, 0.15, 0.218],
[0.15, 0.15, 0.125, 0.15, 0.125],
[0.2, 0.15, 0.125, 0.15, 0.321],
[0.1, 0.15, 0.15, 0.15, 0.142],
[0.15, 0.15, 0.15, 0.15, 0.123],
[0.1, 0.15, 0.175, 0.15, 0.215],
[0.15, 0.15, 0.175, 0.15, 0.115],
[0.1, 0.15, 0.2, 0.15, 0.125],
[0.15, 0.15, 0.2, 0.15, 0.14],
[0.2, 0.15, 0.2, 0.15, 0.339],
[0.25, 0.2, 0.125, 0.15, 0.196],
[0.3, 0.2, 0.125, 0.15, 0.127],
[0.35, 0.2, 0.125, 0.15, 0.304],
[0.4, 0.2, 0.125, 0.15, 0.316],
[0.15, 0.2, 0.15, 0.15, 0.286],
[0.2, 0.2, 0.15, 0.15, 0.169],
[0.25, 0.2, 0.15, 0.15, 0.152],
[0.3, 0.2, 0.15, 0.15, 0.135],
[0.35, 0.2, 0.15, 0.15, 0.217],
[0.4, 0.2, 0.15, 0.15, 0.288],
[0.15, 0.2, 0.175, 0.15, 0.317],
[0.2, 0.2, 0.175, 0.15, 0.208],
[0.25, 0.2, 0.175, 0.15, 0.155],
[0.3, 0.2, 0.175, 0.15, 0.189],
[0.35, 0.2, 0.175, 0.15, 0.283],
[0.4, 0.2, 0.175, 0.15, 0.314],
[0.15, 0.2, 0.2, 0.15, 0.265],
[0.2, 0.2, 0.2, 0.15, 0.081],
[0.25, 0.2, 0.2, 0.15, 0.228],
[0.3, 0.2, 0.2, 0.15, 0.227],
[0.35, 0.2, 0.2, 0.15, 0.305],
[0.1, 0.2, 0.1, 0.15, 0.246],
[0.15, 0.2, 0.1, 0.15, 0.115],
[0.2, 0.2, 0.1, 0.15, 0.153],
[0.15, 0.2, 0.125, 0.15, 0.08],
[0.2, 0.2, 0.125, 0.15, 0.126],
[0.1, 0.2, 0.15, 0.15, 0.282],
[0.15, 0.2, 0.15, 0.15, 0.127],
[0.2, 0.2, 0.15, 0.15, 0.2],
[0.1, 0.2, 0.175, 0.15, 0.334],
[0.15, 0.2, 0.175, 0.15, 0.078],
[0.2, 0.2, 0.175, 0.15, 0.225],
[0.1, 0.2, 0.2, 0.15, 0.32],
[0.15, 0.2, 0.2, 0.15, 0.084],
[0.2, 0.2, 0.2, 0.15, 0.171],
[0.25, 0.25, 0.125, 0.15, 0.273],
[0.3, 0.25, 0.125, 0.15, 0.17],
[0.35, 0.25, 0.125, 0.15, 0.322],
[0.4, 0.25, 0.125, 0.15, 0.323],
[0.2, 0.25, 0.15, 0.15, 0.281],
[0.25, 0.25, 0.15, 0.15, 0.237],
[0.3, 0.25, 0.15, 0.15, 0.099],
[0.35, 0.25, 0.15, 0.15, 0.085],
[0.4, 0.25, 0.15, 0.15, 0.223],
[0.2, 0.25, 0.175, 0.15, 0.316],
[0.25, 0.25, 0.175, 0.15, 0.08],
[0.3, 0.25, 0.175, 0.15, 0.05],
[0.35, 0.25, 0.175, 0.15, 0.158],
[0.4, 0.25, 0.175, 0.15, 0.163],
[0.5, 0.25, 0.175, 0.15, 0.32],
[0.2, 0.25, 0.2, 0.15, 0.215],
[0.25, 0.25, 0.2, 0.15, 0.302],
[0.3, 0.25, 0.2, 0.15, 0.078],
[0.35, 0.25, 0.2, 0.15, 0.12],
[0.4, 0.25, 0.2, 0.15, 0.198],
[0.5, 0.25, 0.2, 0.15, 0.306],
[0.1, 0.25, 0.1, 0.15, 0.332],
[0.15, 0.25, 0.1, 0.15, 0.219],
[0.2, 0.25, 0.1, 0.15, 0.044],
[0.25, 0.25, 0.1, 0.15, 0.296],
[0.3, 0.25, 0.1, 0.15, 0.337],
[0.35, 0.25, 0.1, 0.15, 0.337],
[0.15, 0.25, 0.125, 0.15, 0.189],
[0.2, 0.25, 0.125, 0.15, 0.094],
[0.25, 0.25, 0.125, 0.15, 0.257],
[0.15, 0.25, 0.15, 0.15, 0.235],
[0.2, 0.25, 0.15, 0.15, 0.17],
[0.25, 0.25, 0.15, 0.15, 0.3],
[0.15, 0.25, 0.175, 0.15, 0.185],
[0.2, 0.25, 0.175, 0.15, 0.176],
[0.15, 0.25, 0.2, 0.15, 0.191],
[0.2, 0.25, 0.2, 0.15, 0.169],
[0.25, 0.25, 0.2, 0.15, 0.254],
[0.3, 0.3, 0.125, 0.15, 0.292],
[0.3, 0.3, 0.15, 0.15, 0.227],
[0.35, 0.3, 0.15, 0.15, 0.173],
[0.4, 0.3, 0.15, 0.15, 0.292],
[0.25, 0.3, 0.175, 0.15, 0.206],
[0.3, 0.3, 0.175, 0.15, 0.163],
[0.35, 0.3, 0.175, 0.15, 0.206],
[0.4, 0.3, 0.175, 0.15, 0.177],
[0.5, 0.3, 0.175, 0.15, 0.239],
[0.2, 0.3, 0.2, 0.15, 0.329],
[0.3, 0.3, 0.2, 0.15, 0.111],
[0.35, 0.3, 0.2, 0.15, 0.09],
[0.4, 0.3, 0.2, 0.15, 0.177],
[0.5, 0.3, 0.2, 0.15, 0.238],
[0.15, 0.3, 0.1, 0.15, 0.329],
[0.2, 0.3, 0.1, 0.15, 0.14],
[0.25, 0.3, 0.1, 0.15, 0.11],
[0.15, 0.3, 0.125, 0.15, 0.292],
[0.2, 0.3, 0.125, 0.15, 0.147],
[0.25, 0.3, 0.125, 0.15, 0.056],
[0.3, 0.3, 0.125, 0.15, 0.287],
[0.15, 0.3, 0.15, 0.15, 0.317],
[0.2, 0.3, 0.15, 0.15, 0.142],
[0.25, 0.3, 0.15, 0.15, 0.106],
[0.3, 0.3, 0.15, 0.15, 0.282],
[0.35, 0.3, 0.15, 0.15, 0.337],
[0.15, 0.3, 0.175, 0.15, 0.297],
[0.2, 0.3, 0.175, 0.15, 0.117],
[0.25, 0.3, 0.175, 0.15, 0.202],
[0.3, 0.3, 0.175, 0.15, 0.272],
[0.15, 0.3, 0.2, 0.15, 0.286],
[0.2, 0.3, 0.2, 0.15, 0.176],
[0.25, 0.3, 0.2, 0.15, 0.049],
[0.3, 0.3, 0.2, 0.15, 0.267],
[0.3, 0.35, 0.125, 0.15, 0.315],
[0.3, 0.35, 0.15, 0.15, 0.253],
[0.35, 0.35, 0.15, 0.15, 0.19],
[0.4, 0.35, 0.15, 0.15, 0.286],
[0.25, 0.35, 0.175, 0.15, 0.231],
[0.3, 0.35, 0.175, 0.15, 0.193],
[0.35, 0.35, 0.175, 0.15, 0.201],
[0.4, 0.35, 0.175, 0.15, 0.166],
[0.5, 0.35, 0.175, 0.15, 0.186],
[0.3, 0.35, 0.2, 0.15, 0.182],
[0.35, 0.35, 0.2, 0.15, 0.107],
[0.4, 0.35, 0.2, 0.15, 0.153],
[0.5, 0.35, 0.2, 0.15, 0.189],
[0.2, 0.35, 0.1, 0.15, 0.207],
[0.25, 0.35, 0.1, 0.15, 0.105],
[0.2, 0.35, 0.125, 0.15, 0.219],
[0.25, 0.35, 0.125, 0.15, 0.06],
[0.3, 0.35, 0.125, 0.15, 0.19],
[0.35, 0.35, 0.125, 0.15, 0.259],
[0.2, 0.35, 0.15, 0.15, 0.206],
[0.25, 0.35, 0.15, 0.15, 0.093],
[0.3, 0.35, 0.15, 0.15, 0.198],
[0.35, 0.35, 0.15, 0.15, 0.231],
[0.2, 0.35, 0.175, 0.15, 0.18],
[0.25, 0.35, 0.175, 0.15, 0.183],
[0.3, 0.35, 0.175, 0.15, 0.169],
[0.35, 0.35, 0.175, 0.15, 0.32],
[0.2, 0.35, 0.2, 0.15, 0.24],
[0.25, 0.35, 0.2, 0.15, 0.039],
[0.3, 0.35, 0.2, 0.15, 0.189],
[0.35, 0.35, 0.2, 0.15, 0.323],
[0.3, 0.4, 0.15, 0.15, 0.338],
[0.35, 0.4, 0.15, 0.15, 0.279],
[0.25, 0.4, 0.175, 0.15, 0.319],
[0.3, 0.4, 0.175, 0.15, 0.288],
[0.35, 0.4, 0.175, 0.15, 0.274],
[0.4, 0.4, 0.175, 0.15, 0.241],
[0.5, 0.4, 0.175, 0.15, 0.198],
[0.3, 0.4, 0.2, 0.15, 0.287],
[0.35, 0.4, 0.2, 0.15, 0.212],
[0.4, 0.4, 0.2, 0.15, 0.219],
[0.5, 0.4, 0.2, 0.15, 0.21],
[0.2, 0.4, 0.1, 0.15, 0.232],
[0.25, 0.4, 0.1, 0.15, 0.299],
[0.3, 0.4, 0.1, 0.15, 0.334],
[0.35, 0.4, 0.1, 0.15, 0.284],
[0.2, 0.4, 0.125, 0.15, 0.328],
[0.25, 0.4, 0.125, 0.15, 0.232],
[0.3, 0.4, 0.125, 0.15, 0.34],
[0.5, 0.4, 0.125, 0.15, 0.326],
[0.25, 0.4, 0.15, 0.15, 0.286],
[0.2, 0.4, 0.175, 0.15, 0.335],
[0.3, 0.4, 0.175, 0.15, 0.291],
[0.25, 0.4, 0.2, 0.15, 0.186],
[0.4, 0.5, 0.175, 0.15, 0.319],
[0.5, 0.5, 0.175, 0.15, 0.244],
[0.35, 0.5, 0.2, 0.15, 0.315],
[0.4, 0.5, 0.2, 0.15, 0.293],
[0.5, 0.5, 0.2, 0.15, 0.257],
[0.3, 0.5, 0.1, 0.15, 0.251],
[0.35, 0.5, 0.1, 0.15, 0.182],
[0.5, 0.5, 0.125, 0.15, 0.024],
[0.05, 0.05, 0.1, 0.175, 0.203],
[0.1, 0.05, 0.1, 0.175, 0.299],
[0.05, 0.05, 0.125, 0.175, 0.141],
[0.05, 0.05, 0.15, 0.175, 0.185],
[0.05, 0.05, 0.175, 0.175, 0.11],
[0.05, 0.05, 0.2, 0.175, 0.229],
[0.05, 0.05, 0.1, 0.175, 0.233],
[0.05, 0.05, 0.125, 0.175, 0.279],
[0.05, 0.05, 0.15, 0.175, 0.142],
[0.05, 0.05, 0.175, 0.175, 0.124],
[0.05, 0.05, 0.2, 0.175, 0.18],
[0.1, 0.1, 0.1, 0.175, 0.223],
[0.15, 0.1, 0.1, 0.175, 0.192],
[0.2, 0.1, 0.1, 0.175, 0.253],
[0.1, 0.1, 0.125, 0.175, 0.173],
[0.15, 0.1, 0.125, 0.175, 0.178],
[0.2, 0.1, 0.125, 0.175, 0.242],
[0.1, 0.1, 0.15, 0.175, 0.151],
[0.15, 0.1, 0.15, 0.175, 0.272],
[0.1, 0.1, 0.175, 0.175, 0.149],
[0.15, 0.1, 0.175, 0.175, 0.269],
[0.1, 0.1, 0.2, 0.175, 0.289],
[0.05, 0.1, 0.1, 0.175, 0.206],
[0.1, 0.1, 0.1, 0.175, 0.285],
[0.05, 0.1, 0.125, 0.175, 0.214],
[0.1, 0.1, 0.125, 0.175, 0.135],
[0.05, 0.1, 0.15, 0.175, 0.266],
[0.1, 0.1, 0.15, 0.175, 0.237],
[0.05, 0.1, 0.175, 0.175, 0.306],
[0.1, 0.1, 0.175, 0.175, 0.153],
[0.05, 0.1, 0.2, 0.175, 0.315],
[0.1, 0.1, 0.2, 0.175, 0.323],
[0.15, 0.15, 0.1, 0.175, 0.306],
[0.2, 0.15, 0.1, 0.175, 0.298],
[0.25, 0.15, 0.1, 0.175, 0.297],
[0.3, 0.15, 0.1, 0.175, 0.335],
[0.15, 0.15, 0.125, 0.175, 0.197],
[0.2, 0.15, 0.125, 0.175, 0.258],
[0.25, 0.15, 0.125, 0.175, 0.189],
[0.3, 0.15, 0.125, 0.175, 0.293],
[0.15, 0.15, 0.15, 0.175, 0.114],
[0.2, 0.15, 0.15, 0.175, 0.079],
[0.25, 0.15, 0.15, 0.175, 0.189],
[0.15, 0.15, 0.175, 0.175, 0.168],
[0.2, 0.15, 0.175, 0.175, 0.044],
[0.15, 0.15, 0.2, 0.175, 0.157],
[0.2, 0.15, 0.2, 0.175, 0.165],
[0.25, 0.15, 0.2, 0.175, 0.193],
[0.1, 0.15, 0.1, 0.175, 0.145],
[0.15, 0.15, 0.1, 0.175, 0.128],
[0.1, 0.15, 0.125, 0.175, 0.216],
[0.15, 0.15, 0.125, 0.175, 0.165],
[0.1, 0.15, 0.15, 0.175, 0.147],
[0.15, 0.15, 0.15, 0.175, 0.151],
[0.1, 0.15, 0.175, 0.175, 0.219],
[0.15, 0.15, 0.175, 0.175, 0.16],
[0.1, 0.15, 0.2, 0.175, 0.075],
[0.15, 0.15, 0.2, 0.175, 0.176],
[0.25, 0.2, 0.125, 0.175, 0.2],
[0.3, 0.2, 0.125, 0.175, 0.117],
[0.35, 0.2, 0.125, 0.175, 0.297],
[0.4, 0.2, 0.125, 0.175, 0.308],
[0.15, 0.2, 0.15, 0.175, 0.307],
[0.2, 0.2, 0.15, 0.175, 0.185],
[0.25, 0.2, 0.15, 0.175, 0.157],
[0.3, 0.2, 0.15, 0.175, 0.11],
[0.35, 0.2, 0.15, 0.175, 0.187],
[0.4, 0.2, 0.15, 0.175, 0.266],
[0.15, 0.2, 0.175, 0.175, 0.34],
[0.2, 0.2, 0.175, 0.175, 0.224],
[0.25, 0.2, 0.175, 0.175, 0.127],
[0.3, 0.2, 0.175, 0.175, 0.158],
[0.35, 0.2, 0.175, 0.175, 0.254],
[0.4, 0.2, 0.175, 0.175, 0.283],
[0.15, 0.2, 0.2, 0.175, 0.29],
[0.2, 0.2, 0.2, 0.175, 0.103],
[0.25, 0.2, 0.2, 0.175, 0.233],
[0.3, 0.2, 0.2, 0.175, 0.2],
[0.35, 0.2, 0.2, 0.175, 0.272],
[0.4, 0.2, 0.2, 0.175, 0.33],
[0.1, 0.2, 0.1, 0.175, 0.206],
[0.15, 0.2, 0.1, 0.175, 0.074],
[0.2, 0.2, 0.1, 0.175, 0.204],
[0.1, 0.2, 0.125, 0.175, 0.312],
[0.15, 0.2, 0.125, 0.175, 0.065],
[0.2, 0.2, 0.125, 0.175, 0.184],
[0.1, 0.2, 0.15, 0.175, 0.244],
[0.15, 0.2, 0.15, 0.175, 0.119],
[0.2, 0.2, 0.15, 0.175, 0.259],
[0.1, 0.2, 0.175, 0.175, 0.298],
[0.15, 0.2, 0.175, 0.175, 0.046],
[0.2, 0.2, 0.175, 0.175, 0.287],
[0.1, 0.2, 0.2, 0.175, 0.289],
[0.15, 0.2, 0.2, 0.175, 0.081],
[0.2, 0.2, 0.2, 0.175, 0.224],
[0.25, 0.25, 0.125, 0.175, 0.332],
[0.3, 0.25, 0.125, 0.175, 0.239],
[0.2, 0.25, 0.15, 0.175, 0.321],
[0.25, 0.25, 0.15, 0.175, 0.295],
[0.3, 0.25, 0.15, 0.175, 0.176],
[0.35, 0.25, 0.15, 0.175, 0.156],
[0.4, 0.25, 0.15, 0.175, 0.288],
[0.25, 0.25, 0.175, 0.175, 0.16],
[0.3, 0.25, 0.175, 0.175, 0.128],
[0.35, 0.25, 0.175, 0.175, 0.217],
[0.4, 0.25, 0.175, 0.175, 0.208],
[0.5, 0.25, 0.175, 0.175, 0.332],
[0.2, 0.25, 0.2, 0.175, 0.246],
[0.3, 0.25, 0.2, 0.175, 0.023],
[0.35, 0.25, 0.2, 0.175, 0.137],
[0.4, 0.25, 0.2, 0.175, 0.231],
[0.5, 0.25, 0.2, 0.175, 0.323],
[0.2, 0.25, 0.1, 0.175, 0.164],
[0.25, 0.25, 0.1, 0.175, 0.092],
[0.15, 0.25, 0.125, 0.175, 0.309],
[0.2, 0.25, 0.125, 0.175, 0.168],
[0.25, 0.25, 0.125, 0.175, 0.033],
[0.3, 0.25, 0.125, 0.175, 0.251],
[0.35, 0.25, 0.125, 0.175, 0.332],
[0.15, 0.25, 0.15, 0.175, 0.332],
[0.2, 0.25, 0.15, 0.175, 0.156],
[0.25, 0.25, 0.15, 0.175, 0.085],
[0.3, 0.25, 0.15, 0.175, 0.248],
[0.35, 0.25, 0.15, 0.175, 0.3],
[0.15, 0.25, 0.175, 0.175, 0.315],
[0.2, 0.25, 0.175, 0.175, 0.13],
[0.25, 0.25, 0.175, 0.175, 0.183],
[0.3, 0.25, 0.175, 0.175, 0.236],
[0.15, 0.25, 0.2, 0.175, 0.303],
[0.2, 0.25, 0.2, 0.175, 0.19],
[0.25, 0.25, 0.2, 0.175, 0.022],
[0.3, 0.25, 0.2, 0.175, 0.235],
[0.3, 0.3, 0.125, 0.175, 0.277],
[0.25, 0.3, 0.15, 0.175, 0.343],
[0.3, 0.3, 0.15, 0.175, 0.214],
[0.35, 0.3, 0.15, 0.175, 0.145],
[0.4, 0.3, 0.15, 0.175, 0.241],
[0.25, 0.3, 0.175, 0.175, 0.191],
[0.3, 0.3, 0.175, 0.175, 0.154],
[0.35, 0.3, 0.175, 0.175, 0.152],
[0.4, 0.3, 0.175, 0.175, 0.116],
[0.5, 0.3, 0.175, 0.175, 0.167],
[0.3, 0.3, 0.2, 0.175, 0.173],
[0.35, 0.3, 0.2, 0.175, 0.065],
[0.4, 0.3, 0.2, 0.175, 0.106],
[0.5, 0.3, 0.2, 0.175, 0.166],
[0.2, 0.3, 0.1, 0.175, 0.205],
[0.25, 0.3, 0.1, 0.175, 0.064],
[0.15, 0.3, 0.125, 0.175, 0.331],
[0.2, 0.3, 0.125, 0.175, 0.198],
[0.25, 0.3, 0.125, 0.175, 0.029],
[0.3, 0.3, 0.125, 0.175, 0.192],
[0.35, 0.3, 0.125, 0.175, 0.272],
[0.2, 0.3, 0.15, 0.175, 0.173],
[0.25, 0.3, 0.15, 0.175, 0.053],
[0.3, 0.3, 0.15, 0.175, 0.189],
[0.35, 0.3, 0.15, 0.175, 0.239],
[0.15, 0.3, 0.175, 0.175, 0.338],
[0.2, 0.3, 0.175, 0.175, 0.147],
[0.25, 0.3, 0.175, 0.175, 0.144],
[0.3, 0.3, 0.175, 0.175, 0.181],
[0.35, 0.3, 0.175, 0.175, 0.33],
[0.15, 0.3, 0.2, 0.175, 0.323],
[0.2, 0.3, 0.2, 0.175, 0.207],
[0.25, 0.3, 0.2, 0.175, 0.039],
[0.3, 0.3, 0.2, 0.175, 0.177],
[0.35, 0.3, 0.2, 0.175, 0.34],
[0.3, 0.35, 0.15, 0.175, 0.304],
[0.35, 0.35, 0.15, 0.175, 0.245],
[0.4, 0.35, 0.15, 0.175, 0.331],
[0.25, 0.35, 0.175, 0.175, 0.284],
[0.3, 0.35, 0.175, 0.175, 0.248],
[0.35, 0.35, 0.175, 0.175, 0.252],
[0.4, 0.35, 0.175, 0.175, 0.218],
[0.5, 0.35, 0.175, 0.175, 0.204],
[0.3, 0.35, 0.2, 0.175, 0.23],
[0.35, 0.35, 0.2, 0.175, 0.168],
[0.4, 0.35, 0.2, 0.175, 0.2],
[0.5, 0.35, 0.2, 0.175, 0.212],
[0.2, 0.35, 0.1, 0.175, 0.234],
[0.25, 0.35, 0.1, 0.175, 0.198],
[0.2, 0.35, 0.125, 0.175, 0.284],
[0.25, 0.35, 0.125, 0.175, 0.146],
[0.3, 0.35, 0.125, 0.175, 0.213],
[0.35, 0.35, 0.125, 0.175, 0.248],
[0.4, 0.35, 0.125, 0.175, 0.316],
[0.2, 0.35, 0.15, 0.175, 0.29],
[0.25, 0.35, 0.15, 0.175, 0.185],
[0.3, 0.35, 0.15, 0.175, 0.243],
[0.35, 0.35, 0.15, 0.175, 0.236],
[0.2, 0.35, 0.175, 0.175, 0.264],
[0.25, 0.35, 0.175, 0.175, 0.277],
[0.3, 0.35, 0.175, 0.175, 0.171],
[0.35, 0.35, 0.175, 0.175, 0.314],
[0.2, 0.35, 0.2, 0.175, 0.324],
[0.25, 0.35, 0.2, 0.175, 0.11],
[0.3, 0.35, 0.2, 0.175, 0.243],
[0.35, 0.35, 0.2, 0.175, 0.296],
[0.35, 0.4, 0.15, 0.175, 0.338],
[0.3, 0.4, 0.175, 0.175, 0.342],
[0.35, 0.4, 0.175, 0.175, 0.338],
[0.4, 0.4, 0.175, 0.175, 0.308],
[0.5, 0.4, 0.175, 0.175, 0.265],
[0.3, 0.4, 0.2, 0.175, 0.324],
[0.35, 0.4, 0.2, 0.175, 0.272],
[0.4, 0.4, 0.2, 0.175, 0.288],
[0.5, 0.4, 0.2, 0.175, 0.277],
[0.2, 0.4, 0.1, 0.175, 0.26],
[0.25, 0.4, 0.1, 0.175, 0.193],
[0.2, 0.4, 0.125, 0.175, 0.298],
[0.25, 0.4, 0.125, 0.175, 0.152],
[0.3, 0.4, 0.125, 0.175, 0.171],
[0.35, 0.4, 0.125, 0.175, 0.198],
[0.4, 0.4, 0.125, 0.175, 0.267],
[0.2, 0.4, 0.15, 0.175, 0.297],
[0.25, 0.4, 0.15, 0.175, 0.18],
[0.3, 0.4, 0.15, 0.175, 0.206],
[0.35, 0.4, 0.15, 0.175, 0.189],
[0.4, 0.4, 0.15, 0.175, 0.309],
[0.2, 0.4, 0.175, 0.175, 0.271],
[0.25, 0.4, 0.175, 0.175, 0.26],
[0.3, 0.4, 0.175, 0.175, 0.126],
[0.35, 0.4, 0.175, 0.175, 0.262],
[0.2, 0.4, 0.2, 0.175, 0.33],
[0.25, 0.4, 0.2, 0.175, 0.124],
[0.3, 0.4, 0.2, 0.175, 0.209],
[0.35, 0.4, 0.2, 0.175, 0.244],
[0.25, 0.5, 0.125, 0.175, 0.308],
[0.3, 0.5, 0.125, 0.175, 0.252],
[0.35, 0.5, 0.125, 0.175, 0.199],
[0.4, 0.5, 0.125, 0.175, 0.125],
[0.25, 0.5, 0.15, 0.175, 0.337],
[0.3, 0.5, 0.15, 0.175, 0.308],
[0.35, 0.5, 0.15, 0.175, 0.228],
[0.4, 0.5, 0.15, 0.175, 0.232],
[0.5, 0.5, 0.15, 0.175, 0.307],
[0.3, 0.5, 0.175, 0.175, 0.197],
[0.35, 0.5, 0.175, 0.175, 0.248],
[0.4, 0.5, 0.175, 0.175, 0.326],
[0.25, 0.5, 0.2, 0.175, 0.277],
[0.3, 0.5, 0.2, 0.175, 0.321],
[0.35, 0.5, 0.2, 0.175, 0.179],
[0.05, 0.05, 0.1, 0.2, 0.25],
[0.1, 0.05, 0.1, 0.2, 0.306],
[0.05, 0.05, 0.125, 0.2, 0.179],
[0.1, 0.05, 0.125, 0.2, 0.343],
[0.05, 0.05, 0.15, 0.2, 0.236],
[0.05, 0.05, 0.175, 0.2, 0.153],
[0.1, 0.05, 0.175, 0.2, 0.318],
[0.05, 0.05, 0.2, 0.2, 0.265],
[0.05, 0.05, 0.15, 0.2, 0.303],
[0.05, 0.05, 0.175, 0.2, 0.296],
[0.15, 0.1, 0.1, 0.2, 0.32],
[0.2, 0.1, 0.1, 0.2, 0.342],
[0.1, 0.1, 0.125, 0.2, 0.332],
[0.15, 0.1, 0.125, 0.2, 0.23],
[0.2, 0.1, 0.125, 0.2, 0.31],
[0.1, 0.1, 0.15, 0.2, 0.208],
[0.15, 0.1, 0.15, 0.2, 0.121],
[0.2, 0.1, 0.15, 0.2, 0.317],
[0.1, 0.1, 0.175, 0.2, 0.245],
[0.15, 0.1, 0.175, 0.2, 0.077],
[0.2, 0.1, 0.175, 0.2, 0.269],
[0.1, 0.1, 0.2, 0.2, 0.196],
[0.15, 0.1, 0.2, 0.2, 0.172],
[0.05, 0.1, 0.1, 0.2, 0.273],
[0.1, 0.1, 0.1, 0.2, 0.267],
[0.05, 0.1, 0.125, 0.2, 0.233],
[0.1, 0.1, 0.125, 0.2, 0.206],
[0.05, 0.1, 0.15, 0.2, 0.331],
[0.1, 0.1, 0.15, 0.2, 0.227],
[0.05, 0.1, 0.175, 0.2, 0.326],
[0.1, 0.1, 0.175, 0.2, 0.224],
[0.05, 0.1, 0.2, 0.2, 0.3],
[0.1, 0.1, 0.2, 0.2, 0.205],
[0.15, 0.15, 0.1, 0.2, 0.343],
[0.2, 0.15, 0.1, 0.2, 0.336],
[0.25, 0.15, 0.1, 0.2, 0.323],
[0.15, 0.15, 0.125, 0.2, 0.238],
[0.2, 0.15, 0.125, 0.2, 0.297],
[0.25, 0.15, 0.125, 0.2, 0.189],
[0.3, 0.15, 0.125, 0.2, 0.263],
[0.15, 0.15, 0.15, 0.2, 0.127],
[0.2, 0.15, 0.15, 0.2, 0.078],
[0.25, 0.15, 0.15, 0.2, 0.176],
[0.15, 0.15, 0.175, 0.2, 0.168],
[0.2, 0.15, 0.175, 0.2, 0.076],
[0.15, 0.15, 0.2, 0.2, 0.131],
[0.2, 0.15, 0.2, 0.2, 0.123],
[0.25, 0.15, 0.2, 0.2, 0.202],
[0.1, 0.15, 0.1, 0.2, 0.216],
[0.15, 0.15, 0.1, 0.2, 0.134],
[0.2, 0.15, 0.1, 0.2, 0.32],
[0.1, 0.15, 0.125, 0.2, 0.288],
[0.15, 0.15, 0.125, 0.2, 0.122],
[0.2, 0.15, 0.125, 0.2, 0.264],
[0.1, 0.15, 0.15, 0.2, 0.223],
[0.15, 0.15, 0.15, 0.2, 0.084],
[0.2, 0.15, 0.15, 0.2, 0.304],
[0.1, 0.15, 0.175, 0.2, 0.29],
[0.15, 0.15, 0.175, 0.2, 0.135],
[0.2, 0.15, 0.175, 0.2, 0.336],
[0.1, 0.15, 0.2, 0.2, 0.148],
[0.15, 0.15, 0.2, 0.2, 0.119],
[0.2, 0.15, 0.2, 0.2, 0.26],
[0.25, 0.2, 0.125, 0.2, 0.192],
[0.3, 0.2, 0.125, 0.2, 0.145],
[0.35, 0.2, 0.125, 0.2, 0.314],
[0.4, 0.2, 0.125, 0.2, 0.329],
[0.15, 0.2, 0.15, 0.2, 0.258],
[0.2, 0.2, 0.15, 0.2, 0.15],
[0.25, 0.2, 0.15, 0.2, 0.149],
[0.3, 0.2, 0.15, 0.2, 0.17],
[0.35, 0.2, 0.15, 0.2, 0.258],
[0.4, 0.2, 0.15, 0.2, 0.318],
[0.15, 0.2, 0.175, 0.2, 0.287],
[0.2, 0.2, 0.175, 0.2, 0.187],
[0.25, 0.2, 0.175, 0.2, 0.193],
[0.3, 0.2, 0.175, 0.2, 0.23],
[0.35, 0.2, 0.175, 0.2, 0.323],
[0.15, 0.2, 0.2, 0.2, 0.232],
[0.2, 0.2, 0.2, 0.2, 0.054],
[0.25, 0.2, 0.2, 0.2, 0.224],
[0.3, 0.2, 0.2, 0.2, 0.267],
[0.1, 0.2, 0.1, 0.2, 0.302],
[0.15, 0.2, 0.1, 0.2, 0.178],
[0.2, 0.2, 0.1, 0.2, 0.142],
[0.25, 0.2, 0.1, 0.2, 0.261],
[0.15, 0.2, 0.125, 0.2, 0.124],
[0.2, 0.2, 0.125, 0.2, 0.051],
[0.25, 0.2, 0.125, 0.2, 0.248],
[0.1, 0.2, 0.15, 0.2, 0.328],
[0.15, 0.2, 0.15, 0.2, 0.134],
[0.2, 0.2, 0.15, 0.2, 0.088],
[0.25, 0.2, 0.15, 0.2, 0.271],
[0.15, 0.2, 0.175, 0.2, 0.142],
[0.2, 0.2, 0.175, 0.2, 0.115],
[0.25, 0.2, 0.175, 0.2, 0.323],
[0.1, 0.2, 0.2, 0.2, 0.33],
[0.15, 0.2, 0.2, 0.2, 0.111],
[0.2, 0.2, 0.2, 0.2, 0.056],
[0.25, 0.2, 0.2, 0.2, 0.261],
[0.25, 0.25, 0.125, 0.2, 0.231],
[0.3, 0.25, 0.125, 0.2, 0.122],
[0.35, 0.25, 0.125, 0.2, 0.259],
[0.4, 0.25, 0.125, 0.2, 0.26],
[0.2, 0.25, 0.15, 0.2, 0.272],
[0.25, 0.25, 0.15, 0.2, 0.203],
[0.3, 0.25, 0.15, 0.2, 0.047],
[0.35, 0.25, 0.15, 0.2, 0.04],
[0.4, 0.25, 0.15, 0.2, 0.157],
[0.2, 0.25, 0.175, 0.2, 0.309],
[0.25, 0.25, 0.175, 0.2, 0.022],
[0.3, 0.25, 0.175, 0.2, 0.04],
[0.35, 0.25, 0.175, 0.2, 0.107],
[0.4, 0.25, 0.175, 0.2, 0.13],
[0.5, 0.25, 0.175, 0.2, 0.307],
[0.2, 0.25, 0.2, 0.2, 0.23],
[0.25, 0.25, 0.2, 0.2, 0.258],
[0.3, 0.25, 0.2, 0.2, 0.167],
[0.35, 0.25, 0.2, 0.2, 0.141],
[0.4, 0.25, 0.2, 0.2, 0.175],
[0.5, 0.25, 0.2, 0.2, 0.29],
[0.15, 0.25, 0.1, 0.2, 0.25],
[0.2, 0.25, 0.1, 0.2, 0.107],
[0.3, 0.25, 0.1, 0.2, 0.284],
[0.35, 0.25, 0.1, 0.2, 0.304],
[0.15, 0.25, 0.125, 0.2, 0.24],
[0.2, 0.25, 0.125, 0.2, 0.19],
[0.15, 0.25, 0.15, 0.2, 0.299],
[0.2, 0.25, 0.15, 0.2, 0.273],
[0.15, 0.25, 0.175, 0.2, 0.225],
[0.2, 0.25, 0.175, 0.2, 0.279],
[0.15, 0.25, 0.2, 0.2, 0.25],
[0.2, 0.25, 0.2, 0.2, 0.27],
[0.25, 0.25, 0.2, 0.2, 0.339],
[0.3, 0.3, 0.125, 0.2, 0.285],
[0.3, 0.3, 0.15, 0.2, 0.219],
[0.35, 0.3, 0.15, 0.2, 0.167],
[0.4, 0.3, 0.15, 0.2, 0.29],
[0.25, 0.3, 0.175, 0.2, 0.198],
[0.3, 0.3, 0.175, 0.2, 0.155],
[0.35, 0.3, 0.175, 0.2, 0.204],
[0.4, 0.3, 0.175, 0.2, 0.177],
[0.5, 0.3, 0.175, 0.2, 0.248],
[0.2, 0.3, 0.2, 0.2, 0.319],
[0.3, 0.3, 0.2, 0.2, 0.095],
[0.35, 0.3, 0.2, 0.2, 0.088],
[0.4, 0.3, 0.2, 0.2, 0.18],
[0.5, 0.3, 0.2, 0.2, 0.246],
[0.2, 0.3, 0.1, 0.2, 0.281],
[0.25, 0.3, 0.1, 0.2, 0.072],
[0.15, 0.3, 0.125, 0.2, 0.343],
[0.2, 0.3, 0.125, 0.2, 0.238],
[0.25, 0.3, 0.125, 0.2, 0.12],
[0.3, 0.3, 0.125, 0.2, 0.152],
[0.35, 0.3, 0.125, 0.2, 0.239],
[0.15, 0.3, 0.15, 0.2, 0.34],
[0.2, 0.3, 0.15, 0.2, 0.185],
[0.25, 0.3, 0.15, 0.2, 0.074],
[0.3, 0.3, 0.15, 0.2, 0.122],
[0.35, 0.3, 0.15, 0.2, 0.201],
[0.2, 0.3, 0.175, 0.2, 0.167],
[0.25, 0.3, 0.175, 0.2, 0.024],
[0.3, 0.3, 0.175, 0.2, 0.17],
[0.35, 0.3, 0.175, 0.2, 0.278],
[0.15, 0.3, 0.2, 0.2, 0.331],
[0.2, 0.3, 0.2, 0.2, 0.211],
[0.25, 0.3, 0.2, 0.2, 0.148],
[0.3, 0.3, 0.2, 0.2, 0.103],
[0.35, 0.3, 0.2, 0.2, 0.308],
[0.3, 0.35, 0.15, 0.2, 0.307],
[0.35, 0.35, 0.15, 0.2, 0.244],
[0.4, 0.35, 0.15, 0.2, 0.306],
[0.25, 0.35, 0.175, 0.2, 0.288],
[0.3, 0.35, 0.175, 0.2, 0.257],
[0.35, 0.35, 0.175, 0.2, 0.233],
[0.4, 0.35, 0.175, 0.2, 0.199],
[0.5, 0.35, 0.175, 0.2, 0.157],
[0.3, 0.35, 0.2, 0.2, 0.272],
[0.35, 0.35, 0.2, 0.2, 0.179],
[0.4, 0.35, 0.2, 0.2, 0.174],
[0.5, 0.35, 0.2, 0.2, 0.168],
[0.2, 0.35, 0.1, 0.2, 0.28],
[0.25, 0.35, 0.1, 0.2, 0.122],
[0.2, 0.35, 0.125, 0.2, 0.274],
[0.25, 0.35, 0.125, 0.2, 0.122],
[0.3, 0.35, 0.125, 0.2, 0.081],
[0.35, 0.35, 0.125, 0.2, 0.149],
[0.4, 0.35, 0.125, 0.2, 0.265],
[0.2, 0.35, 0.15, 0.2, 0.246],
[0.25, 0.35, 0.15, 0.2, 0.111],
[0.3, 0.35, 0.15, 0.2, 0.096],
[0.35, 0.35, 0.15, 0.2, 0.12],
[0.4, 0.35, 0.15, 0.2, 0.268],
[0.2, 0.35, 0.175, 0.2, 0.223],
[0.25, 0.35, 0.175, 0.2, 0.151],
[0.3, 0.35, 0.175, 0.2, 0.069],
[0.35, 0.35, 0.175, 0.2, 0.202],
[0.2, 0.35, 0.2, 0.2, 0.277],
[0.25, 0.35, 0.2, 0.2, 0.124],
[0.3, 0.35, 0.2, 0.2, 0.094],
[0.35, 0.35, 0.2, 0.2, 0.213],
[0.4, 0.35, 0.2, 0.2, 0.344],
[0.35, 0.4, 0.15, 0.2, 0.332],
[0.3, 0.4, 0.175, 0.2, 0.331],
[0.35, 0.4, 0.175, 0.2, 0.342],
[0.4, 0.4, 0.175, 0.2, 0.311],
[0.5, 0.4, 0.175, 0.2, 0.285],
[0.3, 0.4, 0.2, 0.2, 0.292],
[0.35, 0.4, 0.2, 0.2, 0.26],
[0.4, 0.4, 0.2, 0.2, 0.295],
[0.5, 0.4, 0.2, 0.2, 0.296],
[0.2, 0.4, 0.1, 0.2, 0.291],
[0.25, 0.4, 0.1, 0.2, 0.226],
[0.2, 0.4, 0.125, 0.2, 0.332],
[0.25, 0.4, 0.125, 0.2, 0.189],
[0.3, 0.4, 0.125, 0.2, 0.168],
[0.35, 0.4, 0.125, 0.2, 0.171],
[0.4, 0.4, 0.125, 0.2, 0.218],
[0.2, 0.4, 0.15, 0.2, 0.332],
[0.25, 0.4, 0.15, 0.2, 0.213],
[0.3, 0.4, 0.15, 0.2, 0.212],
[0.35, 0.4, 0.15, 0.2, 0.172],
[0.4, 0.4, 0.15, 0.2, 0.268],
[0.2, 0.4, 0.175, 0.2, 0.306],
[0.25, 0.4, 0.175, 0.2, 0.288],
[0.3, 0.4, 0.175, 0.2, 0.116],
[0.35, 0.4, 0.175, 0.2, 0.234],
[0.25, 0.4, 0.2, 0.2, 0.161],
[0.3, 0.4, 0.2, 0.2, 0.219],
[0.35, 0.4, 0.2, 0.2, 0.204],
[0.5, 0.5, 0.175, 0.2, 0.298],
[0.35, 0.5, 0.2, 0.2, 0.333],
[0.4, 0.5, 0.2, 0.2, 0.333],
[0.5, 0.5, 0.2, 0.2, 0.31],
[0.25, 0.5, 0.125, 0.2, 0.301],
[0.3, 0.5, 0.125, 0.2, 0.254],
[0.35, 0.5, 0.125, 0.2, 0.206],
[0.4, 0.5, 0.125, 0.2, 0.144],
[0.25, 0.5, 0.15, 0.2, 0.332],
[0.3, 0.5, 0.15, 0.2, 0.309],
[0.35, 0.5, 0.15, 0.2, 0.233],
[0.4, 0.5, 0.15, 0.2, 0.247],
[0.5, 0.5, 0.15, 0.2, 0.329],
[0.3, 0.5, 0.175, 0.2, 0.197],
[0.35, 0.5, 0.175, 0.2, 0.258],
[0.25, 0.5, 0.2, 0.2, 0.268],
[0.3, 0.5, 0.2, 0.2, 0.322],
[0.35, 0.5, 0.2, 0.2, 0.192],
]

list0_10=[
[0.05, 0.05, 0.175, 0.1, 0.2],
[0.15, 0.1, 0.1, 0.1, 0.199],
[0.15, 0.1, 0.125, 0.1, 0.102],
[0.2, 0.1, 0.125, 0.1, 0.194],
[0.15, 0.1, 0.15, 0.1, 0.141],
[0.15, 0.1, 0.175, 0.1, 0.162],
[0.15, 0.15, 0.125, 0.1, 0.169],
[0.15, 0.15, 0.15, 0.1, 0.082],
[0.2, 0.15, 0.15, 0.1, 0.123],
[0.15, 0.15, 0.175, 0.1, 0.141],
[0.2, 0.15, 0.175, 0.1, 0.074],
[0.15, 0.15, 0.2, 0.1, 0.152],
[0.25, 0.2, 0.125, 0.1, 0.065],
[0.3, 0.2, 0.125, 0.1, 0.066],
[0.35, 0.2, 0.125, 0.1, 0.189],
[0.2, 0.2, 0.15, 0.1, 0.151],
[0.25, 0.2, 0.15, 0.1, 0.043],
[0.3, 0.2, 0.15, 0.1, 0.16],
[0.2, 0.2, 0.175, 0.1, 0.188],
[0.25, 0.2, 0.175, 0.1, 0.193],
[0.2, 0.2, 0.2, 0.1, 0.164],
[0.25, 0.2, 0.2, 0.1, 0.098],
[0.3, 0.25, 0.125, 0.1, 0.17],
[0.3, 0.25, 0.15, 0.1, 0.102],
[0.35, 0.25, 0.15, 0.1, 0.096],
[0.25, 0.25, 0.175, 0.1, 0.085],
[0.3, 0.25, 0.175, 0.1, 0.06],
[0.35, 0.25, 0.175, 0.1, 0.168],
[0.4, 0.25, 0.175, 0.1, 0.173],
[0.3, 0.25, 0.2, 0.1, 0.076],
[0.35, 0.25, 0.2, 0.1, 0.13],
[0.3, 0.25, 0.1, 0.1, 0.187],
[0.3, 0.3, 0.1, 0.1, 0.136],
[0.3, 0.3, 0.1, 0.1, 0.125],
[0.3, 0.35, 0.1, 0.1, 0.2],
[0.35, 0.35, 0.125, 0.1, 0.198],
[0.4, 0.35, 0.125, 0.1, 0.183],
[0.5, 0.35, 0.15, 0.1, 0.176],
[0.05, 0.05, 0.15, 0.125, 0.181],
[0.05, 0.05, 0.175, 0.125, 0.163],
[0.1, 0.1, 0.125, 0.125, 0.174],
[0.15, 0.1, 0.125, 0.125, 0.198],
[0.1, 0.1, 0.15, 0.125, 0.124],
[0.1, 0.1, 0.175, 0.125, 0.128],
[0.25, 0.15, 0.125, 0.125, 0.19],
[0.15, 0.15, 0.15, 0.125, 0.117],
[0.2, 0.15, 0.15, 0.125, 0.077],
[0.25, 0.15, 0.15, 0.125, 0.181],
[0.15, 0.15, 0.175, 0.125, 0.164],
[0.2, 0.15, 0.175, 0.125, 0.061],
[0.15, 0.15, 0.2, 0.125, 0.138],
[0.2, 0.15, 0.2, 0.125, 0.141],
[0.25, 0.15, 0.2, 0.125, 0.198],
[0.1, 0.15, 0.1, 0.125, 0.113],
[0.15, 0.15, 0.1, 0.125, 0.195],
[0.1, 0.15, 0.15, 0.125, 0.159],
[0.1, 0.15, 0.175, 0.125, 0.18],
[0.25, 0.2, 0.125, 0.125, 0.133],
[0.3, 0.2, 0.125, 0.125, 0.104],
[0.2, 0.2, 0.15, 0.125, 0.112],
[0.25, 0.2, 0.15, 0.125, 0.087],
[0.3, 0.2, 0.15, 0.125, 0.162],
[0.2, 0.2, 0.175, 0.125, 0.155],
[0.25, 0.2, 0.175, 0.125, 0.192],
[0.2, 0.2, 0.2, 0.125, 0.052],
[0.25, 0.2, 0.2, 0.125, 0.167],
[0.3, 0.25, 0.125, 0.125, 0.131],
[0.3, 0.25, 0.15, 0.125, 0.054],
[0.35, 0.25, 0.15, 0.125, 0.06],
[0.4, 0.25, 0.15, 0.125, 0.188],
[0.25, 0.25, 0.175, 0.125, 0.032],
[0.3, 0.25, 0.175, 0.125, 0.03],
[0.35, 0.25, 0.175, 0.125, 0.135],
[0.4, 0.25, 0.175, 0.125, 0.151],
[0.3, 0.25, 0.2, 0.125, 0.132],
[0.35, 0.25, 0.2, 0.125, 0.139],
[0.4, 0.25, 0.2, 0.125, 0.194],
[0.2, 0.25, 0.1, 0.125, 0.069],
[0.2, 0.25, 0.125, 0.125, 0.148],
[0.35, 0.3, 0.15, 0.125, 0.164],
[0.3, 0.3, 0.175, 0.125, 0.162],
[0.35, 0.3, 0.175, 0.125, 0.185],
[0.4, 0.3, 0.175, 0.125, 0.152],
[0.3, 0.3, 0.2, 0.125, 0.143],
[0.35, 0.3, 0.2, 0.125, 0.076],
[0.4, 0.3, 0.2, 0.125, 0.147],
[0.2, 0.3, 0.1, 0.125, 0.063],
[0.25, 0.3, 0.1, 0.125, 0.194],
[0.2, 0.3, 0.125, 0.125, 0.109],
[0.25, 0.3, 0.125, 0.125, 0.147],
[0.2, 0.3, 0.15, 0.125, 0.145],
[0.25, 0.3, 0.15, 0.125, 0.195],
[0.2, 0.3, 0.175, 0.125, 0.131],
[0.2, 0.3, 0.2, 0.125, 0.168],
[0.25, 0.3, 0.2, 0.125, 0.139],
[0.4, 0.35, 0.175, 0.125, 0.171],
[0.5, 0.35, 0.175, 0.125, 0.167],
[0.35, 0.35, 0.2, 0.125, 0.13],
[0.4, 0.35, 0.2, 0.125, 0.152],
[0.5, 0.35, 0.2, 0.125, 0.173],
[0.3, 0.35, 0.1, 0.125, 0.201],
[0.35, 0.35, 0.1, 0.125, 0.167],
[0.5, 0.4, 0.175, 0.125, 0.114],
[0.4, 0.4, 0.2, 0.125, 0.179],
[0.5, 0.4, 0.2, 0.125, 0.129],
[0.35, 0.4, 0.1, 0.125, 0.159],
[0.5, 0.5, 0.15, 0.125, 0.024],
[0.05, 0.05, 0.1, 0.15, 0.198],
[0.05, 0.05, 0.125, 0.15, 0.168],
[0.05, 0.05, 0.15, 0.15, 0.167],
[0.05, 0.05, 0.175, 0.15, 0.136],
[0.05, 0.05, 0.15, 0.15, 0.201],
[0.15, 0.1, 0.125, 0.15, 0.139],
[0.1, 0.1, 0.15, 0.15, 0.186],
[0.15, 0.1, 0.15, 0.15, 0.165],
[0.15, 0.1, 0.175, 0.15, 0.166],
[0.1, 0.1, 0.125, 0.15, 0.117],
[0.1, 0.1, 0.15, 0.15, 0.199],
[0.1, 0.1, 0.175, 0.15, 0.139],
[0.15, 0.15, 0.15, 0.15, 0.11],
[0.2, 0.15, 0.15, 0.15, 0.114],
[0.15, 0.15, 0.175, 0.15, 0.138],
[0.2, 0.15, 0.175, 0.15, 0.107],
[0.15, 0.15, 0.2, 0.15, 0.089],
[0.2, 0.15, 0.2, 0.15, 0.144],
[0.1, 0.15, 0.1, 0.15, 0.126],
[0.15, 0.15, 0.1, 0.15, 0.078],
[0.15, 0.15, 0.125, 0.15, 0.125],
[0.1, 0.15, 0.15, 0.15, 0.142],
[0.15, 0.15, 0.15, 0.15, 0.123],
[0.15, 0.15, 0.175, 0.15, 0.115],
[0.1, 0.15, 0.2, 0.15, 0.125],
[0.15, 0.15, 0.2, 0.15, 0.14],
[0.25, 0.2, 0.125, 0.15, 0.196],
[0.3, 0.2, 0.125, 0.15, 0.127],
[0.2, 0.2, 0.15, 0.15, 0.169],
[0.25, 0.2, 0.15, 0.15, 0.152],
[0.3, 0.2, 0.15, 0.15, 0.135],
[0.25, 0.2, 0.175, 0.15, 0.155],
[0.3, 0.2, 0.175, 0.15, 0.189],
[0.2, 0.2, 0.2, 0.15, 0.081],
[0.15, 0.2, 0.1, 0.15, 0.115],
[0.2, 0.2, 0.1, 0.15, 0.153],
[0.15, 0.2, 0.125, 0.15, 0.08],
[0.2, 0.2, 0.125, 0.15, 0.126],
[0.15, 0.2, 0.15, 0.15, 0.127],
[0.2, 0.2, 0.15, 0.15, 0.2],
[0.15, 0.2, 0.175, 0.15, 0.078],
[0.15, 0.2, 0.2, 0.15, 0.084],
[0.2, 0.2, 0.2, 0.15, 0.171],
[0.3, 0.25, 0.125, 0.15, 0.17],
[0.3, 0.25, 0.15, 0.15, 0.099],
[0.35, 0.25, 0.15, 0.15, 0.085],
[0.25, 0.25, 0.175, 0.15, 0.08],
[0.3, 0.25, 0.175, 0.15, 0.05],
[0.35, 0.25, 0.175, 0.15, 0.158],
[0.4, 0.25, 0.175, 0.15, 0.163],
[0.3, 0.25, 0.2, 0.15, 0.078],
[0.35, 0.25, 0.2, 0.15, 0.12],
[0.4, 0.25, 0.2, 0.15, 0.198],
[0.2, 0.25, 0.1, 0.15, 0.044],
[0.15, 0.25, 0.125, 0.15, 0.189],
[0.2, 0.25, 0.125, 0.15, 0.094],
[0.2, 0.25, 0.15, 0.15, 0.17],
[0.15, 0.25, 0.175, 0.15, 0.185],
[0.2, 0.25, 0.175, 0.15, 0.176],
[0.15, 0.25, 0.2, 0.15, 0.191],
[0.2, 0.25, 0.2, 0.15, 0.169],
[0.35, 0.3, 0.15, 0.15, 0.173],
[0.3, 0.3, 0.175, 0.15, 0.163],
[0.4, 0.3, 0.175, 0.15, 0.177],
[0.3, 0.3, 0.2, 0.15, 0.111],
[0.35, 0.3, 0.2, 0.15, 0.09],
[0.4, 0.3, 0.2, 0.15, 0.177],
[0.2, 0.3, 0.1, 0.15, 0.14],
[0.25, 0.3, 0.1, 0.15, 0.11],
[0.2, 0.3, 0.125, 0.15, 0.147],
[0.25, 0.3, 0.125, 0.15, 0.056],
[0.2, 0.3, 0.15, 0.15, 0.142],
[0.25, 0.3, 0.15, 0.15, 0.106],
[0.2, 0.3, 0.175, 0.15, 0.117],
[0.2, 0.3, 0.2, 0.15, 0.176],
[0.25, 0.3, 0.2, 0.15, 0.049],
[0.35, 0.35, 0.15, 0.15, 0.19],
[0.3, 0.35, 0.175, 0.15, 0.193],
[0.35, 0.35, 0.175, 0.15, 0.201],
[0.4, 0.35, 0.175, 0.15, 0.166],
[0.5, 0.35, 0.175, 0.15, 0.186],
[0.3, 0.35, 0.2, 0.15, 0.182],
[0.35, 0.35, 0.2, 0.15, 0.107],
[0.4, 0.35, 0.2, 0.15, 0.153],
[0.5, 0.35, 0.2, 0.15, 0.189],
[0.25, 0.35, 0.1, 0.15, 0.105],
[0.25, 0.35, 0.125, 0.15, 0.06],
[0.3, 0.35, 0.125, 0.15, 0.19],
[0.25, 0.35, 0.15, 0.15, 0.093],
[0.3, 0.35, 0.15, 0.15, 0.198],
[0.2, 0.35, 0.175, 0.15, 0.18],
[0.25, 0.35, 0.175, 0.15, 0.183],
[0.3, 0.35, 0.175, 0.15, 0.169],
[0.25, 0.35, 0.2, 0.15, 0.039],
[0.3, 0.35, 0.2, 0.15, 0.189],
[0.5, 0.4, 0.175, 0.15, 0.198],
[0.25, 0.4, 0.2, 0.15, 0.186],
[0.35, 0.5, 0.1, 0.15, 0.182],
[0.5, 0.5, 0.125, 0.15, 0.024],
[0.05, 0.05, 0.125, 0.175, 0.141],
[0.05, 0.05, 0.15, 0.175, 0.185],
[0.05, 0.05, 0.175, 0.175, 0.11],
[0.05, 0.05, 0.15, 0.175, 0.142],
[0.05, 0.05, 0.175, 0.175, 0.124],
[0.05, 0.05, 0.2, 0.175, 0.18],
[0.15, 0.1, 0.1, 0.175, 0.192],
[0.1, 0.1, 0.125, 0.175, 0.173],
[0.15, 0.1, 0.125, 0.175, 0.178],
[0.1, 0.1, 0.15, 0.175, 0.151],
[0.1, 0.1, 0.175, 0.175, 0.149],
[0.1, 0.1, 0.125, 0.175, 0.135],
[0.1, 0.1, 0.175, 0.175, 0.153],
[0.15, 0.15, 0.125, 0.175, 0.197],
[0.25, 0.15, 0.125, 0.175, 0.189],
[0.15, 0.15, 0.15, 0.175, 0.114],
[0.2, 0.15, 0.15, 0.175, 0.079],
[0.25, 0.15, 0.15, 0.175, 0.189],
[0.15, 0.15, 0.175, 0.175, 0.168],
[0.2, 0.15, 0.175, 0.175, 0.044],
[0.15, 0.15, 0.2, 0.175, 0.157],
[0.2, 0.15, 0.2, 0.175, 0.165],
[0.25, 0.15, 0.2, 0.175, 0.193],
[0.1, 0.15, 0.1, 0.175, 0.145],
[0.15, 0.15, 0.1, 0.175, 0.128],
[0.15, 0.15, 0.125, 0.175, 0.165],
[0.1, 0.15, 0.15, 0.175, 0.147],
[0.15, 0.15, 0.15, 0.175, 0.151],
[0.15, 0.15, 0.175, 0.175, 0.16],
[0.1, 0.15, 0.2, 0.175, 0.075],
[0.15, 0.15, 0.2, 0.175, 0.176],
[0.25, 0.2, 0.125, 0.175, 0.2],
[0.3, 0.2, 0.125, 0.175, 0.117],
[0.2, 0.2, 0.15, 0.175, 0.185],
[0.25, 0.2, 0.15, 0.175, 0.157],
[0.3, 0.2, 0.15, 0.175, 0.11],
[0.35, 0.2, 0.15, 0.175, 0.187],
[0.25, 0.2, 0.175, 0.175, 0.127],
[0.3, 0.2, 0.175, 0.175, 0.158],
[0.2, 0.2, 0.2, 0.175, 0.103],
[0.3, 0.2, 0.2, 0.175, 0.2],
[0.15, 0.2, 0.1, 0.175, 0.074],
[0.15, 0.2, 0.125, 0.175, 0.065],
[0.2, 0.2, 0.125, 0.175, 0.184],
[0.15, 0.2, 0.15, 0.175, 0.119],
[0.15, 0.2, 0.175, 0.175, 0.046],
[0.15, 0.2, 0.2, 0.175, 0.081],
[0.3, 0.25, 0.15, 0.175, 0.176],
[0.35, 0.25, 0.15, 0.175, 0.156],
[0.25, 0.25, 0.175, 0.175, 0.16],
[0.3, 0.25, 0.175, 0.175, 0.128],
[0.3, 0.25, 0.2, 0.175, 0.023],
[0.35, 0.25, 0.2, 0.175, 0.137],
[0.2, 0.25, 0.1, 0.175, 0.164],
[0.25, 0.25, 0.1, 0.175, 0.092],
[0.2, 0.25, 0.125, 0.175, 0.168],
[0.25, 0.25, 0.125, 0.175, 0.033],
[0.2, 0.25, 0.15, 0.175, 0.156],
[0.25, 0.25, 0.15, 0.175, 0.085],
[0.2, 0.25, 0.175, 0.175, 0.13],
[0.25, 0.25, 0.175, 0.175, 0.183],
[0.2, 0.25, 0.2, 0.175, 0.19],
[0.25, 0.25, 0.2, 0.175, 0.022],
[0.35, 0.3, 0.15, 0.175, 0.145],
[0.25, 0.3, 0.175, 0.175, 0.191],
[0.3, 0.3, 0.175, 0.175, 0.154],
[0.35, 0.3, 0.175, 0.175, 0.152],
[0.4, 0.3, 0.175, 0.175, 0.116],
[0.5, 0.3, 0.175, 0.175, 0.167],
[0.3, 0.3, 0.2, 0.175, 0.173],
[0.35, 0.3, 0.2, 0.175, 0.065],
[0.4, 0.3, 0.2, 0.175, 0.106],
[0.5, 0.3, 0.2, 0.175, 0.166],
[0.25, 0.3, 0.1, 0.175, 0.064],
[0.2, 0.3, 0.125, 0.175, 0.198],
[0.25, 0.3, 0.125, 0.175, 0.029],
[0.3, 0.3, 0.125, 0.175, 0.192],
[0.2, 0.3, 0.15, 0.175, 0.173],
[0.25, 0.3, 0.15, 0.175, 0.053],
[0.3, 0.3, 0.15, 0.175, 0.189],
[0.2, 0.3, 0.175, 0.175, 0.147],
[0.25, 0.3, 0.175, 0.175, 0.144],
[0.3, 0.3, 0.175, 0.175, 0.181],
[0.25, 0.3, 0.2, 0.175, 0.039],
[0.3, 0.3, 0.2, 0.175, 0.177],
[0.35, 0.35, 0.2, 0.175, 0.168],
[0.4, 0.35, 0.2, 0.175, 0.2],
[0.25, 0.35, 0.1, 0.175, 0.198],
[0.25, 0.35, 0.125, 0.175, 0.146],
[0.25, 0.35, 0.15, 0.175, 0.185],
[0.3, 0.35, 0.175, 0.175, 0.171],
[0.25, 0.35, 0.2, 0.175, 0.11],
[0.25, 0.4, 0.1, 0.175, 0.193],
[0.25, 0.4, 0.125, 0.175, 0.152],
[0.3, 0.4, 0.125, 0.175, 0.171],
[0.35, 0.4, 0.125, 0.175, 0.198],
[0.25, 0.4, 0.15, 0.175, 0.18],
[0.35, 0.4, 0.15, 0.175, 0.189],
[0.3, 0.4, 0.175, 0.175, 0.126],
[0.25, 0.4, 0.2, 0.175, 0.124],
[0.35, 0.5, 0.125, 0.175, 0.199],
[0.4, 0.5, 0.125, 0.175, 0.125],
[0.3, 0.5, 0.175, 0.175, 0.197],
[0.35, 0.5, 0.2, 0.175, 0.179],
[0.05, 0.05, 0.125, 0.2, 0.179],
[0.05, 0.05, 0.175, 0.2, 0.153],
[0.15, 0.1, 0.15, 0.2, 0.121],
[0.15, 0.1, 0.175, 0.2, 0.077],
[0.1, 0.1, 0.2, 0.2, 0.196],
[0.15, 0.1, 0.2, 0.2, 0.172],
[0.25, 0.15, 0.125, 0.2, 0.189],
[0.15, 0.15, 0.15, 0.2, 0.127],
[0.2, 0.15, 0.15, 0.2, 0.078],
[0.25, 0.15, 0.15, 0.2, 0.176],
[0.15, 0.15, 0.175, 0.2, 0.168],
[0.2, 0.15, 0.175, 0.2, 0.076],
[0.15, 0.15, 0.2, 0.2, 0.131],
[0.2, 0.15, 0.2, 0.2, 0.123],
[0.15, 0.15, 0.1, 0.2, 0.134],
[0.15, 0.15, 0.125, 0.2, 0.122],
[0.15, 0.15, 0.15, 0.2, 0.084],
[0.15, 0.15, 0.175, 0.2, 0.135],
[0.1, 0.15, 0.2, 0.2, 0.148],
[0.15, 0.15, 0.2, 0.2, 0.119],
[0.25, 0.2, 0.125, 0.2, 0.192],
[0.3, 0.2, 0.125, 0.2, 0.145],
[0.2, 0.2, 0.15, 0.2, 0.15],
[0.25, 0.2, 0.15, 0.2, 0.149],
[0.3, 0.2, 0.15, 0.2, 0.17],
[0.2, 0.2, 0.175, 0.2, 0.187],
[0.25, 0.2, 0.175, 0.2, 0.193],
[0.2, 0.2, 0.2, 0.2, 0.054],
[0.15, 0.2, 0.1, 0.2, 0.178],
[0.2, 0.2, 0.1, 0.2, 0.142],
[0.15, 0.2, 0.125, 0.2, 0.124],
[0.2, 0.2, 0.125, 0.2, 0.051],
[0.15, 0.2, 0.15, 0.2, 0.134],
[0.2, 0.2, 0.15, 0.2, 0.088],
[0.15, 0.2, 0.175, 0.2, 0.142],
[0.2, 0.2, 0.175, 0.2, 0.115],
[0.15, 0.2, 0.2, 0.2, 0.111],
[0.2, 0.2, 0.2, 0.2, 0.056],
[0.3, 0.25, 0.125, 0.2, 0.122],
[0.3, 0.25, 0.15, 0.2, 0.047],
[0.35, 0.25, 0.15, 0.2, 0.04],
[0.4, 0.25, 0.15, 0.2, 0.157],
[0.25, 0.25, 0.175, 0.2, 0.022],
[0.3, 0.25, 0.175, 0.2, 0.04],
[0.35, 0.25, 0.175, 0.2, 0.107],
[0.4, 0.25, 0.175, 0.2, 0.13],
[0.3, 0.25, 0.2, 0.2, 0.167],
[0.35, 0.25, 0.2, 0.2, 0.141],
[0.4, 0.25, 0.2, 0.2, 0.175],
[0.2, 0.25, 0.1, 0.2, 0.107],
[0.2, 0.25, 0.125, 0.2, 0.19],
[0.35, 0.3, 0.15, 0.2, 0.167],
[0.25, 0.3, 0.175, 0.2, 0.198],
[0.3, 0.3, 0.175, 0.2, 0.155],
[0.4, 0.3, 0.175, 0.2, 0.177],
[0.3, 0.3, 0.2, 0.2, 0.095],
[0.35, 0.3, 0.2, 0.2, 0.088],
[0.4, 0.3, 0.2, 0.2, 0.18],
[0.25, 0.3, 0.1, 0.2, 0.072],
[0.25, 0.3, 0.125, 0.2, 0.12],
[0.3, 0.3, 0.125, 0.2, 0.152],
[0.2, 0.3, 0.15, 0.2, 0.185],
[0.25, 0.3, 0.15, 0.2, 0.074],
[0.3, 0.3, 0.15, 0.2, 0.122],
[0.35, 0.3, 0.15, 0.2, 0.201],
[0.2, 0.3, 0.175, 0.2, 0.167],
[0.25, 0.3, 0.175, 0.2, 0.024],
[0.3, 0.3, 0.175, 0.2, 0.17],
[0.25, 0.3, 0.2, 0.2, 0.148],
[0.3, 0.3, 0.2, 0.2, 0.103],
[0.4, 0.35, 0.175, 0.2, 0.199],
[0.5, 0.35, 0.175, 0.2, 0.157],
[0.35, 0.35, 0.2, 0.2, 0.179],
[0.4, 0.35, 0.2, 0.2, 0.174],
[0.5, 0.35, 0.2, 0.2, 0.168],
[0.25, 0.35, 0.1, 0.2, 0.122],
[0.25, 0.35, 0.125, 0.2, 0.122],
[0.3, 0.35, 0.125, 0.2, 0.081],
[0.35, 0.35, 0.125, 0.2, 0.149],
[0.25, 0.35, 0.15, 0.2, 0.111],
[0.3, 0.35, 0.15, 0.2, 0.096],
[0.35, 0.35, 0.15, 0.2, 0.12],
[0.25, 0.35, 0.175, 0.2, 0.151],
[0.3, 0.35, 0.175, 0.2, 0.069],
[0.25, 0.35, 0.2, 0.2, 0.124],
[0.3, 0.35, 0.2, 0.2, 0.094],
[0.25, 0.4, 0.125, 0.2, 0.189],
[0.3, 0.4, 0.125, 0.2, 0.168],
[0.35, 0.4, 0.125, 0.2, 0.171],
[0.35, 0.4, 0.15, 0.2, 0.172],
[0.3, 0.4, 0.175, 0.2, 0.116],
[0.25, 0.4, 0.2, 0.2, 0.161],
[0.4, 0.5, 0.125, 0.2, 0.144],
[0.3, 0.5, 0.175, 0.2, 0.197],
[0.35, 0.5, 0.2, 0.2, 0.192]
]

list0_025=[
[0.15, 0.1, 0.125, 0.1, 0.102],
[0.15, 0.15, 0.15, 0.1, 0.082],
[0.2, 0.15, 0.175, 0.1, 0.074],
[0.25, 0.2, 0.125, 0.1, 0.065],
[0.3, 0.2, 0.125, 0.1, 0.066],
[0.25, 0.2, 0.15, 0.1, 0.043],
[0.25, 0.2, 0.2, 0.1, 0.098],
[0.3, 0.25, 0.15, 0.1, 0.102],
[0.35, 0.25, 0.15, 0.1, 0.096],
[0.25, 0.25, 0.175, 0.1, 0.085],
[0.3, 0.25, 0.175, 0.1, 0.06],
[0.3, 0.25, 0.2, 0.1, 0.076],
[0.2, 0.15, 0.15, 0.125, 0.077],
[0.2, 0.15, 0.175, 0.125, 0.061],
[0.3, 0.2, 0.125, 0.125, 0.104],
[0.2, 0.2, 0.15, 0.125, 0.112],
[0.25, 0.2, 0.15, 0.125, 0.087],
[0.2, 0.2, 0.2, 0.125, 0.052],
[0.3, 0.25, 0.15, 0.125, 0.054],
[0.35, 0.25, 0.15, 0.125, 0.06],
[0.25, 0.25, 0.175, 0.125, 0.032],
[0.3, 0.25, 0.175, 0.125, 0.03],
[0.2, 0.25, 0.1, 0.125, 0.069],
[0.35, 0.3, 0.2, 0.125, 0.076],
[0.2, 0.3, 0.1, 0.125, 0.063],
[0.2, 0.3, 0.125, 0.125, 0.109],
[0.5, 0.5, 0.15, 0.125, 0.024],
[0.15, 0.15, 0.15, 0.15, 0.11],
[0.2, 0.15, 0.175, 0.15, 0.107],
[0.15, 0.15, 0.2, 0.15, 0.089],
[0.15, 0.15, 0.1, 0.15, 0.078],
[0.2, 0.2, 0.2, 0.15, 0.081],
[0.15, 0.2, 0.125, 0.15, 0.08],
[0.15, 0.2, 0.175, 0.15, 0.078],
[0.15, 0.2, 0.2, 0.15, 0.084],
[0.3, 0.25, 0.15, 0.15, 0.099],
[0.35, 0.25, 0.15, 0.15, 0.085],
[0.25, 0.25, 0.175, 0.15, 0.08],
[0.3, 0.25, 0.175, 0.15, 0.05],
[0.3, 0.25, 0.2, 0.15, 0.078],
[0.2, 0.25, 0.1, 0.15, 0.044],
[0.2, 0.25, 0.125, 0.15, 0.094],
[0.3, 0.3, 0.2, 0.15, 0.111],
[0.35, 0.3, 0.2, 0.15, 0.09],
[0.25, 0.3, 0.1, 0.15, 0.11],
[0.25, 0.3, 0.125, 0.15, 0.056],
[0.25, 0.3, 0.15, 0.15, 0.106],
[0.25, 0.3, 0.2, 0.15, 0.049],
[0.35, 0.35, 0.2, 0.15, 0.107],
[0.25, 0.35, 0.1, 0.15, 0.105],
[0.25, 0.35, 0.125, 0.15, 0.06],
[0.25, 0.35, 0.15, 0.15, 0.093],
[0.25, 0.35, 0.2, 0.15, 0.039],
[0.5, 0.5, 0.125, 0.15, 0.024],
[0.05, 0.05, 0.175, 0.175, 0.11],
[0.2, 0.15, 0.15, 0.175, 0.079],
[0.2, 0.15, 0.175, 0.175, 0.044],
[0.1, 0.15, 0.2, 0.175, 0.075],
[0.3, 0.2, 0.15, 0.175, 0.11],
[0.2, 0.2, 0.2, 0.175, 0.103],
[0.15, 0.2, 0.1, 0.175, 0.074],
[0.15, 0.2, 0.125, 0.175, 0.065],
[0.15, 0.2, 0.175, 0.175, 0.046],
[0.15, 0.2, 0.2, 0.175, 0.081],
[0.3, 0.25, 0.2, 0.175, 0.023],
[0.25, 0.25, 0.1, 0.175, 0.092],
[0.25, 0.25, 0.125, 0.175, 0.033],
[0.25, 0.25, 0.15, 0.175, 0.085],
[0.25, 0.25, 0.2, 0.175, 0.022],
[0.35, 0.3, 0.2, 0.175, 0.065],
[0.4, 0.3, 0.2, 0.175, 0.106],
[0.25, 0.3, 0.1, 0.175, 0.064],
[0.25, 0.3, 0.125, 0.175, 0.029],
[0.25, 0.3, 0.15, 0.175, 0.053],
[0.25, 0.3, 0.2, 0.175, 0.039],
[0.25, 0.35, 0.2, 0.175, 0.11],
[0.15, 0.1, 0.175, 0.2, 0.077],
[0.2, 0.15, 0.15, 0.2, 0.078],
[0.2, 0.15, 0.175, 0.2, 0.076],
[0.15, 0.15, 0.15, 0.2, 0.084],
[0.2, 0.2, 0.2, 0.2, 0.054],
[0.2, 0.2, 0.125, 0.2, 0.051],
[0.2, 0.2, 0.15, 0.2, 0.088],
[0.15, 0.2, 0.2, 0.2, 0.111],
[0.2, 0.2, 0.2, 0.2, 0.056],
[0.3, 0.25, 0.15, 0.2, 0.047],
[0.35, 0.25, 0.15, 0.2, 0.04],
[0.25, 0.25, 0.175, 0.2, 0.022],
[0.3, 0.25, 0.175, 0.2, 0.04],
[0.35, 0.25, 0.175, 0.2, 0.107],
[0.2, 0.25, 0.1, 0.2, 0.107],
[0.3, 0.3, 0.2, 0.2, 0.095],
[0.35, 0.3, 0.2, 0.2, 0.088],
[0.25, 0.3, 0.1, 0.2, 0.072],
[0.25, 0.3, 0.15, 0.2, 0.074],
[0.25, 0.3, 0.175, 0.2, 0.024],
[0.3, 0.3, 0.2, 0.2, 0.103],
[0.3, 0.35, 0.125, 0.2, 0.081],
[0.25, 0.35, 0.15, 0.2, 0.111],
[0.3, 0.35, 0.15, 0.2, 0.096],
[0.3, 0.35, 0.175, 0.2, 0.069],
[0.3, 0.35, 0.2, 0.2, 0.094]
]

#array=np.array(list0_025)
#array=np.array(list0_10)
array=np.array(list0_25)

features = ['B1', 'B2', 'T1', 'T2','D']
df=pd.DataFrame(array,columns=features)

print(df)

#df = px.data.iris()
#print(df)
#features = ['sepal_length', 'sepal_width', 'petal_length', 'petal_width']
X = df[features]
#print(X)

pca = PCA(n_components=2)
#print(pca)
components = pca.fit_transform(X)
#print(components)
#print(pca.explained_variance_)
#print(pca.components_.T)

loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
print(loadings)

#fig = px.scatter(components, x=0, y=1, color=df['D'])
fig = px.scatter(
    components, 
    x=0, 
    y=1, 
    color=df['D']
    )

for i, feature in enumerate(features):
    fig.add_annotation(
        ax=0, ay=0,
        axref="x", ayref="y",
        x=loadings[i, 0],
        y=loadings[i, 1],
        showarrow=True,
        arrowsize=2,
        arrowhead=2,
        xanchor="right",
        yanchor="top"
    )
    fig.add_annotation(
        x=loadings[i, 0],
        y=loadings[i, 1],
        ax=0, ay=0,
        xanchor="center",
        yanchor="bottom",
        text=feature,
        yshift=5,
    )


fig.show()


###

using HDF5

path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\data_old\test\216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_evolution.h5"

c = h5open(path, "r") do file
    read(file)
end

#%%

c = h5open(path, "r") do file
    names(c)
end

#%%

using HDF5
f = h5open(path)
lvl1keys = keys(f)
for key in lvl1keys
    println(key)
    println(f[key])
end

#data = [f[joinpath(key, "rest/of/path")] for key in lvl1keys]

import .GeneralUtilities as GU

dict=GU.load_h5_dict(path)

Keys=keys(dict)

for k in keys(dict)
    println(k)
end

#%%

for i in dict 
    println(i)
end

#%%

include("structure_analysis_modules.jl")
import .GeneralUtilities as GU

#path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\data_old\test\216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_evolution.h5"
path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\multiple_p_6_N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=110.0_GradT=0.1_StepsPerT=0.01_evolution.h5"

function pretty_print(d::Dict, pre=1)
    for (k,v) in d
        if k in ["temperature_vec","total_energy_vec","move_accepted_vec","nr_monte_carlo_steps_per_temperature_vec"]
            println("nothing")
        else
            if typeof(v) <: Dict
                s = "$(repr(k)) => "
                println(join(fill(" ", pre)) * s)
                pretty_print(v, pre+1+length(s))
            else
                println(join(fill(" ", pre)) * "$(repr(k)) => $(repr(v))")
            end
        end
    end
end


dict=GU.load_h5_dict(path)

pretty_print(dict)

dict["theta_ground_state"]


# Merge all files Vincent


# there are some depreciated functions

"""
determine uncertainty on the local volume fraction for given window size
"""
function get_uncertainty_local_volume_fract_variance_random_windows_depreciated(nr_dimensions_data::Int64,
                                                    edge_length_window::Int64,
                                                    volume_fract_tot::Float64,
                                                    sampled_voxels_array::Array{Int64},
                                                    local_volume_fract_variance::Float64)
    
   

    # sometimes the uncertainty becomes negative which can only them from the nr of independent samples
    # therefore print this number for each window size
    println("window length: "*string(edge_length_window)*", nr_independent_samples: "*string(nr_independent_samples))
    println("local_volume_fract_variance: "*string(local_volume_fract_variance))

    # calculate  uncertainty on the local volume fraction
    uncertainty_local_volume_fract_variance = sqrt( 1/( window_volume^3 * nr_independent_samples
                                                        * volume_fract_tot * (1 - volume_fract_tot )^2 )
                                                    + 2 * local_volume_fract_variance 
                                                        / (nr_independent_samples - 1) )

    return uncertainty_local_volume_fract_variance

end



"""
plot window volume times local volume fraction uncertainty as a function of window edge length
to determine whether structure is hyperuniform.
If the plot tends to zero for edge length -> infintity, then the structure is hyperuniform
"""
function plot_hyperuniform_criterion_depreciated(nr_dimensions_data::Int64,
                                    window_edge_length_vec::Vector{Int64}, 
                                    local_volume_fract_variance_vec::Vector{Measurements.Measurement};
                                    title="Hyperuniformity test",
                                    save_plot = false,
                                    save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\\",
                                    save_filename="hyperuniform_test.pdf")

    # get vector of window volumes
    window_volume_vec = window_edge_length_vec .^ nr_dimensions_data

    # plot window volume times local volume fraction variance against window edge length
    hyperuniform_plot = Plots.plot(window_edge_length_vec, 
                Measurements.value.(local_volume_fract_variance_vec) .* window_volume_vec, 
                yerr = Measurements.uncertainty.(local_volume_fract_variance_vec) .* window_volume_vec,
                xlabel="window edge length "*Latex.L"L"*" / voxel edge length",
                ylabel=Fmt.format(Latex.L"\sigma_\tau^2(L) \cdot L^{:d} / ", 
                                    nr_dimensions_data)
                                *"voxel edge length"
                                *Fmt.format(Latex.L"^{:d}", nr_dimensions_data) ,
                legend = false, dpi=250, title=title)

    # if specified by the argument, save the plot
    if save_plot
        Plots.savefig(save_path*save_filename)

    else
        Plots.display(hyperuniform_plot)
    end

    return

end


"""
perform the tensor product of two 3d arrays
"""
function tensor_product_3d(array1::Array, array2::Array)

    # get the sizes of both arrays
    size1 = size(array1)
    size2 = size(array2)

    # initialize output array
    array_out = zeros(size1[1]*size2[1], size1[2]*size2[2], size1[3]*size2[3])

    # loop through all dimensions of both arrays and calculate tensor product element-wise
    for i in eachindex(array1[:,1,1])
        for j in eachindex(array1[1,:,1])
            for k in eachindex(array1[1,1,:])
                for l in eachindex(array2[:,1,1])
                    for m in eachindex(array2[1,:,1])
                        for n in eachindex(array2[1,1,:])

                            array_out[size1[1]*i-1+l, size1[2]*j-1+m, size1[3]*k-1+n] = (array1[i,j,k] 
                                                                                       * array1[l,m,n])
                        end
                    end
                end
            end
        end

    end

    return array_out

end



"""
reduce number of array elements along every dimension by a given factor
"""
function reduce_array(array_to_reduce, reduction_factor)

    # consider only those elements on a grid with the reduction factor as a step size
    reduced_array = array_to_reduce[1:reduction_factor:end, 
                                    1:reduction_factor:end, 
                                    1:reduction_factor:end]

    return reduced_array

end




function plot_hyperuniform_criterion(nr_dimensions_data::Int64,
    window_edge_length_vec::Vector{Int64}, 
    local_volume_fract_variance_vec::Vector{Measurements.Measurement};
    title="Hyperuniformity test",
    save_plot = false,
    save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\\",
    save_filename="hyperuniform_test.pdf")

# get vector of window volumes
window_volume_vec = window_edge_length_vec .^ nr_dimensions_data

# plot window volume times local volume fraction variance against window edge length
hyperuniform_plot = Plots.plot(window_edge_length_vec, 
Measurements.value.(local_volume_fract_variance_vec) .* window_volume_vec, 
ribbon = Measurements.uncertainty.(local_volume_fract_variance_vec) .* window_volume_vec,
xlabel="window edge length "*Latex.L"L"*" / voxel edge length",
ylabel=Fmt.format(Latex.L"\sigma_\tau^2(L) \cdot L^{:d} / ", 
    nr_dimensions_data)
*"voxel edge length"
*Fmt.format(Latex.L"^{:d}", nr_dimensions_data) ,
legend = false, dpi=250, title=title)

# if specified by the argument, save the plot
if save_plot
Plots.savefig(save_path*save_filename)

else
Plots.display(hyperuniform_plot)
end

return

end




"""
get a vector with odd measuring window edge lengths
"""
function get_window_edge_length_vec(mean_edge_length_data::Int64; nr_window_sizes::Int64 = 100 )

    # set maximal window edge length
    max_window_edge_length = Int( round( mean_edge_length_data/4 ) )

    # get vector with odd integer window lengths up to mean_edge_length_data/4
    window_edge_length_vec = ( 2 .* Int.( round.( collect(
                                 LinRange(1, max_window_edge_length/2, nr_window_sizes) 
                                 ) ) )  .+ 1 )

    return window_edge_length_vec

end



"""
determine spectral density as a function of the wavenumber k for isotrope 3d binary
media
"""
function get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data::Int64,
                size_data::Tuple, 
                volume_fract_tot::Float64,
                data_binary::Array{Float64};
                nr_sampling_distances::Int64 = get_nr_sampling_distances(mean_edge_length_data),
                nr_measurements_per_distance::Int64 = 10000,
                nr_wavenumbers::Int64 = 200,
                sampling_distance_vec = get_sampling_distance_vec(mean_edge_length_data; 
                                            nr_sampling_distances = nr_sampling_distances),
                autocovariance_fct_vec = get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                            size_data, 
                                            volume_fract_tot,
                                            data_binary;
                                            nr_sampling_distances = nr_sampling_distances,
                                            nr_measurements_per_distance = nr_measurements_per_distance,
                                            sampling_distance_vec = sampling_distance_vec)[2])
                                     

    # get vector of window edge lengths that will be measured
    wavenumber_vec = get_wavenumber_vec(sampling_distance_vec; nr_wavenumbers = nr_wavenumbers)

    # create vector where for each wavenumber the spectral density and its
    # uncertainty will be stored
    spectral_density_vec = Vector{Complex{Measurements.Measurement}}(undef, nr_wavenumbers)

    # for each wavenumber get spectral density and its uncertainty
    for i in eachindex(wavenumber_vec)

        # get vector of local volume fractions and the number of independent samples
        spectral_density_vec[i] = get_spectral_density_isotrope(autocovariance_fct_vec,
                                                                wavenumber_vec[i];
                                                                sampling_distance_vec=sampling_distance_vec)

        # print current calculation status
        println("wavenumber "*string(wavenumber_vec[i])*" done")

    end

    return[wavenumber_vec, spectral_density_vec]
    
end




"""
Plot absolute value, real and imaginary part of spectral density to 
determine whether structure is hyperuniform.
If the plot tends to zero for k -> 0, then the structure is hyperuniform.
The argument plot_dict_vec contains dictionaries with the following keys:
- wavenumber_vec
- spectral_density_vec
- voxel_edge_length
- label
"""
function plot_spectral_density(plot_dict_vec::Vector,
                                save_path::String;
                                title="Spectral density",
                                save_plot = false,
                                xlims = nothing)

    # set x label which is the same for all plots
    xlabel = "wavenumber " * Latex.L"k / ( 1/ \mathrm{nm} )"

    # create empty plots for real, imaginary parts and absolute value
    re_plot = Plots.plot(xlabel=xlabel,
                                    ylabel=Latex.L"\mathrm{Re}(\tilde{\chi} (k)) ",
                                    legend = true, dpi=250, title=title)

    im_plot = Plots.plot(xlabel=xlabel,
                                    ylabel=Latex.L"\mathrm{Im}(\tilde{\chi} (k)) ",
                                    legend = true, dpi=250, title=title)

    abs_plot = Plots.plot(xlabel=xlabel,
                                    ylabel=Latex.L"\mathrm{Abs}(\tilde{\chi} (k)) ",
                                    legend = true, dpi=250, title=title)

    # if desired, set xlims
    if xlims !== nothing 

        Plots.plot!(re_plot, xlims=xlims)
        Plots.plot!(im_plot, xlims=xlims)
        Plots.plot!(abs_plot, xlims=xlims)

    end

    # loop through plot dictionaries to add plots
    for plot_dict in plot_dict_vec

        # scale wavenumbers with inverse nanometers
        x_vec = plot_dict["wavenumber_vec"] ./ plot_dict["voxel_edge_length"]

        # plot real part of spectral density
        Plots.plot!(re_plot, x_vec, 
                    real.( Measurements.value.(plot_dict["spectral_density_vec"]) ), 
                    ribbon = real.( Measurements.uncertainty.(plot_dict["spectral_density_vec"]) ),
                    label = plot_dict["label"])

        # plot imaginary part of spectral density
        Plots.plot!(im_plot, x_vec, 
                    imag.( Measurements.value.(plot_dict["spectral_density_vec"]) ), 
                    ribbon = imag.( Measurements.uncertainty.(plot_dict["spectral_density_vec"]) ),
                    label = plot_dict["label"])

        # plot absolute value of spectral density
        abs_spectral_density_vec = abs.( plot_dict["spectral_density_vec"]  )

        Plots.plot!(abs_plot, x_vec, 
                    Measurements.value.(abs_spectral_density_vec), 
                    ribbon = Measurements.uncertainty.(abs_spectral_density_vec),
                    label = plot_dict["label"])


    end

    # if specified by the argument, save the plots
    if save_plot && xlims!== nothing
        Plots.savefig(re_plot, save_path*"_spectral_density_re_zoom.png")
        Plots.savefig(im_plot, save_path*"_spectral_density_im_zoom.png")
        Plots.savefig(abs_plot, save_path*"_spectral_density_abs_zoom.png")

    elseif save_plot
        Plots.savefig(re_plot, save_path*"_spectral_density_re.png")
        Plots.savefig(im_plot, save_path*"_spectral_density_im.png")
        Plots.savefig(abs_plot, save_path*"_spectral_density_abs.png")

    # otherwise display the plots
    else
        Plots.display(re_plot)
        Plots.display(im_plot)
        Plots.display(abs_plot)

    end

    return

end


"""
Get array with vectors for which the autocovariance function will be calculated
"""
function get_sampling_vec_array(size_data::Tuple)

    # determine the maximal sampling distances along the three axes
    max_sampling_distances = Int.( floor.( (size_data .- 1) ./ 2 ))

    # determine size of array where sampling vectors will be stored in
    # along one axis (z direction is chosen here) only positive directions are considered,
    # because negative ones would yield redundant information
    # the fourth dimension contains the 3 vector entries
    sampling_vec_array_size = (2*max_sampling_distances[1] + 1 , 
                                2*max_sampling_distances[2] + 1, 
                                max_sampling_distances[3] + 1,
                                3 )

    # initialize array where sampling vectors will be stored in
    sampling_vec_array = Array{Int64}(undef, sampling_vec_array_size...)

    # fill array with sampling vectors
    for i in -max_sampling_distances[1]:max_sampling_distances[1]
        for j in -max_sampling_distances[2]:max_sampling_distances[2]
            for k in 0:max_sampling_distances[3]

                sampling_vec_array[i + max_sampling_distances[1] + 1,
                                    j + max_sampling_distances[2] + 1,
                                    k + 1,
                                    :] = [i,j,k]

            end
        end
    end

    return sampling_vec_array

end




"""
Get vector of vectors of sampled wavenumbers
"""
function get_sampled_wavenumbers_vec_vec(autocovariance_fct_array_values::Array)

    # get vectors of wavenumbers along all three dimensions
    sampled_wavenumbers_vec_vec = collect( FFTW.rfftfreq.( 
                                            size(autocovariance_fct_array_values) ) )

    # convert to Float64
    sampled_wavenumbers_vec_vec_float = []

    for sampled_wavenumbers_vec in sampled_wavenumbers_vec_vec
        push!(sampled_wavenumbers_vec_vec_float, Float64.( sampled_wavenumbers_vec ))

    end

    return  sampled_wavenumbers_vec_vec_float
    
end


"""
Calculate spectral density from autocovariance fct by means of Fast Fourier
Transform
"""
function get_spectral_density(size_data::Tuple, 
                volume_fract_tot::Float64,
                data_binary::Array{Float64};
                nr_measurements_per_direction::Int64 = 1000,
                sampling_vec_array = get_sampling_vec_array(size_data),
                autocovariance_fct_array = get_autocovariance_fct_by_sampling_vec_array(size_data, 
                            volume_fract_tot,
                            data_binary;
                            sampling_vec_array = sampling_vec_array,
                            nr_measurements_per_direction = nr_measurements_per_direction)[3])
    
    # get values of autocovariance function
    autocovariance_fct_array_values = Measurements.value.( autocovariance_fct_array )

    # determine fourier transform of autocovariance function values
    spectral_density_array = FFTW.rfft(autocovariance_fct_array_values)

    # get tuple of vectors of sampled wavenumbers
    sampled_wavenumbers_vec_vec = get_sampled_wavenumbers_vec_vec(
                            autocovariance_fct_array_values)

    # get array of sampled sampled wavevectors 
    sampled_wavevectors_array = get_vector_array(sampled_wavenumbers_vec_vec)

    return[sampled_wavenumbers_vec_vec, sampled_wavevectors_array, spectral_density_array]
    
end



"""
get vector of sampled wavenumbers
"""
function get_sampled_wavenumbers_vec(direction_vec::Vector{Int64}, autocovariance_fct_vec)

    # determine geometrical length of direction vector
    sampling_distance = sqrt(sum( direction_vec .^2 ))

    # determine sampled wavenumbers vec, where the sampling rate is the inverse of the sampling distance
    sampled_wavenumbers_vec = FFTW.rfftfreq( length( autocovariance_fct_vec ), 1/sampling_distance )

    # convert to float before returning
    return Float64.( sampled_wavenumbers_vec )
end


"""
Calculate spectral density from autocovariance fct by means of Fast Fourier
Transform
"""
function get_spectral_density_along_direction(size_data::Tuple, 
                volume_fract_tot::Float64,
                data_binary::Array{Float64};
                direction_vec::Vector{Int64} = [0,0,1],
                nr_measurements_per_direction::Int64 = 1000,
                sampling_vec_array = get_sampling_vec_array(size_data),
                autocovariance_fct_array = get_autocovariance_fct_by_sampling_vec_array(size_data, 
                            volume_fract_tot,
                            data_binary;
                            sampling_vec_array = sampling_vec_array,
                            nr_measurements_per_direction = nr_measurements_per_direction)[3])
    
    # extract autocovariance function vector along given direction
    autocovariance_fct_vec = get_autocovariance_fct_along_direction_vec(direction_vec, autocovariance_fct_array)

    # determine fourier transform of autocovariance function values along direction
    spectral_density_value_vec = FFTW.rfft( Measurements.value.(autocovariance_fct_vec) )

    # determine uncertainty of spectral_density. Since Fourier transform is linear,
    # the output's uncertainty is the Fourier transform of the input's uncertainty
    spectral_density_uncertainty_vec  = FFTW.rfft( Measurements.uncertainty.(autocovariance_fct_vec) )

    # combine values and uncertainties into measurement type
    spectral_density_vec = complex.(Measurements.measurement.( 
                                        real.(spectral_density_value_vec), 
                                        real.(spectral_density_uncertainty_vec)), 
                                    Measurements.measurement.(
                                        imag.(spectral_density_value_vec), 
                                        imag.(spectral_density_uncertainty_vec)))
 
    # get vector of sampled wavenumbers
    sampled_wavenumbers_vec = get_sampled_wavenumbers_vec(direction_vec,
                                            autocovariance_fct_vec)

    return sampled_wavenumbers_vec, spectral_density_vec
    
end



"""
Plot heat map of real part, imaginary part and absolute value of spectral density
by keeping one component of the wavevector fixed.
The fixed wavevector value is given in units of (1/nm)
"""
function plot_spectral_density_heatmap(plot_dict::Dict,
    save_path::String;
    title="Spectral density",
    save_plot = false,
    clims = nothing,
    wavevector_component_to_fix::Int64 = 3,
    wavevector_value_fixed = 0)

    # discriminate between different wavevector components that are fixed
    if wavevector_component_to_fix == 1

        # set vectors of x and y axes
        wavenumber_vec_x = plot_dict["sampled_wavenumbers_vec_vec"][2]
        wavenumber_vec_y = plot_dict["sampled_wavenumbers_vec_vec"][3]

        # find index of fixed wavenumber value
        wavevector_fixed_index = argmin( abs.( 
                                        plot_dict["sampled_wavenumbers_vec_vec"][1]
                                    .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict["spectral_density_array"][wavevector_fixed_index,:,:] 
        
        # set labels and title for the plot
        xlabel = Latex.L"k_y / ( 1/ \mathrm{nm} )" 
        ylabel = Latex.L"k_z / ( 1/ \mathrm{nm} )"
        title = (title*", "
        * Latex.L"k_x = " 
        *Fmt.format(Latex.L"{:.2f}", 
                    plot_dict["sampled_wavenumbers_vec_vec"][1][wavevector_fixed_index] 
                        / plot_dict["voxel_edge_length"] )
        *" "*Latex.L"( 1/ \mathrm{nm} )")


    elseif wavevector_component_to_fix == 2
        
        # set vectors of x and y axes
        wavenumber_vec_x = plot_dict["sampled_wavenumbers_vec_vec"][1]
        wavenumber_vec_y = plot_dict["sampled_wavenumbers_vec_vec"][3]

        # find index of fixed wavenumber value
        wavevector_fixed_value, wavevector_fixed_index = findmin( abs.( 
                                    plot_dict["sampled_wavenumbers_vec_vec"][2]
                                    .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict["spectral_density_array"][:,wavevector_fixed_index,:] 

        
        # set labels and title for the plot
        xlabel = Latex.L"k_x / ( 1/ \mathrm{nm} )" 
        ylabel = Latex.L"k_z / ( 1/ \mathrm{nm} )"
        title = (title*", "
        * Latex.L"k_y = " 
        *Fmt.format(Latex.L"{:.2f}", 
                    plot_dict["sampled_wavenumbers_vec_vec"][2][wavevector_fixed_index] 
                        / plot_dict["voxel_edge_length"] )
        *" "*Latex.L"( 1/ \mathrm{nm} )")

    elseif wavevector_component_to_fix == 3
        
        # set vectors of x and y axes
        wavenumber_vec_x = plot_dict["sampled_wavenumbers_vec_vec"][1]
        wavenumber_vec_y = plot_dict["sampled_wavenumbers_vec_vec"][2]

        # find index of fixed wavenumber value
        wavevector_fixed_value, wavevector_fixed_index = findmin( abs.( 
                                    plot_dict["sampled_wavenumbers_vec_vec"][3] 
                                    .- wavevector_value_fixed ) )

        spectral_density_2d_array = plot_dict["spectral_density_array"][:,:,wavevector_fixed_index] 

        # set labels and title for the plot
        xlabel = Latex.L"k_x / ( 1/ \mathrm{nm} )" 
        ylabel = Latex.L"k_y / ( 1/ \mathrm{nm} )"
        title = (title*", "
        * Latex.L"k_z = " 
        *Fmt.format(Latex.L"{:.2f}", 
                    plot_dict["sampled_wavenumbers_vec_vec"][3][wavevector_fixed_index] 
                        / plot_dict["voxel_edge_length"] )
        *" "*Latex.L"( 1/ \mathrm{nm} )")

    else
        @error ("Wavevector component to fix must 
                be 1, 2 or 3, but is "*string(wavevector_component_to_fix))
    end

    # scale x and y axes in units of 1/nm
    x_axis = FFTW.fftshift( wavenumber_vec_x ) / plot_dict["voxel_edge_length"]
    y_axis = FFTW.fftshift(wavenumber_vec_y ) / plot_dict["voxel_edge_length"]

    # permute dimensions of spectral density array, such that they match the axes
    spectral_density_2d_permuted_array = FFTW.fftshift(permutedims(spectral_density_2d_array) )

    # create plots
    abs_plot = Plots.heatmap(x_axis,
                                y_axis,
                                abs.(spectral_density_2d_permuted_array),
                                xlabel=xlabel,
                                ylabel=ylabel,
                                colorbar_title = Latex.L"\mathrm{Abs}( \tilde{\chi} (\vec{k}) ) " ,
                                legend = true, dpi=250, title=title,
                                c = :bluesreds,
                                aspect_ratio = :equal)

    re_plot = Plots.heatmap(x_axis,
                                y_axis,
                                real.(spectral_density_2d_permuted_array),
                                xlabel=xlabel,
                                ylabel=ylabel,
                                colorbar_title = Latex.L"\mathrm{Re}( \tilde{\chi} (\vec{k}) ) " ,
                                legend = true, dpi=250, title=title,
                                c = :bluesreds,
                                aspect_ratio = :equal)

    im_plot = Plots.heatmap(x_axis,
                                y_axis,
                                imag.(spectral_density_2d_permuted_array),
                                xlabel=xlabel,
                                ylabel=ylabel,
                                colorbar_title = Latex.L"\mathrm{Im}( \tilde{\chi} (\vec{k}) ) " ,
                                legend = true, dpi=250, title=title,
                                c = :bluesreds,
                                aspect_ratio = :equal)

    # set clims if desired
    if clims !== nothing
        Plots.heatmap!(abs_plot, clims = clims)
        Plots.heatmap!(re_plot, clims = clims)
        Plots.heatmap!(im_plot, clims = clims)
    end

    # if specified by the argument, save the plot
    if  save_plot
        Plots.savefig(abs_plot, save_path*"_spectral_density_abs.png")
        Plots.savefig(re_plot, save_path*"_spectral_density_re.png")
        Plots.savefig(im_plot, save_path*"_spectral_density_im.png")

    # otherwise display the plot
    else
        Plots.display(abs_plot)

    end

    return
end



"""
get wavenumber vector such that the discrete Fourier transform can be determined
properly based on the sampled distances
"""
function get_wavenumber_vec(sampling_distance_vec::Vector;
                                nr_wavenumbers::Int64 = 200)

    # determine maximal wavenumber based on the Nyquist–Shannon sampling theorem
    # (Nyquist frequency, https://en.wikipedia.org/wiki/Nyquist%E2%80%93Shannon_sampling_theorem)
    wavenumber_max = 1 / ( 2 * Statistics.mean( sampling_distance_vec[2:end] 
                                                    .- sampling_distance_vec[1:end-1] ) )

    # get minimal wavenumber from maximal sampled distance
    wavenumber_min = 2 * pi / maximum(sampling_distance_vec)

    # get vector of float wavenumbers
    wavenumber_vec = collect(LinRange(wavenumber_min, wavenumber_max, nr_wavenumbers))

    return wavenumber_vec
    
end


"""
get vector of sampled wavenumbers
"""
function get_sampled_wavenumbers_vec(direction_vec::Vector, 
                                    autocovariance_fct_along_direction_vec::Vector)

    # determine geometrical length of direction vector
    sampling_distance = sqrt(sum( direction_vec .^2 ))

    # determine nr of sampling distances
    nr_sampling_distances = length(autocovariance_fct_along_direction_vec)

    # determine fundamental wavenumber
    fundamental_wavenumber = 2*pi/(sampling_distance*nr_sampling_distances)

    # determine sampled wavenumbers vec
    sampled_wavenumbers_vec = collect(1:floor(nr_sampling_distances/2)) .* fundamental_wavenumber

    return sampled_wavenumbers_vec
end


"""
Calculate spectral density from autocovariance fct by means of Fast Fourier
Transform
"""
function get_spectral_density_array_by_fft(size_data::Tuple, 
                volume_fract_tot::Float64,
                data_binary::Array{Float64};
                nr_measurements_per_direction::Int64 = 1000,
                sampling_distance_vec_vec = get_sampling_distance_vec_vec(size_data),
                sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
                autocovariance_fct_array = get_autocovariance_fct_by_sampling_vec_array(size_data, 
                            volume_fract_tot,
                            data_binary;
                            sampling_vec_array = sampling_vec_array,
                            nr_measurements_per_direction = nr_measurements_per_direction)[3],
                save_result = false,
                save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\sample_name",
                voxel_edge_length = 10,
                label = "some structure")

    # save autocovariance_fct array and sampling distances to dictionary
    autocovariance_fct_direction_dict = Dict("sampling_vec_array" => sampling_vec_array,
                                            "sampling_distance_vec_vec" => sampling_distance_vec_vec,
                                            "autocovariance_fct_array" => autocovariance_fct_array,
                                            "voxel_edge_length" => voxel_edge_length,
                                            "label" => label)

    # autocovariance fct was only calculated for half space with z>0 which, due to its mirror
    # symmetry, is sufficient. For the FFT it is however necessary (I think?) to get the values
    # for the entire space, which is calculated here
    (complete_sampling_distance_vec_vec, 
        complete_sampling_vec_array, 
        complete_autocovariance_fct_array) = mirror_autoautocovariance_fct_by_sampling_vec_array(
                                                autocovariance_fct_direction_dict)

    # get values of autocovariance function
    complete_autocovariance_fct_array_values = Measurements.value.( complete_autocovariance_fct_array )

    # determine fourier transform of autocovariance function values
    spectral_density_array = FFTW.rfft(complete_autocovariance_fct_array_values)

    # get tuple of vectors of sampled wavenumbers
    wavenumber_vec_vec = get_wavenumber_vec_vec(complete_autocovariance_fct_array_values)

    # get array of sampled sampled wavevectors 
    wavevector_array = get_vector_array(wavenumber_vec_vec)

    # save results if desired
    if save_result
        
        # create dict to save
        plot_dict = Dict("wavevector_array" => wavevector_array,
                            "wavenumber_vec_vec" => wavenumber_vec_vec,
                            "spectral_density_array" => spectral_density_array,
                            "voxel_edge_length" => voxel_edge_length,
                            "label" => label)

        save_dict_to_h5(copy(plot_dict);
                        save_path=save_path*"_spectral_density_array.h5")

    end

    return[wavenumber_vec_vec, wavevector_array, spectral_density_array]
    
end



"""
load voxel size corrected data from h5 file
"""
function load_binary_data(data_path_h5::String)

    # load data dictionary (this is how h5 files are structured)
    data_binary_dict = FileIO.load(data_path_h5)

    # get data from dictionary. It is expected under the "data" key
    data_binary = data_binary_dict["data_binary"]

    return data_binary

end





"""
Update bond stretching energy for a given bond
"""
function update_local_bond_stretching_energy_keating(graph_dict::Dict, bond::Tuple{Int64, Int64})

    # get bond stretching energy
    bond_stretching_energy = (3/16 * ( 
            graph_dict["spatial_network"][bond...]["distance_squared"] - 1 
                                            )^2 ) 

    # save to dict
    graph_dict["spatial_network"][bond...]["bond_stretching_energy"] = bond_stretching_energy
    
    return graph_dict

end


"""
Update bond bending energy for a given vertex
"""
function update_local_bond_bending_energy_keating!(graph_dict::Dict, vertex_label::Int64)

    # get vector of neighbor labels
    neighbor_label_vec = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex_label))

    # initialize bond bending sum
    bond_bending_sum = 0

    # loop through all bond combinations
    for j in 1:graph_dict["coordination_nr"]

        for k in j+1:graph_dict["coordination_nr"]

            bond_bending_sum += ( 3/8 * graph_dict["bond_bending_const"] 
                * ( LinearAlgebra.dot( sign(neighbor_label_vec[j] - vertex_label) .* 
                            graph_dict["spatial_network"][vertex_label, neighbor_label_vec[j]]["vector"], 
                            sign(neighbor_label_vec[k] - vertex_label) .* 
                            graph_dict["spatial_network"][vertex_label, neighbor_label_vec[k]]["vector"]
                                     ) + 1/3 )^2 )
            
        end

    end

    # calculate bond bending energy
    bond_bending_energy = 3/8 * graph_dict["bond_bending_const"] * bond_bending_sum
    
    # save to dict
    graph_dict["spatial_network"][vertex_label]["bond_bending_energy"] = bond_bending_energy

    return graph_dict

end


"""
Calculate the total energy of a spatial network
"""
function total_energy_keating(graph_dict::Dict)

    total_energy = 0

    # loop through all vertices and sum bond bending energies
    for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        total_energy += graph_dict["spatial_network"][vertex]["bond_bending_energy"]
    end

    # loop through all bonds and sum bond stretching energies
    for edge in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        total_energy += graph_dict["spatial_network"][edge...]["bond_stretching_energy"]
    end

    return total_energy
end



"""
Calculate the local Keating energy for a given vertex from 
the bond bending and stretching energies stored in the dictionary by fully
considering its bonds and not sharing their energy between the two vertices
"""
function local_energy_keating(vertex_label::Int64, graph_dict::Dict)

    # initialize local energy
    local_energy = graph_dict["spatial_network"][vertex_label]["bond_bending_energy"]

    # sum bond stretching energy contributions by considering that each bond
    # is shared by two vertices
    for neighbor in MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex_label)

        local_energy += graph_dict["spatial_network"][vertex_label, neighbor]["bond_stretching_energy"]

    end
    
    return local_energy

end



"""
Calculate the local Keating energy for a given vertex from 
the bond bending and stretching energies stored in the dictionary by
considering the bonds as shared and using half their bond energies
"""
function local_energy_keating_shared_bonds(vertex_label::Int64, graph_dict::Dict)

    # initialize local energy
    local_energy = graph_dict["spatial_network"][vertex]["bond_bending_energy"]

    # sum bond stretching energy contributions by considering that each bond
    # is shared by two vertices
    for neighbor in MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex_label)

        local_energy += 1/2 * graph_dict["spatial_network"][vertex_label, neighbor]["bond_stretching_energy"]

    end
    
    return local_energy
end

    
"""
add bond bending and stretching energies to all vertices and bonds
"""
function add_energies_to_spatial_network!(graph_dict;
    update_local_vertex_energy_fct! = update_local_bond_bending_energy_keating!,
    update_local_bond_energy_fct! = update_local_bond_stretching_energy_keating!,
    total_energy_fct = total_energy_keating)

    # add bond bending energies to vertices
    for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        graph_dict = update_local_vertex_energy_fct!(graph_dict, vertex)

    end

    # add bond stretching energies to bonds
    for edge in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        graph_dict = update_local_bond_energy_fct!(graph_dict, edge)

    end

    # add total energy to graph dict
    graph_dict["total_energy"] = total_energy_fct(graph_dict)

    return graph_dict

end




"""
calculate the total bond bending and stretching Keating
energy of a given network graph
"""
function total_bending_and_stretching_energy_keating(graph_dict::Dict)

    # first calculate bond stretching energy by looping through edges
    bond_stretching_sum = 0

    for edge in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        # get current edge length
        edge_length = graph_dict["spatial_network"][edge...]["distance"]

        # add bond stretching energy of current bond to sum
        bond_stretching_sum += (edge_length^2 - 1 )^2

    end

    # multiply sum with prefactors to determine bond stretching energy
    bond_stretching_energy = 3/16 * bond_stretching_sum

    # calculate bond-bending energy

    
    # get iterator of bond combinations
    bond_combinations_iter = Combinatorics.combinations(
                        collect(1:graph_dict["coordination_nr"]), 2)

    bond_bending_sum = 0

    for current_vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        # get list of neighbors
        neighbors_vec = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], current_vertex))

        # get vectors pointing to neighbors_vec
        bond_vectors_mat = Matrix{Float64}(undef, graph_dict["nr_dimensions"], graph_dict["coordination_nr"])

        for i in eachindex(neighbors_vec)

            # get vector pointing from the current vertex to neighbor
            # where the direction needs to be flipped if the vertex code of
            # the current vertex is inferior to the one of its neighbor
            bond_vectors_mat[:,i] = (sign(neighbors_vec[i] - current_vertex)
                        .* graph_dict["spatial_network"][current_vertex, neighbors_vec[i]]["vector"])

        end

        # loop through all bond combinations to determine corresponding bond-bending energy
        for bond_combination in bond_combinations_iter

            bond_bending_sum += ( LinearAlgebra.dot(bond_vectors_mat[:,bond_combination[1]], 
                                                bond_vectors_mat[:,bond_combination[2]]) 
                                + 1/3 )^2

        end

    end

    # multiply sum with prefactors to determine bond bending energy
    bond_bending_energy = 3/8 * graph_dict["bond_bending_const"]  * bond_bending_sum

    # get total energy
    total_energy = bond_stretching_energy + bond_bending_energy

    return [total_energy, bond_stretching_energy, bond_bending_energy]
end


"""
This function performs a bond switch on a graph.
The argument switched_bond is a tuple of two integers
which is the edge type of the MetaGraphsNext package
"""
function switch_bond_wrong!(graph_dict::Dict, switched_bond::Tuple{Int64, Int64} )

    # break the original bond
    MetaGraphsNext.rem_edge!(graph_dict["spatial_network"], switched_bond...)

    # find the other vertex's neighbors that are the closest to the current vertex
    closest_other_vertices_neighbor_vec = Vector{Int64}(undef, 2)
    vector_to_closest_other_vertices_neighbor_vec = Vector{Vector{Float64}}(undef, 2)
    distance_to_closest_other_vertices_neighbor_vec = Vector{Float64}(undef, 2)

    for i in 1:2

        # get the other vertex's neighbors
        other_vertex_neighbors_vec = MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], switched_bond[(3-i)])

        # get the vertexic position of bond vertex and store it three times
        vertexic_position_vec = graph_dict["spatial_network"][switched_bond[i]]["position"]

        # determine the closest one of the other vertex's neighbors
        closest_distance = Inf

        for neighbor in other_vertex_neighbors_vec

            # determine vector to current neighbor
            if switched_bond[i] < neighbor
                vector_to_neighbor =  get_distance_vector_pbc(
                        vertexic_position_vec,
                        graph_dict["spatial_network"][neighbor]["position"],
                        graph_dict["supercell_edge_length"] )
            else
                vector_to_neighbor =  get_distance_vector_pbc(
                        graph_dict["spatial_network"][neighbor]["position"],
                        vertexic_position_vec,
                        graph_dict["supercell_edge_length"] )
            end

            # determine length of vector to neighbor
            distance_to_neighbor = LinearAlgebra.norm(vector_to_neighbor)

            # store current neighbor if its closer than the previous neighbors
            if distance_to_neighbor < closest_distance
                closest_other_vertices_neighbor_vec[i] = neighbor
                vector_to_closest_other_vertices_neighbor_vec[i] = vector_to_neighbor
                distance_to_closest_other_vertices_neighbor_vec[i] = distance_to_neighbor
            end

        end

    end

    # create the two new bonds
    for i in 1:2

        graph_dict["spatial_network"][switched_bond[i], closest_other_vertices_neighbor_vec[i]] = Dict(
            "vector" => vector_to_closest_other_vertices_neighbor_vec[i], 
            "distance_squared" => distance_to_closest_other_vertices_neighbor_vec[i]^2 )
    end

    return graph_dict

end



"""
Get four vertices in one line around a central bond,
that is one vertex on each side of the bond
"""
function get_four_vertices_around_bond(graph_dict::Dict, 
                                switched_bond::Tuple{Int64, Int64})

    # store four vertices that sit in one line with the switched bond in the center
    vertices_in_line = zeros(Int64, 4)
    vertices_in_line[2:3] = collect(switched_bond)

    # loop through bond vertices
    for i in 1:2

        # pick a random neighbor to which the bond will be cut        
        vertices_in_line[-2+3*i] = setdiff(collect(MetaGraphsNext.neighbor_labels(
                                graph_dict["spatial_network"], switched_bond[i]
                            )) , vertices_in_line)[rand(1:graph_dict["coordination_nr"])]

    end

    return vertices_in_line
end



"""
Approximately relax a single vertex by efficiently determining
the approximated coordinate shift. The corresponding methdo is explained in
10.1142/S0217984987000065
"""
function relax_single_vertex_keating_efficiently!(graph_dict::Dict,
    vertex_to_relax::Int64;
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    update_total_energy::Bool = false)

    # get energy gradient at current vertexic position
    gradient = gradient_keating_efficient(graph_dict, vertex_to_relax)

    # get energy hessian at current vertexic position
    hessian = hessian_keating_efficient(graph_dict, vertex_to_relax)

    # calculate translation vector to approximate energy minimum
    translation_vector = .- LinearAlgebra.inv(hessian)*gradient

    # move vertex 
    graph_dict = move_vertex!(graph_dict, 
                            vertex_to_relax, 
                            translation_vector;
                            update_total_energy = update_total_energy)


    return graph_dict

end


"""
Get mesh from network
"""
function get_mesh_from_network(graph_dict::Dict; bond_radius::Real = 0.05)

    cylinders_merged = nothing

    # loop through bonds
    for bond in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        # get bond's start and target positions and its direction vector
        start_pos = graph_dict["spatial_network"][bond[1]]["position"]
        target_pos = graph_dict["spatial_network"][bond[2]]["position"]
        # direction_vec = graph_dict["spatial_network"][bond...]["vector"]

        # create cylinder surface object
        cylinder_surface = Meshes.CylinderSurface(start_pos, target_pos, bond_radius)
        

        if cylinders_merged == 0
            cylinders_merged = cylinder_surface
        else
            # merge current cylinder surface with previous ones
            cylinders_merged = Meshes.merge(cylinders_merged, cylinder_surface)
        end

    end

    # discretize object of merged cylinders
    cylinder_mesh = Meshes.discretize(cylinders_merged)

    return cylinder_mesh
end


"""
Fully relax a cluster of vertices
"""
function relax_cluster_keating!(graph_dict::Dict,
    cluster_dict::Dict; 
    nr_max_relaxation_cycles::Int64 = 25,
    break_at_relative_cluster_energy_change::Float64 = 0.001,
    reject_during_relaxation_cycle_threshold::Int64 = 5,
    initial_cluster_energy  = cluster_dict["cluster_energy"],
    relax_efficiently::Bool = true,
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    update_total_energy::Bool = false,
    track_cluster_energy::Bool = false)

    # track cluster energy during relaxatio if desired
    if track_cluster_energy
        cluster_energy_vec = [cluster_dict["cluster_energy"]]
    end

    # perform the given number of relaxation cycles
    for cycle_nr in 1:nr_max_relaxation_cycles

        # store previous cluster energy
        previous_cluster_energy = cluster_dict["cluster_energy"]

        graph_dict, cluster_dict = relax_cluster_one_cycle_keating!(graph_dict, 
        cluster_dict;
        relax_efficiently = relax_efficiently,
        relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
        relaxation_optimization_parameter_l = relaxation_optimization_parameter_l,
        update_cluster_energy = true )

        # save cluster energy if desired
        if track_cluster_energy
            push!(cluster_energy_vec, cluster_dict["cluster_energy"]) 
        end

        # break if cluster energy changes less than the given threshold
        relative_cluster_energy_change = (
            abs((previous_cluster_energy - cluster_dict["cluster_energy"])
                    /cluster_dict["cluster_energy"]))

        if (relative_cluster_energy_change < break_at_relative_cluster_energy_change 
                && cycle_nr > reject_during_relaxation_cycle_threshold)
            println("Breaking at cycle nr "*string(cycle_nr))
            break
        end

        # if cycle nr is above the given threshold, check if the relaxation can 
        # be rejected before full relaxation by estimating the final energy
        if cycle_nr > reject_during_relaxation_cycle_threshold

            # to be implemented

        end

    end

    # update total energy if desired
    if update_total_energy
        graph_dict["total_energy"] = (graph_dict["total_energy"] 
                                    + cluster_dict["cluster_energy"]
                                    - initial_cluster_energy)

        graph_dict["total_energy_up_to_date"] = true
    else
        graph_dict["total_energy_up_to_date"] = false
    end
    
    if track_cluster_energy
        return [graph_dict, cluster_energy_vec]
    else
        return graph_dict
    end
end


"""
Pick a random bond that has not been declined since the last accepted move
"""
function get_random_bond(graph_dict::Dict; declined_bonds = [], seed = Nothing)

    # set seed if desired
    if seed !== Nothing
        Random.seed!(seed)
    end

    # determine nr of bonds
    nr_bonds = graph_dict["nr_vertices"] * graph_dict["coordination_nr"] / 2 

    # check if all bonds have been attempted already
    if length(declined_bonds) == nr_bonds
            @warn "All bonds have been attempted without success"
        random_bond = []

    # if the list of declined bonds is already very long
    # pick one of the remaining ones
    elseif length(declined_bonds) > nr_bonds/2
        all_bonds_vec = collect(
                MetaGraphsNext.edge_labels(graph_dict["spatial_network"]))

        random_bond = rand(all_bonds_vec)

    # otherwise get random bond without listing all bonds
    else

        # pick a random vertex
        vertex_1 = rand(1:graph_dict["nr_vertices"])

        # pick a random neighbor
        vertex_2 = collect(MetaGraphsNext.neighbor_labels(
                            graph_dict["spatial_network"], vertex_1)
                            )[rand(1:graph_dict["coordination_nr"])]

        # create bond
        random_bond = Tuple(sort([vertex_1, vertex_2]))

        # find new bond if current one was already declined
        if random_bond in declined_bonds
            random_bond = get_random_bond(graph_dict; declined_bonds = declined_bonds)
        end
    end

    return random_bond
end


"""
This function performs a bond switch on a graph.
The argument switched_bond is a tuple of two integers
which is the edge type of the MetaGraphsNext package
"""
function switch_bond!(graph_dict::Dict,
    switched_bond::Tuple{Int64, Int64} )

    # find the other vertex's neighbors that are the closest to the current vertex
    new_bond_vertex_vec = Vector{Int64}(undef, 2)
    vector_to_new_bond_vertex_vec = Vector{Vector{Float64}}(undef, 2)
    distance_to_new_bond_vertex_vec = Vector{Float64}(undef, 2)

    # get vectors of original neighbors
    original_neighbors_vec_vec = [collect(MetaGraphsNext.neighbor_labels(
        graph_dict["spatial_network"], switched_bond[1]) ),
        collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[2]) )]

    for i in 1:2

        # get the vertex position of bond vertex
        vertex_position_vec = graph_dict["spatial_network"][switched_bond[i]]["position"]

        # get the other bond vertex's neighbors excluding 
        # the bond vertex and the bond vertex's neighbors
        considered_new_bond_vertices_vec = setdiff(original_neighbors_vec_vec[3-i], 
                                            switched_bond[i], 
                                            original_neighbors_vec_vec[i])

        # break if there are no possible new bond vertices
        if considered_new_bond_vertices_vec == []
            new_bond_vec = []
            return [graph_dict, new_bond_vec]
        
        # otherwise, pick a random new bond vertex
        else
            new_bond_vertex_vec[i] = rand(considered_new_bond_vertices_vec)
        end

        # determine vector to new bond vertex
        if switched_bond[i] < new_bond_vertex_vec[i]
            vector_to_new_bond_vertex_vec[i] = get_distance_vector_pbc(
                    vertex_position_vec,
                    graph_dict["spatial_network"][new_bond_vertex_vec[i]]["position"],
                    graph_dict["supercell_edge_length"] )
        else
            vector_to_new_bond_vertex_vec[i] = get_distance_vector_pbc(
                    graph_dict["spatial_network"][new_bond_vertex_vec[i]]["position"],
                    vertex_position_vec,
                    graph_dict["supercell_edge_length"] )
        end

        # determine length of vector to new bond vertex
        distance_to_new_bond_vertex_vec[i] = LinearAlgebra.norm(vector_to_new_bond_vertex_vec[i])

    end

    # create vector to save new bonds
    new_bond_vec = Vector{Tuple{Int64, Int64}}(undef, 2)

    # for each bond vertex, break bond to one neighbor and reconnect to
    # random neighbor of the other vertex
    for i in 1:2

        MetaGraphsNext.rem_edge!(graph_dict["spatial_network"],
            switched_bond[i], new_bond_vertex_vec[3-i])

        new_bond_vec[i] = (switched_bond[i], new_bond_vertex_vec[i])

        graph_dict["spatial_network"][new_bond_vec[i]...] = Dict(
            "vector" => vector_to_new_bond_vertex_vec[i], 
            "distance_squared" => distance_to_new_bond_vertex_vec[i]^2 )
    end

    # note, that total energy is not up to date any more
    graph_dict["total_energy_up_to_date"] = false

    return [graph_dict, new_bond_vec]

end



"""
Get effective hyperuniformity parameter which is the structure factor
at zero momentum normalized by the height of the first peak in the structure factor
as defined in equation 251 in 10.1016/j.physrep.2018.03.001
"""
function get_effective_hyperuniformity_parameter(structure_factor_dict::Dict)

    # locate first peak of structure factor
    pks, vals = Peaks.findmaxima(structure_factor_dict["structure_factor_vec"])

    # cut structure factor data at momentum just above first peak
    structure_factor_cut_vec = structure_factor_dict["structure_factor_vec"][1:pks[1]+1]
    wavenumber_cut_vec = structure_factor_dict["wavenumber_vec"][1:pks[1]+1]

    # set the order of the fitted polynomial
    polynomial_order = 3

    # fit polynomial of given order to cut data
    polynomial_fit = Polynomials.fit(wavenumber_cut_vec, 
                                    structure_factor_cut_vec,
                                    polynomial_order)

    # get extrapolated structure factor at zero momentum
    structure_factor_zero_momentum = polynomial_fit(0)

    # get the two critical momenta where the fitted structure factor is extremal
    critical_momenta = (
    ((-polynomial_fit[2]) 
        .+ [-1, +1 ] .* (sqrt(polynomial_fit[2]^2-3*polynomial_fit[3]*polynomial_fit[1])) )
    ./ (3*polynomial_fit[3]) )

    # get fitted structure factor at (first) peak
    structure_factor_first_peak = maximum( polynomial_fit.(critical_momenta) )

    # get hyperuniformity parameter
    hyperuniformity_parameter = structure_factor_zero_momentum/structure_factor_first_peak

    return [hyperuniformity_parameter, polynomial_fit]
end



"""
Get vector of mean values of q_l (rotationally invariant Steinhardt local 
bond order parameters) for the entire network and for all parameters l up 
to l_max where l is the index of the spherical harmonic Y_{lm}.
"""
function get_q_l_total_network_mean_dict(graph_dict::Dict,
    l_max::Int64)

    # initialize dictionary of q_l averaged over entire network with all values
    # set to 0
    q_l_total_network_mean_dict = Dict{Int64, Float64}()

    for l in 0:l_max
        q_l_total_network_mean_dict[l] = 0.0

    end

    # loop through vertices
    for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        # get vector of steinhardt order parameters for current vertex
        q_l_averaged_single_vertex_dict = (
            get_q_l_averaged_single_vertex_dict(
                graph_dict,
                vertex,
                l_max))

        # for each l, add current vertex' contribution to sum of all vertices
        for l in 0:l_max
            q_l_total_network_mean_dict[l] += (1/graph_dict["nr_vertices"] 
                                        * q_l_averaged_single_vertex_dict[l])
    
        end

    end

    return q_l_total_network_mean_dict
end



"""
Save spatial network to a DOT format file 
"""
function save_spatial_network_to_dot(spatial_network::MetaGraphsNext.MetaGraph,
    filename::String;
    save_path::String 
        = raw"..\structures\random_networks\\")

    # open new file
    open(save_path*filename*".gv", "w") do opened_file

        # write header
        write(opened_file, "graph T {\n")

        # loop through vertices
        for vertex in MetaGraphsNext.labels(spatial_network)

            # write vertex
            write(opened_file, Format.format("    {1} [position = [{2}, {3}, {4}]];\n",
                vertex,
                spatial_network[vertex]["position"][1],
                spatial_network[vertex]["position"][2],
                spatial_network[vertex]["position"][3]))

        end

        # loop through edges
        for edge in MetaGraphsNext.edge_labels(spatial_network)

            # write edge
            write(opened_file, Format.format("    {1} -- {2} [vector = [{3}, {4}, {5}], distance_squared = {6}];\n", 
            edge[1], edge[2], 
            spatial_network[edge...]["vector"][1],
            spatial_network[edge...]["vector"][2],
            spatial_network[edge...]["vector"][3],
            spatial_network[edge...]["distance_squared"]))

        end

        # write footer
        write(opened_file, "}\n")

    end

    return
end 


"""
Save spatial network to a DOT format file and the rest of graph_dict and
evolution_dict to an h5 file
"""
function save_graph_to_h5_and_dot(graph_dict::Dict,
    filename::String;
    evolution_dict = nothing,
    save_path::String 
        = raw"..\structures\random_networks\\")

    # save evolution dict if passed
    if evolution_dict !== nothing
        GU.save_dict_to_h5(evolution_dict;
            save_path=save_path*filename*"_evolution.h5")
    end

    # create copy of graph_dict to not change the original file
    graph_dict_to_save = deepcopy(graph_dict)

    # save graph to dot format
    save_spatial_network_to_dot(graph_dict["spatial_network"], filename; save_path=save_path)

    # remove spatial_network from graph_dict
    delete!(graph_dict_to_save, "spatial_network")

    # save graph dict
    GU.save_dict_to_h5(graph_dict_to_save; save_path=save_path*filename*".h5")

    return
end


"""
Load spatial network from a DOT format file 
"""
function load_spatial_network_from_dot(dict_path::String)

    # create an empty network graph where vertexic positions and edge vectors will be stored
    spatial_network = MetaGraphsNext.MetaGraph(Graphs.Graph(); 
                                        label_type = Int64,
                                        vertex_data_type = Dict{String, Any},
                                        edge_data_type = Dict{String, Any} )

    # read file contents, one line at a time 
    open(dict_path) do opened_file

        # read until end of file
        while ! eof(opened_file) 
        
            # read a new / next line for every iteration		 
            line_string = readline(opened_file)	

            # save edge to graph
            if occursin(" -- ", line_string)

                # get start vertex
                start_vertex_string, rest_string = split(line_string, " -- ")
                start_vertex = parse(Int64, start_vertex_string)

                # get end vertex
                end_vertex_string, rest_string = split(rest_string, " [vector = [")
                end_vertex = parse(Int64, end_vertex_string)

                # get vector and distance squared
                vector_string, rest_string = split(rest_string, "], distance_squared =")
                vector = parse.(Float64, split(vector_string, ", "))

                distance_squared = parse(Float64, rest_string[1:end-3])

                # add edge to graph
                spatial_network[start_vertex, end_vertex] = Dict("vector" => vector, "distance_squared" => distance_squared)
            
            # save vertex to graph
            elseif occursin("position", line_string)

                # get vertex and position
                vertex_string, rest_string = split(line_string, "[position = [")
                vertex = parse(Int64, vertex_string)
                position = parse.(Float64, split( rest_string[1:end-3], ", "))

                # add vertex to graph
                spatial_network[vertex] = Dict("position" => position)
            
            end
        end
    end
    
    return spatial_network
end


"""
Load graph and its properties from a DOT file and a h5 dictionary
"""
function load_graph_from_h5_and_dot(dict_path_without_format::String)

    # load spatial network in MGformat
    spatial_network = load_spatial_network_from_dot(
            dict_path_without_format*".gv")

    # load rest of graph dict
    graph_dict = GU.load_h5_dict(dict_path_without_format*".h5")

    # add spatial network key to graph dict
    graph_dict["spatial_network"] = spatial_network

    return graph_dict
end


"""
Convert a graph in MGformat to a dot format
"""
function convert_MGformat_to_dot(
    filename::String;
    save_path::String 
        = raw"..\structures\random_networks\\")

    # load graph dict
    graph_dict = load_graph_from_h5_and_MGformat(save_path*filename)

    # save graph to dot format
    save_spatial_network_to_dot(graph_dict["spatial_network"], filename; save_path=save_path)

    return
end


"""
Convert all files in a directory in MGformat to dot format
"""
function convert_all_files_in_directory_MGformat_to_dot(directory_path::String)

    # get all files in directory
    filenames = readdir(directory_path)

    # loop through files
    for filename in filenames

        if endswith(filename, ".mg")

            # convert file to dot format
            convert_MGformat_to_dot(filename[1:end-3]; save_path=directory_path)
        end

    end

    return
    
end


"""
Get extremum of a quadratic function given by three points
"""
function get_quadratic_fct_extremum(x_vec::Vector, y_vec::Vector)

    # get coefficients of quadratic function
    a = (x_vec[3] * (-y_vec[1] + y_vec[2]) 
    + x_vec[2] * (y_vec[1] - y_vec[3]) 
    +  x_vec[1] * (-y_vec[2] + y_vec[3])
    )/((x_vec[1] - x_vec[2]) * (x_vec[1] - x_vec[3]) * (x_vec[2] - x_vec[3]))

    b = (x_vec[3]^2 * (y_vec[1] - y_vec[2]) 
    + x_vec[1]^2 * (y_vec[2] - y_vec[3]) 
    + x_vec[2]^2 * (-y_vec[1] + y_vec[3])
    )/((x_vec[1] - x_vec[2]) * (x_vec[1] - x_vec[3]) * (x_vec[2] - x_vec[3]))

    c = (x_vec[1] * x_vec[3] * (-x_vec[1] + x_vec[3]) * y_vec[2] 
    + x_vec[2]^2 * (x_vec[3] * y_vec[1] - x_vec[1] * y_vec[3]) 
    + x_vec[2] * (-x_vec[3]^2 * y_vec[1] + x_vec[1]^2 * y_vec[3])
    )/((x_vec[1] - x_vec[2]) * (x_vec[1] - x_vec[3]) * (x_vec[2] - x_vec[3]))

    # get extremum of quadratic function
    x_extremum = -b/(2*a)
    y_extremum = a*x_extremum^2 + b*x_extremum + c

    return [x_extremum, y_extremum]
end


"""
Get hyperuniformity metric which is the structure factor
at zero momentum normalized by the height of the first peak in the structure factor
as defined in equation 251 in 10.1016/j.physrep.2018.03.001
"""
function get_hyperuniformity_metric(structure_factor_dict::Dict)

    # locate first peak of structure factor
    pks, vals = Peaks.findmaxima(structure_factor_dict["structure_factor_vec"])

    # cut structure factor data at momentum just above first peak
    structure_factor_cut_vec = structure_factor_dict["structure_factor_vec"][1:pks[2]+1]
    wavenumber_cut_vec = structure_factor_dict["wavenumber_vec"][1:pks[2]+1]

    # set the order of the fitted polynomial
    polynomial_order = 5

    # fit polynomial of given order to cut data
    polynomial_fit = Polynomials.fit(wavenumber_cut_vec, 
                                    structure_factor_cut_vec,
                                    polynomial_order)

    # get extrapolated structure factor at zero momentum
    structure_factor_zero_momentum = polynomial_fit(0)

    # get critical momenta which is roots of first derivative of polynomial
    polynomial_derivative = Polynomials.derivative(polynomial_fit)
    critical_momenta = Polynomials.roots(polynomial_derivative)

    # get real critical momenta
    critical_momenta_real = real.(critical_momenta[imag.(critical_momenta) .== 0])

    # get fitted structure factor at highest peak
    structure_factor_first_peak = maximum( polynomial_fit.(critical_momenta_real) )

    # get hyperuniformity metric
    hyperuniformity_metric = structure_factor_zero_momentum/structure_factor_first_peak

    return [hyperuniformity_metric, polynomial_fit]
end



"""
Calculate spectral density from autocovariance fct by means of Fast Fourier
Transform
"""
function get_spectral_density_array_by_fft(complete_autocovariance_fct_direction_dict::Dict;
            save_result = false,
            save_path = raw"..\analysis_data\sample_name",
            voxel_edge_length = nothing,
            label = nothing)

    # get values of autocovariance function
    complete_autocovariance_fct_array_values = Measurements.value.( 
                complete_autocovariance_fct_direction_dict["autocovariance_fct_array"] )

    # determine fourier transform of autocovariance function values
    spectral_density_array_fft_output = FFTW.rfft(complete_autocovariance_fct_array_values)

    # bring spectral densities into an order from negative to positive wavenumbers 
    spectral_density_array = FFTW.fftshift(spectral_density_array_fft_output, [2, 3])

    # get tuple of vectors of sampled wavenumbers
    wavenumber_vec_vec = get_wavenumber_vec_vec(complete_autocovariance_fct_array_values)

    # get array of sampled sampled wavevectors 
    wavevector_array = get_vector_array(wavenumber_vec_vec)

    # create dict to save
    spectral_density_dict = Dict("wavevector_array" => wavevector_array,
        "wavenumber_vec_vec" => wavenumber_vec_vec,
        "spectral_density_array" => spectral_density_array,
        "voxel_edge_length" => complete_autocovariance_fct_direction_dict["voxel_edge_length"],
        "label" => complete_autocovariance_fct_direction_dict["label"])

    # if desired, adjust voxel edge length and label
    spectral_density_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(spectral_density_dict), save_path*"_spectral_density_array.h5")

    end

    return spectral_density_dict
    
end



"""
determine spectral density as a function of wavevector.
Unfortunately in the moment this function takes forever (weeks?).
Instead, I think I need to use the FFT
"""
function get_spectral_density_by_wavevector_array(
                structure_dict::Dict;
                nr_measurements_per_direction::Int64 = 1000,
                sampling_distance_vec_vec = get_sampling_distance_vec_vec(structure_dict["size_data"]),
                sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(structure_dict;
                                            sampling_vec_array = sampling_vec_array,
                                            nr_measurements_per_direction = nr_measurements_per_direction),
                save_result = false,
                save_path = raw"..\analysis_data\sample_name",
                voxel_edge_length = nothing,
                label = nothing)

    # correct the sampling distance vec along the third dimension, where due to the mirror
    # symmetry of the autocovariance fct only positive z values where considered
    corrected_sampling_distance_vec_vec = [
                            autocovariance_fct_direction_dict["sampling_distance_vec_vec"][1], 
                            autocovariance_fct_direction_dict["sampling_distance_vec_vec"][2],
                            vcat(.- reverse(autocovariance_fct_direction_dict["sampling_distance_vec_vec"][3][2:end]),
                                autocovariance_fct_direction_dict["sampling_distance_vec_vec"][3]) ]
    
    # get vector of sampled wavenumbers along the three coordinate directions
    wavenumber_vec_vec = []

    for i in 1:3

        push!(wavenumber_vec_vec, get_wavenumber_vec(corrected_sampling_distance_vec_vec[i]) )

    end

    # create array of wavevectors
    wavevector_array = get_vector_array(wavenumber_vec_vec)

    # determine size of spectral density array
    spectral_density_array_size = size(wavevector_array)[1:3]

    # initialize spectral density array 
    spectral_density_array = Array{Complex{Measurements.Measurement}}(undef, spectral_density_array_size...)

    # for each wavevector, determine the spectral density
    for i in eachindex(wavenumber_vec_vec[1])
        for j in eachindex(wavenumber_vec_vec[2])
            for k in eachindex(wavenumber_vec_vec[3])

                spectral_density_array[i,j,k] = get_spectral_density(
                                                    wavevector_array[i,j,k,:], 
                                                    autocovariance_fct_direction_dict["sampling_distance_vec_vec"], 
                                                    autocovariance_fct_direction_dict["autocovariance_fct_array"])

            end
        end

        println("wavenumber "*string(wavenumber_vec_vec[1][i]*" along x done") )
    end

    # create dict to save
    spectral_density_dict = Dict("wavevector_array" => wavevector_array,
                                "wavenumber_vec_vec" => wavenumber_vec_vec,
                                "spectral_density_array" => spectral_density_array,
                                "nr_measurements_per_direction" 
                                        => autocovariance_fct_direction_dict["nr_measurements_per_direction"],
                                "voxel_edge_length" => autocovariance_fct_direction_dict["voxel_edge_length"],
                                "label" => autocovariance_fct_direction_dict["label"])

    # if desired, adjust voxel edge length and label
    spectral_density_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(spectral_density_dict),
            save_path*"_spectral_density_array.h5")

    end

    return spectral_density_dict

end



"""
Calculate spectral density from autocovariance fct along a given direction.
Somehow this function does not work yet, maybe because the calculation of exponentials for
Measurement type data simply takes too long.
"""
function get_spectral_density_along_direction_by_wavenumber_vec(structure_dict::Dict,
                direction_vec::Vector;
                nr_measurements_per_direction::Int64 = 1000,
                sampling_distance_vec_vec = get_sampling_distance_vec_vec(structure_dict["size_data"]),
                sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(structure_dict;
                                            sampling_vec_array = sampling_vec_array,
                                            nr_measurements_per_direction = nr_measurements_per_direction),
                save_result = false,
                save_path = raw"..\analysis_data\sample_name_direction",
                voxel_edge_length = nothing,
                label = nothing)
    
    # extract autocovariance function vector along given direction
    autocovariance_fct_along_direction_vec = get_autocovariance_fct_along_direction_vec(direction_vec, 
                                                                autocovariance_fct_direction_dict["autocovariance_fct_array"])

    # get sampling distances along given direction
    sampling_distance_along_direction_vec = get_sampling_distance_along_direction_vec(direction_vec, 
                                                                autocovariance_fct_along_direction_vec)

    # get vector of wavenumbers where spectral density will be calculated
    wavenumber_vec = get_wavenumber_vec(sampling_distance_along_direction_vec)

    # initialize spectral density vector
    spectral_density_vec = Vector{Measurements.Measurement}(undef, length(wavenumber_vec))

    # for each wavenumber, determine spectral density
    for i in eachindex(wavenumber_vec)
        spectral_density_vec[i] = get_spectral_density_along_direction(wavenumber_vec[i], 
                                            sampling_distance_along_direction_vec, 
                                            autocovariance_fct_along_direction_vec)

    end

    # create dict to save
    spectral_density_dict = Dict("wavenumber_vec" => wavenumber_vec,
                        "spectral_density_vec" => spectral_density_vec,
                        "nr_measurements_per_direction" => nr_measurements_per_direction,
                        "direction_vec" => direction_vec,
                        "voxel_edge_length" => autocovariance_fct_direction_dict["voxel_edge_length"],
                        "label" => autocovariance_fct_direction_dict["label"])

    # if desired, adjust voxel edge length and label
    spectral_density_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(spectral_density_dict),
            save_path*"_spectral_density_direction.h5")

    end

    return spectral_density_dict
    
end


"""
determine spectral density as a function of wavevector.
Unfortunately this function takes very long (couple of hours for a single wavevector)
"""
function get_spectral_density(wavevector::Vector{Float64}, 
                                sampling_distance_vec_vec::Vector, 
                                autocovariance_fct_array::Array)

    # initialize complex spectral density
    spectral_density = 0.0 + 0.0 * im

    # create an autocovariance function array that is mirrored in the first two dimensions
    autocovariance_fct_array_mirrored = reverse( reverse(autocovariance_fct_array, dims=1), dims=2) 

    for i in eachindex(sampling_distance_vec_vec[1])
        for j in eachindex(sampling_distance_vec_vec[2])

            # initialize complex spectral density sum along z direction
            spectral_density_z_component_sum = 0.0 + 0.0 * im

            for k in eachindex(sampling_distance_vec_vec[3])[2:end]

                spectral_density_z_component_sum += (exp(im*wavevector[3]*sampling_distance_vec_vec[3][k])
                                                            *autocovariance_fct_array_mirrored[i,j,k]
                                                    + exp(-im*wavevector[3]*sampling_distance_vec_vec[3][k])
                                                            *autocovariance_fct_array[i,j,k]   )
                                    
            end

            # add sum along z direction and other terms to sum along x and y directions
            spectral_density += (exp(-im*wavevector[1]*sampling_distance_vec_vec[1][i]) 
                                    * exp(-im*wavevector[2]*sampling_distance_vec_vec[2][j]) 
                                    * ( autocovariance_fct_array[i,j,1] + spectral_density_z_component_sum )
                                    )

            println("sampling distance "*string(sampling_distance_vec_vec[2][j])*" along y done")
        end

        println("sampling distance "*string(sampling_distance_vec_vec[1][i])*" along x done")
    end

    # multiply by the inverse of the sampling volume
    spectral_density *= 1/( (sampling_distance_vec_vec[1][2] - sampling_distance_vec_vec[1][1])
                            * (sampling_distance_vec_vec[2][2] - sampling_distance_vec_vec[2][1])
                            * (sampling_distance_vec_vec[3][2] - sampling_distance_vec_vec[3][1]) )


    return spectral_density
end



"""
get spectral density 
"""
function get_spectral_density_along_direction(wavenumber::Float64, 
                            sampling_distance_vec::Vector, 
                            autocovariance_fct_along_direction_vec::Vector)

    # determine sampling distance
    sampling_distance = sampling_distance_vec[2] - sampling_distance_vec[1]

    # calculate fourier transform for given wavenumber
    spectral_density = (1/sampling_distance 
                            * ( autocovariance_fct_along_direction_vec[1] 
                                + 2 * sum( autocovariance_fct_along_direction_vec[2:end] 
                                                .* cos.( wavenumber .* sampling_distance_vec[2:end] )  ) 
                                )
                        )

    return spectral_density

end


"""
get vector of sampling distances along the given direction
"""
function get_sampling_distance_along_direction_vec(direction_vec, 
                autocovariance_fct_along_direction_vec::Vector)
    
    # determine geometrical length of direction vector
    sampling_distance = sqrt(sum( direction_vec .^2 ))

    # get vector of sampling distances
    sampling_distance_vec = (collect(0:length(autocovariance_fct_along_direction_vec)-1) 
                                .* sampling_distance)

    return sampling_distance_vec

end



"""
function to calculate and save the following measures for a given 3d data set
- local volume fraction variance
- autocovariance function as a function of sampling distance (assuming an isotropic medium)
- spectral density as a function of sampling distance (assuming an isotropic medium)
- autocovariance function as a function of sampling vector (not assuming an isotropic medium)
- spectral density along 6 different directions (not assuming an isotropic medium)
"""
function save_statistical_measures(data_path::String,
                                    save_path::String;
                                    voxel_edge_length = nothing, 
                                    label = nothing,
                                    nr_sampling_distances = nothing,
                                    nr_measurements_per_distance = 10000,
                                    nr_window_sizes = 100,
                                    nr_measurements_per_direction = 1000,
                                    save_autocovariance_fct = true,
                                    save_spectral_density = true,
                                    save_local_volume_fraction_variance = true,
                                    save_autocovariance_fct_direction = true,
                                    save_spectral_density_along_directions = true,
                                    save_complete_autocovariance_fct_direction = true,
                                    save_spectral_density_array = true)

    # load structure dictionary which contains all essential information about the structure
    structure_dict = GU.load_h5_dict(data_path)

    # set voxel edge length and label from structure dict if not specified in the arguments
    if voxel_edge_length === nothing
        voxel_edge_length = structure_dict["voxel_edge_length"]
    end
    if label === nothing
        label = structure_dict["label"]
    end

    # if not specified, determine nr of sampling distances for autocovariance function
    # and spectral density
    if nr_sampling_distances === nothing
        nr_sampling_distances = get_nr_sampling_distances(structure_dict["mean_edge_length_data"] )
    end

    # save autocovariance function if desired
    if save_autocovariance_fct

        # get autocovariance function as a function of sampling distance
        autocovariance_fct_isotrope_dict = get_autocovariance_fct_isotrope_by_sampling_distance_vec(
                                                structure_dict;
                                                nr_sampling_distances = nr_sampling_distances,
                                                nr_measurements_per_distance = nr_measurements_per_distance,
                                                save_result = true,
                                                save_path = save_path,
                                                voxel_edge_length = voxel_edge_length,
                                                label = label)

    end


    # save spectral density if desired
    if save_spectral_density

        # check if autocovariance fct was calculated within this function
        if !save_autocovariance_fct

            # if autocovariance function was not calculated within this function, check if it was
            # calculated before and if so, load the corresponding dictionary
            if isfile(save_path*"_autocovariance_fct.h5")

                # load autocovariance function dict
                autocovariance_fct_isotrope_dict = GU.load_h5_dict(save_path*"_autocovariance_fct.h5")

            else
                autocovariance_fct_isotrope_dict = get_autocovariance_fct_isotrope_by_sampling_distance_vec(
                                                structure_dict;
                                                nr_sampling_distances = nr_sampling_distances,
                                                nr_measurements_per_distance = nr_measurements_per_distance,
                                                save_result = false,
                                                save_path = save_path,
                                                voxel_edge_length = voxel_edge_length,
                                                label = label)
            end
        end
    
        # calculate spectral_density
        spectral_density_isotrope_dict = get_spectral_density_isotrope_by_wavenumber_vec(
                structure_dict;
                nr_sampling_distances = length(autocovariance_fct_isotrope_dict["sampling_distance_vec"]),
                nr_measurements_per_distance = autocovariance_fct_isotrope_dict["nr_measurements_per_distance"],
                sampling_distance_vec = autocovariance_fct_isotrope_dict["sampling_distance_vec"],
                autocovariance_fct_dict = autocovariance_fct_isotrope_dict,
                save_result = true,
                save_path = save_path,
                voxel_edge_length = voxel_edge_length,
                label = label)
                    
    end


    # save local volume fraction if desired
    if save_local_volume_fraction_variance

        # determine local volume fraction variance vector by using measuring windows
        local_volume_fract_variance_dict = get_local_volume_fract_variance_by_window_vec(
                                                            structure_dict;
                                                            nr_window_sizes = nr_window_sizes,
                                                            window_positioning="random",
                                                            window_shape="spherical",
                                                            save_result = true,
                                                            save_path = save_path,
                                                            voxel_edge_length = voxel_edge_length,
                                                            label = label)
    
    end


    # save autocovariance function as a function of sampling direction
    if save_autocovariance_fct_direction

        # get autocovariance function array
        autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                            structure_dict;
                            nr_measurements_per_direction = nr_measurements_per_direction,
                            save_result = true,
                            save_path = save_path,
                            voxel_edge_length = voxel_edge_length,
                            label = label)
    
    end


    # save spectral density along certain directions if desired
    if save_spectral_density_along_directions

        # set vectors along which spectral density will be measured
        direction_vec_vec = [[1,0,0], 
                            [0,1,0],
                            [0,0,1],
                            (1/sqrt(2)) .* [1,-1,0],
                            (1/sqrt(3)) .* [1,1,1],
                            (1/sqrt(6)) .* [1,1,-2]]

        # set vector with names to save the files
        naming_vec = string.( [[1,0,0], 
                                [0,1,0],
                                [0,0,1],
                                [1,-1,0],
                                [1,1,1],
                                [1,1,-2]] )

        # check if autocovariance function was previously calculated in this function
        if !save_autocovariance_fct_direction

            # check if data can be loaded from file
            if isfile( save_path*"_autocovariance_fct_direction.h5" )
            
                # load autocovariance function per direction dict
                autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction.h5")
            
            else
                # get autocovariance function array
                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                    structure_dict;
                    nr_measurements_per_direction = nr_measurements_per_direction,
                    save_result = false,
                    save_path = save_path,
                    voxel_edge_length = voxel_edge_length,
                    label = label)

            end
        end

        for i in eachindex(direction_vec_vec)

            # determine spectral density
            spectral_density_along_direction_dict = get_spectral_density_along_direction_by_wavenumber_vec(
                    structure_dict,
                    direction_vec_vec[i];
                    nr_measurements_per_direction = nr_measurements_per_direction,
                    sampling_distance_vec_vec = autocovariance_fct_direction_dict["sampling_distance_vec_vec"],
                    autocovariance_fct_dict = autocovariance_fct_direction_dict,
                    save_result = true,
                    save_path = save_path*"_"*naming_vec[i],
                    voxel_edge_length = voxel_edge_length,
                    label = label*" "*naming_vec[i])


        end
    end


    # save complete autocovariance function as a function of sampling direction for all spacial directions,
    # not only the half space considered previously
    if save_complete_autocovariance_fct_direction

        # check if autocovariance function was previously calculated in this function
        if !save_autocovariance_fct_direction

            # check if data can be loaded from file
            if isfile( save_path*"_autocovariance_fct_direction.h5" )
            
                # load autocovariance function per direction dict
                autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction.h5")
            
            else
                # get autocovariance function array
                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                    structure_dict;
                    nr_measurements_per_direction = nr_measurements_per_direction,
                    save_result = false,
                    save_path = save_path,
                    voxel_edge_length = voxel_edge_length,
                    label = label)

            end
        end
    
        # calculate complete autocovariance function
        complete_autocovariance_fct_direction_dict = get_complete_autocovariance_fct_by_sampling_vec_array(
            autocovariance_fct_direction_dict;
            save_result = true,
            save_path = save_path)
                    
    end


    # save complete autocovariance function as a function of sampling direction for all spacial directions,
    # not only the half space considered previously
    if save_spectral_density_array

        # check if complete autocovariance function was previously calculated in this function
        if save_complete_autocovariance_fct_direction

            # check if data can be loaded from file
            if isfile( save_path*"_autocovariance_fct_direction_complete.h5" )
            
                # load complete autocovariance function per direction dict
                complete_autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction_complete.h5")
            
            else

                # check if autocovariance function was previously calculated in this function
                if !save_autocovariance_fct_direction
                
                    # check if data can be loaded from file
                    if isfile( save_path*"_autocovariance_fct_direction.h5" )
                    
                        # load autocovariance function per direction dict
                        autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction.h5")
                    
                    else
                        # get autocovariance function array
                        autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                            structure_dict;
                            nr_measurements_per_direction = nr_measurements_per_direction,
                            save_result = false,
                            save_path = save_path,
                            voxel_edge_length = voxel_edge_length,
                            label = label)
                    
                    end
                end

                # get complete autocovariance function array
                complete_autocovariance_fct_direction_dict = get_complete_autocovariance_fct_by_sampling_vec_array(
                    autocovariance_fct_direction_dict;
                    save_result = false,
                    save_path = save_path)

            end
        end
    
        # calculate spectral_density
        spectral_density_array_dict = get_spectral_density_array_by_fft(complete_autocovariance_fct_direction_dict;
                                                                    save_result = true,
                                                                    save_path = save_path,
                                                                    voxel_edge_length = voxel_edge_length,
                                                                    label = label)

        # calculate spectral_density
        spectral_density_array_dict = get_spectral_density_by_wavevector_array_fft(structure_dict;
            nr_measurements_per_direction = nr_measurements_per_direction,
            sampling_distance_vec_vec = get_sampling_distance_vec_vec(structure_dict["size_data"]),
            sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
            autocovariance_fct_direction_dict = autocovariance_fct_direction_dict,
            save_complete_autocovariance_fct_direction = save_complete_autocovariance_fct_direction,
            complete_autocovariance_fct_direction_dict = complete_autocovariance_fct_direction_dict,
            save_result = true,
            save_path = save_path)
                    
    end

    return
    
end


"""
This function plots several statistical measures, that is
- local volume fraction variance
- autocovariance function
- spectral density
"""
function plot_statistical_measures(data_path_vec,
                                save_path::String;
                                voxel_edge_length_vec=nothing,
                                label_vec=nothing,
                                spectral_density_xlims = nothing,
                                autocovariance_fct_direction_x_y_lims = nothing,
                                autocovariance_fct_direction_clims = nothing,
                                spectral_density_heatmaps_clims = (0, 0.1),
                                spectral_density_heatmaps_x_y_lims = nothing,
                                plot_autocovariance_fct_bool = true,
                                plot_spectral_density_bool = true,
                                plot_local_volume_fraction_variance_bool = true,
                                plot_autocovariance_fct_direction_bool = true,
                                plot_spectral_density_along_directions_bool = true,
                                plot_spectral_density_heatmaps_bool = true
                                )

    # if desired plot autocovariance function
    if plot_autocovariance_fct_bool

        # initialize vector for plot dicts
        autocovariance_fct_plot_dict_vec = []

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # load plot dictionary
            autocovariance_fct_plot_dict = GU.load_h5_dict(data_path_vec[i]*"_autocovariance_fct.h5")

            # if desired adjust label and voxel edge length
            if  label_vec !== nothing
                autocovariance_fct_plot_dict["label"] = label_vec[i]
            end

            # if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                autocovariance_fct_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            # add current plot dict to vector
            push!(autocovariance_fct_plot_dict_vec, autocovariance_fct_plot_dict)

        end

        # plot the data
        plot_autocovariance_fct(autocovariance_fct_plot_dict_vec,
                        save_path;
                        title="Autocovariance function",
                        save_plot = true)

        # plot zoom into the data
        plot_autocovariance_fct(autocovariance_fct_plot_dict_vec,
                        save_path;
                        title="Autocovariance function",
                        save_plot = true,
                        ylims = [-0.07, 0.07])

    end

    # if desired plot spectral density
    if plot_spectral_density_bool

        # initialize vector for plot dicts
        spectral_density_plot_dict_vec = []

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # load plot dictionary
            spectral_density_plot_dict = GU.load_h5_dict(data_path_vec[i]*"_spectral_density.h5")

            # if desired adjust label and voxel edge length
            if label_vec !== nothing
                spectral_density_plot_dict["label"] = label_vec[i]
            end

            # if desired adjust label and voxel edge length
            if  voxel_edge_length_vec !== nothing
                spectral_density_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            # add current plot dict to vector
            push!(spectral_density_plot_dict_vec, spectral_density_plot_dict)

        end

        # plot the data
        plot_spectral_density(spectral_density_plot_dict_vec,
                        save_path;
                        save_plot = true)

        # if specified, plot a zoom into the spectral sensity
        if spectral_density_xlims !== nothing

            plot_spectral_density(spectral_density_plot_dict_vec,
                            save_path;
                            save_plot = true,
                            xlims = spectral_density_xlims)
        end

    end


    # if desired plot local volume fraction variance
    if plot_local_volume_fraction_variance_bool

        # initialize vector for plot dicts
        local_volume_fraction_variance_plot_dict_vec = []

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # load plot dictionary
            local_volume_fraction_variance_plot_dict = GU.load_h5_dict(
                                                            data_path_vec[i]*"_volume_fraction_variance.h5")

            # if desired adjust label and voxel edge length
            if  label_vec !== nothing
                local_volume_fraction_variance_plot_dict["label"] = label_vec[i]
            end

            # if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                local_volume_fraction_variance_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            # add current plot dict to vector
            push!(local_volume_fraction_variance_plot_dict_vec, local_volume_fraction_variance_plot_dict)

        end

        # plot the data
        plot_volume_fraction_variance(local_volume_fraction_variance_plot_dict_vec,
                        save_path;
                        save_plot = true)

    end


    # if desired plot autocovariance function heatmaps
    if plot_autocovariance_fct_direction_bool

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # load plot dictionary
            autocovariance_fct_direction_plot_dict = GU.load_h5_dict(
                                                            data_path_vec[i]*"_autocovariance_fct_direction_complete.h5")

            # if desired adjust label and voxel edge length
            if label_vec  !== nothing
                autocovariance_fct_direction_plot_dict["label"] = label_vec[i]
            end

            # if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                autocovariance_fct_direction_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            # plot autocovariance function heatmaps along x-y, x-z and y-z planes
            for j in 1:3

                plot_autocovariance_fct_heatmap(autocovariance_fct_direction_plot_dict,
                        save_path*"_"*autocovariance_fct_direction_plot_dict["label"]*"_"*string(j)*"_fixed";
                        save_plot = true,
                        clims = autocovariance_fct_direction_clims,
                        x_y_lims = autocovariance_fct_direction_x_y_lims,
                        sampling_vector_component_to_fix = j,
                        sampling_vector_value_fixed = 0)
            end

        end
        
    end


    # if desired plot spectral density along three different directions
    if plot_spectral_density_along_directions_bool

        # set vector with names to load the files
        naming_vec = string.( [[1,0,0], 
                                [0,1,0],
                                [0,0,1],
                                [1,-1,0],
                                [1,1,1],
                                [1,1,-2]] )

        # set vector with the two collections of directions that will be plotted in two different plots
        range_vec = [1:3, 4:6]

        # set the names of the two plots
        plot_name_vec = ["_along_axes", "_rotated_axes"]

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # create two plots per sample, one along the axes and one along rotated coordinate axes
            for plot_per_sample_nr in 1:2

                # initialize vector for plot dicts
                spectral_density_direction_plot_dict_vec = []

                spectral_density_direction_plot_dict = Dict()

                # loop through the three directions either along the axes or along rotated coordinate axes
                for j in range_vec[plot_per_sample_nr]

                    # load plot dictionary
                    spectral_density_direction_plot_dict = GU.load_h5_dict(
                                                data_path_vec[i]*"_"*naming_vec[j]*"_spectral_density_direction.h5")

                    # if desired adjust label
                    if label_vec  !== nothing
                        spectral_density_direction_plot_dict["label"] = label_vec[i]*" "*naming_vec[j]
                    end

                    # if desired adjust voxel edge length
                    if voxel_edge_length_vec !== nothing
                        spectral_density_direction_plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
                    end

                    # add current plot dict to vector
                    push!(spectral_density_direction_plot_dict_vec, spectral_density_direction_plot_dict)

                end

                # set saving path by extracting the current sample name from data path
                save_path_specific = ( first( save_path, findlast('\\', save_path) )
                                *SubString(data_path_vec[i], (findlast('\\', data_path_vec[i]) + 1), length(data_path_vec[i]))
                                *"_direction"
                                *plot_name_vec[plot_per_sample_nr]
                                )

                # plot the data
                plot_spectral_density(spectral_density_direction_plot_dict_vec,
                                        save_path_specific;
                                        save_plot = true)

                # if specified, plot a zoom into the spectral sensity
                if spectral_density_xlims !== nothing
                
                    plot_spectral_density(spectral_density_direction_plot_dict_vec,
                    save_path_specific;
                                save_plot = true,
                                xlims = spectral_density_xlims)
                end

            end

        end

    end


    # if desired plot spectral density heatmaps
    if plot_spectral_density_heatmaps_bool

        # loop through data that will be plotted
        for i in eachindex(data_path_vec)

            # load plot dictionary
            spectral_density_array_dict = GU.load_h5_dict(data_path_vec[i]*"_spectral_density_array.h5")

            # if desired adjust label and voxel edge length
            if label_vec  !== nothing
                spectral_density_array_dict["label"] = label_vec[i]
            end

            # if desired adjust label and voxel edge length
            if voxel_edge_length_vec !== nothing
                spectral_density_array_dict["voxel_edge_length"] = voxel_edge_length_vec[i]
            end
                
            # plot spectral density heatmaps along x-y, x-z and y-z planes
            for j in 1:3

                plot_spectral_density_heatmap(spectral_density_array_dict,
                        save_path*"_"*spectral_density_array_dict["label"]*"_"*string(j)*"_fixed";
                        save_plot = true,
                        clims = spectral_density_heatmaps_clims,
                        x_y_lims = spectral_density_heatmaps_x_y_lims,
                        wavevector_component_to_fix = j,
                        wavevector_value_fixed = 0)
            end

        end
        
    end

    return

end



"""
Get the autocovariance function for 3d media with periodic
boundary conditions
"""
function get_autocovariance_fct(sampling_vec::Vector{Int64},
    structure_dict::Dict)

    # initialize vector from which the two point prob. fct. will be calculated later
    two_point_prob_fct_summand_vec = Vector{Float64}(undef, prod(structure_dict["size_data"]) )
    current_index = 1

    # loop through all voxels
    for i in 1:structure_dict["size_data"][1]
        for j in 1:structure_dict["size_data"][2]
            for k in 1:structure_dict["size_data"][3]

                # get indices of current voxel and the one at given sampling vector to it
                x1 = (i,j,k)
                x2 = (mod.(x1 .+ sampling_vec .- 1, structure_dict["size_data"]) .+ 1)

                # calculate the contribution to the two point prob. fct. from these coodinates
                two_point_prob_fct_summand_vec[current_index] = structure_dict["data_binary"][x1...] * structure_dict["data_binary"][x2...]

                current_index += 1

            end
        end
    end

    # calculate 2 point prob. function
    two_point_prob_fct = Statistics.mean( two_point_prob_fct_summand_vec )

    # determine autocovariance function
    autocovariance_fct = two_point_prob_fct - structure_dict["volume_fract_tot"]^2

    return autocovariance_fct
end



"""
Fully relax a cluster of vertices. The cluster energy will always be updated
"""
function relax_cluster_keating!(graph_dict::Dict,
    cluster_dict::Dict, 
    evolution_dict::Dict;
    threshold_cluster_energy = Inf,
    update_total_energy::Bool = false,
    print_progress::Bool = false)

    # make sure that cluster energy is up to date
    if !cluster_dict["cluster_energy_up_to_date"]
        cluster_dict["cluster_energy"] = get_cluster_energy(graph_dict, cluster_dict)
    end

    # store initial cluster energy
    unrelaxed_cluster_energy  = cluster_dict["cluster_energy"]

    # make sure that total energy is up to date if it will be updated later
    if !graph_dict["total_energy_up_to_date"] && update_total_energy
        graph_dict["total_energy"] = get_total_energy_keating(graph_dict)
    end

    # if network is supposed to be relaxed globally, store initial shell nr
    # and set threshold cycle for global relaxation
    if haskey(evolution_dict, "relax_globally_after_threshold_cycle")
        initial_shell_nr = evolution_dict["shell_nr"]
        threshold_cycle_global_relaxation = (evolution_dict["reject_during_relaxation_cycle_threshold"]*2+1)
    else
        threshold_cycle_global_relaxation = evolution_dict["nr_max_relaxation_cycles"] + 1
    end

    # perform the given number of relaxation cycles
    for cycle_nr in 1:evolution_dict["nr_max_relaxation_cycles"]

        # from the threshold cylcle on, relax globally, if desired
        if cylce_nr == threshold_cycle_global_relaxation
            # the following gives shell_nr = 8 for 216 vertices
            evolution_dict["shell_nr"] = Int(ceil(log(graph_dict["nr_vertices"])))+2
        end

        # only update cluster energy, if this is needed to get cluster energy change
        if cycle_nr <= evolution_dict["reject_during_relaxation_cycle_threshold"]-1
            update_cluster_energy = false
        else
            update_cluster_energy = true

            # store previous cluster force and energy
            previous_cluster_force = cluster_dict["cluster_force"]
            previous_cluster_energy = cluster_dict["cluster_energy"]
        end

        # relax cluster for one cycle
        graph_dict, cluster_dict = relax_cluster_one_cycle_keating!(graph_dict, 
        cluster_dict,
        evolution_dict;
        update_total_energy = false,
        update_cluster_energy = update_cluster_energy )

        # if cycle nr is above the given threshold, check if the relaxation can 
        # be breaked when it becomes clear that the total energy will exceed the threshold
        # or because of small relative energy change
        if cycle_nr > evolution_dict["reject_during_relaxation_cycle_threshold"]

            # get vector of last two cluster forces
            cluster_force_vec = [previous_cluster_force, cluster_dict["cluster_force"]]

            # get vector of last two cluster energies
            cluster_energy_vec =[previous_cluster_energy, cluster_dict["cluster_energy"]]

            # estimate relaxed cluster energy
            prefactor_force_squared, relaxed_cluster_energy = get_energy_relaxation_coefficients(
                cluster_force_vec, cluster_energy_vec)

            if print_progress
                println("Prefactor force squared: "*string(prefactor_force_squared))
                println("Relaxed cluster energy: "*string(relaxed_cluster_energy))
            end

            # break if estimated energy change exceeds the given threshold
            if (prefactor_force_squared < 1
                && relaxed_cluster_energy > 1.05*threshold_cluster_energy) 
                
                if print_progress
                    println("Relaxed energy exceeds threshold: breaking at cycle nr "*string(cycle_nr))
                end
                break
            end

            # break if cluster energy changes less than the given threshold
            relative_cluster_energy_change = (
                abs((previous_cluster_energy - cluster_dict["cluster_energy"])
                        /cluster_dict["cluster_energy"]))

            if relative_cluster_energy_change < evolution_dict["break_at_relative_cluster_energy_change"] 
    
                if print_progress
                    println("Negligeable energy change: breaking at cycle nr "*string(cycle_nr))
                end
                break
            end

        end

    end

    # restore initial shell nr in case it was altered during relaxation
    if haskey(evolution_dict, "relax_globally_after_threshold_cycle")
        evolution_dict["shell_nr"] = initial_shell_nr
    end

    # update total energy if desired
    if update_total_energy
        graph_dict["total_energy"] = (graph_dict["total_energy"] 
                                    + cluster_dict["cluster_energy"]
                                    - unrelaxed_cluster_energy)

        graph_dict["total_energy_up_to_date"] = true
    else
        graph_dict["total_energy_up_to_date"] = false
    end
    
    return [graph_dict, cluster_dict]
end


"""
Perform a Monte Carlo move without thermal fluctuations
by switching a bond, relaxing the network and then accepting the network
with Metropolis acceptance probability
"""
function monte_carlo_move!(graph_dict::Dict, 
    evolution_dict::Dict,
    temperature; 
    switched_chain::Tuple{Int64, Int64, Int64, Int64} = get_random_chain(graph_dict),
    print_progress::Bool = false)

    # save original graph dict 
    initial_graph_dict = deepcopy(graph_dict)

    # get initial cluster before bond switch
    initial_cluster_dict = get_cluster_in_shells_dict(
                                    graph_dict, 
                                    switched_chain; 
                                    shell_nr = evolution_dict["shell_nr"])

    # make sure that total energy is up to date
    if !graph_dict["total_energy_up_to_date"]
        graph_dict["total_energy"] = get_total_energy_keating(graph_dict)
        graph_dict["total_energy_up_to_date"] = true
    end

    # initialize cluster weights which will contribute to the acceptance probability
    # when thermal fluctuations are included
    cluster_relaxation_weight = 1
    cluster_excitation_weight = 1

    # set threshold cluster energy for Metropolis acceptance probability
    # if there are no thermal fluctuations considered
    if !evolution_dict["thermal_fluctuations"]
        threshold_cluster_energy = (initial_cluster_dict["cluster_energy"] 
                                    - temperature * log(rand()))
    else
        threshold_cluster_energy = Inf
    end

    # if there are thermal fluctuations, relax cluster first and calculate
    # weights of the corresponding shifts
    if evolution_dict["thermal_fluctuations"]

        # deep copy initial cluster, such that it does not get modified
        cluster_dict = deepcopy(initial_cluster_dict)

        # relax cluster
        graph_dict, cluster_dict = relax_cluster_keating!(graph_dict,
            cluster_dict,
            evolution_dict;
            update_total_energy = false,
            print_progress = print_progress)

        # get cluster weight corresponding to the relaxation translations
        cluster_relaxation_weight = get_cluster_fluctuation_weight(initial_graph_dict, 
                                                            graph_dict, 
                                                            cluster_dict,
                                                            temperature)
    end

    # switch bonds
    graph_dict = switch_chain!(graph_dict, switched_chain)

    # get cluster after bond switch
    cluster_dict = get_cluster_in_shells_dict(
                                    graph_dict, 
                                    switched_chain; 
                                    shell_nr = evolution_dict["shell_nr"])

    # relax cluster around switched chain and only update energy when there won't be
    # thermal fluctuations included afterward
    graph_dict, cluster_dict = relax_cluster_keating!(graph_dict,
        cluster_dict,
        evolution_dict;
        threshold_cluster_energy = threshold_cluster_energy,
        update_total_energy = false)


    # if desired, include thermal fluctuations by randomly shifting all cluster vertices
    if evolution_dict["thermal_fluctuations"]

        # excite cluster with thermal fluctuations and get the corresponding excitation
        # weight
        graph_dict, cluster_dict, cluster_excitation_weight = excite_cluster!(graph_dict,
                                                            cluster_dict,
                                                            temperature;
                                                            update_total_energy = false,
                                                            update_cluster_energy = true)

        # set random threshold energy increase for Metropolis acceptance probability
        threshold_cluster_energy = (initial_cluster_dict["cluster_energy"]
        - temperature 
        * (cluster_excitation_weight/cluster_relaxation_weight) * log(rand()))

    end

    # accept move if energy increase is below threshold
    move_accepted = false

    if cluster_dict["cluster_energy"] <= threshold_cluster_energy
        # update total energy
        graph_dict["total_energy"] = (graph_dict["total_energy"] 
        + cluster_dict["cluster_energy"] - initial_cluster_dict["cluster_energy"])  
        graph_dict["total_energy_up_to_date"] = true   

        move_accepted = true
    else
        graph_dict = initial_graph_dict
    end

    return [graph_dict, move_accepted]
end



"""
Define a anisotropy metric as the normalized variance of the structure factor at the
first peak of the structure factor of the diamond lattice
"""
function get_anisotropy_metric_from_structure_factor(
    structure_factor_angle_averaged_dict::Dict,
    diamond_structure_factor_peak_wavenumber::Float64 = 4.676810,
    diamond_structure_factor_peak_std::Float64 = 8.5573)

    # find the wavenumber that is the closest to the diamond peak wavenumber
    diamond_structure_factor_peak_wavenumber_index = argmin(abs.(
        structure_factor_angle_averaged_dict["wavenumber_vec"] 
        .- diamond_structure_factor_peak_wavenumber))

    # get the structure factor standard deviation around the diamond peak
    peak_structure_factor_std = Measurements.uncertainty(
        structure_factor_angle_averaged_dict["structure_factor_vec"][diamond_structure_factor_peak_wavenumber_index])

    # normalize this standard deviation by the structure factor standard deviation of the diamond peak
    anisotropy_metric_from_structure_factor = peak_structure_factor_std / diamond_structure_factor_peak_std

    return anisotropy_metric_from_structure_factor
end



"""
Define a anisotropy metric as the normalized variance of the spectral density at the
first peak of the spectral density of the diamond lattice
"""
function get_anisotropy_metric_from_spectral_density(
    spectral_density_angle_averaged_dict::Dict,
    diamond_spectral_density_peak_wavenumber::Float64 = 4.680517,
    diamond_spectral_density_peak_std::Float64 = 521.88398)

    # find the wavenumber that is the closest to the diamond peak wavenumber
    diamond_spectral_density_peak_wavenumber_index = argmin(abs.(spectral_density_angle_averaged_dict["wavenumber_vec"] .- diamond_spectral_density_peak_wavenumber))

    # get the spectral density standard deviation around the diamond peak
    peak_spectral_density_std = Measurements.uncertainty(
        spectral_density_angle_averaged_dict["spectral_density_vec"][diamond_spectral_density_peak_wavenumber_index])

    # normalize this standard deviation by the spectral density standard deviation of the diamond peak
    anisotropy_metric_from_spectral_density = peak_spectral_density_std / diamond_spectral_density_peak_std

    return anisotropy_metric_from_spectral_density
end





"""
Calculate the pore size distribution following the method described in
10.1103/PhysRevE.100.053314, modified to work with periodic boundary conditions
and by using random sampling of voxels to speed up the calculation
"""
function get_pore_size_distribution_old(structure_dict::Dict;
    nr_sampled_voxels::Int = 20000,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing,
    print_progress::Bool = false,
    thread_nr::Int64 = 0,
    print_lock = Threads.ReentrantLock())

    # create list of digital spheres with increasing radius
    sphere_pixel_radius_vec = collect(0.5001:0.5:minimum(structure_dict["size_data"]./2))
    digital_sphere_list = [get_digital_sphere(radius) for radius in sphere_pixel_radius_vec]

    # create array with pore radii
    pore_pixel_radius_array = zeros(size(structure_dict["data_binary"])...)

    # initialize counter if progress is printed
    if print_progress
        i = 1
    end

    # sample given number of voxels to speed up calculation
    sampled_coords = StatsBase.sample(CartesianIndices(structure_dict["data_binary"]), 
        nr_sampled_voxels, replace=false)

    # loop through all sampled voxels using cartesian indices
    for coord in sampled_coords

        # check if voxel is in pore
        if !structure_dict["data_binary"][coord]

            # loop through all digital spheres
            for j in eachindex(digital_sphere_list)

                # get digital sphere
                digital_sphere = digital_sphere_list[j]

                # check if voxel is in digital sphere
                if all([!structure_dict["data_binary"][
                    ((coord.I .- [1,1,1] .+ minimum(structure_dict["size_data"]) .+ digital_sphere[k]
                        ) .% structure_dict["size_data"] .+ [1,1,1])...] 
                        for k in eachindex(digital_sphere)])

                    # set pore pixel radius to maximum of current and previous radius
                    for k in eachindex(digital_sphere)

                        pore_pixel_radius_array[
                            ((coord.I .- [1,1,1] .+ minimum(structure_dict["size_data"]) .+ digital_sphere[k]
                        ).% structure_dict["size_data"] .+ [1,1,1])...
                            ] = maximum([
                                pore_pixel_radius_array[
                                    ((coord.I .- [1,1,1] .+ minimum(structure_dict["size_data"]) .+ digital_sphere[k]
                        ).% structure_dict["size_data"] .+ [1,1,1])...],
                                        sphere_pixel_radius_vec[j] 
                                ])
                    end
                
                else
                    break

                end

            end
        end

        # print progress
        if print_progress
            progress_percentage = i/nr_sampled_voxels*100

            lock(print_lock) do
                Format.printfmtln("Pore size distribution calculation progress thread nr {1:d}: {2:.1f} %", 
                    thread_nr, progress_percentage)
            end

            i += 1
        end

    end

    # shape the pore pixel radius array into a vector
    pore_pixel_radius_vec = vec(pore_pixel_radius_array)

    # filter out voxels that are not in a pore
    pore_pixel_radius_filtered_vec = pore_pixel_radius_vec[pore_pixel_radius_vec .> 0.0]

    # create histogram of pore pixel radii
    pixel_radius_histogram = StatsBase.fit(
            StatsBase.Histogram, pore_pixel_radius_filtered_vec, 
            0.2501:0.5:minimum(structure_dict["size_data"]), 
            closed=:left)

    # normalize histogram
    pixel_radius_histogram = LinearAlgebra.normalize(pixel_radius_histogram, mode=:probability)

    # get pore size distribution
    pore_size_distribution = pixel_radius_histogram.weights

    # convert pixel radii to physical radii
    pore_size_vec = sphere_pixel_radius_vec .* structure_dict["voxel_edge_length"]

    

    # create dict to save
    pore_size_distribution_dict = Dict{String, Any}("pore_size_vec" => pore_size_vec,
                                "pore_size_distribution" => pore_size_distribution,
                                "nr_sampled_voxels" => nr_sampled_voxels)

    # add label to dictionary if label is not nothing
    if label !== nothing
        pore_size_distribution_dict["label"] = label
    end

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(pore_size_distribution_dict),
            save_path*"_pore_size_distribution.h5")

    end

    return pore_size_distribution_dict
end



"""
Save spatial network to a GML format file 
"""
function save_spatial_network_to_gml(spatial_network::MetaGraphsNext.MetaGraph,
    filename::String;
    save_path::String 
        = raw"..\structures\random_networks\\")

    # open new file
    open(save_path*filename*".gml", "w") do opened_file

        # write header
        write(opened_file, "graph [ \n")

        # loop through vertices
        for vertex in MetaGraphsNext.labels(spatial_network)

            # write vertex
            write(opened_file, Format.format(
                "  node [\n    id {1}\n    position [ x {2} y {3} z {4} ]\n  ]\n",
                vertex,
                spatial_network[vertex]["position"][1],
                spatial_network[vertex]["position"][2],
                spatial_network[vertex]["position"][3]))

        end

        # loop through edges
        for edge in MetaGraphsNext.edge_labels(spatial_network)

            # write edge
            write(opened_file, Format.format(
                "  edge [\n    source {1}\n    target {2}\n    vector [ x {3} y {4} z {5} ]\n    distance_squared {6}\n  ]\n",
            edge[1],
            edge[2], 
            spatial_network[edge...]["vector"][1],
            spatial_network[edge...]["vector"][2],
            spatial_network[edge...]["vector"][3],
            spatial_network[edge...]["distance_squared"]))

        end

        # write footer
        write(opened_file, "]\n")

    end

    return
end 


"""
Get mesh from network
"""
function save_mesh_from_spatial_network(graph_dict::Dict, filename::String;
    bond_radius::Float64 = 0.3131,
    vector_out_of_supercell_length = 1/2,
    save_path::String = raw"..\structures\random_networks\\",
    duplicate_bonds_close_to_supercell_edge::Bool = true,
    format::String = "obj")

    # create graph dict to plot
    plot_dict = deepcopy(graph_dict)
    
    # cut all bonds that reach out of supercell and replace
    # them by new bonds to duplicated vertices outside of the supercell
    plot_dict = cut_bonds_out_of_supercell!(plot_dict; 
        vector_out_of_supercell_length = vector_out_of_supercell_length)

    # loop through bonds
    for bond in MetaGraphsNext.edge_labels(plot_dict["spatial_network"])

        # get bond's start and target positions and its direction vector
        start_pos = plot_dict["spatial_network"][bond[1]]["position"]
        target_pos = plot_dict["spatial_network"][bond[2]]["position"]
        # direction_vec = plot_dict["spatial_network"][bond...]["vector"]

        # create cylinder object
        current_cylinder = GeometryBasics.Cylinder(
            GeometryBasics.Point( start_pos...),
            GeometryBasics.Point( target_pos...),
            bond_radius)
        
        # mesh cylinder object
        current_cylinder_mesh = GeometryBasics.mesh(current_cylinder)

        # save mesh
        total_path = save_path*filename*"\\"*string(bond[1])*"_"*string(bond[2])*"."*format

        FileIO.save(total_path, current_cylinder_mesh)

        # if one of the two vertices is close to the supercell edge but the vertices
        # are not on opposite sides of the supercell, save another cylinder just outside the
        # supercell on the other side
        if (duplicate_bonds_close_to_supercell_edge
            && (any(start_pos .< bond_radius ) 
            || any(target_pos .< bond_radius ) 
            || any((graph_dict["supercell_edge_length"] .- start_pos) .< bond_radius )
            || any((graph_dict["supercell_edge_length"] .- target_pos) .< bond_radius ) )
            && LinearAlgebra.norm(start_pos .- target_pos) < graph_dict["supercell_edge_length"]/2
            && all(start_pos .< graph_dict["supercell_edge_length"] )
            && all(target_pos .< graph_dict["supercell_edge_length"] )
            && all(start_pos .> 0.0 )
            && all(target_pos .> 0.0 )
            )

            # check on which side of supercell the additional bond should be added and
            # calculate new start and target positions
            if any(start_pos .< bond_radius ) || any(target_pos .< bond_radius ) 
                
                new_start_pos = (start_pos .+ graph_dict["supercell_edge_length"]
                    .* ((start_pos .< bond_radius ) .|| (target_pos .< bond_radius ) ))

                new_target_pos = (target_pos .+ graph_dict["supercell_edge_length"]
                    .* ((start_pos .< bond_radius ) .|| (target_pos .< bond_radius ) ))
            else
                new_start_pos = (start_pos .- graph_dict["supercell_edge_length"]
                    .* (((graph_dict["supercell_edge_length"] .- start_pos) .< bond_radius )
                    .|| ((graph_dict["supercell_edge_length"] .- target_pos) .< bond_radius )))

                new_target_pos = (target_pos .- graph_dict["supercell_edge_length"]
                    .* (((graph_dict["supercell_edge_length"] .- start_pos) .< bond_radius )
                    .|| ((graph_dict["supercell_edge_length"] .- target_pos) .< bond_radius )))
            end

            #println(start_pos, target_pos, new_start_pos, new_target_pos)
                
            # create cylinder object
            current_cylinder = GeometryBasics.Cylinder(
                GeometryBasics.Point( new_start_pos...),
                GeometryBasics.Point( new_target_pos...),
                bond_radius)
            
            # mesh cylinder object
            current_cylinder_mesh = GeometryBasics.mesh(current_cylinder)

            # save mesh
            total_path = save_path*filename*"\\"*string(bond[1])*"_"*string(bond[2])*"_outside_supercell."*format

            FileIO.save(total_path, current_cylinder_mesh)
        end

    end

    return
end


"""
Load spatial network and its properties from a MGformat file and a h5 dictionary
"""
function load_spatial_network_from_h5_and_MGformat(dict_path_without_format::String)

    # load spatial network in MGformat
    spatial_network = MetaGraphsNext.loadgraph(
            dict_path_without_format*".mg", MetaGraphsNext.MGFormat())

    # load rest of spatial network dict
    spatial_network = GU.load_h5_dict(dict_path_without_format*".h5")

    # add spatial network key to graph dict
    spatial_network = spatial_network

    return spatial_network
end




"""
Save spatial network to a GML format file and the rest of spatial_network and
evolution_dict to an h5 file
"""
function save_graph_to_h5_and_gml(spatial_network::MetaGraphsNext.MetaGraph,
    filename::String;
    evolution_dict = nothing,
    save_path::String 
        = raw"..\structures\random_networks\\")

    # save evolution dict if passed
    if evolution_dict !== nothing
        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")
    end

    # create copy of spatial_network to not change the original file
    spatial_network_to_save = deepcopy(spatial_network)

    # save graph to gml format
    save_spatial_network_to_gml(spatial_network, filename; save_path=save_path)

    # remove spatial_network from spatial_network
    delete!(spatial_network_to_save, "spatial_network")

    # save graph dict
    GU.save_dict_to_h5(spatial_network_to_save, save_path*filename*".h5")

    return
end


"""
Load graph and its properties from a GML file and a h5 dictionary
"""
function load_graph_from_h5_and_gml(dict_path_without_format::String)

    # load spatial network in MGformat
    spatial_network = load_spatial_network_from_gml(
            dict_path_without_format*".gml")

    # load rest of graph dict
    spatial_network = GU.load_h5_dict(dict_path_without_format*".h5")

    # add spatial network key to graph dict
    spatial_network = spatial_network

    return spatial_network
end




"""
Save spatial network to an MGformat file and the rest of spatial_network and
evolution_dict to an h5 file
"""
function save_spatial_network_to_h5_and_MGformat(spatial_network::MetaGraphsNext.MetaGraph,
    filename::String;
    evolution_dict = nothing,
    save_path::String 
        = raw"..\structures\random_networks\\")

    # save evolution dict if passed
    if evolution_dict !== nothing
        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")
    end

    # create copy of spatial_network to not change the original file
    spatial_network_to_save = deepcopy(spatial_network)

    # save spatial network to MGformat
    MetaGraphsNext.savegraph(save_path*filename*".mg", spatial_network_to_save["spatial_network"])

    # remove spatial_network from spatial_network
    delete!(spatial_network_to_save, "spatial_network")

    # save spatial network dict
    GU.save_dict_to_h5(spatial_network_to_save, save_path*filename*".h5")

    return
end


"""
Convert a graph in MGformat to a gml format
"""
function convert_MGformat_to_gml(
    filename::String;
    save_path::String 
        = raw"..\structures\random_networks\\")

    # load graph dict
    spatial_network = load_graph_from_h5_and_MGformat(save_path*filename)

    # save graph to gml format
    save_spatial_network_to_gml(spatial_network, filename; save_path=save_path)

    return
end




"""
Convert all files in a directory in MGformat to gml format
"""
function convert_all_files_in_directory_MGformat_to_gml(directory_path::String)

    # get all files in directory
    filenames = readdir(directory_path)

    # loop through files
    for filename in filenames

        if endswith(filename, ".mg")

            # convert file to gml format
            convert_MGformat_to_gml(filename[1:end-3]; save_path=directory_path)
        end

    end

    return
    
end



"""
For each bond in the network, if both its vertices are on the same side of
the supercell but at least one of them lies close to the supercell edge,
duplicate the bond on the other side of the supercell just outside the edge.
This is required when cylinders are assigned to the bonds and it is plotted or
used in an optical simulation.
"""
function duplicate_bonds_close_to_supercell_edge!(
    spatial_network::MetaGraphsNext.MetaGraph;
    bond_radius::Float64 = 0.35)

    # count current vertex
    vertex_count = copy(spatial_network[]["nr_vertices"])

    # loop through bonds
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        # get bond's start and target positions and its direction vector
        start_pos = spatial_network[bond[1]]["position"]
        target_pos = spatial_network[bond[2]]["position"]

        # if one of the two vertices is close to the supercell edge but the
        # vertices are not on opposite sides of the supercell, save another
        # cylinder just outside the supercell on the other side
        if ((any(start_pos .< bond_radius ) 
            || any(target_pos .< bond_radius ) 
            || any((spatial_network[]["supercell_edge_length"] .- start_pos) 
                .< bond_radius )
            || any((spatial_network[]["supercell_edge_length"] .- target_pos) 
                .< bond_radius ) )
            && LinearAlgebra.norm(start_pos .- target_pos) 
                < spatial_network[]["supercell_edge_length"]/2
            && all(start_pos .< spatial_network[]["supercell_edge_length"] )
            && all(target_pos .< spatial_network[]["supercell_edge_length"] )
            && all(start_pos .> 0.0 )
            && all(target_pos .> 0.0 ) )

            # check on which side of supercell the additional bond should be
            # added and calculate new start and target positions
            if (any(start_pos .< bond_radius ) 
                || any(target_pos .< bond_radius ) )
                
                new_start_pos = (
                    start_pos .+ spatial_network[]["supercell_edge_length"]
                    .* ((start_pos .< bond_radius ) 
                    .|| (target_pos .< bond_radius ) ))

                new_target_pos = (target_pos 
                    .+ spatial_network[]["supercell_edge_length"]
                    .* ((start_pos .< bond_radius ) 
                    .|| (target_pos .< bond_radius ) ))
            else
                new_start_pos = (start_pos 
                    .- spatial_network[]["supercell_edge_length"]
                    .* (((spatial_network[]["supercell_edge_length"] 
                        .- start_pos) .< bond_radius )
                    .|| ((spatial_network[]["supercell_edge_length"] 
                        .- target_pos) .< bond_radius )))

                new_target_pos = (target_pos 
                    .- spatial_network[]["supercell_edge_length"]
                    .* (((spatial_network[]["supercell_edge_length"] 
                        .- start_pos) .< bond_radius )
                    .|| ((spatial_network[]["supercell_edge_length"] 
                        .- target_pos) .< bond_radius )))
            end

            # add two new vertices and the bond between them to the spatial
            # network
            spatial_network[vertex_count + 1] = (
                    Dict("position" => new_start_pos) )
            spatial_network[vertex_count + 2] = (
                        Dict("position" => new_target_pos) )

            spatial_network[vertex_count + 1, vertex_count + 2] = (
                Dict("vector" => (new_target_pos .- new_start_pos), 
                    "distance_squared" => (
                LinearAlgebra.norm(new_target_pos .- new_start_pos)^2 )) )

            vertex_count += 2
        end
    end

    spatial_network[]["nr_vertices"] = vertex_count

    return spatial_network
end