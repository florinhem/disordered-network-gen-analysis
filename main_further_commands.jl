


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
#=
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn import datasets
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import plotly.graph_objects as go
=#

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

#features = ['B1', 'B2', 'T1', 'T2','D']
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
#=
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
=#

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






"""
These are the calculations for the pachy weevil from the 10.1002advs.202202145 paper
"""

# set raw data path
data_path_raw = raw"..\structures\pachy_10.1002advs.202202145\3Dvolumes\SD_ff04_a51_box5_vox1.tif"


# load data, correct voxel size anisotropy and save data
data_binary = BDA.get_binary_data_from_colorscale(data_path_raw; 
                                        voxel_size=(1,1,1), 
                                        save_data=true, 
                                        data_path_corrected=data_path_corrected)


# set path to voxel size corrected data
data_path_corrected = raw"..\structures\pachy_10.1002advs.202202145\pachy_blue.h5"
      
# compare hyperuniformity criterion for red and blue patches
data_path_corrected_vec = [raw"..\structures\pachy_10.1002advs.202202145\pachy_blue.h5"]

# set labels for plotting 
label_vec = ["blue patch"] # , "red patch", "simple diamond"

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


# loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )

    # determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 10,
                                                                                                    window_positioning="scanned",
                                                                                                    window_shape="spherical")

    # create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


# plot window volume times local volume fraction uncertainty as a function of window edge length
# to determine whether structure is hyperuniform
BDA.plot_volume_fraction_variance_times_window_volume(nr_dimensions_data,
                            plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of scanned positions",
                            save_filename="pachy_scanned_spherical_windows_blue_red_sd")





# compare hyperuniformity criterion for red and blue patches
data_path_corrected_vec = [raw"..\structures\pachy_10.1002advs.202202145\pachy_blue.h5",
                        raw"..\structures\pachy_10.1002advs.202202145\pachy_red.h5",
                        raw"..\structures\pachy_10.1002advs.202202145\simple_diamond.h5"]

# set labels for plotting 
label_vec = ["blue patch", "red patch", "simple diamond"] # , "red patch", "simple diamond"

# in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
# they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


# loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )

    # determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 100,
                                                                                                    window_positioning="random",
                                                                                                    window_shape="spherical")

    # create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "voxel_edge_length" => voxel_edge_length_vec[i],
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="pachy_volume_fraction_variance_random_spherical_"*label_vec[i])
                                                                                            
end


# plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="pachy_volume_fraction_variance_random_spherical_windows_blue_red_sd")




# compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"..\analysis_data\pachy_volume_fraction_variance_random_spherical_blue patch.h5",
                raw"..\analysis_data\pachy_volume_fraction_variance_random_spherical_red patch.h5",
                raw"..\analysis_data\volume_fraction_variance_random_spherical_D_surface.h5"]

# set labels for plotting 
label_vec = ["blue patch", "red patch", "perfect diamond"] # , "red patch", "simple diamond"

# in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
# they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []

# loop through data that will be analyzed
for i in eachindex(data_path_vec)

    # load dictionary that contains the following keys:
    # "window_edge_length_vec"
    # "local_volume_fract_variance_vec"
    # "voxel_edge_length"
    # "label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    # adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


# plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="pachy_volume_fraction_variance_random_spherical_windows_blue_red_sd")


# perform convergence analysis to determine the nr of measurements per distance 
# when calculating the autocovariance function
BDA.convergence_analysis_autocovariance_fct_nr_measurements_per_distance(size_data,
                            volume_fract_tot,
                            data_binary;
                            save_plot = true )




# loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )


    # get autocovariance function as a function of sampling distance
    sampling_distance_vec, autocovariance_fct_vec = BDA.get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                                                                    size_data, 
                                                                                    volume_fract_tot,
                                                                                    data_binary;
                                                                                    nr_measurements_per_distance = 5000)

    # create dictionary for current plot
    plot_dict = Dict("sampling_distance_vec" => sampling_distance_vec,
                    "autocovariance_fct_vec" => autocovariance_fct_vec,
                    "voxel_edge_length" => voxel_edge_length_vec[i],
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="pachy_autocovariance_fct_"*label_vec[i])
                                                                                            
end


# plot real part, imaginary part and absolute value of spectral 
# density as a function of the wavenumber
BDA.plot_autocovariance_fct(plot_dict_vec;
                        title="Autocovariance function",
                        save_plot = true,
                        save_path=raw"..\plots\\",
                        save_filename="pachy_autocovariance_fct_blue_red_sd")


# compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"..\analysis_data\pachy_autocovariance_fct_blue patch.h5",
                raw"..\analysis_data\pachy_autocovariance_fct_red patch.h5",
                raw"..\analysis_data\autocovariance_fct_D_surface.h5"]

# set labels for plotting 
label_vec = ["blue patch", "red patch", "perfect diamond"] # , "red patch", "simple diamond"

# in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
# they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []

# loop through data that will be analyzed
for i in eachindex(data_path_vec)

    # load dictionary that contains the following keys:
    # "sampling_distance_vec"
    # "autocovariance_fct_vec"
    # "voxel_edge_length"
    # "label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    # adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


# plot local volume fraction variance as a function of window edge length
BDA.plot_autocovariance_fct(plot_dict_vec;
                            title="Autocovariance function",
                            save_plot = true,
                            save_path=raw"..\plots\\",
                            save_filename="pachy_autocovariance_fct_blue_red_sd_zoom")



# loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )


    # get spectral density as a function of the wavenumber
    wavenumber_vec, spectral_density_vec = BDA.get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data,
                                                                                        size_data,
                                                                                        volume_fract_tot,
                                                                                        data_binary;
                                                                                        nr_measurements_per_distance = 5000,
                                                                                        nr_wavenumbers=50) 

    # create dictionary for current plot
    plot_dict = Dict("wavenumber_vec" => wavenumber_vec,
                    "spectral_density_vec" => spectral_density_vec,
                    "voxel_edge_length" => voxel_edge_length_vec[i],
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)
    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="pachy_spectral_density_fct_"*label_vec[i])
                                                                                            
end


# plot real part, imaginary part and absolute value of spectral 
# density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_path=raw"..\plots\\",
                        save_filename="pachy_spectral_density_blue_red_sd")


# compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"..\analysis_data\pachy_spectral_density_fct_blue patch.h5",
                raw"..\analysis_data\pachy_spectral_density_fct_red patch.h5",
                raw"..\analysis_data\spectral_density_D_surface.h5"]

# set labels for plotting 
label_vec = ["blue patch", "red patch", "perfect diamond"] # , "red patch", "simple diamond"

# in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
# they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []

# loop through data that will be analyzed
for i in eachindex(data_path_vec)

    # load dictionary that contains the following keys:
    # "wavenumber_vec"
    # "spectral_density_vec"
    # "voxel_edge_length"
    # "label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    # adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


# plot real part, imaginary part and absolute value of spectral 
# density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_path=raw"..\plots\\",
                        save_filename="pachy_spectral_density_blue_red_sd")


# path of original data
data_path = raw"..\structures\pachy\pachy_blue.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_blue"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "blue", 
                                save_path)


# path of original data
data_path = raw"..\structures\pachy\pachy_red.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_red"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                9, 
                                "red", 
                                save_path)

                    

# set paths where statistical data is stored
data_path_vec = [raw"..\analysis_data\pachy\pachy_blue",
raw"..\analysis_data\pachy\pachy_red",
raw"..\analysis_data\nodal_surfaces\D_surface"]

# set path where plot will be stored
save_path = raw"..\plots\pachy\pachy_blue_red_d"

# plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            voxel_edge_length_vec=[10, 9, 8.5],
            label_vec=["P. c. mirabilis blue", "P. c. mirabilis red", "perfect diamond"]
            )


# path of original data
data_path = raw"..\structures\pachy\pachy_blue.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_blue"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "blue", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)


# path of original data
data_path = raw"..\structures\pachy\pachy_red.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_red"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                9, 
                                "red", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)



# set data path
data_path_vec = [raw"..\analysis_data\pachy\pachy_blue",
                raw"..\analysis_data\pachy\pachy_red"]

# set path to save plot
save_path = raw"..\plots\pachy\\"

BDA.plot_statistical_measures(data_path_vec,
            save_path;
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = true
            )



data_path = raw"..\structures\pachy\pachy_red.h5"


# load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

# set data path
dict_path = raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

# load dict
data_dict = BDA.load_h5_dict(dict_path)

# calculate spectral density
sampled_wavenumbers_vec_vec, sampled_wavevectors_array, spectral_density_array = BDA.get_spectral_density(size_data, 
                                                volume_fract_tot,
                                                data_binary;
                                                nr_measurements_per_direction = 1000,
                                                sampling_vec_array = data_dict["sampling_vec_array"],
                                                autocovariance_fct_array = data_dict["autocovariance_fct_array"])

# path where spectral density is saved
save_path = raw"..\analysis_data\pachy\pachy_red_spectral_density_direction.h5"

# create dict to save
saving_dict = Dict("sampled_wavevectors_array" => sampled_wavevectors_array,
                    "sampled_wavenumbers_vec_vec" => sampled_wavenumbers_vec_vec,
                    "spectral_density_array" => spectral_density_array,
                    "voxel_edge_length" => data_dict["voxel_edge_length"],
                    "label" => data_dict["label"])

BDA.save_dict_to_h5(saving_dict; save_path)



# path where spectral density data is saved
dict_path = raw"..\analysis_data\pachy\pachy_red_spectral_density_direction.h5"


# path where plot is saved
save_path = raw"..\plots\pachy\pachy_red_"

spectral_density_dict = BDA.load_h5_dict(dict_path)

BDA.plot_spectral_density_heatmap(spectral_density_dict,
    save_path;
    save_plot = false,
    clims = (0,200),
    wavevector_component_to_fix = 3,
    wavevector_value_fixed = 0)



# path of original data
data_path = raw"..\structures\pachy\pachy_red.h5"

# load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

# set data path
dict_path = raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

# load dict
data_dict = BDA.load_h5_dict(dict_path)

# set matrix of measured direction vectors (in this case the identity matrix)
direction_vec_mat = [1 0 0; 0 1 0; 0 0 1]

# set vector of labels
label_vec = ["pachy red [1,0,0]", "pachy red [0,1,0]", "pachy red [0,0,1]" ]

for i in 1:3

    direction_vec = direction_vec_mat[:,i]

    # determine spectral density
    sampled_wavenumbers_vec, spectral_density_vec = BDA.get_spectral_density_along_direction(size_data, 
                volume_fract_tot,
                data_binary,
                direction_vec;
                sampling_vec_array = data_dict["sampling_vec_array"],
                autocovariance_fct_array = data_dict["autocovariance_fct_array"])
    
    # path where spectral density is saved
    save_path = raw"..\analysis_data\pachy\\"* label_vec[i] *"_spectral_density_direction.h5"

    # create dict to save
    saving_dict = Dict("wavenumber_vec" => sampled_wavenumbers_vec,
                        "spectral_density_vec" => spectral_density_vec,
                        "voxel_edge_length" => data_dict["voxel_edge_length"],
                        "label" => label_vec[i])

    BDA.save_dict_to_h5(saving_dict; save_path)
end


# path of original data
data_path = raw"..\structures\pachy\pachy_red.h5"

# load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

# set data path
dict_path = raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

# load dict
data_dict = BDA.load_h5_dict(dict_path)

# set matrix of measured direction vectors (in this case the identity matrix)
direction_vec_mat = [1/sqrt(2) 1/sqrt(3) 1/sqrt(6); 
                    -1/sqrt(2) 1/sqrt(3) 1/sqrt(6); 
                    0 1/sqrt(3) -2/sqrt(6)]

# set vector of labels
label_vec = ["pachy red 1/sqrt(2)*[1,-1,0]", "pachy red 1/sqrt(3)*[1,1,1]", "pachy red 1/sqrt(6)*[1,1,-2]" ]


naming_vec = ["pachy red [1,-1,0]", "pachy red [1,1,1]", "pachy red [1,1,-2]" ]

for i in 1:3

    direction_vec = direction_vec_mat[:,i]

    # determine spectral density
    sampled_wavenumbers_vec, spectral_density_vec = BDA.get_spectral_density_along_direction(size_data, 
                volume_fract_tot,
                data_binary,
                direction_vec;
                sampling_vec_array = data_dict["sampling_vec_array"],
                autocovariance_fct_array = data_dict["autocovariance_fct_array"])
    
    # path where spectral density is saved
    save_path = raw"..\analysis_data\pachy\\"* naming_vec[i] *"_spectral_density_direction.h5"

    # create dict to save
    saving_dict = Dict("wavenumber_vec" => sampled_wavenumbers_vec,
                        "spectral_density_vec" => spectral_density_vec,
                        "voxel_edge_length" => data_dict["voxel_edge_length"],
                        "label" => label_vec[i])

    BDA.save_dict_to_h5(saving_dict; save_path)
end


naming_vec = ["pachy red [1,-1,0]", "pachy red [1,1,1]", "pachy red [1,1,-2]" ]


# initialize plot dict vec 
plot_dict_vec = []

for i in 1:3
    
    # path where spectral density is saved
    load_path = raw"..\analysis_data\pachy\\"* naming_vec[i] *"_spectral_density_direction.h5"
    
    # load dict
    data_dict = BDA.load_h5_dict(load_path)

    # add current dict to plot dict vector
    push!(plot_dict_vec, data_dict)
end

# path where plot will be saved
save_path = raw"..\plots\pachy\pachy_red_direction_rotated_axes"

# plot the spectral densities
BDA.plot_spectral_density(plot_dict_vec,
                                save_path,
                                save_plot = true,
                                xlims=[0,0.1])


                    
# path of original data
data_path = raw"..\structures\pachy\pachy_blue.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_blue"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "P. c. mirabilis blue", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)


# path of original data
data_path = raw"..\structures\pachy\pachy_red.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_red"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                9, 
                                "P. c. mirabilis red", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)


# set paths where statistical data is stored
data_path_vec = [raw"..\analysis_data\pachy\pachy_red",
raw"..\analysis_data\pachy\pachy_blue"]

# set path where plot will be stored
save_path = raw"..\plots\pachy\\"

# plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            spectral_density_xlims = [0,0.1],
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = false,
            plot_spectral_density_direction_bool = true
            )


data_path = raw"..\structures\pachy\pachy_red.h5"

# load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

# set data path
dict_path = raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

# load dict
data_dict = BDA.load_h5_dict(dict_path)


# path where spectral density is saved
save_path = raw"..\analysis_data\pachy\pachy_red_spectral_density_direction.h5"

# calculate spectral density
wavenumber_vec_vec, wavevector_array, spectral_density_array = BDA.get_spectral_density_by_wavevector_array(size_data, 
                                                volume_fract_tot,
                                                data_binary;
                                                nr_measurements_per_direction = 1000,
                                                sampling_distance_vec_vec = data_dict["sampling_distance_vec_vec"],
                                                sampling_vec_array = data_dict["sampling_vec_array"],
                                                autocovariance_fct_array = data_dict["autocovariance_fct_array"],
                                                save_result = false,
                save_path = raw"..\analysis_data\pachy\pachy_red",
                voxel_edge_length = 9,
                label = "P. c. mirabilis red")



# set paths for structure dict, autocovariance fct dict and where spectral_density along direction is saved 
structure_path = raw"..\structures\pachy\pachy_blue_structure.h5"
autocovariance_path = raw"..\analysis_data\pachy\pachy_blue_autocovariance_fct.h5"
save_path = raw"..\analysis_data\pachy\pachy_blue"
plot_path = raw"..\plots\pachy\pachy_blue"

# load structure dict
structure_dict = BDA.load_h5_dict(structure_path)

# load autocovariance fct direction dict
autocovariance_dict = BDA.load_h5_dict(autocovariance_path)

spectral_density_dict = BDA.get_spectral_density_isotrope_by_wavenumber_vec(structure_dict;
    nr_sampling_distances = length(autocovariance_dict["sampling_distance_vec"]),
    nr_measurements_per_distance = autocovariance_dict["nr_measurements_per_distance"],
    sampling_distance_vec = autocovariance_dict["sampling_distance_vec"],
    autocovariance_fct_dict = autocovariance_dict,
    sampling_distance_cutoff = 1000,
    save_result = true,
    save_path = save_path)

BDA.plot_spectral_density([spectral_density_dict],
    plot_path;
    title="Spectral density",
    save_plot = true,
    xlims = (0,0.3))



# set paths for structure dict, autocovariance fct dict and where spectral_density along direction is saved 
structure_path = raw"..\structures\pachy\pachy_red_structure.h5"
autocovariance_path = raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"
save_path = raw"..\analysis_data\pachy\pachy_red"
plot_path = raw"..\plots\pachy\pachy_red"

# load structure dict
structure_dict = BDA.load_h5_dict(structure_path)

# load autocovariance fct direction dict
autocovariance_fct_direction_dict = BDA.load_h5_dict(autocovariance_path)

# calculate and save complete autocovariance fct 
complete_autocovariance_dict = BDA.get_complete_autocovariance_fct_by_sampling_vec_array(
    autocovariance_fct_direction_dict;
    save_result = true,
    save_path = save_path)

println("complete dict calculated")

spectral_density_dict = BDA.get_spectral_density_array_by_fft(complete_autocovariance_dict;
        save_result = true,
        save_path = save_path)

println("spectral density dict calculated")

# plot spectral density
BDA.plot_spectral_density_heatmap(spectral_density_dict,
    plot_path;
    title="Spectral density",
    save_plot = true,
    clims = nothing,
    wavevector_component_to_fix = 3,
    wavevector_value_fixed = 0)


data_path_raw = raw"..\structures\pachy\3Dvolumes\Blue_SI.tif"


structure_dict = BDA.get_structure_dict_from_colorscale(data_path_raw; 
    voxel_size=(10,12,10), 
    label = "P. c. mirabilis blue",
    save_result=true, 
    save_path=raw"..\structures\pachy\pachy_blue")


"""
This is where the calculations for nodal surfaces begin
"""


# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

# set properties of generated data
unit_cell_length = 500
nr_unit_cells = 10


# loop through surfaces
for label in label_vec

    # generate 3d binary data for current nodal surface
    data_binary = BDA.get_binary_data_from_nodal_eqn(unit_cell_length, 
                                                nr_unit_cells,
                                                label)


    # save current nodal surface to h5 file
    BDA.save_nodal_surface_data(data_binary,
                            unit_cell_length, 
                            nr_unit_cells,
                            label)
                                                                                            
end



# compare hyperuniformity criterion for red and blue patches
data_path = raw"..\structures\nodal_surfaces\\"

# set surfaces that are analyzed
label_vec = ["D", "G", "P", "I-WP"] 

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


# loop through data that will be analyzed
for label in label_vec

    # determine data path of binary data for current nodal structure
    current_path = data_path*label*"_surface.h5"

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(current_path )

    # determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 100,
                                                                                                    window_positioning="random",
                                                                                                    window_shape="spherical")

    # create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "voxel_edge_length" => 10,
                    "label" => label )

    push!(plot_dict_vec, plot_dict)

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="volume_fraction_variance_random_spherical_"*label*"_surface")
                                                                                            
end


# plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="volume_fraction_variance_random_spherical_nodal_surfaces")



# loop through data that will be analyzed
for label in label_vec

    # determine data path of binary data for current nodal structure
    current_path = data_path*label*"_surface.h5"

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(current_path )

    # get autocovariance function as a function of sampling distance
    sampling_distance_vec, autocovariance_fct_vec = BDA.get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                                                                    size_data, 
                                                                                    volume_fract_tot,
                                                                                    data_binary;
                                                                                    nr_measurements_per_distance = 20000)

    # create dictionary for current plot
    plot_dict = Dict("sampling_distance_vec" => sampling_distance_vec,
                    "autocovariance_fct_vec" => autocovariance_fct_vec,
                    "voxel_edge_length" => 10,
                    "label" => label )

    push!(plot_dict_vec, plot_dict)

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="autocovariance_fct_"*label*"_surface")
                                                                                            
end


# plot local volume fraction variance as a function of window edge length
BDA.plot_autocovariance_fct(plot_dict_vec;
                            title="Autocovariance function",
                            save_plot = true,
                            save_filename="autocovariance_fct_nodal_surfaces")


                            
# create vector for plot_dicts
plot_dict_vec = []

# loop through data that will be analyzed
for label in label_vec

    # determine data path of binary data for current nodal structure
    current_path = data_path*label*"_surface.h5"

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(current_path )

    # load dictionary that contains the following keys:
    # "sampling_distance_vec"
    # "autocovariance_fct_vec"
    # "voxel_edge_length"
    # "label"
    dict = BDA.load_h5_dict("autocovariance_fct_"*label*"_surface")

    # calculate spectral density for loaded autocovariance function
    wavenumber_vec, spectral_density_vec = BDA.get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data,
                size_data, 
                volume_fract_tot,
                data_binary;
                nr_measurements_per_distance = 20000,
                nr_wavenumbers=200,
                sampling_distance_vec = dict["sampling_distance_vec"],
                autocovariance_fct_vec = dict["autocovariance_fct_vec"])

    # create dictionary for current plot
    plot_dict = Dict("wavenumber_vec" => wavenumber_vec,
                    "spectral_density_vec" => spectral_density_vec,
                    "voxel_edge_length" => 10,
                    "label" => label )

    push!(plot_dict_vec, plot_dict)

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
                        save_filename="spectral_density_"*label*"_surface")
                                                                                            
end

# plot real part, imaginary part and absolute value of spectral 
# density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_filename="spectral_density_nodal_surfaces")


# set paths where statistical data is stored
data_path_vec = (raw"..\analysis_data\nodal_surfaces\\"
                    .* ["D", "I-WP", "P", "G"] .* "_surface" )

# set path where plot will be stored
save_path = raw"..\plots\nodal_surfaces\nodal_surfaces"

# plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path)



# get data essentials of stervi data
data_path = raw"..\structures\stervi\stervi_green.h5"
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )
nr_sampling_distances = BDA.get_nr_sampling_distances(mean_edge_length_data)


# now analyze I-WP with the nr of sampling distances of the stervi weevil
# path of original data
data_path = raw"..\structures\nodal_surfaces\I-WP_surface.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\nodal_surfaces\I-WP_surface_fewer_sampling_distances"
    

# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "I-WP", 
                                save_path;
                                nr_sampling_distances = nr_sampling_distances,
                                save_local_volume_fraction_variance=false)



label_vec = ["D", "I-WP", "P", "G"]

for label in label_vec

    # path of original data
    data_path = raw"..\structures\nodal_surfaces\\"*label*"_surface.h5"

    # path where analysis data will be saved
    save_path = raw"..\analysis_data\nodal_surfaces\\"*label*"_surface"

    # calculate and save all statistical measures
    BDA.save_statistical_measures(data_path, 
                                    10, 
                                    label, 
                                    save_path;
                                    nr_sampling_distances = 500,
                                    save_local_volume_fraction_variance=false)

end


# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

# set properties of generated data
unit_cell_length = 500
nr_unit_cells = 1


# loop through surfaces
for label in label_vec

    # generate 3d binary data for current nodal surface
    data_binary = BDA.get_binary_data_from_nodal_eqn(unit_cell_length, 
                                                nr_unit_cells,
                                                label)


    # save current nodal surface to h5 file
    BDA.save_nodal_surface_data(data_binary,
                            unit_cell_length, 
                            nr_unit_cells,
                            "single_unit_cell_"*label)
                                                                                            
end


# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

# set data path
data_path_vec = raw"..\structures\nodal_surfaces\single_unit_cell_" .* label_vec .* "_surface.h5" 


# set path where autocovariance dict will be stored
save_path_vec = raw"..\analysis_data\nodal_surfaces\single_unit_cell_" .* label_vec .* "_surface"


for i in eachindex(data_path_vec)

    # get data and essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_vec[i])

    # get autocovariance function array
    sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(
                        size_data,
                        volume_fract_tot,
                        data_binary,
                        nr_measurements_per_direction=5000,
                        save_result = true,
                        save_path = save_path_vec[i],
                        voxel_edge_length = 10,
                        label = label_vec[i])

    println(label_vec[i]*" done")

end



# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

# set path where autocovariance dict of single unit cell is stored
suc_path_vec = raw"..\analysis_data\nodal_surfaces\single_unit_cell_" .* label_vec .* "_surface_autocovariance_fct_direction.h5"

save_path_vec = raw"..\analysis_data\nodal_surfaces\\" .* label_vec .* "_surface"


for i in eachindex(suc_path_vec)

    # load dict
    suc_dict = BDA.load_h5_dict(suc_path_vec[i])

    # get autocovariance function array
    sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.extrapolate_periodic_data_autocovariance_fct_by_sampling_vec_array(
                suc_dict;
                save_result = true,
                save_path = save_path_vec[i])

    println(label_vec[i]*" done")

end



# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 


# path where analysis data will be saved
data_path_vec = raw"..\structures\nodal_surfaces\\" .* label_vec .* "_surface.h5"

# path where analysis data will be saved
save_path_vec = raw"..\analysis_data\nodal_surfaces\\" .* label_vec .* "_surface"

for i in eachindex(data_path_vec)
    
    # calculate and save all statistical measures
    BDA.save_statistical_measures(data_path_vec[i], 
                                    10, 
                                    label_vec[i], 
                                    save_path_vec[i];
                                    save_autocovariance_fct = false,
                                    save_spectral_density = false,
                                    save_local_volume_fraction_variance = false,
                                    save_autocovariance_fct_direction = false,
                                    save_spectral_density_along_directions = true)


end


# set paths where statistical data is stored
plot_path = raw"..\plots\nodal_surfaces\\"

# plot all statistical measures
BDA.plot_statistical_measures(save_path_vec,
                    plot_path;
                    spectral_density_xlims = [0,0.1],
                    plot_autocovariance_fct_bool = false,
                    plot_spectral_density_bool = false,
                    plot_local_volume_fraction_variance_bool = false,
                    plot_autocovariance_fct_direction_bool = false,
                    plot_spectral_density_direction_bool = true
                    )



# compare hyperuniformity criterion for red and blue patches
data_path = raw"..\structures\nodal_surfaces\\"

# set surfaces that are analyzed
label_vec = ["D", "G", "P", "I-WP"] 

structure_dict_path_vec = data_path .* label_vec .* "_surface"

voxel_edge_length_vec = [10, 10, 10, 10]

for i in eachindex(structure_dict_path_vec)

    # load dictionary
    structure_dict = BDA.load_h5_dict(structure_dict_path_vec[i]* ".h5")

    # get essential information of data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(structure_dict["data"])

    # add this info to dictionary and save it
    structure_dict = Dict("data_binary" => structure_dict["data"], 
                            "volume_fract_tot" => volume_fract_tot, 
                            "size_data" => collect(size_data), 
                            "mean_edge_length_data" => mean_edge_length_data, 
                            "nr_dimensions_data" => nr_dimensions_data,
                            "voxel_edge_length" => voxel_edge_length_vec[i] ,
                            "label" => label_vec[i],
                            "unit_cell_length" => 500,
                            "nr_unit_cells" => 10,
                            "volume_fraction_parameter" => 0 )

    BDA.save_dict_to_h5(structure_dict; save_path=structure_dict_path_vec[i] * "_structure.h5")

end

# compare hyperuniformity criterion for red and blue patches
data_path = raw"..\structures\nodal_surfaces\\single_unit_cell_"

structure_dict_path_vec = data_path .* label_vec .* "_surface"

for i in eachindex(structure_dict_path_vec)

    # load dictionary
    structure_dict = BDA.load_h5_dict(structure_dict_path_vec[i]* ".h5")

    # get essential information of data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(structure_dict["data"])

    # add this info to dictionary and save it
    structure_dict = Dict("data_binary" => structure_dict["data"], 
                            "volume_fract_tot" => volume_fract_tot, 
                            "size_data" => collect(size_data), 
                            "mean_edge_length_data" => mean_edge_length_data, 
                            "nr_dimensions_data" => nr_dimensions_data,
                            "voxel_edge_length" => voxel_edge_length_vec[i] ,
                            "label" => label_vec[i],
                            "unit_cell_length" => 500,
                            "nr_unit_cells" => 1,
                            "volume_fraction_parameter" => 0 )

    BDA.save_dict_to_h5(structure_dict; save_path=structure_dict_path_vec[i] * "_structure.h5")

end


# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"]  # 

# set properties of generated data
unit_cell_length = 500
nr_unit_cells = 10

for label in label_vec

    # generate 3d binary data for current nodal surface
    structure_dict = BDA.get_binary_data_from_nodal_eqn(unit_cell_length, 
                                                nr_unit_cells,
                                                label;
    save_result=true, 
    save_path=raw"..\structures\nodal_surfaces\\"*label*"_surface_structure.h5")

    println(label*" done")
    
end


"""
These are the calculations for the stervi beetle from the 10.1002/adfm.202302720 paper.
The data was sent by Viola and is not directly taken from Zenodo
"""


# set raw data path
data_path_raw_prefix = raw"..\structures\stervi\2d_images_green\slice_"
data_path_raw_suffix = "_max11.tif"

data_path_corrected = raw"..\structures\stervi\stervi_green.h5"

# load data, correct voxel size anisotropy and save data
data_binary = BDA.get_binary_data_from_colorscale_stack(data_path_raw_prefix,
                                        data_path_raw_suffix,
                                        301; 
                                        save_data=true, 
                                        data_path_corrected=data_path_corrected)


                                        
# set raw data path
data_path_raw_prefix = raw"..\structures\stervi\2d_images_blue\slice_"
data_path_raw_suffix = "_max11.tif"

data_path_corrected = raw"..\structures\stervi\stervi_blue.h5"

# load data, correct voxel size anisotropy and save data
data_binary = BDA.get_binary_data_from_colorscale_stack(data_path_raw_prefix,
                                        data_path_raw_suffix,
                                        250; 
                                        save_data=true, 
                                        data_path_corrected=data_path_corrected)


# loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )

    
    # get autocovariance function as a function of sampling distance
    sampling_distance_vec, autocovariance_fct_vec = BDA.get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                                                                    size_data, 
                                                                                    volume_fract_tot,
                                                                                    data_binary;
                                                                                    nr_measurements_per_distance = 10000)

    # create dictionary for current plot
    plot_dict = Dict("sampling_distance_vec" => sampling_distance_vec,
                    "autocovariance_fct_vec" => autocovariance_fct_vec,
                    "voxel_edge_length" => 11,
                    "label" => "green" )

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="stervi_autocovariance_fct_green")


    # calculate spectral density for loaded autocovariance function
    wavenumber_vec, spectral_density_vec = BDA.get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data,
                size_data, 
                volume_fract_tot,
                data_binary;
                nr_measurements_per_distance = 10000,
                nr_wavenumbers=200,
                sampling_distance_vec = sampling_distance_vec,
                autocovariance_fct_vec = autocovariance_fct_vec)

    # create dictionary for current plot
    plot_dict = Dict("wavenumber_vec" => wavenumber_vec,
                    "spectral_density_vec" => spectral_density_vec,
                    "voxel_edge_length" => 11,
                    "label" => "green" )

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
                        save_filename="stervi_spectral_density_green")


    # determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 100,
                                                                                                    window_positioning="random",
                                                                                                    window_shape="spherical")

    # create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "voxel_edge_length" => 11,
                    "label" => "green" )

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="stervi_volume_fraction_variance_random_spherical_green")

                                                                                            
end


# compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"..\analysis_data\stervi_volume_fraction_variance_random_spherical_green.h5",
                raw"..\analysis_data\volume_fraction_variance_random_spherical_I-WP_surface.h5"]

# set surfaces that are analyzed
label_vec = ["S. virescens green", "perfect I-WP"]

# set edge lengths
voxel_edge_length_vec = [11, 6]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


# loop through data that will be analyzed
for i in eachindex(data_path_vec)
    
    # load dictionary that contains the following keys:
    # "window_edge_length_vec"
    # "local_volume_fract_variance_vec"
    # "voxel_edge_length"
    # "label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    # adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end

# plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="stervi_volume_fraction_variance_random_spherical_windows_green_i_wp")



# compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"..\analysis_data\stervi_spectral_density_green.h5",
                raw"..\analysis_data\spectral_density_I-WP_surface.h5"]

# set surfaces that are analyzed
label_vec = ["S. virescens green", "perfect I-WP"]

# set edge lengths
voxel_edge_length_vec = [11, 6]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


# loop through data that will be analyzed
for i in eachindex(data_path_vec)
    
    # load dictionary that contains the following keys:
    # "wavenumber_vec"
    # "spectral_density_vec"
    # "voxel_edge_length"
    # "label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    # adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


# plot real part, imaginary part and absolute value of spectral 
# density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_filename="stervi_spectral_density_green_i_wp")


# path of original data
data_path = raw"..\structures\stervi\stervi_green.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\stervi\stervi_green"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "green", 
                                save_path)

# set paths where statistical data is stored
data_path_vec = [raw"..\analysis_data\stervi\stervi_blue",
raw"..\analysis_data\nodal_surfaces\I-WP_surface",
raw"..\analysis_data\stervi\stervi_green"]

# set path where plot will be stored
save_path = raw"..\plots\stervi\stervi_blue_green_i_wp"

# plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            voxel_edge_length_vec=[11,6,11],
            label_vec=["S. virescens blue", "perfect I-WP", "S. virescens green"]
            )



# path of original data
data_path = raw"..\structures\stervi\stervi_green.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\stervi\stervi_green"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "green", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)


# path of original data
data_path = raw"..\structures\stervi\stervi_blue.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\stervi\stervi_blue"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "blue", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)



# path of original data
data_path = raw"..\structures\stervi\stervi_blue.h5"

save_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"

# get essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

# get autocovariance function array
sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(size_data,
                       volume_fract_tot,
                       data_binary,
                       nr_measurements_per_direction=10)

# create dict to save
saving_dict = Dict("sampling_vec_array" => sampling_vec_array,
                    "autocovariance_fct_array" => autocovariance_fct_array,
                    "voxel_edge_length" => 11,
                    "label" => "stervi blue")

save_dict_to_h5(saving_dict; save_path)


# path where autocovariance fct data is saved
dict_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"

# load dict
data_dict = BDA.load_h5_dict(dict_path)

# get sampling vec array and autocovariance fct array from dict
sampling_vec_array = data_dict["sampling_vec_array"]
autocovariance_fct_array = data_dict["autocovariance_fct"]

# calculate spectral density
sampled_wavevectors_array, spectral_density_array = get_spectral_density(size_data, 
                                                volume_fract_tot,
                                                data_binary;
                                                sampling_vec_array = sampling_vec_array,
                                                autocovariance_fct_array)

# path where spectral density is saved
save_path = raw"..\analysis_data\stervi\stervi_spectral_density_direction_blue_small_sampling.h5"

# create dict to save
saving_dict = Dict("sampled_wavevectors_array" => sampled_wavevectors_array,
                    "spectral_density_array" => spectral_density_array,
                    "voxel_edge_length" => 11,
                    "label" => "stervi blue")

save_dict_to_h5(saving_dict; save_path)



# set data path
data_path = raw"..\structures\stervi\stervi_blue.h5"

# get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)


# path of dict
load_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"

# load dict
data_dict = BDA.load_h5_dict(load_path)

autocovariance_fct_array = data_dict["autocovariance_fct_array"]

# calculate spectral density
sampled_wavenumbers_vec_vec, sampled_wavevectors_array, spectral_density_array = BDA.get_spectral_density(size_data, 
volume_fract_tot,
data_binary;
nr_measurements_per_direction = 50,
autocovariance_fct_array = autocovariance_fct_array)


# create dictionary for current plot
plot_dict = Dict("sampled_wavenumbers_vec_vec" => sampled_wavenumbers_vec_vec,
                "sampled_wavevectors_array" => sampled_wavevectors_array,
                "spectral_density_array" => spectral_density_array,
                "voxel_edge_length" => 11,
                "label" => data_dict["label"])


# path of dict
save_path = raw"..\analysis_data\stervi\stervi_spectral_density_direction_blue_small_sampling.h5"

# save the plot_dict to a H5 file
BDA.save_dict_to_h5(copy(plot_dict); save_path=save_path)



# plot heat map of spectral density in x-y-plane for k_z=0
BDA.plot_spectral_density_heatmap(plot_dict,
                                save_path;
                                save_plot = false)

# wait until key is pressed
readline()




autocovariance_fct_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"
autocovariance_fct_dict = BDA.load_h5_dict(autocovariance_fct_path)

spectral_density_path = raw"..\analysis_data\stervi\stervi_spectral_density_direction_blue_small_sampling.h5"
spectral_density_dict = BDA.load_h5_dict(spectral_density_path)


# set data path
data_path = raw"..\structures\stervi\stervi_blue.h5"

# get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)

sampling_distance_vec_vec = BDA.get_sampling_distance_vec_vec(size_data)

autocovariance_fct_dict["sampling_distance_vec_1"] = sampling_distance_vec_vec[1]
autocovariance_fct_dict["sampling_distance_vec_2"] = sampling_distance_vec_vec[2]
autocovariance_fct_dict["sampling_distance_vec_3"] = sampling_distance_vec_vec[3]

# save the plot_dict to a H5 file
BDA.save_dict_to_h5(copy(autocovariance_fct_dict); save_path=autocovariance_fct_path)



# set data path
data_path = raw"..\structures\stervi\stervi_blue.h5"

# set path where autocovariance dict will be stored
save_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue.h5"

# get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)

# get autocovariance function array
sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(
                    size_data,
                    volume_fract_tot,
                    data_binary,
                    nr_measurements_per_direction=100)

# create dict to save
saving_dict = Dict("sampling_vec_array" => sampling_vec_array,
                    "sampling_distance_vec_vec" => sampling_distance_vec_vec,
                    "autocovariance_fct_array" => autocovariance_fct_array,
                    "voxel_edge_length" => 11,
                    "label" => "stervi blue")

BDA.save_dict_to_h5(saving_dict; save_path)


# set data path
data_path = raw"..\structures\stervi\stervi_blue.h5"

# set path where autocovariance dict will be stored
dict_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue.h5"

# get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)

# load autocovariance array
loaded_dict = BDA.load_h5_dict(dict_path)

plot_path = raw"..\plots\stervi\stervi_direction_blue"


BDA.plot_autocovariance_fct_heatmap(loaded_dict,
    plot_path;
    save_plot = false,
    sampling_vector_component_to_fix = 3,
    sampling_vector_value_fixed = 0 )


    
# set data path
data_path_vec = [raw"..\structures\stervi\stervi_blue.h5",
raw"..\structures\stervi\stervi_green.h5",
raw"..\structures\pachy\pachy_blue.h5",
raw"..\structures\pachy\pachy_red.h5"
]

# set path where autocovariance dict will be stored
save_path_vec = [raw"..\analysis_data\stervi\stervi_blue_autocovariance_fct_direction.h5",
raw"..\analysis_data\stervi\stervi_green_autocovariance_fct_direction.h5",
raw"..\analysis_data\pachy\pachy_blue_autocovariance_fct_direction.h5",
raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"
]

voxel_edge_length_vec = [11,11,10,9]

label_vec = ["stervi blue", "stervi green", "pachy blue","pachy red"]

for i in eachindex(data_path_vec)

    # get data and essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_vec[i])

    # get autocovariance function array
    sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(
                        size_data,
                        volume_fract_tot,
                        data_binary,
                        nr_measurements_per_direction=1000)

    # create dict to save
    saving_dict = Dict("sampling_vec_array" => sampling_vec_array,
                        "sampling_distance_vec_vec" => sampling_distance_vec_vec,
                        "autocovariance_fct_array" => autocovariance_fct_array,
                        "voxel_edge_length" => voxel_edge_length_vec[i],
                        "label" => label_vec[i])

    BDA.save_dict_to_h5(saving_dict; save_path=save_path_vec[i])

    println(label_vec[i]*" done")

end


# set data path
data_path_vec = [raw"..\analysis_data\stervi\stervi_blue",
                raw"..\analysis_data\stervi\stervi_green"]

# set path to save plot
save_path = raw"..\plots\stervi\\"

BDA.plot_statistical_measures(data_path_vec,
            save_path;
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = true
            )

            
# path of original data
data_path = raw"..\structures\stervi\stervi_blue.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\stervi\stervi_blue"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "S. virescens blue", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)

# path of original data
data_path = raw"..\structures\stervi\stervi_green.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\stervi\stervi_green"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "S. virescens green", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)


# set paths where statistical data is stored
data_path_vec = [raw"..\analysis_data\stervi\stervi_green",
                raw"..\analysis_data\stervi\stervi_blue"]

# set path where plot will be stored
save_path = raw"..\plots\stervi\\"

# plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            spectral_density_xlims = [0,0.1],
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = false,
            plot_spectral_density_direction_bool = true
            )


"""
The following commands targeted multiple samples at the same time
"""


data_path_vec = [ raw"..\structures\pachy\pachy_red_structure.h5",
raw"..\structures\pachy\pachy_blue_structure.h5",
raw"..\structures\stervi\stervi_green_structure.h5",
raw"..\structures\stervi\stervi_blue_structure.h5",
raw"..\structures\nodal_surfaces\D_surface_structure.h5",
raw"..\structures\nodal_surfaces\G_surface_structure.h5",
raw"..\structures\nodal_surfaces\P_surface_structure.h5",
raw"..\structures\nodal_surfaces\I-WP_surface_structure.h5"
]

save_path_vec = [ raw"..\analysis_data\pachy\pachy_red",
raw"..\analysis_data\pachy\pachy_blue",
raw"..\analysis_data\stervi\stervi_green",
raw"..\analysis_data\stervi\stervi_blue",
raw"..\analysis_data\nodal_surfaces\D_surface",
raw"..\analysis_data\nodal_surfaces\G_surface",
raw"..\analysis_data\nodal_surfaces\P_surface",
raw"..\analysis_data\nodal_surfaces\I-WP_surface"
]

for i in eachindex(data_path_vec)

    BDA.save_statistical_measures(data_path_vec[i],
        save_path_vec[i];
        save_autocovariance_fct = false,
        save_spectral_density = true,
        save_local_volume_fraction_variance = false,
        save_autocovariance_fct_direction = false,
        save_spectral_density_along_directions = false,
        save_complete_autocovariance_fct_direction = true,
        save_spectral_density_array = true)

    println(string(i)*" done")

end


data_path_vec_vec =  [ [raw"..\analysis_data\pachy\pachy_blue",
    raw"..\analysis_data\pachy\pachy_red",
    raw"..\analysis_data\nodal_surfaces\D_surface"],
[
raw"..\analysis_data\stervi\stervi_blue",
raw"..\analysis_data\nodal_surfaces\I-WP_surface",
raw"..\analysis_data\stervi\stervi_green"],
[ raw"..\analysis_data\nodal_surfaces\P_surface",
    raw"..\analysis_data\nodal_surfaces\I-WP_surface",
    raw"..\analysis_data\nodal_surfaces\D_surface",
raw"..\analysis_data\nodal_surfaces\G_surface"
]
] 

save_path_vec = [ raw"..\plots\pachy\pachy_blue_red_d",
raw"..\plots\stervi\stervi_blue_green_i_wp",
raw"..\plots\nodal_surfaces\nodal_surfaces"
]

voxel_edge_length_vec = [
    [10, 9, 8.5],
    [11,6,11],
    nothing
]

for i in eachindex(data_path_vec_vec)

    BDA.plot_statistical_measures(data_path_vec_vec[i],
        save_path_vec[i];
        voxel_edge_length_vec=voxel_edge_length_vec[i],
        label_vec=nothing,
        spectral_density_xlims = (0,0.1),
        spectral_density_heatmaps_clims = nothing,
        plot_autocovariance_fct_bool = false,
        plot_spectral_density_bool = true,
        plot_local_volume_fraction_variance_bool = false,
        plot_autocovariance_fct_direction_bool = false,
        plot_spectral_density_along_directions_bool = false,
        plot_spectral_density_heatmaps_bool = false
        )

    println(string(i)*" done")

end


# set raw data path
data_path_raw_prefix = raw"..\structures\stervi\2d_images_green\slice_"
data_path_raw_suffix = "_max11.tif"

data_path_corrected = raw"..\structures\stervi\stervi_green.h5"

structure_dict = BDA.get_structure_dict_from_colorscale_stack(data_path_raw_prefix,
    data_path_raw_suffix,
    301; 
    voxel_size=(11,11,11), 
    label = "S. virescens green",
    save_result=true, 
    save_path = raw"..\structures\stervi\stervi_green")



"""
This is where commands for network generation begin
"""


# import my module that contains all functions for the analysis of binary structure data
import .NetworkGeneration as NG
import MetaGraphsNext as MGN


network_dict = NG.get_periodic_network( ; nr_vertices = 150 , 
                            nr_dimensions = 3, 
                            network_type = "diamond")

network_dict["bond_bending_const"] = 0.285

# get list of bonds (edges)
edges_vec = collect(MGN.edge_labels(network_dict["network_graph"]))

# break a random bond
network_dict = NG.switch_bond!(network_dict, edges_vec[4] )

vertex = 1

neighbor_positions_mat = NG.get_neighbor_positions_mat(network_dict, vertex)

local_energy = NG.local_energy_keating(network_dict["network_graph"][vertex], 
        network_dict, neighbor_positions_mat)

gradient = zeros(3)

NG.gradient_keating!(gradient, network_dict["network_graph"][vertex], 
        network_dict, neighbor_positions_mat)

hessian = zeros(3, 3)

NG.hessian_keating!(hessian, network_dict["network_graph"][vertex], 
        network_dict, neighbor_positions_mat)



println(local_energy)
println(gradient)
println(hessian)



graph_dict = NG.get_periodic_network( ; nr_vertices = 150 ,
                            network_type = "diamond")

                            
vertex = 5

# get and print neighbors
vertex_neighbors = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex))

println(vertex_neighbors)

# get position of some vertex
println(graph_dict["spatial_network"][vertex]["position"] )
println(graph_dict["spatial_network"][vertex]["local_energy"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex, 
[1,-0.5,0.3] )

println(graph_dict["spatial_network"][vertex]["position"] )
println(graph_dict["spatial_network"][vertex]["local_energy"] )


# relax vertex
graph_dict = NG.relax_single_vertex_keating!(graph_dict, vertex)

println(graph_dict["spatial_network"][vertex]["position"] )
println(graph_dict["spatial_network"][vertex]["local_energy"] )

graph_dict = NG.get_periodic_network( ; nr_vertices = 1500 ,
                            network_type = "diamond")

                            
vertex_vec = [5,12]

# get and print neighbors
neighbors_in_shells_dict, all_vertices_vec = NG.get_neighbors_in_shells_dict(graph_dict, 
                                    vertex_vec; 
                                    shell_nr = 4)


vertex_vec = [5,12]


# get position of some vertex
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[1], 
[1,-0.5,0.3] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[2], 
[-0.1,0.5,1.5] )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

graph_dict = NG.relax_cluster_one_cycle_keating!(graph_dict, 
vertex_vec; 
shell_nr = 4 )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

graph_dict = NG.relax_cluster_one_cycle_keating!(graph_dict, 
vertex_vec; 
shell_nr = 4 )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )




graph_dict = NG.get_periodic_network( ; nr_vertices = 1500 ,
                            network_type = "diamond")

                            
vertex_vec = [5,12]


# get position of some vertex
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[1], 
[1,-0.5,0.3] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[2], 
[-0.1,0.5,1.5] )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

@time graph_dict = NG.relax_cluster_keating!(graph_dict, vertex_vec; 
nr_cycles = 10,
reject_during_relaxation_cycle_threshold = 5,
shell_nr = 4)

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )


                        
vertex_vec = [5,12]


# get position of some vertex
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[1], 
[1,-0.5,0.3] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[2], 
[-0.1,0.5,1.5] )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )

graph_dict = NG.relax_cluster_one_cycle_keating!(graph_dict, 
vertex_vec; 
shell_nr = 4 )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )

graph_dict = NG.relax_cluster_one_cycle_keating!(graph_dict, 
vertex_vec; 
shell_nr = 4 )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )



vertex_to_relax = 5

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_to_relax, 
[0.2,0.5,-0.9] )

# get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

# get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

# get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

# set energy, gradient and hessian for energy minimization
energy(x) = NG.energy_from_position_keating(x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr )
                                            
gradient!(gradient, x) = NG.gradient_keating!(gradient, x, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)
hessian!(hessian, x) = NG.hessian_keating!(hessian, x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

gradient_an = zeros(3)
gradient!(gradient_an, initial_position)
gradient_num = ForwardDiff.gradient(energy, initial_position)

println(gradient_an)
println(gradient_num)
println("_______")

hessian_an = zeros(3, 3)
hessian!(hessian_an, initial_position)
hessian_num = ForwardDiff.hessian(energy, initial_position)

println(hessian_an)
println(hessian_num)



println(graph_dict["spatial_network"][vertex]["position"])

neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex;
                                    exclude_vertices = [])

next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex)

neighbors_vec = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex))
println(neighbors_vec)

neighbor_nr = 1

next_neighbors_vec = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], neighbors_vec[neighbor_nr]))

for i in 1:4
    println(graph_dict["spatial_network"][next_neighbors_vec[i]]["position"])

end

println("_______")

for i in 1:3
    println(next_neighbor_positions_arr[neighbor_nr,:,i])

end



println("_______")

vertex = 2

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex, 
[1,-0.5,0.3]; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )


# relax vertex
graph_dict = NG.relax_single_vertex_keating!(graph_dict, vertex; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )

println("_______")

vertex = 3

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex, 
[1,-0.5,0.3]; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )


# relax vertex
graph_dict = NG.relax_single_vertex_keating!(graph_dict, vertex; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )



neighbors_vec = collect(MetaGraphsNext.neighbor_labels(
                                        graph_dict["spatial_network"], vertex))

for neighbor in neighbors_vec
    display(graph_dict["spatial_network"][neighbor]["position"])
end
println("_______")

neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex)

display(neighbor_positions_mat)

next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex)

display(next_neighbor_positions_arr[1,:,:])



random_bond = NG.get_random_bond(graph_dict)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][random_bond[1]]["position"] )

temperature = 100.5

graph_dict, move_accepted = NG.monte_carlo_move(graph_dict, 
    temperature; 
    switched_bond = random_bond)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][random_bond[1]]["position"] )


central_vertex = 100

central_vertex_position = graph_dict["spatial_network"][central_vertex]["position"]

    # get central vertices neighbors 
    neighbor_vec = collect(MetaGraphsNext.neighbor_labels(
                                graph_dict["spatial_network"], central_vertex))
# create array to store next to nearest neighbors coordinates
# The first array index labels the number of the direct neighbor
next_neighbor_positions_arr = Array{Float64}(undef, 
                                            graph_dict["coordination_nr"],
                                            graph_dict["nr_dimensions"],
                                            graph_dict["coordination_nr"]-1)

# loop through central vertices neighbors
for i in 1:graph_dict["coordination_nr"]
    current_next_neighbor = 1
    # loop through the current neighbor's neighbors
    for next_neighbor in MetaGraphsNext.neighbor_labels(
                                    graph_dict["spatial_network"], neighbor_vec[i])
        if next_neighbor !== central_vertex
            # get next neighbor's virtual coordinates which might be outside of the 
            # supercell if periodic boundary conditions play a role
            next_neighbor_positions_arr[i,:,current_next_neighbor] = NG.get_virtual_position(
                        central_vertex_position,
                        graph_dict["spatial_network"][next_neighbor]["position"],
                        graph_dict["supercell_edge_length"] )
            current_next_neighbor += 1
        end
    end
end



# get cluster after bond switch
cluster_dict = NG.get_cluster_in_shells_dict(
                                graph_dict, 
                                switched_bond; 
                                shell_nr = shell_nr)

# relax cluster around switched bond
graph_dict = NG.relax_cluster_keating!(graph_dict,
    switched_bond; 
    nr_cycles = nr_cycles,
    reject_during_relaxation_cycle_threshold = reject_during_relaxation_cycle_threshold,
    shell_nr = shell_nr,
    cluster_dict = cluster_dict,
    initial_cluster_energy  = initial_cluster_dict["cluster_energy"],
    update_total_energy = true)

    
switched_bond = (5,12)

display(collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[1]) ))
    
display(collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[2]) ))

println("______")

graph_dict = NG.switch_bond!(graph_dict, switched_bond)

display(collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[1]) ))

display(collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[2]) ))


random_bond = NG.get_random_bond(graph_dict)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][random_bond[1]]["position"] )

temperature = 100.5

graph_dict, move_accepted = NG.monte_carlo_move(graph_dict, 
    temperature; 
    switched_bond = random_bond)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][random_bond[1]]["position"] )


graph_dict = NG.get_periodic_network( ; nr_vertices = 40 , 
nr_dimensions = 3, 
network_type = "diamond",
bond_bending_const = 0.285)

figure = NG.plot_network(graph_dict)

# figure = NG.plot_network(graph_dict)

switched_bond = (5,12) # NG.get_random_bond(graph_dict)

temperature = 100.5

graph_dict, move_accepted, new_bond_vec = NG.monte_carlo_move(graph_dict, 
    temperature; 
    switched_bond = switched_bond)

figure = NG.plot_network(graph_dict; highlight_nodes = switched_bond, highlight_edges = [switched_bond, new_bond_vec...])

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network(graph_dict,
    10, 
    temperature;
    shell_nr = 3,
    print_progress = true)

    
NG.save_graph_to_csv(graph_dict,
"example_graph.csv")


@time graph_dict, cluster_energy_vec = NG.relax_cluster_keating!(graph_dict,
random_bond; 
relax_efficiently = true,
    update_total_energy = true,
    track_cluster_energy = true)

Plots.plot(collect(1:26), cluster_energy_vec)


# compare the gradient in the inefficient and efficient calculation
vertex_to_relax = random_bond[2]

# efficient calculation
gradient_eff = NG.gradient_keating_efficient(graph_dict, vertex_to_relax)

# inefficient calculation

# get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

# get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

# get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

# set energy, gradient and hessian for energy minimization

gradient_ineff = zeros(3)
                                        
NG.gradient_keating!(gradient_ineff, initial_position, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

println("Efficient: "*string(gradient_eff))
println("Inefficient "*string(gradient_ineff))


random_bond = NG.get_random_bond(graph_dict)

graph_dict, new_bond_vec = NG.switch_bond!(graph_dict, random_bond )


# compare the gradient in the inefficient and efficient calculation
vertex_to_relax = random_bond[2]

# efficient calculation
gradient = NG.gradient_keating_efficient(graph_dict, vertex_to_relax)

translation_vector_eff =  NG.get_approximate_translation_vector_keating(gradient, 
        graph_dict["bond_bending_const"];
        relaxation_overshoot_factor_r = 1.5,
        relaxation_optimization_parameter_l = 1)

# inefficient calculation

# get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

# get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

# get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

# set energy, gradient and hessian for energy minimization
energy(x) = NG.energy_from_position_keating(x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr )
                                            
gradient!(gradient, x) = NG.gradient_keating!(gradient, x, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

hessian!(hessian, x) = NG.hessian_keating!(hessian, x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)
# find energy minimum
minimizer_result = Optim.optimize(
                            energy, 
                            gradient!, 
                            hessian!,
                            initial_position, 
                            Optim.Newton())

# get relaxed position and local keating energy
relaxed_position = Optim.minimizer(minimizer_result)

# calculate translation vector for relaxed vertex
translation_vector_ineff = relaxed_position .- initial_position

println("Efficient: "*string(translation_vector_eff))
println("Inefficient "*string(translation_vector_ineff))


@time graph_dict, cluster_energy_vec = NG.relax_cluster_keating!(graph_dict,
random_bond; 
nr_max_relaxation_cycles = 1,
relax_efficiently = true,
    update_total_energy = true,
    track_cluster_energy = true,
    relaxation_optimization_parameter_l=1.3)

figure = NG.plot_network(graph_dict)


random_bond = NG.get_random_bond(graph_dict)

vertex_to_relax = random_bond[2]

graph_dict, new_bond_vec = NG.switch_bond!(graph_dict, random_bond )


# efficient calculation
gradient = NG.gradient_keating_efficient(graph_dict, vertex_to_relax)

hessian = NG.hessian_keating_efficient(graph_dict, vertex_to_relax)

# calculate translation vector to approximate energy minimum
translation_vector_eff = .- LinearAlgebra.inv(hessian)*gradient

# inefficient calculation

# get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

# get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

# get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

# set energy, gradient and hessian for energy minimization
energy(x) = NG.energy_from_position_keating(x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr )
                                            
gradient!(gradient, x) = NG.gradient_keating!(gradient, x, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

hessian!(hessian, x) = NG.hessian_keating!(hessian, x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)
# find energy minimum
minimizer_result = Optim.optimize(
                            energy, 
                            gradient!, 
                            hessian!,
                            initial_position, 
                            Optim.Newton())

# get relaxed position and local keating energy
relaxed_position = Optim.minimizer(minimizer_result)

# calculate translation vector for relaxed vertex
translation_vector_ineff = relaxed_position .- initial_position

println("Efficient: "*string(translation_vector_eff))
println("Inefficient "*string(translation_vector_ineff))



# efficient calculation
neighbor_vec = collect(
        MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], central_vertex))

j=1

# get vector pointing from central vertex to neighbor j
distance_vector_j_eff = (sign(neighbor_vec[j] - central_vertex)
* graph_dict["spatial_network"][central_vertex, neighbor_vec[j]]["vector"])

bond_stretching_term_eff = ( - 3/4 * ( 
            graph_dict["spatial_network"][central_vertex, neighbor_vec[j]]["distance_squared"] - 1 
            ) ) .* distance_vector_j_eff

# inefficient calculation

# get initial position of vertex to relax 
x = graph_dict["spatial_network"][vertex_to_relax]["position"]

# get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

# get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

# set energy, gradient and hessian for energy minimization

# get vector pointing from central vertex to neighbor
distance_vector_j = neighbor_positions_mat[:,j] .- x

# get bond stretching term
bond_stretching_term_ineff = ( - 3/4 * ( LinearAlgebra.norm(distance_vector_j)^2 - 1 ) 
                                ) .* distance_vector_j

println("Efficient: "*string(distance_vector_j_eff))
println("Inefficient "*string(distance_vector_j))



switched_bond = NG.get_random_bond(graph_dict, seed = 1)

graph_dict, new_bond_vec = NG.switch_bond!(graph_dict, switched_bond )

vertexic_position_arr, cluster_energy_arr = NG.compare_relaxation_methods(graph_dict,
switched_bond,
"diamond_disordered_t10_30steps_1" )

temperature = 8

graph_dict, move_accepted, new_bond_vec = NG.monte_carlo_move!(graph_dict, 
    temperature; 
    reject_during_relaxation_cycle_threshold = 10,
        break_at_relative_cluster_energy_change = 0.0001,
        shell_nr = 3,
        relax_efficiently = true,
        thermal_fluctuations = false)

println(graph_dict["total_energy"])
println(NG.get_total_energy_keating(graph_dict))


switched_bond = NG.get_random_bond(graph_dict, seed = 9)
shell_nr = 4

initial_cluster_dict = NG.get_cluster_in_shells_dict(
                graph_dict, 
                switched_bond; 
                shell_nr = shell_nr)


# switch bond
graph_dict, new_bond_vec = NG.switch_bond!(graph_dict, switched_bond )

# get cluster after bond switch
cluster_dict = NG.get_cluster_in_shells_dict(
                graph_dict, 
                switched_bond; 
                shell_nr = shell_nr)

# relax cluster once and update cluster energy
graph_dict, cluster_dict = NG.relax_cluster_keating!(graph_dict,
cluster_dict; 
nr_max_relaxation_cycles = 25,
break_at_relative_cluster_energy_change = 0.001,
reject_during_relaxation_cycle_threshold = 10,
relax_efficiently = true,
update_total_energy = true)

# calculate new total energy and compare to actual total energy
smart_total_energy = graph_dict["total_energy"] 

actual_total_energy = NG.get_total_energy_keating(graph_dict)

println(string(smart_total_energy))
println(string(actual_total_energy))


temperature = 5
shell_nr = 3


graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network(graph_dict,
    20, 
    temperature; 
    nr_max_relaxation_cycles = 25,
        break_at_relative_cluster_energy_change = 0.001,
        reject_during_relaxation_cycle_threshold = 10,
        relax_efficiently = true,
        shell_nr = shell_nr,
    print_progress = true,
    random_evolution_seed = 3,
    thermal_fluctuations = false)

# calculate new total energy and compare to actual total energy
smart_total_energy = graph_dict["total_energy"] 

actual_total_energy = NG.get_total_energy_keating(graph_dict)

println(string(smart_total_energy))
println(string(actual_total_energy))


structure_factor_dict = NA.get_structure_factor_isotrope_by_wavenumber_vec(
        graph_dict)

hyperuniformity_parameter = NA.get_effective_hyperuniformity_parameter(structure_factor_dict)
println("hyperuniformity parameter: "*string(hyperuniformity_parameter))

local_nr_variance_dict = NA.get_local_nr_variance_by_window_radius_vec(graph_dict;
structure_factor_dict = structure_factor_dict)

Plots.plot(local_nr_variance_dict["window_radius_vec"], 
local_nr_variance_dict["local_nr_variance_vec"] ./ local_nr_variance_dict["window_radius_vec"].^3)

NG.plot_network(graph_dict)


graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence(graph_dict,
        evolution_dict;
        print_progress = true,
        print_every_nr_attempted_bond_switches = 50)

evolution_dict["temperature_vec"] = [0]
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = [30]

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence(graph_dict,
        evolution_dict;
        move_accepted_vec = move_accepted_vec,
        total_energy_vec = total_energy_vec,
        print_progress = true,
        print_every_nr_attempted_bond_switches = 50)

Plots.plot(collect(1:length(total_energy_vec)), total_energy_vec)

NG.save_graph_to_h5_and_MGformat(deepcopy(graph_dict), "example")


evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [0.1],
nr_monte_carlo_steps_per_temperature_vec = [1])

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 20)

Plots.plot(collect(1:length(total_energy_vec)), total_energy_vec)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "64_vertices_T_0.1"

save_path = raw"..\structures\random_networks\\"

NG.save_mesh_from_network(graph_dict_to_save, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict_to_save,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)
        


network_names = ["64_vertices_T_0.1", "1000_vertices_T_1_quenched", "1000_vertices_T_2_quenched", "1000_vertices_T_4_quenched"]

save_path = raw"..\structures\random_networks\\"

for name in network_names
    graph_dict_to_save = NG.load_graph_from_h5_and_MGformat(save_path*name)

    NG.save_mesh_from_network(graph_dict_to_save, name; save_path = save_path)

end


evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 ,temperature_vec = [0.25, 0],
nr_monte_carlo_steps_per_temperature_vec = [3, 50])

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
        move_accepted_vec=move_accepted_vec,
        total_energy_vec=total_energy_vec,
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

Plots.plot(collect(1:length(total_energy_vec)), total_energy_vec)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "1000_vertices_T_0.25_quenched"

save_path = raw"..\structures\random_networks\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)

        
nr_samples = 10

evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 )

evolution_dict["temperature_vec"] = zeros(nr_samples) .+ 1
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = ones(Int64, nr_samples)
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][1] = 3 

graph_dict = NG.get_periodic_network(evolution_dict)


graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200,
    save_network_after_each_step = true,
    filename = "1000_vertices_T_1",)

nr_samples = 5

evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000)

evolution_dict["temperature_vec"] = zeros(nr_samples) .+ 0.5
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = ones(Int64, nr_samples)
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][1] = 3 

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200,
    save_network_after_each_step = true,
    filename = "1000_vertices_T_0.5",)

evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000)

evolution_dict["temperature_vec"] = zeros(nr_samples) .+ 0.015625
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = ones(Int64, nr_samples)
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][1] = 3 

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200,
    save_network_after_each_step = true,
    filename = "1000_vertices_T_0.015625",)


load_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

load_name = "1000_vertices_T_1_quenched"

# load dictionary with 1000 vertices which was heated to T=1 and then quenched
graph_dict = NG.load_graph_from_h5_and_MGformat(load_path*load_name)

save_path = raw"..\structures\random_networks\\"

NG.save_mesh_from_network(graph_dict, load_name*"_thick_bonds"; save_path = load_path, bond_radius = 0.3131)

structure_factor_dict = NA.get_structure_factor_isotrope_by_wavenumber_vec(
        graph_dict)

network_names = [ "1000_vertices_T_1_quenched", "1000_vertices_T_4_quenched"]

load_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

graph_dict_1 = NG.load_graph_from_h5_and_MGformat(load_path*network_names[1])

structure_factor_dict_1 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict_1;
    sampling_distance_step_length = 0.025,
    maximal_sampling_distance = 4*graph_dict_1["supercell_edge_length"],
    save_result = false,
    save_path = raw"..\analysis_data\random_networks\1000_vertices_T_1_quenched_high_sampling_rate",
    label = nothing)


graph_dict_4 = NG.load_graph_from_h5_and_MGformat(load_path*network_names[2])

structure_factor_dict_4 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict_4;
    sampling_distance_step_length = 0.025,
    maximal_sampling_distance = 4*graph_dict_4["supercell_edge_length"],
    save_result = false,
    save_path = raw"..\analysis_data\random_networks\1000_vertices_T_4_quenched_high_sampling_rate",
    label = nothing)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [1, 0],
nr_monte_carlo_steps_per_temperature_vec = [1, 50], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "216_vertices_T_1_quenched"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 512 ,temperature_vec = [1, 0],
nr_monte_carlo_steps_per_temperature_vec = [1, 50], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "512_vertices_T_1_quenched"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)


Plots.plot(collect(1:length(total_energy_vec)) ./ (512*18), total_energy_vec, xlabel = "steps", ylabel = "total energy", label = "ylabel")


# path where structures are stored
load_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

# load graph with 512 vertices
graph_dict_5 = NG.load_graph_from_h5_and_MGformat(load_path*"512_vertices_T_1_quenched")

# determine structure factor
structure_factor_dict_5 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict_5;
    sampling_distance_step_length = 0.025,
    maximal_sampling_distance = 4*graph_dict_5["supercell_edge_length"],
    save_result = true,
    save_path = raw"..\analysis_data\random_networks\512_vertices_T_1_quenched_high_sampling_rate",
    label = "512_vertices_T_1_quenched_high_sampling_rate")


# load graph with 216 vertices
graph_dict_2 = NG.load_graph_from_h5_and_MGformat(load_path*"216_vertices_T_1_quenched")

# determine structure factor
structure_factor_dict_2 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict_2;
    sampling_distance_step_length = 0.025,
    maximal_sampling_distance = 4*graph_dict_2["supercell_edge_length"],
    save_result = true,
    save_path = raw"..\analysis_data\random_networks\216_vertices_T_1_quenched_high_sampling_rate",
    label = "216_vertices_T_1_quenched_high_sampling_rate")


# load structure factor for 1000 vertices
dict_path_1 = raw"..\analysis_data\random_networks\1000_vertices_T_1_quenched_high_sampling_rate_structure_factor_isotrope.h5"

structure_factor_dict_1 = GU.load_h5_dict(dict_path_1)

dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.125, 0.25, 0.5, 1, 2, 4, 6, 8]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [2, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end




dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 6, 8]

for i in eachindex(temperatures)

    graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*"216_vertices_T_"*string(temperatures[i])*"_quenched")

    structure_factor_dict = NA.get_structure_factor_bartlett_isotrope_by_wavenumber_vec(
        graph_dict;
        sampling_distance_step_length = 0.05,
        maximal_sampling_distance = 4*graph_dict["supercell_edge_length"],
        save_result = true,
        save_path = raw"..\analysis_data\random_networks\216_vertices_T_"*string(temperatures[i])*"_quenched",
        label = "T = "*string(temperatures[i]))

end


evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [1],
nr_monte_carlo_steps_per_temperature_vec = [0.01], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "64_vertices_slight_disorder"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path, bond_radius = 0.3131)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [1],
nr_monte_carlo_steps_per_temperature_vec = [0.01], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

filename = "64_vertices_perfect_diamond_thick_bonds"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path, bond_radius = 0.3131)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [0.5, 0],
nr_monte_carlo_steps_per_temperature_vec = [0.1, 30], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
        total_energy_vec = total_energy_vec,
        move_accepted_vec = move_accepted_vec,
        print_progress = true,
    print_every_nr_attempted_bond_switches = 50)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "64_vertices_T_0.5_quenched_thick_bonds"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path, bond_radius = 0.3131)

NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [1],
nr_monte_carlo_steps_per_temperature_vec = [0.03], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
        total_energy_vec = total_energy_vec,
        move_accepted_vec = move_accepted_vec,
        print_progress = true,
    print_every_nr_attempted_bond_switches = 50)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "64_vertices_slightly_more_disorder_thick_bonds"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path, bond_radius = 0.3131)

NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.125, 0.25, 0.5, 1, 2, 4, 8]

evolution_dict = GU.load_h5_dict(dict_path*"216_vertices_T_"*string(temperatures[7])*"_quenched_evolution.h5")

Plots.plot(collect(1:length(evolution_dict["total_energy_vec"])) ./ (216*18), evolution_dict["total_energy_vec"], xlabel = "steps", ylabel = "total energy", label = "ylabel")




dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.1, 0.15, 0.2, 0.3, 0.4]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [2, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end


temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_cooling_gradient(temperature;
    temperature_decrease_per_monte_carlo_step = 0.1,
    quench = true )

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_cool_0.1_per_mc_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
    temperature_gradient = 0.1,
    quench = true )

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_cool_0.1_per_mc_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [0.01, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_heated_for_0.01_steps_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end


temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [0.05, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_heated_for_0.05_steps_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [5, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_heated_for_5.0_steps_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_heated_for_10.0_steps_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end


central_vertex = 10

l_max = 12

evolution_dict = NA.get_evolution_dict(nr_vertices = 64, network_type = "diamond")
graph_dict_diamond = NG.get_periodic_network(evolution_dict)

single_vertex_q_l_diamond = NA.get_q_l_averaged_single_vertex_dict(graph_dict_diamond,
central_vertex,
l_max)

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, single_vertex_q_l_diamond[i])
    println()
end

evolution_dict = NA.get_evolution_dict(nr_vertices = 16, network_type = "bcc")
graph_dict_bcc = NG.get_periodic_network(evolution_dict)

single_vertex_q_l_bcc = NA.get_q_l_averaged_single_vertex_dict(graph_dict_bcc,
central_vertex,
l_max)

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, single_vertex_q_l_bcc[i])
    println()
end

evolution_dict = NA.get_evolution_dict(nr_vertices = 16, network_type = "fcc")
graph_dict_fcc = NG.get_periodic_network(evolution_dict)

single_vertex_q_l_fcc = NA.get_q_l_averaged_single_vertex_dict(graph_dict_fcc,
central_vertex,
l_max)

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, single_vertex_q_l_fcc[i])
    println()
end

evolution_dict = NA.get_evolution_dict(nr_vertices = 16, network_type = "primitive cubic")
graph_dict_primitive_cubic = NG.get_periodic_network(evolution_dict)

single_vertex_q_l_primitive_cubic = NA.get_q_l_averaged_single_vertex_dict(graph_dict_primitive_cubic,
central_vertex,
l_max)

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, single_vertex_q_l_primitive_cubic[i])
    println()
end


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.75, 1.0, 2.0, 4.0]

for temperature in temperatures

    for nr_heating_mc_steps in [0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 5.0, 10.0]

        evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
        nr_monte_carlo_steps_per_temperature_vec = [nr_heating_mc_steps, 50], min_ring_size = 3)

        graph_dict = NG.get_periodic_network(evolution_dict)

        graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
                evolution_dict; 
            print_progress = true,
            print_every_nr_attempted_bond_switches = 500)

        evolution_dict["total_energy_vec"] = total_energy_vec
        evolution_dict["move_accepted_vec"] = move_accepted_vec

        NG.plot_network(graph_dict)

        filename = "216_vertices_T_"*string(temperature)*"_heated_for_"*string(nr_heating_mc_steps)*"_steps_quenched"

        save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

        NG.save_graph_to_h5_and_MGformat(graph_dict,
            filename;
            evolution_dict = evolution_dict,
            save_path 
                = save_path)

    end

end


temperatures = [0.75, 1.0, 2.0, 4.0]

for temperature in temperatures

    for temperature_gradient in [0.025, 0.05, 0.1, 0.2]

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_gradient = temperature_gradient, 
        nr_monte_carlo_steps_per_temperature = 0.01,
        quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
        nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

        graph_dict = NG.get_periodic_network(evolution_dict)

        graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
                evolution_dict; 
            print_progress = true,
            print_every_nr_attempted_bond_switches = 500)

        evolution_dict["total_energy_vec"] = total_energy_vec
        evolution_dict["move_accepted_vec"] = move_accepted_vec

        NG.plot_network(graph_dict)

        filename = "216_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

        save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

        NG.save_graph_to_h5_and_MGformat(graph_dict,
            filename;
            evolution_dict = evolution_dict,
            save_path 
                = save_path)

    end

end


temperatures = [0.75, 1.0, 2.0, 4.0]

for temperature in temperatures

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_cooling_gradient(temperature;
    temperature_decrease_per_monte_carlo_step = 0.1,
    quench = true )

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_cool_0.1_per_mc_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [1.0, 2.0, 4.0, 0.75]


for nr_heating_mc_steps in [0.5, 1.0, 5.0, 10.0 ]
    for temperature in temperatures

        evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
        nr_monte_carlo_steps_per_temperature_vec = [nr_heating_mc_steps, 50], min_ring_size = 3)

        graph_dict = NG.get_periodic_network(evolution_dict)

        graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
                evolution_dict; 
            print_progress = true,
            print_every_nr_attempted_bond_switches = 500)

        evolution_dict["total_energy_vec"] = total_energy_vec
        evolution_dict["move_accepted_vec"] = move_accepted_vec

        NG.plot_network(graph_dict)

        filename = "216_vertices_T_"*string(temperature)*"_heated_for_"*string(nr_heating_mc_steps)*"_steps_quenched"

        save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

        NG.save_graph_to_h5_and_MGformat(graph_dict,
            filename;
            evolution_dict = evolution_dict,
            save_path 
                = save_path)

    end

end


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

graph_dict_low_t = NG.load_graph_from_h5_and_MGformat(dict_path*"216_vertices_T_0.1_heated_for_0.5_steps_quenched")

graph_dict_high_t = NG.load_graph_from_h5_and_MGformat(dict_path*"216_vertices_T_1.0_heated_for_0.5_steps_quenched")

l_max = 12

q_l_low_t = NA.get_q_l_total_network_mean_dict(graph_dict_low_t, l_max)

q_l_high_t = NA.get_q_l_total_network_mean_dict(graph_dict_high_t, l_max)

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, q_l_low_t[i])
    println()
end

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, q_l_high_t[i])
    println()
end


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.1, 1.0, 2.0]


for nr_heating_mc_steps in [0.5]
    for temperature in temperatures

        evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 ,temperature_vec = [temperature, 0],
        nr_monte_carlo_steps_per_temperature_vec = [nr_heating_mc_steps, 50], min_ring_size = 3)

        graph_dict = NG.get_periodic_network(evolution_dict)

        graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
                evolution_dict; 
            print_progress = true,
            print_every_nr_attempted_bond_switches = 1000)

        evolution_dict["total_energy_vec"] = total_energy_vec
        evolution_dict["move_accepted_vec"] = move_accepted_vec

        filename = "1000_vertices_T_"*string(temperature)*"_heated_for_"*string(nr_heating_mc_steps)*"_steps_quenched"

        save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

        NG.save_graph_to_h5_and_MGformat(graph_dict,
            filename;
            evolution_dict = evolution_dict,
            save_path 
                = save_path)

    end

end


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.1, 1.0, 2.0]

for i in eachindex(temperatures)

    graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*
    "1000_vertices_T_"
    *string(temperatures[i])*"_heated_for_0.5_steps_quenched")

    structure_factor_dict = NA.get_structure_factor_bartlett_isotrope_by_wavenumber_vec(
        graph_dict;
        sampling_distance_step_length = 0.05,
        maximal_sampling_distance = 4*graph_dict["supercell_edge_length"],
        save_result = true,
        save_path = raw"..\analysis_data\random_networks\1000_vertices_T_"
        *string(temperatures[i])*"_heated_for_0.5_steps_quenched",
        print_progress = true,
        label = "T = "*string(temperatures[i]))

end


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

filename = "216_vertices_T_0.2_heated_for_1.0_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*filename)

NG.save_graph_to_h5_and_gml(graph_dict,
"my_graph")

dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\my_graph"

loaded_graph_dict = NG.load_graph_from_h5_and_gml(dict_path)

NG.plot_network(graph_dict)


directory_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\without_ring_size_limitation\\"

NG.convert_all_files_in_directory_MGformat_to_gml(directory_path)


# load some network
dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"
filename = "216_vertices_T_0.2_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)
evolution_dict = GU.load_h5_dict(dict_path*filename*"_evolution.h5")
evolution_dict["reject_during_relaxation_cycle_threshold"] = 5

# perform a bond switch
switched_chain = NG.get_random_chain(graph_dict; 
                                min_ring_size = evolution_dict["min_ring_size"], seed=15)

initial_cluster_dict = NG.get_cluster_in_shells_dict(
                                    graph_dict, 
                                    switched_chain; 
                                    shell_nr = evolution_dict["shell_nr"])

NG.switch_chain!(graph_dict,
    switched_chain )

# fully relax cluster

cluster_dict = NG.get_cluster_in_shells_dict(
                                    graph_dict, 
                                    switched_chain; 
                                    shell_nr = evolution_dict["shell_nr"])

temperature = 50
threshold_cluster_energy = initial_cluster_dict["cluster_energy"] - temperature * log(rand())

graph_dict, new_cluster_dict = NG.relax_cluster_keating!(graph_dict,
    cluster_dict, 
    evolution_dict;
    threshold_cluster_energy = threshold_cluster_energy,
    update_total_energy = false,
    print_progress = true)


    
evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 ,temperature_vec = [0.5],
    nr_monte_carlo_steps_per_temperature_vec = [1], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)
@time graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)


@ProfileView.profview NG.monte_carlo_move!(graph_dict, 
evolution_dict,
temperature;
print_progress = false)

@ProfileView.profview NG.monte_carlo_move!(graph_dict, 
evolution_dict,
temperature;
print_progress = false)

@time graph_dict, move_accepted = NG.monte_carlo_move!(graph_dict, 
evolution_dict,
temperature;
print_progress = false)

@time graph_dict, move_accepted = NG.monte_carlo_move!(graph_dict, 
evolution_dict,
temperature;
print_progress = false)

@time graph_dict, move_accepted = NG.monte_carlo_move!(graph_dict, 
evolution_dict,
temperature;
print_progress = false)


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_multiple_runs/test_networks_run_1/"
save_path = "../structures/random_networks/216_vertices_multiple_runs/test_networks_run_2/"

NG.generate_graphs_from_evolution_dicts_in_directory(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 100,
    print_progress = true)


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_multiple_runs/216_vertices_run_1/"
save_path = "../structures/random_networks/216_vertices_multiple_runs/216_vertices_run_4/"

println("Starting network generation")

NG.generate_graphs_from_evolution_dicts_in_directory(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 500,
    print_progress = true)



evolution_dicts_directory_path = "../structures/random_networks/216_vertices_multiple_runs/216_vertices_run_1/"

for i in 3:5

    save_path = "../structures/random_networks/216_vertices_multiple_runs/216_vertices_run_"*string(i)*"/"

    println("Starting run "*string(i))

    NG.generate_graphs_from_evolution_dicts_in_directory(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 500,
    print_progress = true)

end
    

dict_path = raw"..\analysis_data\random_networks\\"
filename = "216_vertices_T_0.2_heated_for_0.5_steps_quenched"

structure_factor_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_array.h5")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path= dict_path*filename)


dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\216_vertices_run_4\\"
filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

NG.plot_network(graph_dict)


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"
filename = "1000_vertices_T_0.1_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(graph_dict;
save_result = true,
save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=0.1, heated for 0.5 steps, quenched")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=0.1, heated for 0.5 steps, quenched")


filename = "1000_vertices_T_1.0_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(graph_dict;
save_result = true,
save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")


filename = "1000_vertices_perfect_diamond"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(graph_dict;
save_result = true,
save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")

dict_path = raw"..\analysis_data\random_networks\\"
filename = "1000_vertices_T_0.1_heated_for_0.5_steps_quenched"

structure_factor_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_array.h5")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
gaussian_filter_sigma_x = 2*pi/20, 
gaussian_filter_filtered_data_x_step_length = 2*pi/20,
label = "1000 vertices, T=0.1, heated for 0.5 steps, quenched")


filename = "1000_vertices_T_1.0_heated_for_0.5_steps_quenched"


structure_factor_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_array.h5")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
gaussian_filter_sigma_x = 2*pi/20, 
gaussian_filter_filtered_data_x_step_length = 2*pi/20,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")


filename = "1000_vertices_perfect_diamond"


structure_factor_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_array.h5")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
gaussian_filter_sigma_x = 2*pi/20, 
gaussian_filter_filtered_data_x_step_length = 2*pi/20,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_anneal_quench_multiple_runs\evolution_dicts\\"

randomization_temperature_vec = [2., 4., 8.]
randomization_nr_monte_carlo_steps_vec = [2., 4., 8.]
annealing_temperature_vec = [0.36, 0.5, 1.]

for randomization_temperature in randomization_temperature_vec
    for randomization_nr_monte_carlo_steps in randomization_nr_monte_carlo_steps_vec
        for annealing_temperature in annealing_temperature_vec

            temperature_vec = [randomization_temperature]
            nr_monte_carlo_steps_per_temperature_vec = [randomization_nr_monte_carlo_steps]

            for i in 1:5
                append!(temperature_vec, [annealing_temperature, 0])
                append!(nr_monte_carlo_steps_per_temperature_vec, [randomization_nr_monte_carlo_steps, 50])
            end

            evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

            filename = "216_vertices_randomization_T_"*string(randomization_temperature)*"_randomization_nr_MC_steps_"*string(randomization_nr_monte_carlo_steps)*"_annealing_T_"*string(annealing_temperature)*"_quenched_evolution.h5"

            GU.save_dict_to_h5(evolution_dict;
                        save_path=save_path*filename)

        end
    end
end


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_anneal_quench_multiple_runs/evolution_dicts/"

for i in 1:5

    save_path = "../structures/random_networks/216_vertices_anneal_quench_multiple_runs/run_"*string(i)*"/"

    println("Starting run "*string(i))

    NG.generate_graphs_from_evolution_dicts_in_directory(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 500,
    print_progress = true)

end


dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\216_vertices_run_4\\"
filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_dict = BDA.get_binary_data_from_spatial_network(graph_dict;
    bond_radius = 0.35,
    filename = filename,
    save_result=true)

save_path = raw"..\analysis_data\random_networks\\"*filename

# load dict
complete_autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction_complete.h5")

spectral_density_dict = BDA.get_spectral_density_by_wavevector_array_fft(structure_dict;
    save_complete_autocovariance_fct_direction_dict = false,
    save_result = true,
    save_path = save_path,
    complete_autocovariance_fct_direction_dict = complete_autocovariance_fct_direction_dict)


# loop through folders
for i in 1:5

    graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\run_"*string(i)*"\\"

    structure_dict_save_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\run_"*string(i)*"\\"

    # loop through all files in folder
    for filename in readdir(graph_dict_path)

        # check if file is a h5 file
        if endswith(filename, ".gml")

            # load graph dict
            graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename[1:end-4])

            # create and save voxelized data array
            structure_dict = BDA.get_binary_data_from_spatial_network(graph_dict;
            bond_radius = 0.35,
            voxel_edge_length = 0.2,
            save_path = structure_dict_save_path,
            filename = filename[1:end-4],
            save_result=true)
#
        end
    end

end


graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\216_vertices_run_4\\"
filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"

#graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename)
structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\run_4\\"
structure_dict = GU.load_h5_dict(structure_dict_path*filename*"_structure.h5")
save_path = raw"..\analysis_data\random_networks\\"*filename

# get autocovariance fct for the structure
autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction_pbc.h5")

spectral_density_dict = BDA.get_spectral_density_by_wavevector_array_fft_pbc(structure_dict;
    save_autocovariance_fct_direction_dict = false,
    save_result = true,
    save_path = save_path,
    autocovariance_fct_direction_dict = autocovariance_fct_direction_dict
    )


filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"

#graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename)
structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\run_4\\"
structure_dict = GU.load_h5_dict(structure_dict_path*filename*"_structure.h5")
save_path = raw"..\analysis_data\random_networks\\"*filename

# get autocovariance fct for the structure
autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction_pbc.h5")

# reorganize the autocovariance fct array by shifting all entries by half the length of the array
autocovariance_fct_array = cat(autocovariance_fct_direction_dict["autocovariance_fct_array"][Int(ceil(size(autocovariance_fct_direction_dict["autocovariance_fct_array"])[1]/2)):end,:,:], autocovariance_fct_direction_dict["autocovariance_fct_array"][1:Int(ceil(size(autocovariance_fct_direction_dict["autocovariance_fct_array"])[1]/2))-1,:,:], dims = 1)

autocovariance_fct_array = cat(autocovariance_fct_array[:,Int(ceil(size(autocovariance_fct_array)[2]/2)):end,:], autocovariance_fct_array[:,1:Int(ceil(size(autocovariance_fct_array)[2]/2))-1,:], dims = 2)

autocovariance_fct_array = cat(autocovariance_fct_array[:,:,Int(ceil(size(autocovariance_fct_array)[3]/2)):end], autocovariance_fct_array[:,:,1:Int(ceil(size(autocovariance_fct_array)[3]/2))-1], dims = 3)

autocovariance_fct_direction_dict["autocovariance_fct_array"] = autocovariance_fct_array

spectral_density_dict = BDA.get_spectral_density_by_wavevector_array_fft_pbc(structure_dict;
    save_autocovariance_fct_direction_dict = false,
    save_result = true,
    save_path = save_path,
    autocovariance_fct_direction_dict = autocovariance_fct_direction_dict
    )


dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\216_vertices_run_1\\"
filename = "216_vertices_T_2.0_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_dict = NA.get_binary_data_from_spatial_network(graph_dict;
    bond_radius = 0.35,
    voxel_edge_length = 1/2, 
    filename = filename,
    save_result=false)

@ProfileView.profview autocovariance_fct_direction_dict = NA.get_autocovariance_fct_by_sampling_indices_array(structure_dict;
save_result = false,
print_progress = true)

@ProfileView.profview autocovariance_fct_direction_dict = NA.get_autocovariance_fct_by_sampling_indices_array(structure_dict;
save_result = false,
print_progress = true)



for i in 1:5

    structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_multiple_runs/run_"*string(i)*"/"
    save_path = "../analysis_data/random_networks/216_vertices_multiple_runs/run_"*string(i)*"/"

    NA.get_autocovariance_fct_direction_from_filenames_multithreading(structure_dicts_path;
    print_progress = true,
    save_path = save_path)

end


print_lock = Threads.ReentrantLock()
graph_dicts_path = "../structures/random_networks/216_vertices_multiple_runs/"

structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_multiple_runs/"

autocovariance_dicts_path = "../analysis_data/random_networks/216_vertices_multiple_runs/"

NA.get_binary_data_from_spatial_network_multithreading(graph_dicts_path,
structure_dicts_path;
print_progress = true,
print_lock = print_lock,
nr_runs = 5,
bond_radius = 0.35,
voxel_edge_length = 0.1)



print_lock = Threads.ReentrantLock()

structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_multiple_runs/"

save_path = "../analysis_data/random_networks/216_vertices_multiple_runs/"

NA.get_autocovariance_fct_direction_from_filenames_multithreading(structure_dicts_path;
print_progress = true,
save_path = save_path,
print_lock = print_lock)


evolution_dicts_directory_path = raw"..\structures\random_networks\216_vertices_multiple_runs\\"

nr_monte_carlo_steps_for_quenching_vec = NA.get_monte_carlo_steps_for_quenching_vec(evolution_dicts_directory_path)


print_lock = Threads.ReentrantLock()
evolution_dicts_directory_path = "../structures/random_networks/216_vertices_anneal_quench_multiple_runs/evolution_dicts/"

i = 1

save_path = "../structures/random_networks/216_vertices_anneal_quench_multiple_runs/run_"*string(i)*"/"
println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 500,
print_progress = true,
print_lock = print_lock)



dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

for i in 1:5

    current_dict_path = dict_path*"run_"*string(i)*"\\"

    filenames = readdir(current_dict_path)
    
    for filename in filenames
        loaded_dict = GU.load_h5_dict(current_dict_path*filename)

        new_filename = filename[1:Int((length(filename)-13)/2)]*"_structure.h5"

        GU.save_dict_to_h5(loaded_dict, current_dict_path*new_filename)

    end

end


function get_spectral_densities(structure_dict_path, analysis_data_path)

    for i in 1:5

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames
            structure_dict = GU.load_h5_dict(current_structure_dict_path*structure_dict_filename)
    
            autocovariance_fct_direction_dict = GU.load_h5_dict(current_analysis_data_path*structure_dict_filename[1:end-13]*"_autocovariance_fct_direction.h5")
    
            spectral_density_dict = NA.get_spectral_density_by_wavevector_array_fft(structure_dict;
                save_autocovariance_fct_direction_dict = false,
                save_result = true,
                save_path = current_analysis_data_path*structure_dict_filename[1:end-13],
                autocovariance_fct_direction_dict = autocovariance_fct_direction_dict
                )
    
            spectral_density_angle_averaged_dict = NA.get_spectral_density_angle_averaged(spectral_density_dict;
                gaussian_filter = true,
                gaussian_filter_sigma_x = 2*pi/25, 
                gaussian_filter_filtered_data_x_step_length = 2*pi/25,
                save_result = true,
            save_path = current_analysis_data_path*structure_dict_filename[1:end-13])
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\\"

get_spectral_densities(structure_dict_path, analysis_data_path)


randomization_temperature = 8.0
randomization_nr_monte_carlo_steps = 2.0
annealing_temperature = 1.0

i=1

temperature_vec = [randomization_temperature]
nr_monte_carlo_steps_per_temperature_vec = [randomization_nr_monte_carlo_steps]

for i in 1:5
    append!(temperature_vec, [annealing_temperature, 0])
    append!(nr_monte_carlo_steps_per_temperature_vec, [randomization_nr_monte_carlo_steps, 50])
end

evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)
graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec


filename = "216_vertices_randomization_T_"*string(randomization_temperature)*"_randomization_nr_MC_steps_"*string(randomization_nr_monte_carlo_steps)*"_annealing_T_"*string(annealing_temperature)*"_quenched_evolution.h5"

save_path = "../structures/random_networks/216_vertices_anneal_quench_multiple_runs/run_"*string(i)*"/"

NG.save_graph_to_h5_and_MGformat(graph_dict,
            filename;
            evolution_dict = evolution_dict,
            save_path 
                = save_path)



filename_disorder = "216_vertices_T_2.0_heated_for_0.5_steps_quenched"
filename_order = "216_vertices_T_0.4_heated_for_0.25_steps_quenched"

filename_other = "216_vertices_T_0.25_heated_for_0.05_steps_quenched"

data_path_order = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"*filename_order
data_path_disorder = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"*filename_disorder

data_path_other = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"*filename_other

spectral_density_angle_averaged_dict = GU.load_h5_dict(data_path_other*"_spectral_density_angle_averaged.h5")

function two_gaussians(x, p)
    a2, b2, c2, a3, b3, c3 = p
    return a2 * exp.(-b2 .* (x .- c2).^2) .+ a3 * exp.(-b3 .* (x .- c3).^2)
end

function exp_gaussian(x, p)
    a1, b1, a2, b2, c2 = p
    return a1 * exp.(-b1 .* x) .+ a2 * exp.(-b2 .* (x .- c2).^2)
end

function exp_decay(x, p)
    a1, b1 = p
    return a1 * exp.(-b1 .* x)
end


# define function to estimate the fit parameters from data
function estimate_fit_parameters(x_vec, y_vec)

    # get first parameters by assuming that the exponential decay is the dominant feature
    # and linear close to x=0
    a1 = ( x_vec[2]*y_vec[1] - x_vec[1]*y_vec[2] )/(x_vec[2] - x_vec[1])
    b1 = 3*(y_vec[1] - y_vec[2])/( x_vec[2]*y_vec[1] - x_vec[1]*y_vec[2] )

    # locate peaks of spectral density
    pks, vals = Peaks.findmaxima( y_vec )

    # get parameters for two gaussians
    a2 = vals[1]
    b2 = 2.
    c2 = x_vec[pks[1]]
    a3 = vals[2]
    b3 = 2.
    c3 = x_vec[pks[2]]

    # if the exponential decay is not the dominant feature, get parameters for gaussians
    if a1 < 0
        return [a2, b2, c2, a3, b3, c3]
    else
        # if there is no dominant peak, use only exponential decay
        if a2 < 0.5
            return [a1, b1]
        else
            return [a1, b1, a2, b2, c2 ]
        end
    end

end

# get wavenumber vector and structure factor vector
x_vec = spectral_density_angle_averaged_dict["wavenumber_vec"]
y_vec = Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"])

# estimate fit parameters
p0 = estimate_fit_parameters(x_vec, y_vec)

# get function to fit
if length(p0) == 6
    fit_function = two_gaussians
elseif length(p0) == 5
    fit_function = exp_decay_gaussian
else
    fit_function = exp_decay
end

# fit function to data
fit_result = LsqFit.curve_fit(fit_function, x_vec, y_vec, p0)

# get fit parameters and their uncertainties
fit_param = Measurements.measurement.(fit_result.param, 
    sqrt.(LinearAlgebra.diag(LsqFit.estimate_covar(fit_result))) )

Plots.plot(x_vec, y_vec, label="data")
Plots.plot!(x_vec, fit_function(x_vec, p0), label="initial guess", ls=:dot)
Plots.plot!(x_vec, fit_function(x_vec, Measurements.value.(fit_param)), label="fit", ls=:dash)
Plots.ylims!(0,100)


graph_dict_path = raw"..\structure_analysis\structures\random_networks\216_vertices_anneal_quench_multiple_runs\run_1\\"


all_filenames = readdir(graph_dict_path)

filenames_filtered = filter(filename -> endswith(filename, ".gml"), all_filenames)

filenames = [filename[1:end-4] for filename in filenames_filtered]

for filename in filenames
    println(filename)

    # Use a regular expression to match all numbers (both integers and floating point numbers)
    regex = r"\d+\.?\d*"

    # Find all matches in the string
    matches = eachmatch(regex, filename)

    # Convert the matches to numbers (Float64)
    numbers = [parse(Float64, match.match) for match in matches]

    graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename)

    new_filename = "216_vertices_rand_T_"*string(numbers[2])*"_rand_nr_MC_steps_"*string(numbers[3])*"_anneal_T_"*string(numbers[4])*"_quenched"

    NG.save_graph_to_h5_and_gml(graph_dict,
    new_filename;
    evolution_dict = GU.load_h5_dict(graph_dict_path*filename*"_evolution.h5"),
    save_path  = graph_dict_path)
end


graph_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_anneal_quench_multiple_runs\run_1\\"

all_filenames = readdir(graph_dict_path)
filenames = filter(filename -> endswith(filename, "_evolution.h5"), all_filenames)
final_energy_vec = Float64[]

for filename in filenames
    println(filename)

    evolution_dict = GU.load_h5_dict(graph_dict_path*filename)

    push!(final_energy_vec, evolution_dict["total_energy_vec"][end])
end

filenames_sorted = filenames[sortperm(final_energy_vec)   ]
sort!(final_energy_vec)

graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filenames_sorted[26][1:end-13])
NG.plot_spatial_network(graph_dict)



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_anneal_quench_multiple_runs\evolution_dicts_2\\"

randomization_temperature_vec = [2., 4., 8.]
randomization_nr_monte_carlo_steps_vec = [2., 4., 8.]
annealing_temperature_vec = [0.36, 0.5, 1.]

for randomization_temperature in randomization_temperature_vec
    for randomization_nr_monte_carlo_steps in randomization_nr_monte_carlo_steps_vec
        for annealing_temperature in annealing_temperature_vec

            temperature_vec = [randomization_temperature]
            nr_monte_carlo_steps_per_temperature_vec = [randomization_nr_monte_carlo_steps]

            for i in 1:5
                append!(temperature_vec, [annealing_temperature, 0])
                append!(nr_monte_carlo_steps_per_temperature_vec, [6, 50])
            end

            evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

            filename = "216_vertices_rand_T_"*string(randomization_temperature)*"_rand_nr_MC_steps_"*string(randomization_nr_monte_carlo_steps)*"_anneal_T_"*string(annealing_temperature)*"_quenched_evolution.h5"

            GU.save_dict_to_h5(evolution_dict, save_path*filename)

        end
    end
end



print_lock = Threads.ReentrantLock()
evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts/"

i = 2

save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
print_lock = print_lock)


graph_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\run_2\\"

filename = "216_vertices_rand_T_8.0_rand_nr_MC_steps_4.0_anneal_T_0.36_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename)
NG.plot_spatial_network(graph_dict)

evolution_dict = GU.load_h5_dict(graph_dict_path*filename*"_evolution.h5")

Plots.plot(collect(1:length(evolution_dict["total_energy_vec"]))./(216*18), evolution_dict["total_energy_vec"]./216, xlabel="MC steps", ylabel="Energy", label="energy per vertex")


temperature_vec = [2,0,2,0]
nr_monte_carlo_steps_per_temperature_vec = [0.01, 50, 0.01, 50]

evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

graph_dict = NG.get_periodic_network(evolution_dict)
graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 50)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

Plots.plot(collect(1:length(evolution_dict["total_energy_vec"]))./(64*18), evolution_dict["total_energy_vec"]./64, xlabel="MC steps", ylabel="Energy", label="energy per vertex")



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_3\\"

randomization_temperature = 1.45

randomization_nr_monte_carlo_steps_vec = [4., 6., 8., 10.]
cooling_nr_monte_carlo_steps_vec = [0.5, 1., 2., 4.]

temperature_vec = vcat(collect(1.45:-0.1:0.35), [0])

for randomization_nr_monte_carlo_steps in randomization_nr_monte_carlo_steps_vec
    for cooling_nr_monte_carlo_steps in cooling_nr_monte_carlo_steps_vec
        
        nr_monte_carlo_steps_per_temperature_vec = vcat([randomization_nr_monte_carlo_steps], cooling_nr_monte_carlo_steps .* ones(11), [50] )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
        nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

        filename = "216_vertices_rand_nr_MC_steps_"*string(randomization_nr_monte_carlo_steps)*"_cool_nr_MC_steps_"*string(cooling_nr_monte_carlo_steps)*"_quenched_evolution.h5"

        GU.save_dict_to_h5(evolution_dict, save_path*filename)
    end
end

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_3\\"

randomization_temperature_vec = [1.45]
randomization_nr_monte_carlo_steps_vec = [4., 6., 8., 10.]
annealing_temperature_vec = [0.36]

for randomization_temperature in randomization_temperature_vec
    for randomization_nr_monte_carlo_steps in randomization_nr_monte_carlo_steps_vec
        for annealing_temperature in annealing_temperature_vec

            temperature_vec = [randomization_temperature]
            nr_monte_carlo_steps_per_temperature_vec = [randomization_nr_monte_carlo_steps]

            for i in 1:5
                append!(temperature_vec, [annealing_temperature, 0])
                append!(nr_monte_carlo_steps_per_temperature_vec, [6, 50])
            end

            evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

            filename = "216_vertices_rand_T_"*string(randomization_temperature)*"_rand_nr_MC_steps_"*string(randomization_nr_monte_carlo_steps)*"_anneal_T_"*string(annealing_temperature)*"_quenched_evolution.h5"

            GU.save_dict_to_h5(evolution_dict, save_path*filename)

        end
    end
end


print_lock = Threads.ReentrantLock()

i = 3

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_3/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"



println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
print_lock = print_lock)


save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_4\\"

evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [1.45],
nr_monte_carlo_steps_per_temperature_vec = [10], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(10.0)*"_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [1.45, 0],
nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(10.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = collect(1.45:-0.05:0)
nr_monte_carlo_steps_per_temperature_vec = vcat([10], 0.5 .* ones(length(temperature_vec)-2), [50])

evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(10.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = collect(1.45:-0.05:0)
nr_monte_carlo_steps_per_temperature_vec = vcat([10], 1.0 .* ones(length(temperature_vec)-2), [50])

evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(10.0)*"_cool_nr_MC_steps_"*string(1.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


print_lock = Threads.ReentrantLock()

i = 4

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_4/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"


println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
print_lock = print_lock)



graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_multiple_runs\run_4\\"

evolution_dicts_directory_path = raw"..\structures\random_networks\216_vertices_multiple_runs\evolution_dicts_check\\"

filenames = ["216_vertices_T_0.1_cool_0.1_per_mc_quenched", "216_vertices_T_0.5_heated_for_1.0_steps_quenched", "216_vertices_T_0.25_heated_for_5.0_steps_quenched"]

for filename in filenames

    # load evolution dict
    evolution_dict = GU.load_h5_dict(graph_path * filename * "_evolution.h5")

    # add missing keys to the evolution dict
    evolution_dict["relax_globally_after_threshold_cycle"] = true
    evolution_dict["mean_nr_monte_carlo_steps_for_quenching"] = 13.7

    # save the evolution dict to new folder
    GU.save_dict_to_h5(evolution_dict, evolution_dicts_directory_path * filename * "_evolution.h5")

end


print_lock = Threads.ReentrantLock()


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_multiple_runs/evolution_dicts_check/"
save_path = "../structures/random_networks/216_vertices_multiple_runs/run_check/"


NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
print_lock = print_lock)


graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_multiple_runs\run_2\\"

filenames = ["216_vertices_T_0.1_heated_for_0.5_steps_quenched",
"216_vertices_T_0.2_heated_for_0.5_steps_quenched",
"216_vertices_T_0.4_heated_for_0.5_steps_quenched",
"216_vertices_T_0.5_heated_for_0.5_steps_quenched",
"216_vertices_T_1.0_heated_for_0.5_steps_quenched",
"216_vertices_T_2.0_heated_for_0.5_steps_quenched",
"216_vertices_T_0.1_heated_for_1.0_steps_quenched",
"216_vertices_T_0.1_heated_for_10.0_steps_quenched",
"216_vertices_T_0.4_heated_for_0.1_steps_quenched",
"216_vertices_T_0.4_heated_for_0.25_steps_quenched",
"216_vertices_T_0.5_heated_for_10.0_steps_quenched",
]


evolution_dicts_directory_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\1728_vertices\evolution_dicts\\"

for filename in filenames

    # load evolution dict
    evolution_dict = GU.load_h5_dict(graph_path * filename * "_evolution.h5")

    # add missing keys to the evolution dict
    evolution_dict["relax_globally_after_threshold_cycle"] = true
    evolution_dict["mean_nr_monte_carlo_steps_for_quenching"] = 13.7
    evolution_dict["nr_vertices"] = 1728

    # save the evolution dict to new folder
    GU.save_dict_to_h5(evolution_dict, evolution_dicts_directory_path * "1728" *filename[4:end] * "_evolution.h5")

end


evolution_dicts_directory_path = raw"../structures/random_networks/1728_vertices/evolution_dicts/"

print_lock = Threads.ReentrantLock()

save_path = "../structures/random_networks/1728_vertices/run_cubic/"


NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = true,
print_lock = print_lock)


function get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)

    for i in 1:5

        current_graph_dict_path = graph_dict_path*"run_"*string(i)*"\\"

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames

            structure_dict = GU.load_h5_dict(current_structure_dict_path*structure_dict_filename)
    
            autocovariance_fct_direction_dict = GU.load_h5_dict(current_analysis_data_path*structure_dict_filename[1:end-13]*"_autocovariance_fct_direction.h5")
    
            volume_fract_variance_dict = NA.get_volume_fract_variance(autocovariance_fct_direction_dict;
            save_result = true,
            save_path = current_analysis_data_path*structure_dict_filename[1:end-13])
            
            graph_dict = NG.load_graph_from_h5_and_gml(current_graph_dict_path*structure_dict_filename[1:end-13])

            structure_factor_dict = NA.get_structure_factor_by_wavevector_array(graph_dict;
                save_result = true,
                save_path = current_analysis_data_path*structure_dict_filename[1:end-13])

            structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict;
                gaussian_filter = true,
                gaussian_filter_sigma_x = 2*pi/25, 
                gaussian_filter_filtered_data_x_step_length = 2*pi/25,
                save_result = true,
                save_path = current_analysis_data_path*structure_dict_filename[1:end-13])

            #local_nr_variance_dict = NA.get_local_nr_variance_by_window_radius_vec(
            #    graph_dict;
            #    structure_factor_dict = structure_factor_angle_averaged_dict,
            #    window_radius_step_length = 0.2,
            #    save_result = true,
            #    save_path = current_analysis_data_path*structure_dict_filename[1:end-13])
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\\"

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\\"

get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)


graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_multiple_runs\run_2\\"

data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"

filename = "216_vertices_T_0.2_heated_for_0.05_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(graph_path*filename)

structure_factor_angle_averaged_dict = GU.load_h5_dict(data_path*filename*"_structure_factor_angle_averaged.h5")

Plots.plot(structure_factor_angle_averaged_dict["wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict["structure_factor_vec"]) , xlims = (0,10))


graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_multiple_runs\run_2\\"

data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"

filename = "216_vertices_T_0.1_heated_for_5.0_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(graph_path*filename)

spectral_density_angle_averaged_dict = GU.load_h5_dict(data_path*filename*"_spectral_density_angle_averaged.h5")

anisotropy_metric_from_spectral_density = NA.get_anisotropy_metric_from_spectral_density(spectral_density_angle_averaged_dict)


function get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)

    for i in 1:5

        current_graph_dict_path = graph_dict_path*"run_"*string(i)*"\\"

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames
            
            graph_dict = NG.load_graph_from_h5_and_gml(current_graph_dict_path*structure_dict_filename[1:end-13])

            correlation_functions_dict = NA.get_correlation_functions(graph_dict;
                distance_histogram_bin_width = 0.02,
                save_result = true,
                save_path = current_analysis_data_path*structure_dict_filename[1:end-13])
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\\"

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\\"

get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)



graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\run_2\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"

filename = "216_vertices_T_0.2_heated_for_0.01_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename)

small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(filename,
    graph_path,
    analysis_data_path;
    save_result = false,
    )


function get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)

    for i in 1:5

        current_graph_dict_path = graph_dict_path*"run_"*string(i)*"\\"

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames

            small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(structure_dict_filename[1:end-13],
            current_graph_dict_path,
            current_analysis_data_path,
            save_result = true)
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\\"

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\\"

get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)



analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"

all_filenames = readdir(analysis_data_path)

# get filenames of small scale order anisotropy_metric_from_spectral_density
order_metrics_filenames = [filename for filename in all_filenames if occursin("small_scale_order_metrics", filename)]

# initialize vectors for all order metrics
total_keating_energy_vec = Vector{Float64}(undef, length(order_metrics_filenames))
bond_length_std_vec = Vector{Float64}(undef, length(order_metrics_filenames))
bond_angle_std_vec = Vector{Float64}(undef, length(order_metrics_filenames))
dihedral_angle_std_vec = Vector{Float64}(undef, length(order_metrics_filenames))
q_l_vec_vec = Vector{Vector{Measurements.Measurement{Float64}}}(undef, length(order_metrics_filenames))
cluster_metric_vec = Vector{Float64}(undef, length(order_metrics_filenames))
anisotropy_metric_from_structure_factor_vec = Vector{Float64}(undef, length(order_metrics_filenames))
anisotropy_metric_from_spectral_density_vec = Vector{Float64}(undef, length(order_metrics_filenames))


# loop through order metric filenames
for i in eachindex(order_metrics_filenames)

    # load order metrics
    order_metrics_dict = GU.load_h5_dict(analysis_data_path*order_metrics_filenames[i])

    # get all order metrics and save them to the corresponding vectors
    total_keating_energy_vec[i] = order_metrics_dict["total_keating_energy"]
    bond_length_std_vec[i] = order_metrics_dict["bond_length_std"]
    bond_angle_std_vec[i] = order_metrics_dict["bond_angle_std"]
    dihedral_angle_std_vec[i] = order_metrics_dict["dihedral_angle_std"]
    q_l_vec_vec[i] = order_metrics_dict["q_l_vec"]
    cluster_metric_vec[i] = order_metrics_dict["cluster_metric"]
    anisotropy_metric_from_structure_factor_vec[i] = order_metrics_dict["anisotropy_metric_from_structure_factor"]
    anisotropy_metric_from_spectral_density_vec[i] = order_metrics_dict["anisotropy_metric_from_spectral_density"]

end

# sort all vectors with respect to the keating energy
sort!(total_keating_energy_vec)

order_metrics_filenames = order_metrics_filenames[sortperm(total_keating_energy_vec)]

bond_length_std_vec = bond_length_std_vec[sortperm(total_keating_energy_vec)]
bond_angle_std_vec = bond_angle_std_vec[sortperm(total_keating_energy_vec)]
dihedral_angle_std_vec = dihedral_angle_std_vec[sortperm(total_keating_energy_vec)]
q_l_vec_vec = q_l_vec_vec[sortperm(total_keating_energy_vec)]
cluster_metric_vec = cluster_metric_vec[sortperm(total_keating_energy_vec)]
anisotropy_metric_from_structure_factor_vec = anisotropy_metric_from_structure_factor_vec[sortperm(total_keating_energy_vec)]
anisotropy_metric_from_spectral_density_vec = anisotropy_metric_from_spectral_density_vec[sortperm(total_keating_energy_vec)]

Plots.scatter(bond_angle_std_vec[1:end-2], bond_length_std_vec[1:end-2])
Plots.scatter(anisotropy_metric_from_structure_factor_vec, anisotropy_metric_from_spectral_density_vec)
Plots.scatter(total_keating_energy_vec, anisotropy_metric_from_spectral_density_vec)

Plots.scatter(total_keating_energy_vec[1:end-2], cluster_metric_vec[1:end-2])


filename = order_metrics_filenames[80][1:end-29]
graph_dict = NG.load_graph_from_h5_and_gml(graph_path*filename)

evolution_dict = GU.load_h5_dict(graph_path*filename*"_evolution.h5")
evolution_dict["relax_globally_after_threshold_cycle"] = true
evolution_dict["reject_during_relaxation_cycle_threshold"]  = 5
random_chain = NG.get_random_chain(graph_dict)

total_energy = NG.get_total_energy_keating(graph_dict)

original_graph_dict = deepcopy(graph_dict)

graph_dict = NG.relax_network_keating!(graph_dict,
random_chain,
evolution_dict;
threshold_total_energy = Inf,
update_total_energy = true,
print_progress = true)

println(graph_dict["total_energy"])


function relax_all_networks_globally(graph_path)

    for i in 1:5
        # get current path
        current_path = graph_path * "run_" * string(i) * "\\"

        # get all files in directory
        filenames = readdir(current_path)

        # get filenames of all evolution dicts
        filenames_evolution_dicts = filter(filename -> endswith(filename, "_evolution.h5"), filenames)

        for filename in filenames_evolution_dicts

            println(filename)

            # load evolution dict
            evolution_dict = GU.load_h5_dict(current_path * filename)

            # load graph
            graph_dict = NG.load_graph_from_h5_and_gml(current_path * filename[1:end-13])

            evolution_dict["relax_globally_after_threshold_cycle"] = true
            evolution_dict["reject_during_relaxation_cycle_threshold"]  = 5
            evolution_dict["mean_nr_monte_carlo_steps_for_quenching"]  = 13.7

            # relax network
            graph_dict = NG.relax_network_keating!(graph_dict,
            NG.get_random_chain(graph_dict),
            evolution_dict;
            threshold_total_energy = Inf,
            update_total_energy = true,
            print_progress = false)

            # save graph
            NG.save_graph_to_h5_and_gml(graph_dict, filename[1:end-13]; evolution_dict = evolution_dict,
            save_path = current_path)
        end

    end

    return
end

graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_globally_relaxed\\"

relax_all_networks_globally(graph_path)



graph_dicts_path = "../structures/random_networks/216_vertices_globally_relaxed/"
structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_globally_relaxed/"
analysis_data_path = "../analysis_data/random_networks/216_vertices_globally_relaxed/"

print_lock = Threads.ReentrantLock()

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
structure_dicts_path,
analysis_data_path,
print_progress = true,
runs_vec = [2],
print_lock = print_lock)


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_2\\"

order_metrics_dict = NA.get_small_length_scale_order_metrics_all_files(analysis_data_path;
    l_max_steinhardt_q_l = 12,
    save_result = true,)


for i in [1,3,4,5]
    analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_"*string(i)*"\\"

    order_metrics_dict = NA.get_small_length_scale_order_metrics_all_files(analysis_data_path;
        l_max_steinhardt_q_l = 12,
        save_result = true,)

    println("run $(i) done")

end



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_5\\"

temperature_vec = vcat(0.07 .* ones(20), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(20), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(80.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.07 .* ones(20), collect(0.07:-0.01:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(20), ones(length(collect(0.07:-0.01:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [0.07, 0],
nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(80.0)*"_cool_nr_MC_steps_"*string(1.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.07 .* ones(20), collect(0.07:-0.01:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(20), 0.5 .* ones(length(collect(0.07:-0.01:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [0.07, 0],
nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(80.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


temperature_vec = vcat(0.07 .* ones(10), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(10), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(40.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.07 .* ones(10), collect(0.07:-0.01:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(10), ones(length(collect(0.07:-0.01:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [0.07, 0],
nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(40.0)*"_cool_nr_MC_steps_"*string(1.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.07 .* ones(10), collect(0.07:-0.01:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(10), 0.5 .* ones(length(collect(0.07:-0.01:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [0.07, 0],
nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(40.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


print_lock = Threads.ReentrantLock()

i = 5

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_"*string(i)*"/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = true,
print_lock = print_lock)


graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_multiple_runs\run_2\\"

filenames = ["216_vertices_T_0.125_heat_cool_0.1_per_mc_quenched",
"216_vertices_T_0.2_heat_cool_0.1_per_mc_quenched",
"216_vertices_T_0.3_heat_cool_0.1_per_mc_quenched",
"216_vertices_T_0.1_cool_0.1_per_mc_quenched",
"216_vertices_T_0.15_cool_0.1_per_mc_quenched",
"216_vertices_T_0.2_cool_0.1_per_mc_quenched",
]

evolution_dicts_directory_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\1728_vertices\evolution_dicts_2\\"

for filename in filenames

    # load evolution dict
    evolution_dict = GU.load_h5_dict(graph_path * filename * "_evolution.h5")

    # add missing keys to the evolution dict
    evolution_dict["relax_globally_after_threshold_cycle"] = true
    evolution_dict["mean_nr_monte_carlo_steps_for_quenching"] = 13.7
    evolution_dict["reject_during_relaxation_cycle_threshold"] = 5
    evolution_dict["nr_vertices"] = 1728

    # save the evolution dict to new folder
    GU.save_dict_to_h5(evolution_dict, evolution_dicts_directory_path * "1728" *filename[4:end] * "_evolution.h5")

end


save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_6\\"


temperature_vec = vcat(0.065 .* ones(10), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(10), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_T_"*string(0.065)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.062 .* ones(10), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(10), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_T_"*string(0.062)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


# 2 threads
print_lock = Threads.ReentrantLock()

i = 6

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_"*string(i)*"/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = true,
print_lock = print_lock)


# 6 threads
evolution_dicts_directory_path = raw"../structures/random_networks/1728_vertices/evolution_dicts_2/"

print_lock = Threads.ReentrantLock()

save_path = "../structures/random_networks/1728_vertices/run_2/"


NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = true,
print_lock = print_lock)


analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_globally_relaxed\run_2\\"

all_filenames = readdir(analysis_data_path)

# get filenames of small scale order anisotropy_metric_from_spectral_density
order_metrics_filenames = [filename for filename in all_filenames if occursin("small_scale_order_metrics", filename)]

# initialize vectors for all order metrics
anisotropy_metric_from_structure_factor_vec = Vector{Float64}(undef, length(order_metrics_filenames))
anisotropy_metric_from_spectral_density_vec = Vector{Float64}(undef, length(order_metrics_filenames))


# loop through order metric filenames
for i in eachindex(order_metrics_filenames)

    # load order metrics
    order_metrics_dict = GU.load_h5_dict(analysis_data_path*order_metrics_filenames[i])

    # get all order metrics and save them to the corresponding vectors
    total_keating_energy_vec[i] = order_metrics_dict["total_keating_energy"]
    anisotropy_metric_from_structure_factor_vec[i] = order_metrics_dict["anisotropy_metric_from_structure_factor"]
    anisotropy_metric_from_spectral_density_vec[i] = order_metrics_dict["anisotropy_metric_from_spectral_density"]

end

# sort all vectors with respect to the keating energy
sort!(total_keating_energy_vec)

anisotropy_metric_from_structure_factor_vec = anisotropy_metric_from_structure_factor_vec[sortperm(total_keating_energy_vec)]
anisotropy_metric_from_spectral_density_vec = anisotropy_metric_from_spectral_density_vec[sortperm(total_keating_energy_vec)]

Plots.scatter(bond_angle_std_vec[1:end-2], bond_length_std_vec[1:end-2])
Plots.scatter(anisotropy_metric_from_structure_factor_vec, anisotropy_metric_from_spectral_density_vec)
Plots.scatter(total_keating_energy_vec, anisotropy_metric_from_spectral_density_vec)

Plots.scatter(total_keating_energy_vec[1:end-2], cluster_metric_vec[1:end-2])



filename_diamond = "216_vertices_T_0.1_heated_for_0.01_steps_quenched"

graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_globally_relaxed\run_2\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_2\\"

graph_dict = NG.load_graph_from_h5_and_gml(graph_path*filename_diamond) 
NG.plot_spatial_network(graph_dict)

structure_factor_angle_averaged_dict_diamond = GU.load_h5_dict(analysis_data_path*filename_diamond*"_structure_factor_angle_averaged.h5")

spectral_density_angle_averaged_dict_diamond = GU.load_h5_dict(analysis_data_path*filename_diamond*"_spectral_density_angle_averaged.h5")

wavenumbers_to_check_vec = 2*pi*collect(0.2:0.1:1.0)

# get the wavenumbers that lie the clostest to the wavenumbers_to_check_vec
index_vec = [argmin(abs.(spectral_density_angle_averaged_dict_diamond["wavenumber_vec"] .- wavenumbers_to_check_vec[i])) for i in eachindex(wavenumbers_to_check_vec)]

summed_spectral_density_to_check = sum( spectral_density_angle_averaged_dict_diamond["spectral_density_vec"][index_vec] )

std_value_ratio = (Measurements.uncertainty(summed_spectral_density_to_check) / Measurements.value(summed_spectral_density_to_check))



function get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)

    for i in 1:5

        current_graph_dict_path = graph_dict_path*"run_"*string(i)*"\\"

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames
            small_scale_order_metrics_dict = GU.load_h5_dict(current_analysis_data_path*structure_dict_filename[1:end-13]*"_small_scale_order_metrics.h5")

            structure_factor_angle_averaged_dict = GU.load_h5_dict(current_analysis_data_path*structure_dict_filename[1:end-13]*"_structure_factor_angle_averaged.h5")

            anisotropy_metric_from_structure_factor = NA.get_anisotropy_metric_from_structure_factor(
    structure_factor_angle_averaged_dict)

            spectral_density_angle_averaged_dict = GU.load_h5_dict(current_analysis_data_path*structure_dict_filename[1:end-13]*"_spectral_density_angle_averaged.h5")

            anisotropy_metric_from_spectral_density = NA.get_anisotropy_metric_from_spectral_density(
                spectral_density_angle_averaged_dict)

            small_scale_order_metrics_dict["anisotropy_metric_from_structure_factor"] = anisotropy_metric_from_structure_factor
            small_scale_order_metrics_dict["anisotropy_metric_from_spectral_density"] = anisotropy_metric_from_spectral_density

            GU.save_dict_to_h5(small_scale_order_metrics_dict, current_analysis_data_path*structure_dict_filename[1:end-13]*"_small_scale_order_metrics.h5")
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

graph_dict_path = raw"..\structures\random_networks\216_vertices_globally_relaxed\\"

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_globally_relaxed\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_globally_relaxed\\"

get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)


function my_function()

    for i in 1:5
        analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_"*string(i)*"\\"

        order_metrics_dict = NA.get_small_length_scale_order_metrics_all_files(analysis_data_path;
            l_max_steinhardt_q_l = 12,
            save_result = true,)

        println("run $(i) done")

    end

    return
end

my_function()


# 6 threads
evolution_dicts_directory_path = raw"../structures/random_networks/1728_vertices/evolution_dicts_2/"

print_lock = Threads.ReentrantLock()

save_path = "../structures/random_networks/1728_vertices/run_2/"


NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = true,
further_evolve_previous_networks = true,
print_lock = print_lock)


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\\"

order_metrics_dict = Dict()

# loop through folders and append all order metrics to the order_metrics_dict
for i in 1:5
    
    current_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_"*string(i)*"\\"

    current_order_metrics_dict = GU.load_h5_dict(current_analysis_data_path*"all_order_metrics.h5")

    for (key, value) in current_order_metrics_dict
        if haskey(order_metrics_dict, key)
            order_metrics_dict[key] = vcat(order_metrics_dict[key], value)
        else
            order_metrics_dict[key] = value
        end
    end
    
end

order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec", "cluster_metric_vec"]

# sort all vectors in order of the total keating energy
for order_metric_name in order_metrics_names
    order_metrics_dict[order_metric_name] = order_metrics_dict[order_metric_name][sortperm(order_metrics_dict["total_keating_energy_vec"])]
end
order_metrics_dict["filenames_vec"] = order_metrics_dict["filenames_vec"][sortperm(order_metrics_dict["total_keating_energy_vec"])]
sort!(order_metrics_dict["total_keating_energy_vec"])

mask_vec = [contains.(order_metrics_dict["filenames_vec"], "heated_for_0.1_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_0.25_steps"),
    contains.(order_metrics_dict["filenames_vec"], "heated_for_0.5_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_1.0_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_5.0_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_10.0_steps"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.025"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.05"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.1"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.2"),
( contains.(order_metrics_dict["filenames_vec"], "cool_0.1")
    .& .!(contains.(order_metrics_dict["filenames_vec"], "heat")) ),
]

filename_vec = ["heated_for_0.1_steps", "heated_for_0.25_steps", "heated_for_0.5_steps", "heated_for_1.0_steps", "heated_for_5.0_steps", "heated_for_10.0_steps", "heat_cool_0.025", "heat_cool_0.05", "heat_cool_0.1", "heat_cool_0.2", "cool_0.1"]

title_vec = ["heated for 0.1 steps", "heated for 0.25 steps", "heated for 0.5 steps", "heated for 1.0 steps", "heated for 5.0 steps", "heated for 10.0 steps", "heat and cool 0.025/MC step", "heat and cool 0.05/MC step", "heat and cool 0.1/MC step", "heat and cool 0.2/MC step", "cool 0.1/MC step"]

i=9

mask = mask_vec[i]
filtered_filenames_vec = order_metrics_dict["filenames_vec"][mask]
filtered_total_keating_energy_vec = order_metrics_dict["total_keating_energy_vec"][mask]
filtered_bond_length_std_vec = order_metrics_dict["bond_length_std_vec"][mask]
filtered_bond_angle_std_vec = order_metrics_dict["bond_angle_std_vec"][mask]
filtered_dihedral_angle_std_vec = order_metrics_dict["dihedral_angle_std_vec"][mask]
filtered_anisotropy_metric_from_structure_factor_vec = order_metrics_dict["anisotropy_metric_from_structure_factor_vec"][mask]
filtered_anisotropy_metric_from_spectral_density_vec = order_metrics_dict["anisotropy_metric_from_spectral_density_vec"][mask]
filtered_cluster_metric_vec = order_metrics_dict["cluster_metric_vec"][mask]

pattern = r"T_([0-9\.]+)"
extracted_numbers = [match(pattern, s).captures[1] for s in filtered_filenames_vec]
temperatures = parse.(Float64, extracted_numbers)
min_temp = minimum(temperatures)
max_temp = maximum(temperatures)
normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

graph_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_globally_relaxed\run_2\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_globally_relaxed\run_2_thick_bonds\\"

run_2_index_vec = []

for temperature in [0.1, 0.125, 0.15,  0.2, 0.25, 0.3, 0.4, 0.5]

    current_filename = "216_vertices_T_"*string(temperature)*"_heat_cool_0.1_per_mc_quenched"

    graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*current_filename)

    run_2_index = argmin(abs.(graph_dict["total_energy"] .- filtered_total_keating_energy_vec))
    push!(run_2_index_vec, run_2_index)

    NG.save_mesh_from_spatial_network(graph_dict, current_filename;
    bond_radius = 0.3131,
    save_path=save_path)
end

println("run_2")
println(run_2_index_vec)
filtered_total_keating_energy_vec[run_2_index_vec] ./216



graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\run_6\\"

filename = "216_vertices_rand_T_0.062_quenched_3"

graph_dict = NG.load_graph_from_h5_and_gml(graph_path * filename)
evolution_dict = GU.load_h5_dict(graph_path * filename * "_evolution.h5")
NG.plot_spatial_network(graph_dict)


save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_7\\"


temperature_vec = vcat(0.062 .* ones(2), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(2), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(8.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.062 .* ones(3), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(3), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(12.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.062 .* ones(4), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(4), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(16.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.062 .* ones(2), collect(0.06:-0.02:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(2), 0.5 .* ones(length(collect(0.06:-0.02:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(8.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat( 0.062 .* ones(3), collect(0.06:-0.02:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(3), 0.5 .* ones(length(collect(0.06:-0.02:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(12.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


temperature_vec = vcat(0.062 .* ones(4), collect(0.06:-0.02:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(4), 0.5 .* ones(length(collect(0.06:-0.02:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(16.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


i = 7

print_lock = Threads.ReentrantLock()

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_"*string(i)*"/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = false,
print_lock = print_lock)



dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.285\run_1\\"

save_path_1 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\evolution_dicts\\"

save_path_2 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.36\evolution_dicts\\"

# load all files in dict path that end with "evolution.h5"
# and save them to save_path

# get all files in dict_path
files = readdir(dict_path)

# filter files that end with "evolution.h5"
files = filter(x -> occursin("evolution.h5", x), files)

# load all files and save them to save_path
for file in files

    # load the file
    evolution_dict = GU.load_h5_dict(dict_path*file)

    evolution_dict["bond_bending_const"] = 0.21
    
    # save the file
    GU.save_dict_to_h5(evolution_dict, save_path_1*file)

    evolution_dict["bond_bending_const"] = 0.36

    # save the file
    GU.save_dict_to_h5(evolution_dict, save_path_2*file)
end



i = 5

print_lock = Threads.ReentrantLock()

evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.21/evolution_dicts/"

save_path = "../structures/random_networks/216_vertices_bond_bending_0.21/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = false,
print_lock = print_lock)



graph_dict_path = "../structures/random_networks/1728_vertices/run_2/"

structure_dict_path = "../structures/random_networks/binary_structures/1728_vertices/"

analysis_data_path = "../analysis_data/random_networks/1728_vertices/run_2/"

filename = "1728_vertices_T_0.2_heat_cool_0.1_per_mc_quenched"


NA.get_all_dicts_from_graph_single_file(filename,
    graph_dict_path,
    structure_dict_path,
    analysis_data_path;
    print_progress = true)



graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\1728_vertices\run_2\\"
analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1728_vertices\run_2\\"

filename = "1728_vertices_T_0.2_heat_cool_0.1_per_mc_quenched"

small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(filename,
    graph_path,
    analysis_data_path;
    save_result = true,
    )

filename = "1728_vertices_T_0.125_heat_cool_0.1_per_mc_quenched"

small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(filename,
        graph_path,
        analysis_data_path;
        save_result = true,
        )



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\1000_vertices_bond_bending_0.285\evolution_dicts\\"

temperatures = [0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24]

temperature_gradient = 0.1

for temperature in temperatures

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
    temperature_gradient = temperature_gradient, 
        nr_monte_carlo_steps_per_temperature = 0.01,
        quench = true )

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 ,temperature_vec = temperature_vec,
        nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

    filename = "1000_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

    GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

end


i = 1

print_lock = Threads.ReentrantLock()

evolution_dicts_directory_path = "../structures/random_networks/1000_vertices_bond_bending_0.285/evolution_dicts/"

save_path = "../structures/random_networks/1000_vertices_bond_bending_0.285/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = false,
print_lock = print_lock)



temperatures = [0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24]

temperature_gradient = 0.1

nr_vertices = 512

bond_bending_const_vec = [0.21, 0.28, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_gradient = temperature_gradient, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end


nr_vertices = 1000

bond_bending_const_vec = [0.21, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_gradient = temperature_gradient, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end



nr_vertices = 216

bond_bending_const_vec = [0.21, 0.285, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"_heat_cool\\evolution_dicts\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_gradient = temperature_gradient, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.36_heat_cool/evolution_dicts/"
save_path = "../structures/random_networks/216_vertices_bond_bending_0.36_heat_cool/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


evolution_dicts_directory_path = "../structures/random_networks/512_vertices_bond_bending_0.285/evolution_dicts/"
save_path = "../structures/random_networks/512_vertices_bond_bending_0.285/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:2),
    print_lock = Threads.ReentrantLock())



temperatures = [0.11, 0.13, 0.15, 0.17]

temperature_gradient = 0.1

nr_vertices = 512

bond_bending_const_vec = [0.21, 0.285, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts_2\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_gradient = temperature_gradient, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end


nr_vertices = 1000

bond_bending_const_vec = [0.21, 0.285, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts_2\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_gradient = temperature_gradient, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end



nr_vertices = 216

bond_bending_const_vec = [0.21, 0.285, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"_heat_cool\\evolution_dicts_2\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_gradient = temperature_gradient, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end



evolution_dicts_directory_path = "../structures/random_networks/512_vertices_bond_bending_0.285/evolution_dicts_2/"
save_path = "../structures/random_networks/512_vertices_bond_bending_0.285/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())
    

evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.285_heat_cool/evolution_dicts_2/"
save_path = "../structures/random_networks/216_vertices_bond_bending_0.285_heat_cool/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


evolution_dicts_directory_path = "../structures/random_networks/1000_vertices_bond_bending_0.21/evolution_dicts/"
save_path = "../structures/random_networks/1000_vertices_bond_bending_0.21/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.21_heat_cool/evolution_dicts/"
save_path = "../structures/random_networks/216_vertices_bond_bending_0.21_heat_cool/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


graph_dict_path = "../structures/random_networks/1000_vertices_bond_bending_0.285/run_1/"

structure_dict_path = "../structures/random_networks/binary_structures/1000_vertices_bond_bending_0.285/run_1/"

analysis_data_path = "../analysis_data/random_networks/1000_vertices_bond_bending_0.285/run_1/"

# loop through all files in folder
for filename_with_format in readdir(graph_dict_path)

    # check if file is a h5 file
    if endswith(filename_with_format, ".gml")

        filename = filename_with_format[1:end-4]

        NA.get_all_dicts_from_graph_single_file(filename,
            graph_dict_path,
            structure_dict_path,
            analysis_data_path;
            print_progress = true)

    end
end


network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285\run_2\216_vertices_T_0.15_heat_cool_0.1_per_mc_quenched_structure.h5"

structure_dict_network = GU.load_h5_dict(network_path)

# load pachy weevil data

pachy_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\pachy_blue_structure.h5"

structure_dict_pachy = GU.load_h5_dict(pachy_path)

# adjust pachy weevil data to use same keys as random network

structure_dict_pachy["data_binary"] = structure_dict_pachy["data_binary"][:,1:308,1:308]
structure_dict_pachy["size_data"] = (308,308,308)


structure_dict_pachy["volume_fract_tot"] = sum(structure_dict_pachy["data_binary"])/structure_dict_pachy["size_data"][1]^3

# voxel size = 10 nm
# mean bond length =~ 160 nm
# voxel edge length in units of bond length = 0.0625

structure_dict_pachy["voxel_edge_length"] = 0.0625

structure_dict_pachy["mean_edge_length_data"] = structure_dict_pachy["size_data"][1] * structure_dict_pachy["voxel_edge_length"]

# save pachy weevil data

GU.save_dict_to_h5(structure_dict_pachy, raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\pachy_blue_structure_adjusted.h5")


pachy_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\pachy_red_structure.h5"

structure_dict_pachy = GU.load_h5_dict(pachy_path)


structure_dict_pachy["data_binary"] = structure_dict_pachy["data_binary"][:,1:303,1:303]
structure_dict_pachy["size_data"] = (303,303,303)


structure_dict_pachy["volume_fract_tot"] = sum(structure_dict_pachy["data_binary"])/structure_dict_pachy["size_data"][1]^3

# voxel size = 9 nm
# mean bond length =~ 185 nm
# voxel edge length in units of bond length = 0.04865

structure_dict_pachy["voxel_edge_length"] = 0.04865

structure_dict_pachy["mean_edge_length_data"] = structure_dict_pachy["size_data"][1] * structure_dict_pachy["voxel_edge_length"]

# save pachy weevil data

GU.save_dict_to_h5(structure_dict_pachy, raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\pachy_red_structure_adjusted.h5")



structure_dict_path = "../structures/biological/"

analysis_data_path = "../analysis_data/biological/"

filename = "pachy_blue"

NA.get_all_dicts_from_voxelized_structure(filename,
    structure_dict_path,
    analysis_data_path;
    print_progress = true,
    print_lock = Threads.ReentrantLock())

filename = "pachy_red"

NA.get_all_dicts_from_voxelized_structure(filename,
    structure_dict_path,
    analysis_data_path;
    print_progress = true,
    print_lock = Threads.ReentrantLock())


network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285\run_2\216_vertices_T_0.15_heat_cool_0.1_per_mc_quenched_structure.h5"

structure_dict_network = GU.load_h5_dict(network_path)



evolution_dicts_directory_path = "../structures/random_networks/512_vertices_bond_bending_0.21/evolution_dicts/"
save_path = "../structures/random_networks/512_vertices_bond_bending_0.21/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:2),
    print_lock = Threads.ReentrantLock())


network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285\run_2\216_vertices_T_0.4_heat_cool_0.1_per_mc_quenched_structure.h5"

structure_dict_network = GU.load_h5_dict(network_path)

pore_size_distribution_dict = NA.get_pore_size_distribution(structure_dict_network)

pore_size_distribution_second_moment = NA.get_pore_size_distribution_second_moment(pore_size_distribution_dict)

function get_pore_size_distributions()
    file_count = 0

    structure_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285\run_"

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.285\run_"

    for i in 1:5

        current_structure_path = structure_path * string(i) * "\\"
        current_save_path = save_path * string(i) * "\\"

        # read all files in the current structure path
        structure_files = readdir(current_structure_path)

        for file in structure_files
            # load structure
            structure_dict = GU.load_h5_dict(current_structure_path * file)

            filename = file[1:end-13]

            # get pore size distribution
            pore_size_distribution_dict = NA.get_pore_size_distribution(structure_dict;
            save_result = true,
            save_path = current_save_path * filename)

            # get pore size distribution second moment
            pore_size_distribution_second_moment = NA.get_pore_size_distribution_second_moment(pore_size_distribution_dict)

            # load small scale order metrics dict 
            small_scale_order_metrics_dict = GU.load_h5_dict(current_save_path * filename * "_small_scale_order_metrics.h5")

            # save small scale order metrics dict with pore size distribution second moment
            small_scale_order_metrics_dict["pore_size_distribution_second_moment"] = pore_size_distribution_second_moment

            GU.save_dict_to_h5(small_scale_order_metrics_dict,
            current_save_path * filename*"_small_scale_order_metrics.h5")

            file_count += 1
            println("File ", file_count, " done.")
            
        end
        
    end

    return

end

get_pore_size_distributions()



# load pachy weevil data

pachy_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\pachy_red_autocovariance_fct_direction_complete_unadjusted.h5"

autocovariance_dict_pachy = GU.load_h5_dict(pachy_path)

# cut all arrays to cubic shape

previous_size = size(autocovariance_dict_pachy["autocovariance_fct_array"])[1:3]
central_indices = Int.( (previous_size .+ 1) ./ 2)
half_cubic_edge_length = Int( ( minimum(previous_size[1:3] ) - 1) /2)


autocovariance_dict_pachy["autocovariance_fct_array"] = autocovariance_dict_pachy["autocovariance_fct_array"][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length, central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length, central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length, :]

autocovariance_dict_pachy["autocovariance_fct_array_uncertainty"] = Measurements.uncertainty.(autocovariance_dict_pachy["autocovariance_fct_array"])
autocovariance_dict_pachy["autocovariance_fct_array"] = Measurements.value.(autocovariance_dict_pachy["autocovariance_fct_array"])

autocovariance_dict_pachy["sampling_index_array"] = autocovariance_dict_pachy["sampling_vec_array"][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length, central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length, central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length, :]

autocovariance_dict_pachy["sampling_index_vec_vec"] = [autocovariance_dict_pachy["sampling_distance_vec_vec"][1][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length], autocovariance_dict_pachy["sampling_distance_vec_vec"][2][central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length], autocovariance_dict_pachy["sampling_distance_vec_vec"][3][central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length]]

# voxel size = 9 nm
# mean bond length =~ 185 nm
# voxel edge length in units of bond length = 0.04865
autocovariance_dict_pachy["voxel_edge_length"] = 0.04865
autocovariance_dict_pachy["sampling_distance_vec_vec"] = autocovariance_dict_pachy["sampling_index_vec_vec"] .* autocovariance_dict_pachy["voxel_edge_length"]
autocovariance_dict_pachy["sampling_distance_array"] = autocovariance_dict_pachy["sampling_index_array"] .* autocovariance_dict_pachy["voxel_edge_length"]

# save pachy weevil data

GU.save_dict_to_h5(autocovariance_dict_pachy, raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\pachy_red_autocovariance_fct_direction.h5")



# load pachy weevil data

pachy_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\pachy_blue_autocovariance_fct_direction_complete_unadjusted.h5"

autocovariance_dict_pachy = GU.load_h5_dict(pachy_path)

# cut all arrays to cubic shape

previous_size = size(autocovariance_dict_pachy["autocovariance_fct_array"])[1:3]
central_indices = Int.( (previous_size .+ 1) ./ 2)
half_cubic_edge_length = Int( ( minimum(previous_size[1:3] ) - 1) /2)


autocovariance_dict_pachy["autocovariance_fct_array"] = autocovariance_dict_pachy["autocovariance_fct_array"][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length, central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length, central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length, :]

autocovariance_dict_pachy["autocovariance_fct_array_uncertainty"] = Measurements.uncertainty.(autocovariance_dict_pachy["autocovariance_fct_array"])
autocovariance_dict_pachy["autocovariance_fct_array"] = Measurements.value.(autocovariance_dict_pachy["autocovariance_fct_array"])

autocovariance_dict_pachy["sampling_index_array"] = autocovariance_dict_pachy["sampling_vec_array"][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length, central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length, central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length, :]

autocovariance_dict_pachy["sampling_index_vec_vec"] = [autocovariance_dict_pachy["sampling_distance_vec_vec"][1][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length], autocovariance_dict_pachy["sampling_distance_vec_vec"][2][central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length], autocovariance_dict_pachy["sampling_distance_vec_vec"][3][central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length]]

# voxel size = 10 nm
# mean bond length =~ 160 nm
# voxel edge length in units of bond length = 0.0625

autocovariance_dict_pachy["voxel_edge_length"] = 0.0625
autocovariance_dict_pachy["sampling_distance_vec_vec"] = autocovariance_dict_pachy["sampling_index_vec_vec"] .* autocovariance_dict_pachy["voxel_edge_length"]
autocovariance_dict_pachy["sampling_distance_array"] = autocovariance_dict_pachy["sampling_index_array"] .* autocovariance_dict_pachy["voxel_edge_length"]

# save pachy weevil data

GU.save_dict_to_h5(autocovariance_dict_pachy, raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\pachy_blue_autocovariance_fct_direction.h5")


function get_spectral_density_from_autocovariance_fct(filename)

    # set path to analysis data
    analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\\"

    # load structure dict
    pachy_structure_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\\"*filename*"_structure.h5"

    # load autocovariance function dict
    pachy_auto_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\\"*filename*"_autocovariance_fct_direction.h5"

    # load structure dict
    structure_dict = GU.load_h5_dict(pachy_structure_path)

    # load autocovariance function dict
    autocovariance_fct_direction_dict = GU.load_h5_dict(pachy_auto_path)

    # get spectral density by wavevector array
    spectral_density_dict = NA.get_spectral_density_by_wavevector_array_fft(structure_dict;
    save_autocovariance_fct_direction_dict = false,
    save_result = true,
    save_path = analysis_data_path*filename,
    autocovariance_fct_direction_dict = autocovariance_fct_direction_dict)

    # get angle averaged spectral density
    spectral_density_angle_averaged_dict = NA.get_spectral_density_angle_averaged(spectral_density_dict;
    gaussian_filter = true,
    gaussian_filter_sigma_x = 2*pi/25, 
    gaussian_filter_filtered_data_x_step_length = 2*pi/25,
    save_result = true,
    save_path = analysis_data_path*filename)

    return

end

get_spectral_density_from_autocovariance_fct("pachy_blue")

get_spectral_density_from_autocovariance_fct("pachy_red")



evolution_dicts_directory_path = "../structures/random_networks/512_vertices_bond_bending_0.21/evolution_dicts/"
save_path = "../structures/random_networks/512_vertices_bond_bending_0.21/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(3:5),
    print_lock = Threads.ReentrantLock())


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"

filename = "512_vertices_perfect_diamond"

# load the network
evolution_dict = NA.get_evolution_dict(;nr_vertices = 512)

graph_dict = NG.get_periodic_network(evolution_dict)

NG.save_graph_to_h5_and_gml(graph_dict, filename; 
            save_path = path)

graph_dict = NG.load_graph_from_h5_and_gml(path * filename)

NG.plot_spatial_network(graph_dict)

NG.save_mesh_from_spatial_network(graph_dict, filename; save_path = path, bond_radius = 0.3131)


graph_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"

structure_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\diamonds\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

filename = "216_vertices_perfect_diamond"

NA.get_all_dicts_from_graph_single_file(filename,
    graph_dict_path,
    structure_dict_path,
    analysis_data_path;
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    print_lock = Threads.ReentrantLock())



graph_dicts_path = "../structures/random_networks/1000_vertices_bond_bending_0.21/"
structure_dicts_path = "../structures/random_networks/binary_structures/1000_vertices_bond_bending_0.21/"
analysis_data_path = "../analysis_data/random_networks/1000_vertices_bond_bending_0.21/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
    structure_dicts_path,
    analysis_data_path::String;
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


graph_dicts_path = "../structures/random_networks/512_vertices_bond_bending_0.21/"
structure_dicts_path = "../structures/random_networks/binary_structures/512_vertices_bond_bending_0.21/"
analysis_data_path = "../analysis_data/random_networks/512_vertices_bond_bending_0.21/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
        structure_dicts_path,
        analysis_data_path::String;
        structure_factor_diamond_std_value_ratio = 1,
        spectral_density_diamond_std_value_ratio = 1,
        pore_size_distribution_nr_sampled_voxels = 20000,
        print_progress = true,
        runs_vec = collect(1:5),
        print_lock = Threads.ReentrantLock())


graph_dicts_path = "../structures/random_networks/216_vertices_bond_bending_0.21_heat_cool/"
structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_bond_bending_0.21_heat_cool/"
analysis_data_path = "../analysis_data/random_networks/216_vertices_bond_bending_0.21_heat_cool/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
                structure_dicts_path,
                analysis_data_path::String;
                structure_factor_diamond_std_value_ratio = 1,
                spectral_density_diamond_std_value_ratio = 1,
                pore_size_distribution_nr_sampled_voxels = 20000,
                print_progress = true,
                runs_vec = collect(1:5),
                print_lock = Threads.ReentrantLock())


graph_dicts_path = "../structures/random_networks/1000_vertices_bond_bending_0.285/"
structure_dicts_path = "../structures/random_networks/binary_structures/1000_vertices_bond_bending_0.285/"
analysis_data_path = "../analysis_data/random_networks/1000_vertices_bond_bending_0.285/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
    structure_dicts_path,
    analysis_data_path::String;
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


graph_dicts_path = "../structures/random_networks/512_vertices_bond_bending_0.285/"
structure_dicts_path = "../structures/random_networks/binary_structures/512_vertices_bond_bending_0.285/"
analysis_data_path = "../analysis_data/random_networks/512_vertices_bond_bending_0.285/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
        structure_dicts_path,
        analysis_data_path::String;
        structure_factor_diamond_std_value_ratio = 1,
        spectral_density_diamond_std_value_ratio = 1,
        pore_size_distribution_nr_sampled_voxels = 20000,
        print_progress = true,
        runs_vec = collect(1:5),
        print_lock = Threads.ReentrantLock())


graph_dicts_path = "../structures/random_networks/216_vertices_bond_bending_0.285_heat_cool/"
structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_bond_bending_0.285_heat_cool/"
analysis_data_path = "../analysis_data/random_networks/216_vertices_bond_bending_0.285_heat_cool/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
                structure_dicts_path,
                analysis_data_path::String;
                structure_factor_diamond_std_value_ratio = 1,
                spectral_density_diamond_std_value_ratio = 1,
                pore_size_distribution_nr_sampled_voxels = 20000,
                print_progress = true,
                runs_vec = collect(1:5),
                print_lock = Threads.ReentrantLock())


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_1\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\geometrical_models\216_vertices_bond_bending_0.21\run_1\\"

filename = "216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(path * filename)

#NG.plot_spatial_network(graph_dict)

NG.save_mesh_from_spatial_network(graph_dict, filename; save_path = save_path, bond_radius = 0.35,
vector_out_of_supercell_length = 1)


paths = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285_heat_cool\run_1\216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_structure.h5",
    raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\512_vertices_bond_bending_0.285\run_1\512_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_structure.h5",
    raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\1000_vertices_bond_bending_0.285\run_1\1000_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_structure.h5"]


for path in paths

    structure_dict = GU.load_h5_dict(path)

    NA.get_digital_sphere_mask_dict(structure_dict["size_data"];
        save_result = true,
        save_path = raw"..\analysis_data\random_networks\digital_sphere_masks\\")

    println("done")

end



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

files = ["216_vertices_perfect_diamond_small_scale_order_metrics.h5", "512_vertices_perfect_diamond_small_scale_order_metrics.h5", "1000_vertices_perfect_diamond_small_scale_order_metrics.h5"]

anisotropy_metric_from_structure_factor_vec = Vector{Float64}(undef, 3)
anisotropy_metric_from_spectral_density_vec = Vector{Float64}(undef, 3)

for i in eachindex(files)

    current_dict = GU.load_h5_dict(path * files[i])

    anisotropy_metric_from_structure_factor_vec[i] = current_dict["anisotropy_metric_from_structure_factor"]
    anisotropy_metric_from_spectral_density_vec[i] = current_dict["anisotropy_metric_from_spectral_density"]

    println("Anisotropy metric from structure factor for ", files[i], ": ", anisotropy_metric_from_structure_factor_vec[i])
    println("Anisotropy metric from spectral density for ", files[i], ": ", anisotropy_metric_from_spectral_density_vec[i])
end


analysis_data_paths = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.21_heat_cool\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.285_heat_cool\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.285\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.285\\"]

for analysis_data_path in analysis_data_paths
    for i in 1:5
        NA.get_small_length_scale_order_metrics_all_files(analysis_data_path*"run_$i\\";
        l_max_steinhardt_q_l = 12,
        save_result = true,)
    end
    println("done with $analysis_data_path")
end


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.36/evolution_dicts/"

save_path = "../structures/random_networks/216_vertices_bond_bending_0.36/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\geometrical_models\diamonds\\"

nr_vertices_vec = [216, 512, 1000]

bond_radius_vec = [0.26, 0.35]

for nr_vertices_vec in nr_vertices_vec
    filename = string(nr_vertices_vec, "_vertices_perfect_diamond")
    graph_dict = NG.load_graph_from_h5_and_gml(path*filename)

    for bond_radius in bond_radius_vec

        save_filename = string(filename, "_bond_radius_", bond_radius)
        
        NG.save_mesh_from_spatial_network(graph_dict, save_filename; save_path = save_path, bond_radius = bond_radius,
        vector_out_of_supercell_length = 1)
    end
end


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\1000_vertices_bond_bending_0.285\run_2\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\geometrical_models\1000_vertices_bond_bending_0.285\run_2\\"

temperatures = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))

for temperature in temperatures
    filename = string("1000_vertices_T_", temperature, "_heat_cool_0.1_per_mc_quenched")

    graph_dict = NG.load_graph_from_h5_and_gml(path*filename)

    save_filename = string(filename, "_br_0.26")
        
    NG.save_mesh_from_spatial_network(graph_dict, save_filename; save_path = save_path, bond_radius = 0.26,
        vector_out_of_supercell_length = 1)
end



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.285_heat_cool\run_2\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"


filename = "216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(path*filename)

NG.save_spatial_network_to_gml(graph_dict["spatial_network"], filename, save_path=save_path)


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"

load_path = raw"C:\\Users\\HemmannF\\OneDrive - Université de Fribourg\\structure_analysis\\structures\\random_networks\\216_vertices_bond_bending_0.285\\run_1\\216_vertices_T_0.25_heated_for_0.1_steps_quenched.gml"
path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"

spatial_network = NG.load_spatial_network_from_gml(load_path)

NG.save_spatial_network_to_gml(spatial_network,
    "216_vertices_T_0.25_heated_for_0.1_steps_quenched";
    save_path=path)

path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"
function myfunc()
    count = 0
    filepath_list = []
    # Iterate through all directories and subdirectories
    for (root, dirs, files) in walkdir(path)
        for file in files
            if endswith(file, ".gml")

                joined_path = joinpath(root, file)
                spatial_network = NG.load_spatial_network_from_gml(joined_path)
                NG.save_spatial_network_to_gml(spatial_network,
                file[1:end-4];
                    save_path=root*"\\")

                #println(joinpath(root, file))  # Print the full path to the .gml file
                push!(filepath_list, joinpath(root, file))

                #print every 50th file
                count += 1
                if count % 50 == 0
                    println(count, " ", joined_path)
                end
            end
        end
    end

    return
end

myfunc()


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.285\run_2\\"

filename = "216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched"

# load the network
graph_dict = NG.load_graph_from_h5_and_gml(path*filename)

my_graph_dict = MetaGraphsNext.MetaGraph(
    Graphs.Graph();  
    label_type=Int64, 
    vertex_data_type=Dict{String, Any},  
    edge_data_type=Dict{String, Any}, 
    graph_data=Dict{String, Any}("coordination_nr"   => 4)
)



function my_func()

    nr_vertices_vec = [216, 512, 1000]
    bond_bending_vec = [0.21, 0.285, 0.36]

    # loop through each number of vertices
    for nr_vertices in nr_vertices_vec

        # loop through each bond bending
        for bond_bending in bond_bending_vec

            path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"*string(nr_vertices, "_vertices_bond_bending_", bond_bending, "\\")

            for run in 1:5
                current_path = path*"run_"*string(run)*"\\"

                # get all files in the current path
                files = readdir(current_path)

                # get vector of all gml files
                gml_files = [file for file in files if endswith(file, ".gml")]

                # loop through each file that is a gml file
                for gml_file in gml_files
                    filename  = gml_file[1:end-4]

                    spatial_network = NG.load_spatial_network_from_h5_and_gml(current_path*filename)

                    # save spatial network
                    NG.save_spatial_network_to_gml(
                        spatial_network,
                        filename;
                        save_path = current_path)

                    # delete h5 file
                    rm(current_path*filename*".h5")
                end
            end
        end
    end
end

my_func()



nr_vertices = 512

spatial_network_path = "../structures/random_networks/"*string(nr_vertices)*"_vertices_bond_bending_0.36/"
structure_dict_path = "../structures/random_networks/binary_structures/"*string(nr_vertices)*"_vertices_bond_bending_0.36/"
analysis_data_path = "../analysis_data/random_networks/"*string(nr_vertices)*"_vertices_bond_bending_0.36/"

NA.get_all_dicts_from_networks_multithreading(
    spatial_network_path,
    structure_dict_path,
    analysis_data_path;
    bond_radius = 0.35,
    voxel_edge_length = 0.1,
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())



spatial_network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.36\run_1\\"

filename = "216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched"

structure_dict_path = raw"C:\Users\HemmannF\Downloads\\"

analysis_data_path = structure_dict_path

NA.get_all_dicts_from_network_single_file(
    filename,
    spatial_network_path,
    structure_dict_path,
    analysis_data_path;
    pore_size_distribution_nr_sampled_voxels = 1000,
    print_progress = true)



analysis_data_paths = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.36\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.36\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.36\\"]

for analysis_data_path in analysis_data_paths
    for i in 1:5
        NA.get_small_length_scale_order_metrics_all_files(analysis_data_path*"run_$i\\";
        l_max_steinhardt_q_l = 12,
        save_result = true,)
    end
    println("done with $analysis_data_path")
end

analysis_data_paths = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.285\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.285\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.285\\"]

for analysis_data_path in analysis_data_paths
    for i in 1:5
        NA.get_small_length_scale_order_metrics_all_files(analysis_data_path*"run_$i\\";
        l_max_steinhardt_q_l = 12,
        save_result = true,)
    end
    println("done with $analysis_data_path")
end


function my_func()

    nr_vertices_vec = [512, 1000]
    bond_bending_vec = [0.21, 0.285, 0.36]

    # loop through each number of vertices
    for nr_vertices in nr_vertices_vec

        # loop through each bond bending
        for bond_bending in bond_bending_vec

            path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"*string(nr_vertices, "_vertices_bond_bending_", bond_bending, "\\")

            analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\\"*string(nr_vertices, "_vertices_bond_bending_", bond_bending, "\\")

            for run in 1:5

                println("Starting run ", run, " with ", nr_vertices, " vertices and bond bending ", bond_bending)

                current_path = path*"run_"*string(run)*"\\"
                current_analysis_data_path = analysis_data_path*"run_"*string(run)*"\\"

                # get all files in the current path
                files = readdir(current_path)

                # get vector of all gml files
                gml_files = [file for file in files if endswith(file, ".gml")]

                # loop through each file that is a gml file
                for gml_file in gml_files
                    filename  = gml_file[1:end-4]

                    spatial_network = NG.load_spatial_network_from_gml(current_path*filename*".gml")

                    # get correlation functions
                    correlation_functions_dict = NA.get_correlation_functions(
                        spatial_network;
                        distance_histogram_bin_width = 0.02,
                        save_result = true,
                        save_path = current_analysis_data_path*filename,
                        label = nothing)

                    # get all order metrics that contain information about small length scales
                    small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(
                        filename,
                        current_path,
                        current_analysis_data_path;
                        l_max_steinhardt_q_l = 12,
                        structure_factor_diamond_std_value_ratio 
                            = 1,
                        spectral_density_diamond_std_value_ratio 
                            = 1,
                        save_result = true,
                        )
                end
            end
        end
    end
end

my_func()


current_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"
current_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

# get all files in the current path
files = readdir(current_path)

# get vector of all gml files
gml_files = [file for file in files if endswith(file, ".gml")]

# loop through each file that is a gml file
for gml_file in gml_files
    filename  = gml_file[1:end-4]
    spatial_network = NG.load_spatial_network_from_h5_and_gml(current_path*filename)

    NG.save_spatial_network_to_gml(
        spatial_network,
        filename;
        save_path = current_path)

    # delete h5 file
    rm(current_path*filename*".h5")
end


function my_func()
    
    current_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"
    current_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

    # get all files in the current path
    files = readdir(current_path)
    # get vector of all gml files
    gml_files = [file for file in files if endswith(file, ".gml")]
    
    # loop through each file that is a gml file
    for gml_file in gml_files
        filename  = gml_file[1:end-4]
        spatial_network = NG.load_spatial_network_from_gml(current_path*filename*".gml")

        # get correlation functions
        correlation_functions_dict = NA.get_correlation_functions(
            spatial_network;
            distance_histogram_bin_width = 0.02,
            save_result = true,
            save_path = current_analysis_data_path*filename,
            label = nothing)
        # get all order metrics that contain information about small length scales
        small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(
            filename,
            current_path,
            current_analysis_data_path;
            l_max_steinhardt_q_l = 12,
            structure_factor_diamond_std_value_ratio 
                = 1,
            spectral_density_diamond_std_value_ratio 
                = 1,
            save_result = true,
            )
    end
end

my_func()


nr_vertices = 512

spatial_network_path = "../structures/random_networks/"*string(nr_vertices)*"_vertices_bond_bending_0.21/"
structure_dict_path = "../structures/random_networks/binary_structures/"*string(nr_vertices)*"_vertices_bond_bending_0.21/"
analysis_data_path = "../analysis_data/random_networks/"*string(nr_vertices)*"_vertices_bond_bending_0.21/"

NA.get_all_dicts_from_networks_multithreading(
    spatial_network_path,
    structure_dict_path,
    analysis_data_path;
    bond_radius = 0.35,
    voxel_edge_length = 0.1,
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


function my_func()
    
    path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"

    nr_vertices_vec = [216, 512, 1000]
    bond_bending_vec = [0.21, 0.285, 0.36]

    for nr_vertices in nr_vertices_vec
        for bond_bending in bond_bending_vec
            for run in 1:5
                current_path = path*string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\run_"*string(run)*"\\"

                save_path = path*"networks_for_simulation\\"*string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\run_"*string(run)*"\\"

                # get all files in the current path
                files = readdir(current_path)
                # get vector of all gml files
                gml_files = [file for file in files if endswith(file, ".gml")]

                # loop through each file that is a gml file
                for gml_file in gml_files
                    filename  = gml_file[1:end-4]
                    spatial_network = NG.load_spatial_network_from_gml(current_path*filename*".gml")

                    spatial_network_for_simulation = NG.get_spatial_network_for_simulation!(
                        spatial_network;
                        vector_out_of_supercell_length = 1,
                        duplicate_bonds_close_to_supercell_edge = true,
                        save_result = true,
                        filename = filename,
                        save_path = save_path)
                end
            end

            println("finished run for ", nr_vertices, " vertices and bond bending ", bond_bending)
        end
    end

    current_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"
    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\networks_for_simulation\diamonds\\"

    # get all files in the current path
    files = readdir(current_path)
    # get vector of all gml files
    gml_files = [file for file in files if endswith(file, ".gml")]

    # loop through each file that is a gml file
    for gml_file in gml_files
        filename  = gml_file[1:end-4]
        spatial_network = NG.load_spatial_network_from_gml(current_path*filename*".gml")

        spatial_network_for_simulation = NG.get_spatial_network_for_simulation!(
            spatial_network;
            vector_out_of_supercell_length = 1,
            duplicate_bonds_close_to_supercell_edge = true,
            save_result = true,
            filename = filename,
            save_path = save_path)
    end

end

my_func()



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_1\\"

filename = "216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched"

spatial_network = NG.load_spatial_network_from_gml(path * filename * ".gml")

spatial_network_for_sim = NG.get_spatial_network_for_simulation!(
    spatial_network;
    vector_out_of_supercell_length = 1,
    duplicate_bonds_close_to_supercell_edge = true,
    save_result = false)

NG.plot_spatial_network(spatial_network_for_sim)


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.285\run_1\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\geometrical_models\216_vertices_bond_bending_0.285\run_1\\"

temperatures = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))

for temperature in temperatures
    filename = string("216_vertices_T_", temperature, "_heat_cool_0.1_per_mc_quenched")

    spatial_network = NG.load_spatial_network_from_gml(path*filename*".gml")

    save_filename = string(filename, "_br_0.26")
        
    NG.save_mesh_from_spatial_network(spatial_network, save_filename; save_path = save_path, bond_radius = 0.26,
        vector_out_of_supercell_length = 1)
end



evolution_dicts_directory_path = "../structures/random_networks_evolution/216_vertices_bond_bending_0.285/evolution_dicts/"

temperatures = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.2, 0.22, 0.24]

nr_vertices = 216
bond_bending = 0.285
temperature_increase_per_monte_carlo_step = 0.1
nr_monte_carlo_steps_per_temperature = 0.1

for temperature in temperatures
    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_increase_per_monte_carlo_step = temperature_increase_per_monte_carlo_step, 
        nr_monte_carlo_steps_per_temperature = nr_monte_carlo_steps_per_temperature,
        quench = false )

    # extend the temperature vec by 50 elements of 0 temperature and the corresponding number of monte carlo steps
    temperature_vec = vcat(temperature_vec, zeros(Int(50/nr_monte_carlo_steps_per_temperature)))
    nr_monte_carlo_steps_per_temperature_vec = vcat(nr_monte_carlo_steps_per_temperature_vec, ones(Int(50/nr_monte_carlo_steps_per_temperature)) .* nr_monte_carlo_steps_per_temperature)

    evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
        nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
        bond_bending_const = bond_bending, thermal_fluctuations = false)

    filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_increase_per_monte_carlo_step)*"_per_mc_quenched"
    GU.save_dict_to_h5(evolution_dict, evolution_dicts_directory_path*filename*"_evolution.h5")
end

save_path = "../structures/random_networks_evolution/216_vertices_bond_bending_0.285/"

NG.generate_graphs_from_evolution_dicts_in_directory(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = true,
    further_evolve_previous_networks = false,
    print_lock = Threads.ReentrantLock())


function get_order_metric_evolution()

    # set the path to the directory where the networks are stored
    load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks_th_fluct\216_vertices_bond_bending_0.285\\"

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks_th_fluct\216_vertices_bond_bending_0.285\\"

    # get all filenames in the load_path directory that contain the string "T_0.1_heat_cool_0.1_per_mc_quenched"
    temperature_vec = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.2, 0.22, 0.24]

    bond_angle_std_vec_vec = []
    bond_length_std_vec_vec = []

    for temperature in temperature_vec
        # get the filenames for the current temperature which contain the current temperature
        # and the term "gml"
        filenames = [f for f in readdir(load_path) if occursin("T_"*string(temperature)*"_heat_cool_0.1_per_mc_quenched_", f)]
        filenames = [f for f in filenames if endswith(f, "0.gml")]

        # get the integers in the end of all filenames and sort the filenames accordingly
        for f in filenames
            if match(r"quenched_(\d+)\.gml$", f) == nothing
                println("Error: no match found for filename: $(f)")
            end
        end

        filename_integers = [parse(Int, match(r"quenched_(\d+)\.gml$", f).captures[1]) for f in filenames]

        # sort filenames and filename_integers according to filename_integers
        order = sortperm(filename_integers)

        sorted_filenames = filenames[order]
        filename_integers = filename_integers[order]
        
        # for each file, load the network and calculate the bond angle and the bond length
        # standard deviation
        bond_angle_std_vec = Vector{Float64}()
        bond_length_std_vec = Vector{Float64}()

        filename_counter = 0

        for filename in sorted_filenames
            # load the network
            network = NG.load_spatial_network_from_gml(load_path*filename)

            # calculate the bond angle and bond length standard deviation
            bond_angle_std, bond_angle_vec = NA.get_bond_angle_std(network)
            bond_length_std, bond_length_vec = NA.get_bond_length_std(network)

            if !isa(bond_angle_std, Float64) || !isa(bond_length_std, Float64)
                println("Error: bond_angle_std or bond_length_std is not a float")
                println("bond_angle_std: $(bond_angle_std), bond_length_std: $(bond_length_std)")
                println("filename: $(filename)")
                continue
            end

            push!(bond_angle_std_vec, float(bond_angle_std))
            push!(bond_length_std_vec, float(bond_length_std))

            # every 100th network, print Progress
            if filename_counter % 100 == 0
                println("Temperature: $(temperature), progress: $(filename_counter)/$(length(sorted_filenames))")
                println("Bond angle std: $(bond_angle_std), bond length std: $(bond_length_std)")
            end

            filename_counter += 1

        end

        # create a dictionary with the bond angle and bond length standard deviations
        order_metric_evolution_dict = Dict{String, Any}(
            "temperature" => temperature,
            "filename_integers" => filename_integers,
            "bond_angle_std_vec" => bond_angle_std_vec,
            "bond_length_std_vec"  => bond_length_std_vec)

        GU.save_dict_to_h5(order_metric_evolution_dict, save_path*"T_$(temperature)_order_metric_evolution.h5")
        
    end

end

get_order_metric_evolution()



temperature_vec = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.2, 0.22, 0.24]

load_path = load_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks_evolution\216_vertices_bond_bending_0.285\\"

# load network
network = NG.load_spatial_network_from_gml(load_path*"216_vertices_T_0.2_heat_cool_0.1_per_mc_quenched_3000.gml")

evolution_dict = GU.load_h5_dict(load_path*"216_vertices_T_0.2_heat_cool_0.1_per_mc_quenched_3000_evolution.h5")

println(sum(evolution_dict["move_accepted_vec"]))

# plot network
NG.plot_spatial_network(network)



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_comparison\\"

bond_bending_vec = [0.21, 0.285, 0.36]

run_vec = [1, 2, 3, 4, 5]

# loop through all bond bending values and runs. For each run,
# load the dicitionary all_order_metrics.h5. Save each vector in the dictionary in a column of
# DataFrame all_data_dict

all_data_dict = Dict()

for bond_bending in bond_bending_vec
    for run in run_vec
        filename = path * string(bond_bending) * raw"\run_" * string(run) * raw"\all_order_metrics.h5"
        
        if bond_bending == 0.21 && run == 1
            current_order_metrics_dict = GU.load_h5_dict(filename)

            for key in keys(current_order_metrics_dict)
                all_data_dict[key] = current_order_metrics_dict[key]
            end

            all_data_dict["run"] = fill(run, length(current_order_metrics_dict["filenames_vec"]))
            all_data_dict["bond_bending_const"] = fill(bond_bending, length(current_order_metrics_dict["filenames_vec"]))
        else
            current_order_metrics_dict = GU.load_h5_dict(filename)

            for key in keys(current_order_metrics_dict)
                all_data_dict[key] = vcat(all_data_dict[key], current_order_metrics_dict[key])
            end

            all_data_dict["run"] = vcat(all_data_dict["run"], fill(run, length(current_order_metrics_dict["filenames_vec"])) )
            all_data_dict["bond_bending_const"] = vcat(all_data_dict["bond_bending_const"], fill(bond_bending, length(current_order_metrics_dict["filenames_vec"])) ) 
        end
    end
end

all_data_frame = DataFrames.DataFrame(all_data_dict)
CSV.write(save_path * "order_metrics_all_networks.csv", all_data_frame)



path_vec = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_5\\",
            raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.36\run_1\\"]

save_path_vec = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\geometrical_models\216_vertices_bond_bending_0.21\run_5\\",
                 raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\geometrical_models\216_vertices_bond_bending_0.36\run_1\\"]

filename_vec = ["216_vertices_T_0.11_heat_cool_0.1_per_mc_quenched", "216_vertices_T_0.11_heat_cool_0.1_per_mc_quenched"]

for i in eachindex(path_vec)

    spatial_network = NG.load_spatial_network_from_gml(path_vec[i]*filename_vec[i]*".gml")

    save_filename = string(filename_vec[i], "_br_0.26")
        
    NG.save_mesh_from_spatial_network(spatial_network, save_filename; save_path = save_path_vec[i], bond_radius = 0.26,
        vector_out_of_supercell_length = 1)
end



nr_vertices = 216
temperature_gradient = 0.1
bond_bending_const_vec = [0.135, 0.435]

temperatures_vec = [collect(0.07:0.01:0.17), vcat(collect(0.12:0.01:0.2), collect(0.22:0.02:0.26))]

#for bond_bending in bond_bending_const_vec
for i in eachindex(bond_bending_const_vec)
    bond_bending = bond_bending_const_vec[i]
    temperatures = temperatures_vec[i]

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_gradient = temperature_gradient, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end


nr_vertices = 216
temperature_gradient = 0.1
bond_bending_const_vec = [0.06, 0.51]

temperatures_vec = [collect(0.06:0.01:0.16), vcat(collect(0.13:0.01:0.21), collect(0.23:0.02:0.27))]

#for bond_bending in bond_bending_const_vec
for i in eachindex(bond_bending_const_vec)
    bond_bending = bond_bending_const_vec[i]
    temperatures = temperatures_vec[i]

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_gradient = temperature_gradient, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.135/evolution_dicts/"
save_path = "../structures/random_networks/216_vertices_bond_bending_0.135/"

NG.generate_spatial_networks_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())



bond_bending = 0.06

spatial_network_path = "../structures/random_networks/216_vertices_bond_bending_"*string(bond_bending)*"/"
structure_dict_path = "../structures/random_networks/binary_structures/216_vertices_bond_bending_"*string(bond_bending)*"/"
analysis_data_path = "../analysis_data/random_networks/216_vertices_bond_bending_"*string(bond_bending)*"/"

NA.get_all_dicts_from_networks_multithreading(
    spatial_network_path,
    structure_dict_path,
    analysis_data_path;
    bond_radius = 0.35,
    voxel_edge_length = 0.1,
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())



analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_"

bond_bending_vec = [0.06, 0.135,  0.435, 0.51] #0.21, 0.285, 0.36,

for bond_bending in bond_bending_vec

    for run in 1:5
        NA.get_small_length_scale_order_metrics_all_files(
            analysis_data_path*string(bond_bending)*"\\run_"*string(run)*"\\",;
            l_max_steinhardt_q_l = 12,
            save_result = true,)
    end

    println("Finished bond_bending = ", bond_bending)
end



analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_"

bond_bending_vec = [0.06, 0.135, 0.21, 0.285, 0.36, 0.435, 0.51]

for bond_bending in bond_bending_vec

    for run in 1:5
        # for the current folder, get all filenames that contain the string "pore_size_distribution"
        all_filenames = readdir(analysis_data_path*string(bond_bending)*"\\run_"*string(run)*"\\")
        pore_size_distribution_filenames = [filename for filename in all_filenames
            if occursin("pore_size_distribution", filename)]

        # get the second moment of the pore size distribution for each file
        for pore_size_distribution_filename in pore_size_distribution_filenames
            pore_size_distribution_dict = GU.load_h5_dict(analysis_data_path*string(bond_bending)*"\\run_"*string(run)*"\\"*pore_size_distribution_filename)

            small_scale_filename = pore_size_distribution_filename[1:end-26]*"_small_scale_order_metrics.h5"
            small_scale_order_metrics_dict = GU.load_h5_dict(analysis_data_path*string(bond_bending)*"\\run_"*string(run)*"\\"*small_scale_filename)

            second_moment = NA.get_pore_size_distribution_second_moment(pore_size_distribution_dict)
            small_scale_order_metrics_dict["pore_size_distribution_second_moment"] = second_moment

            GU.save_dict_to_h5(small_scale_order_metrics_dict, analysis_data_path*string(bond_bending)*"\\run_"*string(run)*"\\"*small_scale_filename)

        end

        NA.get_small_length_scale_order_metrics_all_files(
            analysis_data_path*string(bond_bending)*"\\run_"*string(run)*"\\",;
            l_max_steinhardt_q_l = 12,
            save_result = true,)
    end

    println("Finished bond_bending = ", bond_bending)
end


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_comparison\\"

bond_bending_vec = [0.06, 0.135, 0.21, 0.285, 0.36, 0.435, 0.51]

run_vec = [1, 2, 3, 4, 5]

# loop through all bond bending values and runs. For each run,
# load the dicitionary all_order_metrics.h5. Save each vector in the dictionary in a column of
# DataFrame all_data_dict

all_data_dict = Dict()

for bond_bending in bond_bending_vec
    for run in run_vec
        filename = path * string(bond_bending) * raw"\run_" * string(run) * raw"\all_order_metrics.h5"
        
        if bond_bending == 0.06 && run == 1
            current_order_metrics_dict = GU.load_h5_dict(filename)

            for key in keys(current_order_metrics_dict)
                all_data_dict[key] = current_order_metrics_dict[key]
            end

            all_data_dict["run"] = fill(run, length(current_order_metrics_dict["filenames_vec"]))
            all_data_dict["bond_bending_const"] = fill(bond_bending, length(current_order_metrics_dict["filenames_vec"]))
        else
            current_order_metrics_dict = GU.load_h5_dict(filename)

            for key in keys(current_order_metrics_dict)
                all_data_dict[key] = vcat(all_data_dict[key], current_order_metrics_dict[key])
            end

            all_data_dict["run"] = vcat(all_data_dict["run"], fill(run, length(current_order_metrics_dict["filenames_vec"])) )
            all_data_dict["bond_bending_const"] = vcat(all_data_dict["bond_bending_const"], fill(bond_bending, length(current_order_metrics_dict["filenames_vec"])) ) 
        end
    end
end

all_data_frame = DataFrames.DataFrame(all_data_dict)
CSV.write(save_path * "more_order_metrics_all_networks.csv", all_data_frame)



function my_func()
    
    path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"

    nr_vertices_vec = [216]
    bond_bending_vec = [0.06, 0.135, 0.435, 0.51]

    for nr_vertices in nr_vertices_vec
        for bond_bending in bond_bending_vec
            for run in 1:5
                current_path = path*string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\run_"*string(run)*"\\"

                save_path = path*"networks_for_simulation\\"*string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\run_"*string(run)*"\\"

                # get all files in the current path
                files = readdir(current_path)
                # get vector of all gml files
                gml_files = [file for file in files if endswith(file, ".gml")]

                # loop through each file that is a gml file
                for gml_file in gml_files
                    filename  = gml_file[1:end-4]
                    spatial_network = NG.load_spatial_network_from_gml(current_path*filename*".gml")

                    spatial_network_for_simulation = NG.get_spatial_network_for_simulation!(
                        spatial_network;
                        vector_out_of_supercell_length = 1,
                        duplicate_bonds_close_to_supercell_edge = true,
                        save_result = true,
                        filename = filename,
                        save_path = save_path)
                end
            end

            println("finished run for ", nr_vertices, " vertices and bond bending ", bond_bending)
        end
    end

end

my_func()