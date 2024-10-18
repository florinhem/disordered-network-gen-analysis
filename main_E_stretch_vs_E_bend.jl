
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
                filename = ("multiple_gml_6"
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

                println(x,x_err)
                println(y,y_err)

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
                    show=true)

            end
        end
    end

    minimum_temperature=minimum(maximal_temperature_array)
    maximum_temperature=maximum(maximal_temperature_array)
    
    minimum_nr_vertices=minimum(nr_vertices_array)
    maximum_nr_vertices=maximum(nr_vertices_array)

    plot_save_path = raw".\my_networks\stretch_vs_bend\\"
    plot_filename = ("stretch_vs_bend_4"
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
    maximal_temperature_array=[0.0,0.5,1.0,1.5],
    nr_trials_per_temperature=1,
    bond_bending_const=0.285,
    shape_array=[:circle,:rect]
)