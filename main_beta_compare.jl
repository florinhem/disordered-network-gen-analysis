
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

