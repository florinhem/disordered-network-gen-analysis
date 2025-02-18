# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

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

                            for i in eachindex(nr_trials_per_temperature_array)

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

                                    
                                    for vertex in MetaGraphsNext.labels(spatial_network)
                                        println(spatial_network[vertex]["coordination_nr"])
                                    end
                                    

                                    println(NA.get_dihedral_angle_vec(spatial_network))
                                    println(length(NA.get_dihedral_angle_vec(spatial_network)))
                                    println(NA.get_dihedral_angle_ratio_peak_to_avg(spatial_network,10/360*2*pi))

            

                                    println(NG.relax_single_vertex_keating_efficiently!(spatial_network,1))

                                    println(NG.gradient_keating_efficient(spatial_network,1))

                                    println(NG.hessian_keating_efficient(spatial_network,1))

                                    println(length(NG.get_all_chains(spatial_network)))

                                    

                                    #=println(spatial_network)
                                    println(spatial_network[])
                                    println(spatial_network[1])
                                    #println(spatial_network[]["coordination_nr_vec"])
                                    println(spatial_network[1]["coordination_nr"])=#


                                    
                                    evolution_dict=GU.load_h5_dict(total_path*"_evolution.h5")
                                    #println(evolution_dict)
                                    println(NG.get_poisson_random_network(evolution_dict))

                                    #println(NG.load_spatial_network_from_gml(total_path*".gml"))

                                    println(NG.get_neighbor_positions_mat(spatial_network,1;exclude_vertices=[]))

                                    display(NG.get_next_neighbor_positions_arr(spatial_network,1))
                                    
                                    display(NG.get_next_neighbor_positions_arr(spatial_network,2))

                                    println(NG.get_incorrectly_coordinated_vertices(spatial_network))
                                    
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
end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[64],             #[216],
    maximal_temperature_array=[0.2],    #[0.1,0.125,0.15,0.175,0.2],
    bond_bending_const_array=[0.3],     #[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5],
    temperature_gradient_array=[0.1],     #[0.1],   
    nr_monte_carlo_steps_per_temperature_array=[0.01],     #[0.01], 
    theta_ground_state_array=[110.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_a2_CN",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_rad_ma_png_2_"
)

