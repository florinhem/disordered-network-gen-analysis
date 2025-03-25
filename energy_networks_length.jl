
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")    #*#

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
Plots.plotlyjs()
import .Threads
import Statistics
import LinearAlgebra

function save_multiple_N_T_trials_beta_gml(
    ;
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature_array,
    edge_length_unit_cell_array,
    network_type
    )
    
    println(Threads.nthreads())

    bond_stretching_energy=[]
    nr_vertices=nr_vertices_array[1]
    maximal_temperature=maximal_temperature_array[1]
    bond_bending_const=bond_bending_const_array[1]
    temperature_gradient=temperature_gradient_array[1]
    nr_monte_carlo_steps_per_temperature=nr_monte_carlo_steps_per_temperature_array[1]
    theta_ground_state=theta_ground_state_array[1]
    nr_trials_per_temperature=nr_trials_per_temperature_array[1]

    for edge_length_unit_cell in edge_length_unit_cell_array
            
        println("$nr_vertices"*", "*
		"$maximal_temperature"*", "*
		"$bond_bending_const"*", "*
		"$temperature_gradient"*", "*
        "$nr_monte_carlo_steps_per_temperature"*", "*
        "$theta_ground_state"*", "*
        "$nr_trials_per_temperature"*", "*
        "$edge_length_unit_cell" )
    
        evolution_dict = NA.get_evolution_dict(;
            nr_vertices = nr_vertices, 
            network_type=network_type, 
            bond_bending_const=bond_bending_const, 
            min_ring_size=3,
            theta_ground_state=theta_ground_state
        )

        spatial_network = NG.get_periodic_network(evolution_dict, edge_length_unit_cell)

        println("sigma_L, $((NA.get_bond_length_std(spatial_network))[1])")
        println("sigma_A, $((NA.get_bond_angle_std(spatial_network))[1])")

        spatial_network_copy=deepcopy(spatial_network)

        E=NG.get_stretching_energy_keating(spatial_network_copy)
        append!(bond_stretching_energy,E)
        #println("l=$(edge_length_unit_cell), E: $(E)")
                    
    end
    println("edge_length_unit_cell_array, $edge_length_unit_cell_array")
    println("bond_stretching_energy, $bond_stretching_energy")

    p=Plots.plot(
        edge_length_unit_cell_array,
        bond_stretching_energy,
        xlabel="Edge Length Unit Cell",
        ylabel="Bond Stretching Energy",
        label="Bond Stretching Energy"
    )
    
    # srd
    if cmp(network_type , "srd") == 0
        l1=2
        l2=4/sqrt(2)
        n1=6
        n2=4

        l_weighted_average=(l1*n1+l2*n2)/(n1+n2)

    elseif cmp(network_type , "ctn") == 0
        x=0.2082
        V1=[x, x, x]
        V2=[3/8, 0.0, 1/4]

        # the difference between V1 and V2 gives us the length
        D1=V1 .- V2
        L1=LinearAlgebra.norm(D1)

        l_weighted_average=1/L1
    end
    Plots.vline!(p, [l_weighted_average], color=:red, linewidth=2,label="Weighted Average at, $(round(l_weighted_average,digits=3))")

end

#network_type="dia"     #diamond
#network_type="srs"      #gyroid
#network_type="srd"
network_type="ctn"

save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array=[216],
    maximal_temperature_array=[0.05],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[100.0],
    nr_trials_per_temperature_array=[1],
    edge_length_unit_cell_array=[3.6,3.65,3.675,3.703279044748403,3.725,3.75,3.8],
    network_type=network_type
)