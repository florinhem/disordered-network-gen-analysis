
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
    network_type_array,
    nr_vertices_array,
    bond_bending_const_array,
    theta_ground_state_array,
    p_array
    )
    
    println(Threads.nthreads())

    # Pair nr_vertices and network_type using zip
    vertex_network_pairs = zip(network_type_array, nr_vertices_array)

    # Create the iterator for all combinations
    Iter = collect(Iterators.product(
        vertex_network_pairs,
        bond_bending_const_array,
        theta_ground_state_array,
        p_array
    ))

    data=[]

    data_lock = Threads.ReentrantLock()

    Threads.@threads for (
        (network_type, nr_vertices), 
        bond_bending_const, 
        theta_ground_state, 
        p) in Iter

        println(
            "$network_type"*", "*
            "$nr_vertices"*", "*
            "$bond_bending_const"*", "*
            "$theta_ground_state"*", "*
            "$p" )

        evolution_dict = NA.get_evolution_dict(;
            nr_vertices = nr_vertices, 
            network_type=network_type, 
            bond_bending_const=bond_bending_const, 
            min_ring_size=3,
            theta_ground_state=theta_ground_state
        )

        spatial_network = NG.get_periodic_network(evolution_dict)

        #println("sigma_L, $((NA.get_bond_length_std(spatial_network))[1])")
        #println("sigma_A, $((NA.get_bond_angle_std(spatial_network))[1])")

        spatial_network_copy=deepcopy(spatial_network)

        # E start
        E_start=NG.get_total_energy_keating(spatial_network_copy)
        T_melt_array=[]

        for i in collect(range(0, 0, length=nr_vertices)) #We choose nr_vertices so that on average we get once every vertex
            spatial_network=deepcopy(spatial_network_copy)

            # E end
            switched_chain = NG.get_random_chain(
                spatial_network; 
                declined_chains = [],
                remaining_chains = [],
                min_ring_size = evolution_dict["min_ring_size"])

            # switch bonds
            spatial_network = NG.switch_chain!(spatial_network, switched_chain)

            spatial_network = NG.relax_network_keating!(
                spatial_network,
                switched_chain,
                evolution_dict)

            E_end=NG.get_total_energy_keating(spatial_network)
            dE=E_end-E_start
            T_melt=-dE/log(p)
            append!(T_melt_array,T_melt)
        end

        #println("T_melt_array, $T_melt_array")
        #println("T_melt_min, $(minimum(T_melt_array))")

        #println("[\"$network_type\", $nr_vertices, $bond_bending_const, $theta_ground_state, $p, $(minimum(T_melt_array))],")
        result=(network_type, nr_vertices, bond_bending_const, theta_ground_state, p, minimum(T_melt_array))
        
        lock(data_lock)
        try
            push!(data, result)
        finally
            # Always unlock the lock, even if an error occurs
            unlock(data_lock)
        end

    end

    sort!(data, by = x -> (x[1], x[2], x[3], x[4], x[5], x[6])) 


    for (nt, nv, b, t, p, T) in data
        println("[\"$nt\", $nv, $b, $t, $p, $T ],")
    end
end



save_multiple_N_T_trials_beta_gml(;
    network_type_array=["dia"],         #["dia", "srs", "srd", "ctn", "pto", "lcs"],
    nr_vertices_array=[8] .* 2 .^ 3,     #[8, 8, 10, 28, 14, 24] .* 2 .^ 3,
    bond_bending_const_array=[0.0],     #0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0
    theta_ground_state_array=[180.0],   #[110.0, 180.0],
    p_array=[0.06]                             #0.06
)
