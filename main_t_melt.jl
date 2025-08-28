
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

#import MetaGraphsNext
#import Graphs
#import Plots
#Plots.plotlyjs()
#import .Threads
#import Statistics
#import LinearAlgebra


"""
Print the melting temperatures for various network configurations.
"""
function print_melting_temperatures(
    ;
    network_type_vec,
    nr_vertices_vec,
    bond_bending_const_vec,
    theta_ground_state_vec,
    acceptance_probability_vec
    )

    println(Threads.nthreads())
    # Pair nr_vertices and network_type using zip
    vertex_network_pairs = zip(network_type_vec, nr_vertices_vec)

    # Create the iterator for all combinations
    Iter = collect(Iterators.product(
        vertex_network_pairs,
        bond_bending_const_vec,
        theta_ground_state_vec,
        acceptance_probability_vec
    ))
    data=[]
    data_lock = Threads.ReentrantLock()
    Threads.@threads for (
        (network_type, nr_vertices),
        bond_bending_const,
        theta_ground_state,
        acceptance_probability) in Iter
        println(
            "$network_type"*", "*
            "$nr_vertices"*", "*
            "$bond_bending_const"*", "*
            "$theta_ground_state"*", "*
            "$acceptance_probability" )
        evolution_dict = NA.get_evolution_dict(;
            nr_vertices = nr_vertices,
            network_type=network_type,
            bond_bending_const=bond_bending_const,
            min_ring_size=3,
            theta_ground_state=theta_ground_state
        )
        spatial_network = NG.get_periodic_network(evolution_dict)
        spatial_network_copy=deepcopy(spatial_network)

        # E start
        E_start=NG.get_total_energy_keating(spatial_network_copy)
        T_melt_vec=[]

        #We choose nr_vertices so that on average we get once every vertex
        for i in collect(range(0, 0, length=nr_vertices)) 
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
            T_melt=-dE/log(acceptance_probability)
            append!(T_melt_vec,T_melt)
        end
        
        result=(network_type, nr_vertices, bond_bending_const, 
            theta_ground_state, acceptance_probability, minimum(T_melt_vec))
        lock(data_lock)
        try
            push!(data, result)
        finally
            # Always unlock the lock, even if an error occurs
            unlock(data_lock)
        end
    end
    
    sort!(data, by = x -> (x[1], x[2], x[3], x[4], x[5], x[6]))
    for (nt, nv, b, t, acceptance_probability, T) in data
        println("[\"$nt\", $nv, $b, $t, $acceptance_probability, $T ],")
    end

    return
end


network_type_vec = ["pto", "srd"]
nr_vertices_vec = [112, 270]
bond_bending_const_vec = [0.0, 0.25, 0.5, 0.75, 1.0]
theta_ground_state_vec = [180.0]
acceptance_probability_vec = [0.001]

save_multiple_N_T_trials_beta_gml(
    ;
    network_type_vec,
    nr_vertices_vec,
    bond_bending_const_vec,
    theta_ground_state_vec,
    acceptance_probability_vec
    )