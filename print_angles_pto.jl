
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
    network_type,
    save_path,
    filename_start
    )
    
    println(Threads.nthreads())


    nr_vertices=nr_vertices_array[1]
    maximal_temperature=maximal_temperature_array[1]
    bond_bending_const=bond_bending_const_array[1]
    temperature_gradient=temperature_gradient_array[1]
    nr_monte_carlo_steps_per_temperature=nr_monte_carlo_steps_per_temperature_array[1]
    theta_ground_state=theta_ground_state_array[1]
    trial=nr_trials_per_temperature_array[1]

    println("$nr_vertices"*", "*
    "$maximal_temperature"*", "*
    "$bond_bending_const"*", "*
    "$temperature_gradient"*", "*
    "$nr_monte_carlo_steps_per_temperature"*", "*
    "$theta_ground_state"*", "*
    "$trial" )

    evolution_dict = NA.get_evolution_dict(;
        nr_vertices = nr_vertices, 
        network_type=network_type, 
        bond_bending_const=bond_bending_const, 
        min_ring_size=3,
        theta_ground_state=theta_ground_state
    )

    spatial_network = NG.get_periodic_network(evolution_dict)

    for vertex_label in MetaGraphsNext.labels(spatial_network)
   
        neighbor_label_vec::Vector{Int64} = collect(MetaGraphsNext.neighbor_labels(
            spatial_network, 
            vertex_label))

        vertex_coordination_nr=length(neighbor_label_vec)

        for j in 1:(vertex_coordination_nr-1)
            sign1::Int64=sign(neighbor_label_vec[j] - vertex_label)
            vector_j::Vector{Float64}=(sign1*
                spatial_network[vertex_label, neighbor_label_vec[j]]["vector"])
            
            for k in j+1:vertex_coordination_nr
                sign2::Int64=sign1*sign(neighbor_label_vec[k] - vertex_label)
                vector_k::Vector{Float64}=(sign2*
                    spatial_network[vertex_label, neighbor_label_vec[k]]["vector"])
                
                dot_product = LinearAlgebra.dot(vector_j, vector_k)
                norm_product = LinearAlgebra.norm(vector_j) * LinearAlgebra.norm(vector_k)
                angle_rad = LinearAlgebra.acos(clamp(dot_product / norm_product, -1.0, 1.0))
                angle_deg = LinearAlgebra.rad2deg(angle_rad)
                angle_deg_rounded = round(angle_deg; digits=2)
                println("Angle between bonds: $angle_deg_rounded degrees")
                
                
            end
        end
    end
    
    return
end


println("A")

save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array=[10*1^3],
    maximal_temperature_array=[0.0001],
    bond_bending_const_array=[0.3],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    network_type="srd",
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/srd_simulation/",
    filename_start="angles"
)


println("B")

save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array=[28*1^3],
    maximal_temperature_array=[0.0001],
    bond_bending_const_array=[0.3],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    network_type="ctn",
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/srd_simulation/",
    filename_start="angles"
)


println("C")

save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array=[14*1^3],
    maximal_temperature_array=[0.0001],
    bond_bending_const_array=[0.3],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    network_type="pto",
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/srd_simulation/",
    filename_start="angles"
)

println("D")

save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array=[24*1^3],
    maximal_temperature_array=[0.0001],
    bond_bending_const_array=[0.3],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    network_type="lcs",
    save_path ="C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/srd_simulation/",
    filename_start="angles"
)