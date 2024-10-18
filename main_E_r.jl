
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
            temperature_gradient = 10.0, 
            nr_monte_carlo_steps_per_temperature = 2/(18*216),
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
    
    
    Plots.scatter!(r,E,xlabel="bond length [1/d]",ylabel="streching energy [1/(α d^2)]",label="Measured")

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

    # plot the theoretical and taylor function around the equilibrium
    theta_theoretical=collect(0:0.05:pi)
    theta_equilibrium=acos(-1/3)
    r_norm=1
    E_bend=3/8 * (r_norm.^2*cos.(theta_theoretical) .+1/3).^2
    second_taylor_constant=3/8 * bond_bending_const * 12 * (sin(theta_equilibrium)^2-cos(theta_equilibrium)^2-1/3*cos(theta_equilibrium)) 
    E_taylor=1/2 * second_taylor_constant *(theta_theoretical .- theta_equilibrium).^2

    Plots.plot(theta_theoretical,E_bend,label="E_bend",legend=:topright)
    Plots.plot!(theta_theoretical,E_taylor,label="E_taylor")
    Plots.plot!([theta_equilibrium], seriestype="vline", label=false, color=:blue)

    
    # create the network
    evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices, network_type="diamond", bond_bending_const=bond_bending_const, min_ring_size=3)
    spatial_network = NG.get_periodic_network(evolution_dict)


    # heat and cool down
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
        print_every_nr_attempted_bond_switches = 10)


    # prepare for scatter plotting
    θ=[]
    E=[]

    for vertex_label in MetaGraphsNext.labels(spatial_network)

        neighbor_label_vec::Vector{Int64} = collect(MetaGraphsNext.neighbor_labels(
        spatial_network, 
        vertex_label))

        for j in 1:(spatial_network[]["coordination_nr"]-1)
            s1::Int64=sign(neighbor_label_vec[j] - vertex_label)
            v1::Vector{Float64}=spatial_network[vertex_label, neighbor_label_vec[j]][
                "vector"]
            a=s1.*v1
    
            for k in j+1:spatial_network[]["coordination_nr"]
                s2::Int64=sign(neighbor_label_vec[k] - vertex_label)
                v2::Vector{Float64}=spatial_network[vertex_label, neighbor_label_vec[k]]["vector"]
                b=s2.*v2
                append!(θ,acos(LinearAlgebra.dot(a,b)/(LinearAlgebra.norm(a)*LinearAlgebra.norm(b))))
                append!(E,NG.local_bond_bending_energy_keating(spatial_network, vertex_label))
                
            end
        end
    end

    println(θ)
    println(E)


    # prepare for std around equilibrium
    θ_std=Statistics.std(θ)

    Plots.plot!([theta_equilibrium-θ_std], seriestype="vline", label=false, color=:red)
    Plots.plot!([theta_equilibrium+θ_std], seriestype="vline", label=false, color=:red)

    Plots.scatter!(θ,E,xlabel="bond angle [rad]",ylabel="bending energy [1/(α d^2)]",label="Measured")


    # save picture
    save_path = raw".\my_networks\E_bend\\"
    filename = ("E_bend_17_"
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
    maximal_temperature=0.8,
    bond_bending_const=0.285)
