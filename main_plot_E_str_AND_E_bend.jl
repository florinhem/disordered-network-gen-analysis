
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
