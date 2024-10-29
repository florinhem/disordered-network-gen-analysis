
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