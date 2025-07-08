


# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
Plots.plotlyjs()
import .Threads
import Glob

function plot_single_network(;
    nr_vertices,
    maximal_temperature,
    bond_bending_const,
    temperature_gradient,
    nr_monte_carlo_steps_per_temperature,
    theta_ground_state,
    nr_trials_per_temperature,
    save_path,
    network_type,
    filename_start,
    plot_appendix,
    plot_folder
)

    filename = (filename_start
        *"_NW="*"$network_type"
        *"_N="*"$nr_vertices"
        *"_T="*"$maximal_temperature"
        *"_Beta="*"$bond_bending_const"
        *"_GradT="*"$temperature_gradient"
        *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
        *"_Theta_GS="*"$theta_ground_state"
        *"_Trial="*"$nr_trials_per_temperature"
        )

    total_path=save_path*filename

    path_array=Glob.glob(filename_start*"*",save_path)
    #println("path_array, $path_array")

    if(total_path*".gml" in path_array)
        #println("scatter done")
        println("yes, $filename")
        spatial_network=NG.load_spatial_network_from_gml(total_path*".gml")
        A=NG.plot_spatial_network_2(spatial_network)

        plot_name = filename*plot_appendix
        plot_path = plot_folder*plot_name
        Plots.savefig(A,plot_path)
    else
        println("no, $filename")
    end
end



plot_single_network(;
    nr_vertices=8*3^3,
    maximal_temperature=1.5,
    bond_bending_const=0.0,
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state=180.0,
    nr_trials_per_temperature=1,
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",
    network_type="dia",
    filename_start="mts_5",
    plot_appendix="_1.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)

#=
plot_single_network(;
    nr_vertices=216,
    maximal_temperature=0.0001,
    bond_bending_const=0.1,
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state=180.0,
    nr_trials_per_temperature=1,
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/ft_1/",
    network_type="lcs",
    filename_start="test_1",
    plot_appendix="_2.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
    =#

#=
plot_single_network(;
    nr_vertices=216,
    maximal_temperature=0.2,
    bond_bending_const=0.1,
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state=180.0,
    nr_trials_per_temperature=1,
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/ft_1/",
    network_type="srd",
    filename_start="ft_1",
    plot_appendix="_2.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
    =#

#=
plot_single_network(;
    nr_vertices=216,
    maximal_temperature=0.00001,
    bond_bending_const=0.1,
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state=180.0,
    nr_trials_per_temperature=1,
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/networks_dia_srs_srd_ctn/",
    network_type="srd",
    filename_start="ft_1",
    plot_appendix="_1.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)

=#

#=

network_type="diamond"

plot_single_network(;
    nr_vertices=216,
    maximal_temperature=0.35,
    bond_bending_const=0.5,
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state=100.0,
    nr_trials_per_temperature=1,
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/networks_noBr/",
    network_type=network_type,
    filename_start="m_$(network_type)_noBr_1",
    plot_appendix="_1.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
    =#

#=
network_type="ctn"

plot_single_network(;
    nr_vertices=216,
    maximal_temperature=0.1,
    bond_bending_const=0.2,
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state=100.0,
    nr_trials_per_temperature=1,
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/networks_noBr/",
    network_type=network_type,
    filename_start="m_$(network_type)_noBr_1",
    plot_appendix="_4.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
=#

#=
network_type="srd"

plot_single_network(;
    nr_vertices=216,
    maximal_temperature=0.4,
    bond_bending_const=1.0,
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state=180.0,
    nr_trials_per_temperature=1,
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/networks_noBr/",
    network_type=network_type,
    filename_start="m_$(network_type)_noBr_1",
    plot_appendix="_3.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
    =#

    #=
    network_type="dia"

plot_single_network(;
    nr_vertices=216,
    maximal_temperature=0.00000000001,
    bond_bending_const=700.0,
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state=180.0,
    nr_trials_per_temperature=1,
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/networks_dia_srd_ctn/",
    network_type=network_type,
    filename_start="m_$(network_type)_2",
    plot_appendix="_3.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
    =#