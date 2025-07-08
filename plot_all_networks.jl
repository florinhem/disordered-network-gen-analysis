


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
    network_type_array,
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature_array,
    save_path,
    filename_start,
    plot_appendix,
    plot_folder
    )

    println(Threads.nthreads())

    Iter=collect(Iterators.product(
        network_type_array,
        nr_vertices_array,
        maximal_temperature_array,
        bond_bending_const_array,
        temperature_gradient_array,
        nr_monte_carlo_steps_per_temperature_array,
        theta_ground_state_array,
        nr_trials_per_temperature_array
        ))

    Threads.@threads for (
        network_type,
        nr_vertices,
        maximal_temperature,
        bond_bending_const,
        temperature_gradient,
        nr_monte_carlo_steps_per_temperature,
        theta_ground_state,
        trial) in Iter

        filename = (filename_start
            *"_NW="*"$network_type"
            *"_N="*"$nr_vertices"
            *"_T="*"$maximal_temperature"
            *"_Beta="*"$bond_bending_const"
            *"_GradT="*"$temperature_gradient"
            *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
            *"_Theta_GS="*"$theta_ground_state"
            *"_Trial="*"$trial"
            )

        total_path=save_path*filename

        path_array=Glob.glob(filename_start*"*",save_path)
        #println("path_array, $path_array")

        if(total_path*".gml" in path_array)
            #println("scatter done")
            println("yes, $filename")
            spatial_network=NG.load_spatial_network_from_gml(total_path*".gml")
            A=NG.plot_spatial_network_2(spatial_network)
            display(A)
            return
            plot_name = filename*plot_appendix
            plot_path = plot_folder*plot_name
            Plots.savefig(A,plot_path)
        else
            println("no, $filename")
        end
    end
end



plot_single_network(;
    network_type_array=["dia"],
    nr_vertices_array=[8*3^3],
    maximal_temperature_array=[0.5,1.0,1.5],
    bond_bending_const_array=[0.0,0.25,0.5,0.75,1.0],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/melting_temp_search/",
    filename_start="mts_5",
    plot_appendix="_1.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
#=
plot_single_network(;
    network_type_array=["srd_one_unitcell"],
    nr_vertices_array=[19*1^3],
    maximal_temperature_array=[0.3],
    bond_bending_const_array=[0.3],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/ft_1/",
    filename_start="MC=1_Q=No_1",
    plot_appendix="_1.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
=#
#=
plot_single_network(;
    network_type_array=["lcs"],
    nr_vertices_array=[24*2^3],
    maximal_temperature_array=[0.1],
    bond_bending_const_array=[0.2],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/ft_1/",
    filename_start="test_1",
    plot_appendix="_5.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
=#
#=
plot_single_network(;
    network_type_array=["pto"],
    nr_vertices_array=[14*2^3],
    maximal_temperature_array=[110.0],
    bond_bending_const_array=[0.1],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/ft_1/",
    filename_start="test_1",
    plot_appendix="_4.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
    =#

#=
plot_single_network(;
    network_type_array=["pto"],
    nr_vertices_array=[14*3^3],
    maximal_temperature_array=[0.4],
    bond_bending_const_array=[0.1],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/ft_1/",
    filename_start="test_1",
    plot_appendix="_3.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
    =#

#=
plot_single_network(;
    network_type_array=["dia", "srs", "ctn"],
    nr_vertices_array=[216],
    maximal_temperature_array=[0.1,0.15,0.2],
    bond_bending_const_array=[0.1,0.05,0.025],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/ft_1/",
    filename_start="ft_1",
    plot_appendix="_2.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)=#