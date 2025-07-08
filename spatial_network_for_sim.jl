
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

function spatial_network_for_simulation(;
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
    plot_folder,
    sim_appendix,
    sim_folder
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

            spatial_network=NG.get_spatial_network_for_simulation!(
                spatial_network;
                vector_out_of_supercell_length = 1,
                duplicate_bonds_close_to_supercell_edge = true,
                save_result = true,
                filename = filename*sim_appendix,
                save_path = sim_folder)


            A=NG.plot_spatial_network_2(spatial_network)

            plot_name = filename*plot_appendix
            plot_path = plot_folder*plot_name
            Plots.savefig(A,plot_path)
        else
            println("no, $filename")
        end
    end
end



spatial_network_for_simulation(;
    network_type_array=["pto"],
    nr_vertices_array=[14*2^3],
    maximal_temperature_array=[0.35],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/",
    filename_start="thesis_1",
    plot_appendix="_1.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/",
    sim_appendix="-",
    sim_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/for_sim/"
)

#=
spatial_network_for_simulation(;
    network_type_array=["ctn"],
    nr_vertices_array=[28*3^3],
    maximal_temperature_array=[1.15],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/",
    filename_start="thesis_1",
    plot_appendix="_2.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/",
    sim_appendix="-",
    sim_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/for_sim/"
)
=#

#=
spatial_network_for_simulation(;
    network_type_array=["srd"],
    nr_vertices_array=[10*2^3],
    maximal_temperature_array=[0.22],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/",
    filename_start="thesis_1",
    plot_appendix="_2.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/",
    sim_appendix="-",
    sim_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/for_sim/"
)=#

#=
spatial_network_for_simulation(;
    network_type_array=["pto"],
    nr_vertices_array=[14*1^3],
    maximal_temperature_array=[0.08],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/",
    filename_start="thesis_1",
    plot_appendix="_2.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/",
    sim_appendix="-",
    sim_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/for_sim/"
)
    =#

#=
spatial_network_for_simulation(;
    network_type_array=["ctn"],
    nr_vertices_array=[28*1^3],
    maximal_temperature_array=[0.13],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/",
    filename_start="thesis_1",
    plot_appendix="_2.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/",
    sim_appendix="-",
    sim_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/for_sim/"
)
    =#


#=
spatial_network_for_simulation(;
    network_type_array=["pto"],
    nr_vertices_array=[14*3^3],
    maximal_temperature_array=[1.17],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/",
    filename_start="thesis_1",
    plot_appendix="_2.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/",
    sim_appendix="-",
    sim_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/for_sim/"
)=#





#=
spatial_network_for_simulation(;
    network_type_array=["srd"],
    nr_vertices_array=[10*3^3],
    maximal_temperature_array=[0.33],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/",
    filename_start="thesis_1",
    plot_appendix="_2.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/",
    sim_appendix="-",
    sim_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/for_sim/"
)
    =#

#=
spatial_network_for_simulation(;
    network_type_array=["ctn"],
    nr_vertices_array=[28*1^3],
    maximal_temperature_array=[1.15],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/",
    filename_start="thesis_1",
    plot_appendix="_2.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/",
    sim_appendix="-",
    sim_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/thesis/for_sim/"
)
    =#

#=
spatial_network_for_simulation(;
    network_type_array=["ctn"],
    nr_vertices_array=[28*3^3],
    maximal_temperature_array=[1.15],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/disordered_pto_ctn/",
    filename_start="dpc_1",
    plot_appendix="_1.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/",
    sim_appendix="_1",
    sim_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/disordered_pto_ctn/for_sim/"
)
=#

#=
spatial_network_for_simulation(;
    network_type_array=["pto"],
    nr_vertices_array=[14*3^3],
    maximal_temperature_array=[1.56],
    bond_bending_const_array=[0.25],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[180.0],
    nr_trials_per_temperature_array=[1],
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/disordered_pto_ctn/",
    filename_start="dpc_1",
    plot_appendix="_1.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/",
    sim_appendix="_1",
    sim_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/disordered_pto_ctn/for_sim/"
)=#