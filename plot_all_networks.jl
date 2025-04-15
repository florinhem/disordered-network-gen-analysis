


# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")    #*#

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import .Threads
import Glob
import TestImages
import Plots
Plots.plotlyjs()

#=p = Plots.plot(
        axis = nothing,
        layout = @layout([a [b [c [d;_]; _]; _]]),
        size = (800,400)
    )

    for i=1:4
        Plots.plot!(p[i], img, ratio=1)
    end=#

function combine_plots_grid(plots #=::Vector{Plots.Subplot{Plots.PlotlyJSBackend}} =#, m::Int, n::Int; size=(1200, 800))
    # Ensure the number of plots matches the grid size
    if length(plots) > m * n
        error("Number of plots exceeds the grid size (m x n).")
    end

    # Create a layout for the grid
    grid_layout = Plots.@layout([Plots.grid(m, n)])

    # Combine the plots into the grid
    combined_plot = Plots.plot(plots...; layout=grid_layout, size=size)

    # Display the combined plot
    display(combined_plot)
end

function plot_single_network(;
    nr_vertices,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient,
    nr_monte_carlo_steps_per_temperature,
    theta_ground_state,
    nr_trials_per_temperature,
    save_path,
    filename_start,
    plot_appendix,
    plot_folder
    )

#=
    img = TestImages.testimage("blobs")
    p = Plots.plot(
        axis = nothing,
        
        layout = @layout([a b; c d]),
        size = (800,400)
    )

    for i=1:4
        Plots.plot!(p[i], img, ratio=1,xlabel="x_$i",)
    end
    return p
    =#

    #img = TestImages.testimage("blobs")
    p = Plots.plot(
        axis = nothing,
        layout = Plots.@layout(
            [a a; 
            a a 
            ]),
        size = (800, 400)
    )

    #=
                [a a a a a; 
            a a a a a;
            a a a a a;
            a a a a a;
            a a a a a
            ]),
    =#

    m=2
    n=2
    #=
    for i in 1:m*n
        Plots.plot!(p[i], img, ratio=1, xlabel="x_$i")
    end

    return p
    =#

    # Try to save the picture and not plot it
    # Change the A plots to the p plots.
    # Move the p[i] plots into the 2 for loops
    PlotsList::Vector{Plots.Subplot{Plots.PlotlyJSBackend}}=[]
    filename="test"

    for i in eachindex(maximal_temperature_array)

        maximal_temperature=maximal_temperature_array[i]

        for j in eachindex(bond_bending_const_array)

            bond_bending_const=bond_bending_const_array[j]

            filename = (filename_start
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
                Plots.xlabel!("x")
                Plots.ylabel!("y")
                Plots.zlabel!("z")
                
                append!(PlotsList,A)
                #position=(i-1)*m+(j-1)+1
                #Plots.plot!(p[position], A, ratio=1, xlabel="x_$position")
            else
                println("no, $filename")
            end
        end
    end

    combine_plots_grid(PlotsList, m, n)

    plot_name = filename*plot_appendix
    plot_path = plot_folder*plot_name
    #Plots.savefig(p,plot_path)
end

network_type="ctn"

plot_single_network(;
    nr_vertices=216,
    maximal_temperature_array=[0.1,0.2],
    bond_bending_const_array=[0.1,0.2],
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state=100.0,
    nr_trials_per_temperature=1,
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/networks_noBr/",
    filename_start="m_$(network_type)_noBr_1",
    plot_appendix="_4.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)


#=





# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")    #*#

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import .Threads
import Glob
import TestImages
import Plots
Plots.plotlyjs()

#=p = Plots.plot(
        axis = nothing,
        layout = @layout([a [b [c [d;_]; _]; _]]),
        size = (800,400)
    )

    for i=1:4
        Plots.plot!(p[i], img, ratio=1)
    end=#

function plot_single_network(;
    nr_vertices,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient,
    nr_monte_carlo_steps_per_temperature,
    theta_ground_state,
    nr_trials_per_temperature,
    save_path,
    filename_start,
    plot_appendix,
    plot_folder
    )

#=
    img = TestImages.testimage("blobs")
    p = Plots.plot(
        axis = nothing,
        
        layout = @layout([a b; c d]),
        size = (800,400)
    )

    for i=1:4
        Plots.plot!(p[i], img, ratio=1,xlabel="x_$i",)
    end
    return p
    =#
    img = TestImages.testimage("blobs")
    p = Plots.plot(
        axis = nothing,
        layout = Plots.@layout(
            [a a a a a; 
            a a a a a;
            a a a a a;
            a a a a a;
            a a a a a
            ]),
        size = (800, 400)
    )

    m=5
    n=5
    for i in 1:m*n
        Plots.plot!(p[i], img, ratio=1, xlabel="x_$i")
    end

    return p

    # Try to save the picture and not plot it
    # Change the A plots to the p plots.
    # Move the p[i] plots into the 2 for loops

    

    for i in eachindex(maximal_temperature_array)

        maximal_temperature=maximal_temperature_array[i]

        for j in eachindex(bond_bending_const_array)

            bond_bending_const=bond_bending_const_array[j]

            filename = (filename_start
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
    end
end

network_type="ctn"

plot_single_network(;
    nr_vertices=216,
    maximal_temperature_array=[0.1],
    bond_bending_const_array=[0.2],
    temperature_gradient=0.1,
    nr_monte_carlo_steps_per_temperature=0.01,
    theta_ground_state=100.0,
    nr_trials_per_temperature=1,
    save_path = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/networks_noBr/",
    filename_start="m_$(network_type)_noBr_1",
    plot_appendix="_4.png",
    plot_folder = raw"C:/Users/GlauserV/OneDrive - Université de Fribourg/Anlagen/AMI/Projekt/GitFlorin/code_photonic_structures/simulations/analysis_plot/"
)
=#