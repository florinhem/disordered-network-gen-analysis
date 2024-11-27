
include("structure_analysis_modules.jl")
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
using MetaGraphsNext


path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\\"
filename=raw"m_BTMC_q_t__N=216_T=0.205_Beta=0.0_GradT=0.1_StepsPerT=0.01_Theta_GS=30.0_Trial=1.gml"
spatial_network = NG.load_spatial_network_from_gml(path*filename)

#=
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216, network_type="diamond", bond_bending_const=0.285, min_ring_size=3)
spatial_network = NG.get_periodic_network(evolution_dict)
=#
function plot_SN(;
    filename::String,
    spatial_network::MetaGraphsNext.MetaGraph
    )

    P=Plots.plot()

    for vertex in MetaGraphsNext.labels(spatial_network)
        x=spatial_network[vertex]["position"][1]
        y=spatial_network[vertex]["position"][2]
        z=spatial_network[vertex]["position"][3]
        println(x)
        Plots.plot!([x],[y],[z],seriestype=:scatter,legend=false,camera = (50, 40))
    end

    plot_save_path = raw".\simulations\analysis_plot\\"
    plot_filename_start = "Plot_NW_"

    plot_total_path=(plot_save_path
        *plot_filename_start
        *filename
        *".png")

    Plots.savefig(P,plot_total_path)

end

plot_SN(
    filename=filename,
    spatial_network=spatial_network
)