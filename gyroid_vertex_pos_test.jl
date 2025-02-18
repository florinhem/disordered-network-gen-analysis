# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
import Colors
import Glob
import DataFrames
import LaTeXStrings
using StatsPlots
import GraphMakie

function main_test()
    #=diamond=NG.get_diamond_network(8)
    println("(typeof(diamond)), $(typeof(diamond))")
    println("diamond, $diamond")
    SN=NG.convert_original_graph_to_spatial_network(diamond)
    f=NG.plot_spatial_network_2(SN)=#

    #=
    gyroid=NG.get_gyroid_network(8*8)
    println("gyroid, $gyroid")
    SN=NG.convert_original_graph_to_spatial_network(gyroid)
    f=NG.plot_spatial_network_2(SN)
    =#

    #=
    gyroid=NG.get_gyroid_network(8*8*8)
    println("gyroid, $gyroid")
    SN=NG.convert_original_graph_to_spatial_network(gyroid)
    f=NG.plot_spatial_network_2(SN)
    =#
end

main_test()