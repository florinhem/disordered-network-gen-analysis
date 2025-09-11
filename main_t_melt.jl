
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

#import MetaGraphsNext
#import Graphs
#import Plots
#Plots.plotlyjs()
#import .Threads
#import Statistics
#import LinearAlgebra

import ProfileView


network_type_vec = ["dia"]
nr_vertices_vec = [216] # 1000, 1728
bond_bending_const_vec = [0.25]
theta_ground_state_vec = [109.5, 180.0]
acceptance_probability_vec = [0.001]

NA.print_melting_temperatures(
    ;
    network_type_vec,
    nr_vertices_vec,
    bond_bending_const_vec,
    theta_ground_state_vec,
    acceptance_probability_vec
    )
