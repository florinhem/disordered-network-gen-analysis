"""
Copy this code and paste it into the terminal with julia
"""


# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

# import the benchmark module for using @btime
using BenchmarkTools

# import the plot module
using Plots

# prepare short tests
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216, network_type="diamond", bond_bending_const=0.285, min_ring_size=3)
spatial_network = NG.get_periodic_network(evolution_dict)

@benchmark NG.get_total_energy_keating2(spatial_network)

@benchmark NG.get_total_energy_keating(spatial_network)

#@benchmark NG.get_total_energy_keating2(spatial_network)

println("end")