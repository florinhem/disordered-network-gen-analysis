"""
This module contains all functions with tests of functions or modules
"""
module MyTests

# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

# import the text module for using @test
using Test

#import the module MetaGraphsNext
using MetaGraphsNext


# prepare short tests
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216, network_type="diamond", bond_bending_const=0.285, min_ring_size=3)
spatial_network = NG.get_periodic_network(evolution_dict)

# testing with @test
@test true
@test spatial_network[]["supercell_edge_length"]===6.9282032302755105








#Energy test
for vertex in MetaGraphsNext.labels(spatial_network)
    println(vertex)
    println(NG.local_bond_bending_energy_keating(spatial_network, vertex))
    println(NG.local_bond_bending_energy_keating2(spatial_network, vertex))
    @test NG.local_bond_bending_energy_keating(spatial_network, vertex)===
          NG.local_bond_bending_energy_keating2(spatial_network, vertex)
    println("yes test has passed") 
end








# testing multiple tests with table (passed, failed, total, time)
@testset "All tests" begin

    @testset "First" begin
        @test false | true
        @test false | true
    end

    @testset "Second" begin
        # is run, ok
        @test true 
    end
end




end