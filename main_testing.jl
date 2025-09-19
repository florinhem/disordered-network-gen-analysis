
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import GeometryBasics
import Polylabel
import Random
import Distributions


network_type_vec = ["pcu_cn_4_5_6"]
nr_vertices_vec = [216] 
bond_bending_const_vec = [0.0, 0.25, 0.5, 0.75, 1.0]
theta_ground_state_vec = [180.0]
acceptance_probability_vec = [0.001]
relax_globally_after_threshold_cycle_vec = [true]
shell_nr_vec = [4]

NA.print_melting_temperatures(
    ;
    network_type_vec,
    nr_vertices_vec,
    bond_bending_const_vec,
    theta_ground_state_vec,
    acceptance_probability_vec,
    relax_globally_after_threshold_cycle_vec,
    shell_nr_vec
    )