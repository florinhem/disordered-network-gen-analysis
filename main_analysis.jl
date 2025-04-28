
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU


import Plots
import Graphs
import MetaGraphsNext
import LinearAlgebra
import Polynomials
#import Statistics
#import Measurements

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

save_path =  raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_bond_bending_0.21\run_1\\"

diamond_bonds = GU.load_h5_dict(save_path*"diamond_bonds_structure_factor_angle_averaged.h5")
disorder_bonds = GU.load_h5_dict(save_path*"disorder_bonds_structure_factor_angle_averaged.h5")
disorder_2_bonds = GU.load_h5_dict(save_path*"disorder_2_bonds_structure_factor_angle_averaged.h5")
disorder_3_bonds = GU.load_h5_dict(save_path*"disorder_3_bonds_structure_factor_angle_averaged.h5")

slope_measurement_time_time, alpha = NA.get_hyperuniformity_alpha(diamond_bonds)
slope_measurement_time_time_disorder, alpha_disorder = NA.get_hyperuniformity_alpha(disorder_bonds)
slope_measurement_time_time_disorder_2, alpha_disorder_2 = NA.get_hyperuniformity_alpha(disorder_2_bonds)
slope_measurement_time_time_disorder_3, alpha_disorder_3 = NA.get_hyperuniformity_alpha(disorder_3_bonds)

println("diamond_bonds: slope_measurement_time_time = ", slope_measurement_time_time, ", alpha = ", alpha)
println("disorder_bonds: slope_measurement_time_time = ", slope_measurement_time_time_disorder, ", alpha = ", alpha_disorder)
println("disorder_2_bonds: slope_measurement_time_time = ", slope_measurement_time_time_disorder_2, ", alpha = ", alpha_disorder_2)
println("disorder_3_bonds: slope_measurement_time_time = ", slope_measurement_time_time_disorder_3, ", alpha = ", alpha_disorder_3)
