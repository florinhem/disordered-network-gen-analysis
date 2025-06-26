
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU


import Plots
import Graphs
import MetaGraphsNext
import LinearAlgebra
import Polynomials
import Format
#import Statistics
#import Measurements

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)



nr_vertices = 224
network_type = "ctn"
theta_ground_state = 180.0
beta = 0.6
t_max = 5.0
t_gradient = 1.0

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\neural_network_targeted\ctn\evolution_dict_2\\"

temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(t_max,
    temperature_gradient = t_gradient, 
    nr_monte_carlo_steps_per_temperature = 0.01,
    quench = true )
evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,     
    temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
    bond_bending_const = beta, network_type = network_type,
    theta_ground_state = theta_ground_state,)
filename = Format.format("ctn_beta_{1:.4f}_t_max_{2:.4f}_t_gradient_{3:.4f}", beta, t_max, t_gradient)
GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")