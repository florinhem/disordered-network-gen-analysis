
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


nr_samples = 1000

# choose random beta values for the samples between 0 and 1
beta_vec = rand(nr_samples)

# get the melting temperature for the beta values
t_melt_vec = [NA.get_melting_temperature("lcs", beta) for beta in beta_vec]

# get random values of t_max between t_melt/3 and 3*t_melt 
t_max_vec = t_melt_vec .* (1/3 .+ 8/3 .* rand(nr_samples))

# get random values of the heating/cooling gradient between t_melt/10 and 
# t_melt
t_gradient_vec = t_melt_vec .* (1/10 .+ 9/10 .* rand(nr_samples))


nr_vertices = 192
network_type = "lcs"
theta_ground_state = 180.0

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\neural_network_networks\lcs\evolution_dicts_3\\"

for i in 1:nr_samples
    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(t_max_vec[i],
        temperature_gradient = t_gradient_vec[i], 
        nr_monte_carlo_steps_per_temperature = 0.01,
        quench = true )


    evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,     
        temperature_vec = temperature_vec,
        nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
        bond_bending_const = beta_vec[i], network_type = network_type,
        theta_ground_state = theta_ground_state,)

    filename = Format.format("lcs_beta_{1:.4f}_t_max_{2:.4f}_t_gradient_{3:.4f}", beta_vec[i], t_max_vec[i], t_gradient_vec[i])
    GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")
end