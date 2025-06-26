
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

evolution_dicts_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\neural_network_targeted\test_networks\evolution_dicts\\"

nr_samples = 20

network_types = ["ctn", "dia", "lcs", "srs"]
beta_vec = [1/4, 2/4, 3/4] 
t_max_over_t_melt_vec = [1/2, 3/2]
t_gradient_over_t_melt_vec = [2/5, 4/5]

nr_vertices_vec = [216, 216, 192, 216] # diamond, diamond, lcs, srs

theta_ground_state = 180.0

for i in eachindex(network_types)
    for beta in beta_vec
        for t_max_over_t_melt in t_max_over_t_melt_vec
            for t_gradient_over_t_melt in t_gradient_over_t_melt_vec

                nr_vertices = nr_vertices_vec[i]
                network_type = network_types[i]
                

                # get the melting temperature for the beta values
                t_melt = NA.get_melting_temperature(network_type, beta)

                # get random values of t_max and t_gradient
                t_max = t_melt * t_max_over_t_melt
                t_gradient = t_melt * t_gradient_over_t_melt

                temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(t_max,
                    temperature_gradient = t_gradient, 
                    nr_monte_carlo_steps_per_temperature = 0.01,
                    quench = true )

                evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,
                    temperature_vec = temperature_vec,
                    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
                    bond_bending_const = beta, network_type = network_type,
                    theta_ground_state = theta_ground_state,)

                filename = Format.format("$(network_type)_beta_{1:.4f}_t_max_{2:.4f}_t_gradient_{3:.4f}", beta, t_max, t_gradient)

                GU.save_dict_to_h5(evolution_dict, evolution_dicts_path*filename*"_evolution.h5")
            end
        end
    end
end
