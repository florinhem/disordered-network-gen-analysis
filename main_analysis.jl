
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import DataFrames
import CSV

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

nr_vertices = 216
temperature_gradient = 0.1
bond_bending_const_vec = [0.085, 0.11, 0.165, 0.6, 0.8, 1.0 ]

temperatures_vec = [collect(0.06:0.01:0.17), 
collect(0.06:0.01:0.17), 
collect(0.08:0.01:0.2), 
vcat(collect(0.13:0.01:0.19), collect(0.21:0.02:0.27)),
vcat(collect(0.14:0.01:0.2), collect(0.22:0.02:0.28)),
vcat(collect(0.16:0.01:0.23), collect(0.25:0.02:0.35))]

#for bond_bending in bond_bending_const_vec
for i in eachindex(bond_bending_const_vec)
    bond_bending = bond_bending_const_vec[i]
    temperatures = temperatures_vec[i]

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_gradient = temperature_gradient, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_gradient)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end


