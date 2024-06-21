
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .BinaryDataAnalysis as BDA
import .GeneralUtilities as GU

import Plots

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

# julia --threads 23


i = 5

print_lock = Threads.ReentrantLock()

evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.21/evolution_dicts/"

save_path = "../structures/random_networks/216_vertices_bond_bending_0.21/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = false,
print_lock = print_lock)
