
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .BinaryDataAnalysis as BDA
import .GeneralUtilities as GU

import Peaks
import Measurements
import Plots
import LsqFit

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007

# julia --threads 23

print_lock = Threads.ReentrantLock()

i = 4

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_4/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"


println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
print_lock = print_lock)