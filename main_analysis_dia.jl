
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

print_lock = Threads.ReentrantLock()

spatial_networks_path = "../structures/local_relaxation/random/dia/"
analysis_data_path = "../analysis_data/local_relaxation/random/dia/"

NA.get_all_dicts_from_networks_multithreading(
spatial_networks_path,
analysis_data_path;
print_progress = true,
runs_vec = collect(1:2),
print_lock = print_lock)
