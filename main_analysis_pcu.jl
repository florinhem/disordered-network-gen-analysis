
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

sleep(46*3600)

print_lock = Threads.ReentrantLock()

spatial_networks_path = "../structures/local_relaxation/random/pcu_cn_4_5_6/"
analysis_data_path = "../analysis_data/local_relaxation/random/pcu_cn_4_5_6/"

NA.get_all_dicts_from_networks_multithreading(
spatial_networks_path,
analysis_data_path;
print_progress = true,
runs_vec = [2,3],
print_lock = print_lock)
