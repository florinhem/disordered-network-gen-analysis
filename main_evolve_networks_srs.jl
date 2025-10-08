
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

evolution_dicts_directory_path = "../structures/local_relaxation/random/srs/evolution_dicts_2/"

save_path = "../structures/local_relaxation/random/srs/" 

print_every_nr_attempted_bond_switches = 200
print_progress = true
save_network_after_each_temperature = false
further_evolve_previous_networks = false
runs_vec = [2] #collect(1:2)
random_evolution_seed = -1
print_lock = Threads.ReentrantLock()


NG.generate_spatial_networks_from_evolution_dicts_in_directory_multiple_runs(
        evolution_dicts_directory_path,
        save_path;
        print_every_nr_attempted_bond_switches=print_every_nr_attempted_bond_switches,
        print_progress=print_progress,
        random_evolution_seed=random_evolution_seed,
        save_network_after_each_temperature=save_network_after_each_temperature,
        further_evolve_previous_networks=further_evolve_previous_networks,
        runs_vec=runs_vec,
        print_lock=print_lock)
