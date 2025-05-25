
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU


# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

sleep(10)
println("Sleep seems to work")

sleep(24*3600) # sleep for 24 hours

save_path = "../structures/neural_network_networks/lcs/"

evolution_dicts_directory_path = "../structures/neural_network_networks/lcs/evolution_dicts_1/"

print_every_nr_attempted_bond_switches = 200
print_progress = true
save_network_after_each_temperature = false
further_evolve_previous_networks = false
runs_vec = [1] #collect(1:2)
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