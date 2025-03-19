
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

bond_bending = 1.0

evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_"*string(bond_bending)*"/evolution_dicts/"
save_path = "../structures/random_networks/216_vertices_bond_bending_"*string(bond_bending)*"/"

print_every_nr_attempted_bond_switches = 200
print_progress = true
save_network_after_each_temperature = false
further_evolve_previous_networks = false
runs_vec = collect(1:5)
random_evolution_seed = -1
print_lock = Threads.ReentrantLock()


# get all filenames by reading the evolution dict directory
filenames = readdir(evolution_dicts_directory_path)
filenames_evolution_dicts = filter(
    filename -> endswith(filename, "_evolution.h5"), filenames)

# get vector of filenames and save paths for multiple runs
save_path_all_runs_vec = Vector{String}(undef, 0)
filenames_evolution_dicts_all_runs_vec = Vector{String}(undef, 0)

for run in runs_vec
    append!(save_path_all_runs_vec, (save_path .* "run_" 
        .* string.( Int.( ones(length(filenames_evolution_dicts)) 
        .* run ) ) .* "/" ) )
    append!(filenames_evolution_dicts_all_runs_vec, 
        filenames_evolution_dicts)
end

# create tuples out of the elements of both vectors
save_path_filename_tuple_vec = Vector{Tuple{String, String}}(undef, 
    length(save_path_all_runs_vec))

for i in eachindex(save_path_filename_tuple_vec)
    save_path_filename_tuple_vec[i] = (save_path_all_runs_vec[i],
        filenames_evolution_dicts_all_runs_vec[i] )
end

# split filenames and save_paths into chunks for multi-threading
save_path_filename_tuple_chunks = Iterators.partition(
    save_path_filename_tuple_vec, 
    length(save_path_filename_tuple_vec) ÷ Threads.nthreads())

println(Threads.nthreads())

map(save_path_filename_tuple_chunks) do save_path_filename_tuple_chunk
    Threads.@spawn NG.generate_spatial_networks_from_evolution_dicts_single_thread_multiple_runs(
        save_path_filename_tuple_chunk,
        evolution_dicts_directory_path;
        print_every_nr_attempted_bond_switches 
            = print_every_nr_attempted_bond_switches,
        print_progress = print_progress,
        random_evolution_seed = random_evolution_seed,
        save_network_after_each_temperature 
            = save_network_after_each_temperature,
        further_evolve_previous_networks 
            = further_evolve_previous_networks,
        print_lock = print_lock)
end
