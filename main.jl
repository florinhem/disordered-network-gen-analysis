
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU
import Profile
import ProfileView

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007


evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 ,temperature_vec = [0.5],
    nr_monte_carlo_steps_per_temperature_vec = [0.05], min_ring_size = 3)

total_energy_vec = Vector{Float64}()
move_accepted_vec = Vector{Bool}()

graph_dict = NG.get_periodic_network(evolution_dict)
@time graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
        total_energy_vec = total_energy_vec,
        move_accepted_vec = move_accepted_vec,
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200)

@time graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
    evolution_dict; 
    total_energy_vec = total_energy_vec,
    move_accepted_vec = move_accepted_vec,
print_progress = true,
print_every_nr_attempted_bond_switches = 200)

@ProfileView.profview NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200)

@ProfileView.profview NG.evolve_network_temperature_sequence!(graph_dict,
    evolution_dict;
print_progress = true,
print_every_nr_attempted_bond_switches = 200)
