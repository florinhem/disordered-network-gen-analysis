
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\neural_networks\predictions\local_relaxation\ctn\\"

nr_layers = 4
nr_neurons = 75

filename = "ctn_predictions_nr_layers_$(nr_layers)_nr_neurons_$(nr_neurons)_full.h5"

save_filename = "ctn_predictions_nr_layers_$(nr_layers)_nr_neurons_$(nr_neurons)_full_minimal_loss_order_metrics.h5"

data_dict = GU.load_h5_dict(analysis_data_path*filename)

predictions_array = data_dict["predictions_array"]
loss_array = data_dict["loss_array"]

# permute dims of the arrays to have the shape (bond_bending_const, t_max, t_gradient )
predictions_array = permutedims(predictions_array, (4, 3, 2, 1))
loss_array = permutedims(loss_array, (3, 2, 1))

bond_bending_const_vec = data_dict["bond_bending_const_vec"]
t_max_vec = data_dict["t_max_vec"]
t_gradient_vec = data_dict["t_gradient_vec"]

target_bond_bonding_const = 4.8
target_t_max = 5.25
target_t_gradient = 1.5

# find the index of the closest values in the vectors
bond_bending_index = argmin(abs.(bond_bending_const_vec .- target_bond_bonding_const))
t_max_index = argmin(abs.(t_max_vec .- target_t_max))
t_gradient_index = argmin(abs.(t_gradient_vec .- target_t_gradient))

# print the target values and the closest values found
println("Target bond bending constant: $target_bond_bonding_const, closest value found: $(bond_bending_const_vec[bond_bending_index])")
println("Target t_max: $target_t_max, closest value found: $(t_max_vec[t_max_index])")
println("Target t_gradient: $target_t_gradient, closest value found: $(t_gradient_vec[t_gradient_index])")

# find the prediction at these indices 
order_metric_vec_predicted = predictions_array[bond_bending_index, t_max_index, t_gradient_index, :]

# First fixed part
part1 = [
    "bond_length_std",
    "bond_angle_std",
    "dihedral_angle_entropy",
    "bond_orientation_entropy",
    "coordination_nr_mean",
    "coordination_nr_std"
]

# q_l_value_0 ... q_l_value_12
part2 = ["q_l_value_$(i)" for i in 0:12]

# q_l_uncertainty_0 ... q_l_uncertainty_12
part3 = ["q_l_uncertainty_$(i)" for i in 0:12]

# Last fixed part
part4 = [
    "vertex_homogeneity_metric",
    "uncoordinated_neighbor_distance",
    "ring_size_mean",
    "ring_size_std",
    "ring_radius_mean",
    "ring_radius_std",
    "critical_pore_radius",
    "anisotropy_metric_from_structure_factor",
    "anisotropy_metric_from_structure_factor_bonds",
    "hyperuniformity_alpha_value",
]


# save the order metric vec as a dictionary to a .h5 file with the q_l_lavues 
# and the q_l_uncertainties as vectors 
order_metric_dict = Dict{String, Any}()
for i in 1:length(part1)
    order_metric_dict[part1[i]] = order_metric_vec_predicted[i]
end

order_metric_dict["q_l_vec_values"] = order_metric_vec_predicted[length(part1)+1 : length(part1)+length(part2)]
order_metric_dict["q_l_vec_uncertainties"] = order_metric_vec_predicted[length(part1)+length(part2)+1 : length(part1)+length(part2)+length(part3)]

for i in 1:length(part4)
    order_metric_dict[part4[i]] = order_metric_vec_predicted[length(part1)+length(part2)+length(part3)+i]
end

# print all keys and values of the order_metric_dict
for (key, value) in order_metric_dict
    println("$key => $value")
end 


#GU.save_dict_to_h5(order_metric_dict, 
#            analysis_data_path*save_filename)