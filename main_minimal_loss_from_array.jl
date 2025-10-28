
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Statistics

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\neural_networks\predictions\local_relaxation\ctn\\"

nr_layers = 4
nr_neurons = 67
nr_pca_components = 10

filename = "ctn_predictions_nr_layers_$(nr_layers)_nr_neurons_$(nr_neurons)_full_pca_$(nr_pca_components)_weights_very_focussed.h5"

save_filename = "ctn_predictions_nr_layers_$(nr_layers)_nr_neurons_$(nr_neurons)_full_pca_$(nr_pca_components)_minimal_loss_weights_very_focussed_order_metrics.h5"

data_dict = GU.load_h5_dict(analysis_data_path*filename)

predictions_array = data_dict["predictions_array"]
loss_array = data_dict["loss_array"]

# permute dims of the arrays to have the shape (bond_bending_const, t_max, t_gradient )
predictions_array = permutedims(predictions_array, (4, 3, 2, 1))
loss_array = permutedims(loss_array, (3, 2, 1))

bond_bending_const_vec = data_dict["bond_bending_const_vec"]
t_max_vec = data_dict["t_max_vec"]
t_gradient_vec = data_dict["t_gradient_vec"]

# find the 3d window of size (3,3,3) in the loss array where the average loss
# is the smallest
function get_window_smallest_average(loss_array; window_size=3)
    min_loss_value = Inf
    min_i = 0
    min_j = 0
    min_k = 0

    half_window = div(window_size, 2)

    for i in 1+half_window:(size(loss_array, 1)-half_window)
        for j in 1+half_window:(size(loss_array, 2)-half_window)
            for k in 1+half_window:(size(loss_array, 3)-half_window)
                window = loss_array[(i-half_window):(i+half_window), (j-half_window):(j+half_window), (k-half_window):(k+half_window)]
                window_mean = Statistics.mean(window)
                if window_mean < min_loss_value
                    min_loss_value = window_mean
                    min_i = i
                    min_j = j
                    min_k = k
                end
            end
        end
    end
    return (min_i, min_j, min_k), min_loss_value
end

min_loss_index, min_loss_value = get_window_smallest_average(loss_array, window_size=7)


# print the values of the parameters at this index
println("Minimum positive loss: $min_loss_value")
println("At bond_bending_const = $(bond_bending_const_vec[min_loss_index[1]])")
println("At t_max = $(t_max_vec[min_loss_index[2]])")
println("At t_gradient = $(t_gradient_vec[min_loss_index[3]])")

max_loss_value = maximum(loss_array[loss_array .< 999.0])
max_loss_index = findfirst(x -> x == max_loss_value, loss_array)

# print the values of the parameters at this index
println("Maximum loss: $max_loss_value")
println("At bond_bending_const = $(bond_bending_const_vec[max_loss_index[1]])")
println("At t_max = $(t_max_vec[max_loss_index[2]])")
println("At t_gradient = $(t_gradient_vec[max_loss_index[3]])")

# plot a heatmap of the predictions for fixed bond_bending_const = 5.0 with a
# logarithmic color scale
fixed_bond_bending_const = bond_bending_const_vec[min_loss_index[1]]
bond_bending_const_index = min_loss_index[1]

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
    "ring_size_mean",
    "ring_size_std",
    "ring_radius_mean",
    "ring_radius_std",
    "critical_pore_radius",
    "anisotropy_metric_from_structure_factor",
    "anisotropy_metric_from_structure_factor_bonds",
    "hyperuniformity_alpha_value",
    "hyperuniformity_alpha_uncertainty"
]

order_metric_vec_at_minimal_loss = predictions_array[min_loss_index[1], min_loss_index[2], min_loss_index[3], :]

# save the order metric vec as a dictionary to a .h5 file with the q_l_lavues 
# and the q_l_uncertainties as vectors 
order_metric_dict = Dict{String, Any}()
for i in 1:length(part1)
    order_metric_dict[part1[i]] = order_metric_vec_at_minimal_loss[i]
end

order_metric_dict["q_l_vec_values"] = order_metric_vec_at_minimal_loss[length(part1)+1 : length(part1)+length(part2)]
order_metric_dict["q_l_vec_uncertainties"] = order_metric_vec_at_minimal_loss[length(part1)+length(part2)+1 : length(part1)+length(part2)+length(part3)]

for i in 1:length(part4)
    order_metric_dict[part4[i]] = order_metric_vec_at_minimal_loss[length(part1)+length(part2)+length(part3)+i]
end

# print all keys and values of the order_metric_dict
for (key, value) in order_metric_dict
    println("$key => $value")
end 


GU.save_dict_to_h5(order_metric_dict, 
            analysis_data_path*save_filename)