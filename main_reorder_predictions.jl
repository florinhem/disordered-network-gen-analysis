
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import MetaGraphsNext

predictions_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\neural_networks\predictions\local_relaxation\ctn\\"

prediction_filename = "ctn_predictions_nr_layers_4_nr_neurons_72_full_pca_10.h5"

predictions_dict = GU.load_h5_dict(predictions_path*prediction_filename)


# reverse all dimensions of the arrays in predictions_dict to correct the
# python row-major ordering
for (k, v) in predictions_dict
    predictions_dict[k] = permutedims(v, reverse(1:ndims(v)))
end


# print shape of all keys in the predictions_dict
for (k, v) in predictions_dict
    println("Key: $k, Shape: ", size(v))
end

GU.save_dict_to_h5(predictions_dict, predictions_path*"ctn_predictions_nr_layers_4_nr_neurons_72_full_pca_10_reordered.h5")