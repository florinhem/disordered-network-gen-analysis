
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Random

# initial guess
p0_ctn = [0.28, 0.7, 0.25, 0.2, 0.34, 0.05]
p0_dia = [0.28, 0.7, 0.25, 0.2, 0.4, 0.05]
p0_lcs = [0.28, 0.7, 0.25, 0.2, 0.38, 0.05]
p0_srs = [0.28, 0.7, 0.25, 0.2, 0.3, 0.05]


network_type_vec = ["ctn", "dia", "lcs", "srs"]
p0_list = [p0_ctn, p0_dia, p0_lcs, p0_srs]

i = 4
network_type = network_type_vec[i]


# load the order metrics dict
analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\neural_network_networks\\" * network_type * raw"\\"

# get the r_t dict path
r_t_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\photonics\tidy3d\simulation_data\neural_network_networks\\" * network_type * raw"\run_1_2_r_t_low_n\\"

files = readdir(r_t_dict_path, join=true)
r_t_files = filter(x -> endswith(x, "r_t_only.hdf5"), files)

# shuffle files
Random.shuffle!(r_t_files)

#file_nr = 1
for file_nr in 1:10
    r_t_dict = GU.load_h5_dict(r_t_files[file_nr])

    # get the base filename
    base_filename = basename(r_t_files[file_nr])[1:end-28]

    NA.get_reflection_peak_height_to_width(r_t_dict;
        p0=[0.28, 0.7, 0.25, 0.2, 0.3, 0.05],
        save_plot = true,
        save_path = raw"..\..\photonics\tidy3d\plots\neural_network_networks\\" * network_type * raw"\run_1_2_r_t_low_n\\"*base_filename)

end