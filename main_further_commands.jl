
"""
These are the calculations for the pachy weevil from the 10.1002advs.202202145 paper
"""

# set raw data path
data_path_raw = raw"..\structures\pachy_10.1002advs.202202145\3Dvolumes\SD_ff04_a51_box5_vox1.tif"


# load data, correct voxel size anisotropy and save data
data_binary = BDA.get_binary_data_from_colorscale(data_path_raw; 
                                        voxel_size=(1,1,1), 
                                        save_data=true, 
                                        data_path_corrected=data_path_corrected)


# set path to voxel size corrected data
data_path_corrected = raw"..\structures\pachy_10.1002advs.202202145\pachy_blue.h5"
      
# compare hyperuniformity criterion for red and blue patches
data_path_corrected_vec = [raw"..\structures\pachy_10.1002advs.202202145\pachy_blue.h5"]

# set labels for plotting 
label_vec = ["blue patch"] # , "red patch", "simple diamond"

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


# loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )

    # determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 10,
                                                                                                    window_positioning="scanned",
                                                                                                    window_shape="spherical")

    # create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


# plot window volume times local volume fraction uncertainty as a function of window edge length
# to determine whether structure is hyperuniform
BDA.plot_volume_fraction_variance_times_window_volume(nr_dimensions_data,
                            plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of scanned positions",
                            save_filename="pachy_scanned_spherical_windows_blue_red_sd")





# compare hyperuniformity criterion for red and blue patches
data_path_corrected_vec = [raw"..\structures\pachy_10.1002advs.202202145\pachy_blue.h5",
                        raw"..\structures\pachy_10.1002advs.202202145\pachy_red.h5",
                        raw"..\structures\pachy_10.1002advs.202202145\simple_diamond.h5"]

# set labels for plotting 
label_vec = ["blue patch", "red patch", "simple diamond"] # , "red patch", "simple diamond"

# in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
# they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


# loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )

    # determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 100,
                                                                                                    window_positioning="random",
                                                                                                    window_shape="spherical")

    # create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "voxel_edge_length" => voxel_edge_length_vec[i],
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="pachy_volume_fraction_variance_random_spherical_"*label_vec[i])
                                                                                            
end


# plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="pachy_volume_fraction_variance_random_spherical_windows_blue_red_sd")




# compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"..\analysis_data\pachy_volume_fraction_variance_random_spherical_blue patch.h5",
                raw"..\analysis_data\pachy_volume_fraction_variance_random_spherical_red patch.h5",
                raw"..\analysis_data\volume_fraction_variance_random_spherical_D_surface.h5"]

# set labels for plotting 
label_vec = ["blue patch", "red patch", "perfect diamond"] # , "red patch", "simple diamond"

# in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
# they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []

# loop through data that will be analyzed
for i in eachindex(data_path_vec)

    # load dictionary that contains the following keys:
    # "window_edge_length_vec"
    # "local_volume_fract_variance_vec"
    # "voxel_edge_length"
    # "label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    # adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


# plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="pachy_volume_fraction_variance_random_spherical_windows_blue_red_sd")


# perform convergence analysis to determine the nr of measurements per distance 
# when calculating the autocovariance function
BDA.convergence_analysis_autocovariance_fct_nr_measurements_per_distance(size_data,
                            volume_fract_tot,
                            data_binary;
                            save_plot = true )




# loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )


    # get autocovariance function as a function of sampling distance
    sampling_distance_vec, autocovariance_fct_vec = BDA.get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                                                                    size_data, 
                                                                                    volume_fract_tot,
                                                                                    data_binary;
                                                                                    nr_measurements_per_distance = 5000)

    # create dictionary for current plot
    plot_dict = Dict("sampling_distance_vec" => sampling_distance_vec,
                    "autocovariance_fct_vec" => autocovariance_fct_vec,
                    "voxel_edge_length" => voxel_edge_length_vec[i],
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="pachy_autocovariance_fct_"*label_vec[i])
                                                                                            
end


# plot real part, imaginary part and absolute value of spectral 
# density as a function of the wavenumber
BDA.plot_autocovariance_fct(plot_dict_vec;
                        title="Autocovariance function",
                        save_plot = true,
                        save_path=raw"..\plots\\",
                        save_filename="pachy_autocovariance_fct_blue_red_sd")


# compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"..\analysis_data\pachy_autocovariance_fct_blue patch.h5",
                raw"..\analysis_data\pachy_autocovariance_fct_red patch.h5",
                raw"..\analysis_data\autocovariance_fct_D_surface.h5"]

# set labels for plotting 
label_vec = ["blue patch", "red patch", "perfect diamond"] # , "red patch", "simple diamond"

# in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
# they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []

# loop through data that will be analyzed
for i in eachindex(data_path_vec)

    # load dictionary that contains the following keys:
    # "sampling_distance_vec"
    # "autocovariance_fct_vec"
    # "voxel_edge_length"
    # "label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    # adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


# plot local volume fraction variance as a function of window edge length
BDA.plot_autocovariance_fct(plot_dict_vec;
                            title="Autocovariance function",
                            save_plot = true,
                            save_path=raw"..\plots\\",
                            save_filename="pachy_autocovariance_fct_blue_red_sd_zoom")



# loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )


    # get spectral density as a function of the wavenumber
    wavenumber_vec, spectral_density_vec = BDA.get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data,
                                                                                        size_data,
                                                                                        volume_fract_tot,
                                                                                        data_binary;
                                                                                        nr_measurements_per_distance = 5000,
                                                                                        nr_wavenumbers=50) 

    # create dictionary for current plot
    plot_dict = Dict("wavenumber_vec" => wavenumber_vec,
                    "spectral_density_vec" => spectral_density_vec,
                    "voxel_edge_length" => voxel_edge_length_vec[i],
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)
    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="pachy_spectral_density_fct_"*label_vec[i])
                                                                                            
end


# plot real part, imaginary part and absolute value of spectral 
# density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_path=raw"..\plots\\",
                        save_filename="pachy_spectral_density_blue_red_sd")


# compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"..\analysis_data\pachy_spectral_density_fct_blue patch.h5",
                raw"..\analysis_data\pachy_spectral_density_fct_red patch.h5",
                raw"..\analysis_data\spectral_density_D_surface.h5"]

# set labels for plotting 
label_vec = ["blue patch", "red patch", "perfect diamond"] # , "red patch", "simple diamond"

# in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
# they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []

# loop through data that will be analyzed
for i in eachindex(data_path_vec)

    # load dictionary that contains the following keys:
    # "wavenumber_vec"
    # "spectral_density_vec"
    # "voxel_edge_length"
    # "label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    # adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


# plot real part, imaginary part and absolute value of spectral 
# density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_path=raw"..\plots\\",
                        save_filename="pachy_spectral_density_blue_red_sd")


# path of original data
data_path = raw"..\structures\pachy\pachy_blue.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_blue"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "blue", 
                                save_path)


# path of original data
data_path = raw"..\structures\pachy\pachy_red.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_red"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                9, 
                                "red", 
                                save_path)

                    

# set paths where statistical data is stored
data_path_vec = [raw"..\analysis_data\pachy\pachy_blue",
raw"..\analysis_data\pachy\pachy_red",
raw"..\analysis_data\nodal_surfaces\D_surface"]

# set path where plot will be stored
save_path = raw"..\plots\pachy\pachy_blue_red_d"

# plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            voxel_edge_length_vec=[10, 9, 8.5],
            label_vec=["P. c. mirabilis blue", "P. c. mirabilis red", "perfect diamond"]
            )


# path of original data
data_path = raw"..\structures\pachy\pachy_blue.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_blue"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "blue", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)


# path of original data
data_path = raw"..\structures\pachy\pachy_red.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_red"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                9, 
                                "red", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)



# set data path
data_path_vec = [raw"..\analysis_data\pachy\pachy_blue",
                raw"..\analysis_data\pachy\pachy_red"]

# set path to save plot
save_path = raw"..\plots\pachy\\"

BDA.plot_statistical_measures(data_path_vec,
            save_path;
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = true
            )



data_path = raw"..\structures\pachy\pachy_red.h5"


# load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

# set data path
dict_path = raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

# load dict
data_dict = BDA.load_h5_dict(dict_path)

# calculate spectral density
sampled_wavenumbers_vec_vec, sampled_wavevectors_array, spectral_density_array = BDA.get_spectral_density(size_data, 
                                                volume_fract_tot,
                                                data_binary;
                                                nr_measurements_per_direction = 1000,
                                                sampling_vec_array = data_dict["sampling_vec_array"],
                                                autocovariance_fct_array = data_dict["autocovariance_fct_array"])

# path where spectral density is saved
save_path = raw"..\analysis_data\pachy\pachy_red_spectral_density_direction.h5"

# create dict to save
saving_dict = Dict("sampled_wavevectors_array" => sampled_wavevectors_array,
                    "sampled_wavenumbers_vec_vec" => sampled_wavenumbers_vec_vec,
                    "spectral_density_array" => spectral_density_array,
                    "voxel_edge_length" => data_dict["voxel_edge_length"],
                    "label" => data_dict["label"])

BDA.save_dict_to_h5(saving_dict; save_path)



# path where spectral density data is saved
dict_path = raw"..\analysis_data\pachy\pachy_red_spectral_density_direction.h5"


# path where plot is saved
save_path = raw"..\plots\pachy\pachy_red_"

spectral_density_dict = BDA.load_h5_dict(dict_path)

BDA.plot_spectral_density_heatmap(spectral_density_dict,
    save_path;
    save_plot = false,
    clims = (0,200),
    wavevector_component_to_fix = 3,
    wavevector_value_fixed = 0)



# path of original data
data_path = raw"..\structures\pachy\pachy_red.h5"

# load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

# set data path
dict_path = raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

# load dict
data_dict = BDA.load_h5_dict(dict_path)

# set matrix of measured direction vectors (in this case the identity matrix)
direction_vec_mat = [1 0 0; 0 1 0; 0 0 1]

# set vector of labels
label_vec = ["pachy red [1,0,0]", "pachy red [0,1,0]", "pachy red [0,0,1]" ]

for i in 1:3

    direction_vec = direction_vec_mat[:,i]

    # determine spectral density
    sampled_wavenumbers_vec, spectral_density_vec = BDA.get_spectral_density_along_direction(size_data, 
                volume_fract_tot,
                data_binary,
                direction_vec;
                sampling_vec_array = data_dict["sampling_vec_array"],
                autocovariance_fct_array = data_dict["autocovariance_fct_array"])
    
    # path where spectral density is saved
    save_path = raw"..\analysis_data\pachy\\"* label_vec[i] *"_spectral_density_direction.h5"

    # create dict to save
    saving_dict = Dict("wavenumber_vec" => sampled_wavenumbers_vec,
                        "spectral_density_vec" => spectral_density_vec,
                        "voxel_edge_length" => data_dict["voxel_edge_length"],
                        "label" => label_vec[i])

    BDA.save_dict_to_h5(saving_dict; save_path)
end


# path of original data
data_path = raw"..\structures\pachy\pachy_red.h5"

# load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

# set data path
dict_path = raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

# load dict
data_dict = BDA.load_h5_dict(dict_path)

# set matrix of measured direction vectors (in this case the identity matrix)
direction_vec_mat = [1/sqrt(2) 1/sqrt(3) 1/sqrt(6); 
                    -1/sqrt(2) 1/sqrt(3) 1/sqrt(6); 
                    0 1/sqrt(3) -2/sqrt(6)]

# set vector of labels
label_vec = ["pachy red 1/sqrt(2)*[1,-1,0]", "pachy red 1/sqrt(3)*[1,1,1]", "pachy red 1/sqrt(6)*[1,1,-2]" ]


naming_vec = ["pachy red [1,-1,0]", "pachy red [1,1,1]", "pachy red [1,1,-2]" ]

for i in 1:3

    direction_vec = direction_vec_mat[:,i]

    # determine spectral density
    sampled_wavenumbers_vec, spectral_density_vec = BDA.get_spectral_density_along_direction(size_data, 
                volume_fract_tot,
                data_binary,
                direction_vec;
                sampling_vec_array = data_dict["sampling_vec_array"],
                autocovariance_fct_array = data_dict["autocovariance_fct_array"])
    
    # path where spectral density is saved
    save_path = raw"..\analysis_data\pachy\\"* naming_vec[i] *"_spectral_density_direction.h5"

    # create dict to save
    saving_dict = Dict("wavenumber_vec" => sampled_wavenumbers_vec,
                        "spectral_density_vec" => spectral_density_vec,
                        "voxel_edge_length" => data_dict["voxel_edge_length"],
                        "label" => label_vec[i])

    BDA.save_dict_to_h5(saving_dict; save_path)
end


naming_vec = ["pachy red [1,-1,0]", "pachy red [1,1,1]", "pachy red [1,1,-2]" ]


# initialize plot dict vec 
plot_dict_vec = []

for i in 1:3
    
    # path where spectral density is saved
    load_path = raw"..\analysis_data\pachy\\"* naming_vec[i] *"_spectral_density_direction.h5"
    
    # load dict
    data_dict = BDA.load_h5_dict(load_path)

    # add current dict to plot dict vector
    push!(plot_dict_vec, data_dict)
end

# path where plot will be saved
save_path = raw"..\plots\pachy\pachy_red_direction_rotated_axes"

# plot the spectral densities
BDA.plot_spectral_density(plot_dict_vec,
                                save_path,
                                save_plot = true,
                                xlims=[0,0.1])


                    
# path of original data
data_path = raw"..\structures\pachy\pachy_blue.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_blue"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "P. c. mirabilis blue", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)


# path of original data
data_path = raw"..\structures\pachy\pachy_red.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\pachy\pachy_red"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                9, 
                                "P. c. mirabilis red", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)


# set paths where statistical data is stored
data_path_vec = [raw"..\analysis_data\pachy\pachy_red",
raw"..\analysis_data\pachy\pachy_blue"]

# set path where plot will be stored
save_path = raw"..\plots\pachy\\"

# plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            spectral_density_xlims = [0,0.1],
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = false,
            plot_spectral_density_direction_bool = true
            )


data_path = raw"..\structures\pachy\pachy_red.h5"

# load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

# set data path
dict_path = raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

# load dict
data_dict = BDA.load_h5_dict(dict_path)


# path where spectral density is saved
save_path = raw"..\analysis_data\pachy\pachy_red_spectral_density_direction.h5"

# calculate spectral density
wavenumber_vec_vec, wavevector_array, spectral_density_array = BDA.get_spectral_density_by_wavevector_array(size_data, 
                                                volume_fract_tot,
                                                data_binary;
                                                nr_measurements_per_direction = 1000,
                                                sampling_distance_vec_vec = data_dict["sampling_distance_vec_vec"],
                                                sampling_vec_array = data_dict["sampling_vec_array"],
                                                autocovariance_fct_array = data_dict["autocovariance_fct_array"],
                                                save_result = false,
                save_path = raw"..\analysis_data\pachy\pachy_red",
                voxel_edge_length = 9,
                label = "P. c. mirabilis red")



# set paths for structure dict, autocovariance fct dict and where spectral_density along direction is saved 
structure_path = raw"..\structures\pachy\pachy_blue_structure.h5"
autocovariance_path = raw"..\analysis_data\pachy\pachy_blue_autocovariance_fct.h5"
save_path = raw"..\analysis_data\pachy\pachy_blue"
plot_path = raw"..\plots\pachy\pachy_blue"

# load structure dict
structure_dict = BDA.load_h5_dict(structure_path)

# load autocovariance fct direction dict
autocovariance_dict = BDA.load_h5_dict(autocovariance_path)

spectral_density_dict = BDA.get_spectral_density_isotrope_by_wavenumber_vec(structure_dict;
    nr_sampling_distances = length(autocovariance_dict["sampling_distance_vec"]),
    nr_measurements_per_distance = autocovariance_dict["nr_measurements_per_distance"],
    sampling_distance_vec = autocovariance_dict["sampling_distance_vec"],
    autocovariance_fct_dict = autocovariance_dict,
    sampling_distance_cutoff = 1000,
    save_result = true,
    save_path = save_path)

BDA.plot_spectral_density([spectral_density_dict],
    plot_path;
    title="Spectral density",
    save_plot = true,
    xlims = (0,0.3))



# set paths for structure dict, autocovariance fct dict and where spectral_density along direction is saved 
structure_path = raw"..\structures\pachy\pachy_red_structure.h5"
autocovariance_path = raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"
save_path = raw"..\analysis_data\pachy\pachy_red"
plot_path = raw"..\plots\pachy\pachy_red"

# load structure dict
structure_dict = BDA.load_h5_dict(structure_path)

# load autocovariance fct direction dict
autocovariance_fct_direction_dict = BDA.load_h5_dict(autocovariance_path)

# calculate and save complete autocovariance fct 
complete_autocovariance_dict = BDA.get_complete_autocovariance_fct_by_sampling_vec_array(
    autocovariance_fct_direction_dict;
    save_result = true,
    save_path = save_path)

println("complete dict calculated")

spectral_density_dict = BDA.get_spectral_density_array_by_fft(complete_autocovariance_dict;
        save_result = true,
        save_path = save_path)

println("spectral density dict calculated")

# plot spectral density
BDA.plot_spectral_density_heatmap(spectral_density_dict,
    plot_path;
    title="Spectral density",
    save_plot = true,
    clims = nothing,
    wavevector_component_to_fix = 3,
    wavevector_value_fixed = 0)


data_path_raw = raw"..\structures\pachy\3Dvolumes\Blue_SI.tif"


structure_dict = BDA.get_structure_dict_from_colorscale(data_path_raw; 
    voxel_size=(10,12,10), 
    label = "P. c. mirabilis blue",
    save_result=true, 
    save_path=raw"..\structures\pachy\pachy_blue")


"""
This is where the calculations for nodal surfaces begin
"""


# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

# set properties of generated data
unit_cell_length = 500
nr_unit_cells = 10


# loop through surfaces
for label in label_vec

    # generate 3d binary data for current nodal surface
    data_binary = BDA.get_binary_data_from_nodal_eqn(unit_cell_length, 
                                                nr_unit_cells,
                                                label)


    # save current nodal surface to h5 file
    BDA.save_nodal_surface_data(data_binary,
                            unit_cell_length, 
                            nr_unit_cells,
                            label)
                                                                                            
end



# compare hyperuniformity criterion for red and blue patches
data_path = raw"..\structures\nodal_surfaces\\"

# set surfaces that are analyzed
label_vec = ["D", "G", "P", "I-WP"] 

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


# loop through data that will be analyzed
for label in label_vec

    # determine data path of binary data for current nodal structure
    current_path = data_path*label*"_surface.h5"

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(current_path )

    # determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 100,
                                                                                                    window_positioning="random",
                                                                                                    window_shape="spherical")

    # create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "voxel_edge_length" => 10,
                    "label" => label )

    push!(plot_dict_vec, plot_dict)

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="volume_fraction_variance_random_spherical_"*label*"_surface")
                                                                                            
end


# plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="volume_fraction_variance_random_spherical_nodal_surfaces")



# loop through data that will be analyzed
for label in label_vec

    # determine data path of binary data for current nodal structure
    current_path = data_path*label*"_surface.h5"

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(current_path )

    # get autocovariance function as a function of sampling distance
    sampling_distance_vec, autocovariance_fct_vec = BDA.get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                                                                    size_data, 
                                                                                    volume_fract_tot,
                                                                                    data_binary;
                                                                                    nr_measurements_per_distance = 20000)

    # create dictionary for current plot
    plot_dict = Dict("sampling_distance_vec" => sampling_distance_vec,
                    "autocovariance_fct_vec" => autocovariance_fct_vec,
                    "voxel_edge_length" => 10,
                    "label" => label )

    push!(plot_dict_vec, plot_dict)

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="autocovariance_fct_"*label*"_surface")
                                                                                            
end


# plot local volume fraction variance as a function of window edge length
BDA.plot_autocovariance_fct(plot_dict_vec;
                            title="Autocovariance function",
                            save_plot = true,
                            save_filename="autocovariance_fct_nodal_surfaces")


                            
# create vector for plot_dicts
plot_dict_vec = []

# loop through data that will be analyzed
for label in label_vec

    # determine data path of binary data for current nodal structure
    current_path = data_path*label*"_surface.h5"

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(current_path )

    # load dictionary that contains the following keys:
    # "sampling_distance_vec"
    # "autocovariance_fct_vec"
    # "voxel_edge_length"
    # "label"
    dict = BDA.load_h5_dict("autocovariance_fct_"*label*"_surface")

    # calculate spectral density for loaded autocovariance function
    wavenumber_vec, spectral_density_vec = BDA.get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data,
                size_data, 
                volume_fract_tot,
                data_binary;
                nr_measurements_per_distance = 20000,
                nr_wavenumbers=200,
                sampling_distance_vec = dict["sampling_distance_vec"],
                autocovariance_fct_vec = dict["autocovariance_fct_vec"])

    # create dictionary for current plot
    plot_dict = Dict("wavenumber_vec" => wavenumber_vec,
                    "spectral_density_vec" => spectral_density_vec,
                    "voxel_edge_length" => 10,
                    "label" => label )

    push!(plot_dict_vec, plot_dict)

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
                        save_filename="spectral_density_"*label*"_surface")
                                                                                            
end

# plot real part, imaginary part and absolute value of spectral 
# density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_filename="spectral_density_nodal_surfaces")


# set paths where statistical data is stored
data_path_vec = (raw"..\analysis_data\nodal_surfaces\\"
                    .* ["D", "I-WP", "P", "G"] .* "_surface" )

# set path where plot will be stored
save_path = raw"..\plots\nodal_surfaces\nodal_surfaces"

# plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path)



# get data essentials of stervi data
data_path = raw"..\structures\stervi\stervi_green.h5"
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )
nr_sampling_distances = BDA.get_nr_sampling_distances(mean_edge_length_data)


# now analyze I-WP with the nr of sampling distances of the stervi weevil
# path of original data
data_path = raw"..\structures\nodal_surfaces\I-WP_surface.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\nodal_surfaces\I-WP_surface_fewer_sampling_distances"
    

# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "I-WP", 
                                save_path;
                                nr_sampling_distances = nr_sampling_distances,
                                save_local_volume_fraction_variance=false)



label_vec = ["D", "I-WP", "P", "G"]

for label in label_vec

    # path of original data
    data_path = raw"..\structures\nodal_surfaces\\"*label*"_surface.h5"

    # path where analysis data will be saved
    save_path = raw"..\analysis_data\nodal_surfaces\\"*label*"_surface"

    # calculate and save all statistical measures
    BDA.save_statistical_measures(data_path, 
                                    10, 
                                    label, 
                                    save_path;
                                    nr_sampling_distances = 500,
                                    save_local_volume_fraction_variance=false)

end


# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

# set properties of generated data
unit_cell_length = 500
nr_unit_cells = 1


# loop through surfaces
for label in label_vec

    # generate 3d binary data for current nodal surface
    data_binary = BDA.get_binary_data_from_nodal_eqn(unit_cell_length, 
                                                nr_unit_cells,
                                                label)


    # save current nodal surface to h5 file
    BDA.save_nodal_surface_data(data_binary,
                            unit_cell_length, 
                            nr_unit_cells,
                            "single_unit_cell_"*label)
                                                                                            
end


# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

# set data path
data_path_vec = raw"..\structures\nodal_surfaces\single_unit_cell_" .* label_vec .* "_surface.h5" 


# set path where autocovariance dict will be stored
save_path_vec = raw"..\analysis_data\nodal_surfaces\single_unit_cell_" .* label_vec .* "_surface"


for i in eachindex(data_path_vec)

    # get data and essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_vec[i])

    # get autocovariance function array
    sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(
                        size_data,
                        volume_fract_tot,
                        data_binary,
                        nr_measurements_per_direction=5000,
                        save_result = true,
                        save_path = save_path_vec[i],
                        voxel_edge_length = 10,
                        label = label_vec[i])

    println(label_vec[i]*" done")

end



# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

# set path where autocovariance dict of single unit cell is stored
suc_path_vec = raw"..\analysis_data\nodal_surfaces\single_unit_cell_" .* label_vec .* "_surface_autocovariance_fct_direction.h5"

save_path_vec = raw"..\analysis_data\nodal_surfaces\\" .* label_vec .* "_surface"


for i in eachindex(suc_path_vec)

    # load dict
    suc_dict = BDA.load_h5_dict(suc_path_vec[i])

    # get autocovariance function array
    sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.extrapolate_periodic_data_autocovariance_fct_by_sampling_vec_array(
                suc_dict;
                save_result = true,
                save_path = save_path_vec[i])

    println(label_vec[i]*" done")

end



# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 


# path where analysis data will be saved
data_path_vec = raw"..\structures\nodal_surfaces\\" .* label_vec .* "_surface.h5"

# path where analysis data will be saved
save_path_vec = raw"..\analysis_data\nodal_surfaces\\" .* label_vec .* "_surface"

for i in eachindex(data_path_vec)
    
    # calculate and save all statistical measures
    BDA.save_statistical_measures(data_path_vec[i], 
                                    10, 
                                    label_vec[i], 
                                    save_path_vec[i];
                                    save_autocovariance_fct = false,
                                    save_spectral_density = false,
                                    save_local_volume_fraction_variance = false,
                                    save_autocovariance_fct_direction = false,
                                    save_spectral_density_along_directions = true)


end


# set paths where statistical data is stored
plot_path = raw"..\plots\nodal_surfaces\\"

# plot all statistical measures
BDA.plot_statistical_measures(save_path_vec,
                    plot_path;
                    spectral_density_xlims = [0,0.1],
                    plot_autocovariance_fct_bool = false,
                    plot_spectral_density_bool = false,
                    plot_local_volume_fraction_variance_bool = false,
                    plot_autocovariance_fct_direction_bool = false,
                    plot_spectral_density_direction_bool = true
                    )



# compare hyperuniformity criterion for red and blue patches
data_path = raw"..\structures\nodal_surfaces\\"

# set surfaces that are analyzed
label_vec = ["D", "G", "P", "I-WP"] 

structure_dict_path_vec = data_path .* label_vec .* "_surface"

voxel_edge_length_vec = [10, 10, 10, 10]

for i in eachindex(structure_dict_path_vec)

    # load dictionary
    structure_dict = BDA.load_h5_dict(structure_dict_path_vec[i]* ".h5")

    # get essential information of data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(structure_dict["data"])

    # add this info to dictionary and save it
    structure_dict = Dict("data_binary" => structure_dict["data"], 
                            "volume_fract_tot" => volume_fract_tot, 
                            "size_data" => collect(size_data), 
                            "mean_edge_length_data" => mean_edge_length_data, 
                            "nr_dimensions_data" => nr_dimensions_data,
                            "voxel_edge_length" => voxel_edge_length_vec[i] ,
                            "label" => label_vec[i],
                            "unit_cell_length" => 500,
                            "nr_unit_cells" => 10,
                            "volume_fraction_parameter" => 0 )

    BDA.save_dict_to_h5(structure_dict; save_path=structure_dict_path_vec[i] * "_structure.h5")

end

# compare hyperuniformity criterion for red and blue patches
data_path = raw"..\structures\nodal_surfaces\\single_unit_cell_"

structure_dict_path_vec = data_path .* label_vec .* "_surface"

for i in eachindex(structure_dict_path_vec)

    # load dictionary
    structure_dict = BDA.load_h5_dict(structure_dict_path_vec[i]* ".h5")

    # get essential information of data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(structure_dict["data"])

    # add this info to dictionary and save it
    structure_dict = Dict("data_binary" => structure_dict["data"], 
                            "volume_fract_tot" => volume_fract_tot, 
                            "size_data" => collect(size_data), 
                            "mean_edge_length_data" => mean_edge_length_data, 
                            "nr_dimensions_data" => nr_dimensions_data,
                            "voxel_edge_length" => voxel_edge_length_vec[i] ,
                            "label" => label_vec[i],
                            "unit_cell_length" => 500,
                            "nr_unit_cells" => 1,
                            "volume_fraction_parameter" => 0 )

    BDA.save_dict_to_h5(structure_dict; save_path=structure_dict_path_vec[i] * "_structure.h5")

end


# set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"]  # 

# set properties of generated data
unit_cell_length = 500
nr_unit_cells = 10

for label in label_vec

    # generate 3d binary data for current nodal surface
    structure_dict = BDA.get_binary_data_from_nodal_eqn(unit_cell_length, 
                                                nr_unit_cells,
                                                label;
    save_result=true, 
    save_path=raw"..\structures\nodal_surfaces\\"*label*"_surface_structure.h5")

    println(label*" done")
    
end


"""
These are the calculations for the stervi beetle from the 10.1002/adfm.202302720 paper.
The data was sent by Viola and is not directly taken from Zenodo
"""


# set raw data path
data_path_raw_prefix = raw"..\structures\stervi\2d_images_green\slice_"
data_path_raw_suffix = "_max11.tif"

data_path_corrected = raw"..\structures\stervi\stervi_green.h5"

# load data, correct voxel size anisotropy and save data
data_binary = BDA.get_binary_data_from_colorscale_stack(data_path_raw_prefix,
                                        data_path_raw_suffix,
                                        301; 
                                        save_data=true, 
                                        data_path_corrected=data_path_corrected)


                                        
# set raw data path
data_path_raw_prefix = raw"..\structures\stervi\2d_images_blue\slice_"
data_path_raw_suffix = "_max11.tif"

data_path_corrected = raw"..\structures\stervi\stervi_blue.h5"

# load data, correct voxel size anisotropy and save data
data_binary = BDA.get_binary_data_from_colorscale_stack(data_path_raw_prefix,
                                        data_path_raw_suffix,
                                        250; 
                                        save_data=true, 
                                        data_path_corrected=data_path_corrected)


# loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    # load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )

    
    # get autocovariance function as a function of sampling distance
    sampling_distance_vec, autocovariance_fct_vec = BDA.get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                                                                    size_data, 
                                                                                    volume_fract_tot,
                                                                                    data_binary;
                                                                                    nr_measurements_per_distance = 10000)

    # create dictionary for current plot
    plot_dict = Dict("sampling_distance_vec" => sampling_distance_vec,
                    "autocovariance_fct_vec" => autocovariance_fct_vec,
                    "voxel_edge_length" => 11,
                    "label" => "green" )

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="stervi_autocovariance_fct_green")


    # calculate spectral density for loaded autocovariance function
    wavenumber_vec, spectral_density_vec = BDA.get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data,
                size_data, 
                volume_fract_tot,
                data_binary;
                nr_measurements_per_distance = 10000,
                nr_wavenumbers=200,
                sampling_distance_vec = sampling_distance_vec,
                autocovariance_fct_vec = autocovariance_fct_vec)

    # create dictionary for current plot
    plot_dict = Dict("wavenumber_vec" => wavenumber_vec,
                    "spectral_density_vec" => spectral_density_vec,
                    "voxel_edge_length" => 11,
                    "label" => "green" )

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
                        save_filename="stervi_spectral_density_green")


    # determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 100,
                                                                                                    window_positioning="random",
                                                                                                    window_shape="spherical")

    # create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "voxel_edge_length" => 11,
                    "label" => "green" )

    
    # save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="stervi_volume_fraction_variance_random_spherical_green")

                                                                                            
end


# compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"..\analysis_data\stervi_volume_fraction_variance_random_spherical_green.h5",
                raw"..\analysis_data\volume_fraction_variance_random_spherical_I-WP_surface.h5"]

# set surfaces that are analyzed
label_vec = ["S. virescens green", "perfect I-WP"]

# set edge lengths
voxel_edge_length_vec = [11, 6]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


# loop through data that will be analyzed
for i in eachindex(data_path_vec)
    
    # load dictionary that contains the following keys:
    # "window_edge_length_vec"
    # "local_volume_fract_variance_vec"
    # "voxel_edge_length"
    # "label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    # adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end

# plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="stervi_volume_fraction_variance_random_spherical_windows_green_i_wp")



# compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"..\analysis_data\stervi_spectral_density_green.h5",
                raw"..\analysis_data\spectral_density_I-WP_surface.h5"]

# set surfaces that are analyzed
label_vec = ["S. virescens green", "perfect I-WP"]

# set edge lengths
voxel_edge_length_vec = [11, 6]

# create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


# loop through data that will be analyzed
for i in eachindex(data_path_vec)
    
    # load dictionary that contains the following keys:
    # "wavenumber_vec"
    # "spectral_density_vec"
    # "voxel_edge_length"
    # "label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    # adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


# plot real part, imaginary part and absolute value of spectral 
# density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_filename="stervi_spectral_density_green_i_wp")


# path of original data
data_path = raw"..\structures\stervi\stervi_green.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\stervi\stervi_green"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "green", 
                                save_path)

# set paths where statistical data is stored
data_path_vec = [raw"..\analysis_data\stervi\stervi_blue",
raw"..\analysis_data\nodal_surfaces\I-WP_surface",
raw"..\analysis_data\stervi\stervi_green"]

# set path where plot will be stored
save_path = raw"..\plots\stervi\stervi_blue_green_i_wp"

# plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            voxel_edge_length_vec=[11,6,11],
            label_vec=["S. virescens blue", "perfect I-WP", "S. virescens green"]
            )



# path of original data
data_path = raw"..\structures\stervi\stervi_green.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\stervi\stervi_green"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "green", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)


# path of original data
data_path = raw"..\structures\stervi\stervi_blue.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\stervi\stervi_blue"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "blue", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)



# path of original data
data_path = raw"..\structures\stervi\stervi_blue.h5"

save_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"

# get essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

# get autocovariance function array
sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(size_data,
                       volume_fract_tot,
                       data_binary,
                       nr_measurements_per_direction=10)

# create dict to save
saving_dict = Dict("sampling_vec_array" => sampling_vec_array,
                    "autocovariance_fct_array" => autocovariance_fct_array,
                    "voxel_edge_length" => 11,
                    "label" => "stervi blue")

save_dict_to_h5(saving_dict; save_path)


# path where autocovariance fct data is saved
dict_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"

# load dict
data_dict = BDA.load_h5_dict(dict_path)

# get sampling vec array and autocovariance fct array from dict
sampling_vec_array = data_dict["sampling_vec_array"]
autocovariance_fct_array = data_dict["autocovariance_fct"]

# calculate spectral density
sampled_wavevectors_array, spectral_density_array = get_spectral_density(size_data, 
                                                volume_fract_tot,
                                                data_binary;
                                                sampling_vec_array = sampling_vec_array,
                                                autocovariance_fct_array)

# path where spectral density is saved
save_path = raw"..\analysis_data\stervi\stervi_spectral_density_direction_blue_small_sampling.h5"

# create dict to save
saving_dict = Dict("sampled_wavevectors_array" => sampled_wavevectors_array,
                    "spectral_density_array" => spectral_density_array,
                    "voxel_edge_length" => 11,
                    "label" => "stervi blue")

save_dict_to_h5(saving_dict; save_path)



# set data path
data_path = raw"..\structures\stervi\stervi_blue.h5"

# get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)


# path of dict
load_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"

# load dict
data_dict = BDA.load_h5_dict(load_path)

autocovariance_fct_array = data_dict["autocovariance_fct_array"]

# calculate spectral density
sampled_wavenumbers_vec_vec, sampled_wavevectors_array, spectral_density_array = BDA.get_spectral_density(size_data, 
volume_fract_tot,
data_binary;
nr_measurements_per_direction = 50,
autocovariance_fct_array = autocovariance_fct_array)


# create dictionary for current plot
plot_dict = Dict("sampled_wavenumbers_vec_vec" => sampled_wavenumbers_vec_vec,
                "sampled_wavevectors_array" => sampled_wavevectors_array,
                "spectral_density_array" => spectral_density_array,
                "voxel_edge_length" => 11,
                "label" => data_dict["label"])


# path of dict
save_path = raw"..\analysis_data\stervi\stervi_spectral_density_direction_blue_small_sampling.h5"

# save the plot_dict to a H5 file
BDA.save_dict_to_h5(copy(plot_dict); save_path=save_path)



# plot heat map of spectral density in x-y-plane for k_z=0
BDA.plot_spectral_density_heatmap(plot_dict,
                                save_path;
                                save_plot = false)

# wait until key is pressed
readline()




autocovariance_fct_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"
autocovariance_fct_dict = BDA.load_h5_dict(autocovariance_fct_path)

spectral_density_path = raw"..\analysis_data\stervi\stervi_spectral_density_direction_blue_small_sampling.h5"
spectral_density_dict = BDA.load_h5_dict(spectral_density_path)


# set data path
data_path = raw"..\structures\stervi\stervi_blue.h5"

# get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)

sampling_distance_vec_vec = BDA.get_sampling_distance_vec_vec(size_data)

autocovariance_fct_dict["sampling_distance_vec_1"] = sampling_distance_vec_vec[1]
autocovariance_fct_dict["sampling_distance_vec_2"] = sampling_distance_vec_vec[2]
autocovariance_fct_dict["sampling_distance_vec_3"] = sampling_distance_vec_vec[3]

# save the plot_dict to a H5 file
BDA.save_dict_to_h5(copy(autocovariance_fct_dict); save_path=autocovariance_fct_path)



# set data path
data_path = raw"..\structures\stervi\stervi_blue.h5"

# set path where autocovariance dict will be stored
save_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue.h5"

# get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)

# get autocovariance function array
sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(
                    size_data,
                    volume_fract_tot,
                    data_binary,
                    nr_measurements_per_direction=100)

# create dict to save
saving_dict = Dict("sampling_vec_array" => sampling_vec_array,
                    "sampling_distance_vec_vec" => sampling_distance_vec_vec,
                    "autocovariance_fct_array" => autocovariance_fct_array,
                    "voxel_edge_length" => 11,
                    "label" => "stervi blue")

BDA.save_dict_to_h5(saving_dict; save_path)


# set data path
data_path = raw"..\structures\stervi\stervi_blue.h5"

# set path where autocovariance dict will be stored
dict_path = raw"..\analysis_data\stervi\stervi_autocovariance_fct_direction_blue.h5"

# get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)

# load autocovariance array
loaded_dict = BDA.load_h5_dict(dict_path)

plot_path = raw"..\plots\stervi\stervi_direction_blue"


BDA.plot_autocovariance_fct_heatmap(loaded_dict,
    plot_path;
    save_plot = false,
    sampling_vector_component_to_fix = 3,
    sampling_vector_value_fixed = 0 )


    
# set data path
data_path_vec = [raw"..\structures\stervi\stervi_blue.h5",
raw"..\structures\stervi\stervi_green.h5",
raw"..\structures\pachy\pachy_blue.h5",
raw"..\structures\pachy\pachy_red.h5"
]

# set path where autocovariance dict will be stored
save_path_vec = [raw"..\analysis_data\stervi\stervi_blue_autocovariance_fct_direction.h5",
raw"..\analysis_data\stervi\stervi_green_autocovariance_fct_direction.h5",
raw"..\analysis_data\pachy\pachy_blue_autocovariance_fct_direction.h5",
raw"..\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"
]

voxel_edge_length_vec = [11,11,10,9]

label_vec = ["stervi blue", "stervi green", "pachy blue","pachy red"]

for i in eachindex(data_path_vec)

    # get data and essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_vec[i])

    # get autocovariance function array
    sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(
                        size_data,
                        volume_fract_tot,
                        data_binary,
                        nr_measurements_per_direction=1000)

    # create dict to save
    saving_dict = Dict("sampling_vec_array" => sampling_vec_array,
                        "sampling_distance_vec_vec" => sampling_distance_vec_vec,
                        "autocovariance_fct_array" => autocovariance_fct_array,
                        "voxel_edge_length" => voxel_edge_length_vec[i],
                        "label" => label_vec[i])

    BDA.save_dict_to_h5(saving_dict; save_path=save_path_vec[i])

    println(label_vec[i]*" done")

end


# set data path
data_path_vec = [raw"..\analysis_data\stervi\stervi_blue",
                raw"..\analysis_data\stervi\stervi_green"]

# set path to save plot
save_path = raw"..\plots\stervi\\"

BDA.plot_statistical_measures(data_path_vec,
            save_path;
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = true
            )

            
# path of original data
data_path = raw"..\structures\stervi\stervi_blue.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\stervi\stervi_blue"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "S. virescens blue", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)

# path of original data
data_path = raw"..\structures\stervi\stervi_green.h5"

# path where analysis data will be saved
save_path = raw"..\analysis_data\stervi\stervi_green"
    
# calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "S. virescens green", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)


# set paths where statistical data is stored
data_path_vec = [raw"..\analysis_data\stervi\stervi_green",
                raw"..\analysis_data\stervi\stervi_blue"]

# set path where plot will be stored
save_path = raw"..\plots\stervi\\"

# plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            spectral_density_xlims = [0,0.1],
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = false,
            plot_spectral_density_direction_bool = true
            )


"""
The following commands targeted multiple samples at the same time
"""


data_path_vec = [ raw"..\structures\pachy\pachy_red_structure.h5",
raw"..\structures\pachy\pachy_blue_structure.h5",
raw"..\structures\stervi\stervi_green_structure.h5",
raw"..\structures\stervi\stervi_blue_structure.h5",
raw"..\structures\nodal_surfaces\D_surface_structure.h5",
raw"..\structures\nodal_surfaces\G_surface_structure.h5",
raw"..\structures\nodal_surfaces\P_surface_structure.h5",
raw"..\structures\nodal_surfaces\I-WP_surface_structure.h5"
]

save_path_vec = [ raw"..\analysis_data\pachy\pachy_red",
raw"..\analysis_data\pachy\pachy_blue",
raw"..\analysis_data\stervi\stervi_green",
raw"..\analysis_data\stervi\stervi_blue",
raw"..\analysis_data\nodal_surfaces\D_surface",
raw"..\analysis_data\nodal_surfaces\G_surface",
raw"..\analysis_data\nodal_surfaces\P_surface",
raw"..\analysis_data\nodal_surfaces\I-WP_surface"
]

for i in eachindex(data_path_vec)

    BDA.save_statistical_measures(data_path_vec[i],
        save_path_vec[i];
        save_autocovariance_fct = false,
        save_spectral_density = true,
        save_local_volume_fraction_variance = false,
        save_autocovariance_fct_direction = false,
        save_spectral_density_along_directions = false,
        save_complete_autocovariance_fct_direction = true,
        save_spectral_density_array = true)

    println(string(i)*" done")

end


data_path_vec_vec =  [ [raw"..\analysis_data\pachy\pachy_blue",
    raw"..\analysis_data\pachy\pachy_red",
    raw"..\analysis_data\nodal_surfaces\D_surface"],
[
raw"..\analysis_data\stervi\stervi_blue",
raw"..\analysis_data\nodal_surfaces\I-WP_surface",
raw"..\analysis_data\stervi\stervi_green"],
[ raw"..\analysis_data\nodal_surfaces\P_surface",
    raw"..\analysis_data\nodal_surfaces\I-WP_surface",
    raw"..\analysis_data\nodal_surfaces\D_surface",
raw"..\analysis_data\nodal_surfaces\G_surface"
]
] 

save_path_vec = [ raw"..\plots\pachy\pachy_blue_red_d",
raw"..\plots\stervi\stervi_blue_green_i_wp",
raw"..\plots\nodal_surfaces\nodal_surfaces"
]

voxel_edge_length_vec = [
    [10, 9, 8.5],
    [11,6,11],
    nothing
]

for i in eachindex(data_path_vec_vec)

    BDA.plot_statistical_measures(data_path_vec_vec[i],
        save_path_vec[i];
        voxel_edge_length_vec=voxel_edge_length_vec[i],
        label_vec=nothing,
        spectral_density_xlims = (0,0.1),
        spectral_density_heatmaps_clims = nothing,
        plot_autocovariance_fct_bool = false,
        plot_spectral_density_bool = true,
        plot_local_volume_fraction_variance_bool = false,
        plot_autocovariance_fct_direction_bool = false,
        plot_spectral_density_along_directions_bool = false,
        plot_spectral_density_heatmaps_bool = false
        )

    println(string(i)*" done")

end


# set raw data path
data_path_raw_prefix = raw"..\structures\stervi\2d_images_green\slice_"
data_path_raw_suffix = "_max11.tif"

data_path_corrected = raw"..\structures\stervi\stervi_green.h5"

structure_dict = BDA.get_structure_dict_from_colorscale_stack(data_path_raw_prefix,
    data_path_raw_suffix,
    301; 
    voxel_size=(11,11,11), 
    label = "S. virescens green",
    save_result=true, 
    save_path = raw"..\structures\stervi\stervi_green")



"""
This is where commands for network generation begin
"""


# import my module that contains all functions for the analysis of binary structure data
import .NetworkGeneration as NG
import MetaGraphsNext as MGN


network_dict = NG.get_periodic_network( ; nr_vertices = 150 , 
                            nr_dimensions = 3, 
                            network_type = "diamond")

network_dict["bond_bending_const"] = 0.285

# get list of bonds (edges)
edges_vec = collect(MGN.edge_labels(network_dict["network_graph"]))

# break a random bond
network_dict = NG.switch_bond!(network_dict, edges_vec[4] )

vertex = 1

neighbor_positions_mat = NG.get_neighbor_positions_mat(network_dict, vertex)

local_energy = NG.local_energy_keating(network_dict["network_graph"][vertex], 
        network_dict, neighbor_positions_mat)

gradient = zeros(3)

NG.gradient_keating!(gradient, network_dict["network_graph"][vertex], 
        network_dict, neighbor_positions_mat)

hessian = zeros(3, 3)

NG.hessian_keating!(hessian, network_dict["network_graph"][vertex], 
        network_dict, neighbor_positions_mat)



println(local_energy)
println(gradient)
println(hessian)



graph_dict = NG.get_periodic_network( ; nr_vertices = 150 ,
                            network_type = "diamond")

                            
vertex = 5

# get and print neighbors
vertex_neighbors = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex))

println(vertex_neighbors)

# get position of some vertex
println(graph_dict["spatial_network"][vertex]["position"] )
println(graph_dict["spatial_network"][vertex]["local_energy"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex, 
[1,-0.5,0.3] )

println(graph_dict["spatial_network"][vertex]["position"] )
println(graph_dict["spatial_network"][vertex]["local_energy"] )


# relax vertex
graph_dict = NG.relax_single_vertex_keating!(graph_dict, vertex)

println(graph_dict["spatial_network"][vertex]["position"] )
println(graph_dict["spatial_network"][vertex]["local_energy"] )

graph_dict = NG.get_periodic_network( ; nr_vertices = 1500 ,
                            network_type = "diamond")

                            
vertex_vec = [5,12]

# get and print neighbors
neighbors_in_shells_dict, all_vertices_vec = NG.get_neighbors_in_shells_dict(graph_dict, 
                                    vertex_vec; 
                                    shell_nr = 4)


vertex_vec = [5,12]


# get position of some vertex
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[1], 
[1,-0.5,0.3] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[2], 
[-0.1,0.5,1.5] )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

graph_dict = NG.relax_cluster_one_cycle_keating!(graph_dict, 
vertex_vec; 
shell_nr = 4 )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

graph_dict = NG.relax_cluster_one_cycle_keating!(graph_dict, 
vertex_vec; 
shell_nr = 4 )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )




graph_dict = NG.get_periodic_network( ; nr_vertices = 1500 ,
                            network_type = "diamond")

                            
vertex_vec = [5,12]


# get position of some vertex
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[1], 
[1,-0.5,0.3] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[2], 
[-0.1,0.5,1.5] )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

@time graph_dict = NG.relax_cluster_keating!(graph_dict, vertex_vec; 
nr_cycles = 10,
reject_during_relaxation_cycle_threshold = 5,
shell_nr = 4)

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )


                        
vertex_vec = [5,12]


# get position of some vertex
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[1], 
[1,-0.5,0.3] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[2], 
[-0.1,0.5,1.5] )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )

graph_dict = NG.relax_cluster_one_cycle_keating!(graph_dict, 
vertex_vec; 
shell_nr = 4 )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )

graph_dict = NG.relax_cluster_one_cycle_keating!(graph_dict, 
vertex_vec; 
shell_nr = 4 )

println("___")
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )



vertex_to_relax = 5

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_to_relax, 
[0.2,0.5,-0.9] )

# get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

# get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

# get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

# set energy, gradient and hessian for energy minimization
energy(x) = NG.energy_from_position_keating(x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr )
                                            
gradient!(gradient, x) = NG.gradient_keating!(gradient, x, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)
hessian!(hessian, x) = NG.hessian_keating!(hessian, x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

gradient_an = zeros(3)
gradient!(gradient_an, initial_position)
gradient_num = ForwardDiff.gradient(energy, initial_position)

println(gradient_an)
println(gradient_num)
println("_______")

hessian_an = zeros(3, 3)
hessian!(hessian_an, initial_position)
hessian_num = ForwardDiff.hessian(energy, initial_position)

println(hessian_an)
println(hessian_num)



println(graph_dict["spatial_network"][vertex]["position"])

neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex;
                                    exclude_vertices = [])

next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex)

neighbors_vec = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex))
println(neighbors_vec)

neighbor_nr = 1

next_neighbors_vec = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], neighbors_vec[neighbor_nr]))

for i in 1:4
    println(graph_dict["spatial_network"][next_neighbors_vec[i]]["position"])

end

println("_______")

for i in 1:3
    println(next_neighbor_positions_arr[neighbor_nr,:,i])

end



println("_______")

vertex = 2

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex, 
[1,-0.5,0.3]; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )


# relax vertex
graph_dict = NG.relax_single_vertex_keating!(graph_dict, vertex; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )

println("_______")

vertex = 3

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )

# move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex, 
[1,-0.5,0.3]; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )


# relax vertex
graph_dict = NG.relax_single_vertex_keating!(graph_dict, vertex; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )



neighbors_vec = collect(MetaGraphsNext.neighbor_labels(
                                        graph_dict["spatial_network"], vertex))

for neighbor in neighbors_vec
    display(graph_dict["spatial_network"][neighbor]["position"])
end
println("_______")

neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex)

display(neighbor_positions_mat)

next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex)

display(next_neighbor_positions_arr[1,:,:])



random_bond = NG.get_random_bond(graph_dict)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][random_bond[1]]["position"] )

temperature = 100.5

graph_dict, move_accepted = NG.monte_carlo_move(graph_dict, 
    temperature; 
    switched_bond = random_bond)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][random_bond[1]]["position"] )


central_vertex = 100

central_vertex_position = graph_dict["spatial_network"][central_vertex]["position"]

    # get central vertices neighbors 
    neighbor_vec = collect(MetaGraphsNext.neighbor_labels(
                                graph_dict["spatial_network"], central_vertex))
# create array to store next to nearest neighbors coordinates
# The first array index labels the number of the direct neighbor
next_neighbor_positions_arr = Array{Float64}(undef, 
                                            graph_dict["coordination_nr"],
                                            graph_dict["nr_dimensions"],
                                            graph_dict["coordination_nr"]-1)

# loop through central vertices neighbors
for i in 1:graph_dict["coordination_nr"]
    current_next_neighbor = 1
    # loop through the current neighbor's neighbors
    for next_neighbor in MetaGraphsNext.neighbor_labels(
                                    graph_dict["spatial_network"], neighbor_vec[i])
        if next_neighbor !== central_vertex
            # get next neighbor's virtual coordinates which might be outside of the 
            # supercell if periodic boundary conditions play a role
            next_neighbor_positions_arr[i,:,current_next_neighbor] = NG.get_virtual_position(
                        central_vertex_position,
                        graph_dict["spatial_network"][next_neighbor]["position"],
                        graph_dict["supercell_edge_length"] )
            current_next_neighbor += 1
        end
    end
end



# get cluster after bond switch
cluster_dict = NG.get_cluster_in_shells_dict(
                                graph_dict, 
                                switched_bond; 
                                shell_nr = shell_nr)

# relax cluster around switched bond
graph_dict = NG.relax_cluster_keating!(graph_dict,
    switched_bond; 
    nr_cycles = nr_cycles,
    reject_during_relaxation_cycle_threshold = reject_during_relaxation_cycle_threshold,
    shell_nr = shell_nr,
    cluster_dict = cluster_dict,
    initial_cluster_energy  = initial_cluster_dict["cluster_energy"],
    update_total_energy = true)

    
switched_bond = (5,12)

display(collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[1]) ))
    
display(collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[2]) ))

println("______")

graph_dict = NG.switch_bond!(graph_dict, switched_bond)

display(collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[1]) ))

display(collect(MetaGraphsNext.neighbor_labels(
            graph_dict["spatial_network"], switched_bond[2]) ))


random_bond = NG.get_random_bond(graph_dict)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][random_bond[1]]["position"] )

temperature = 100.5

graph_dict, move_accepted = NG.monte_carlo_move(graph_dict, 
    temperature; 
    switched_bond = random_bond)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][random_bond[1]]["position"] )


graph_dict = NG.get_periodic_network( ; nr_vertices = 40 , 
nr_dimensions = 3, 
network_type = "diamond",
bond_bending_const = 0.285)

figure = NG.plot_network(graph_dict)

# figure = NG.plot_network(graph_dict)

switched_bond = (5,12) # NG.get_random_bond(graph_dict)

temperature = 100.5

graph_dict, move_accepted, new_bond_vec = NG.monte_carlo_move(graph_dict, 
    temperature; 
    switched_bond = switched_bond)

figure = NG.plot_network(graph_dict; highlight_nodes = switched_bond, highlight_edges = [switched_bond, new_bond_vec...])

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network(graph_dict,
    10, 
    temperature;
    shell_nr = 3,
    print_progress = true)

    
NG.save_graph_to_csv(graph_dict,
"example_graph.csv")


@time graph_dict, cluster_energy_vec = NG.relax_cluster_keating!(graph_dict,
random_bond; 
relax_efficiently = true,
    update_total_energy = true,
    track_cluster_energy = true)

Plots.plot(collect(1:26), cluster_energy_vec)


# compare the gradient in the inefficient and efficient calculation
vertex_to_relax = random_bond[2]

# efficient calculation
gradient_eff = NG.gradient_keating_efficient(graph_dict, vertex_to_relax)

# inefficient calculation

# get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

# get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

# get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

# set energy, gradient and hessian for energy minimization

gradient_ineff = zeros(3)
                                        
NG.gradient_keating!(gradient_ineff, initial_position, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

println("Efficient: "*string(gradient_eff))
println("Inefficient "*string(gradient_ineff))


random_bond = NG.get_random_bond(graph_dict)

graph_dict, new_bond_vec = NG.switch_bond!(graph_dict, random_bond )


# compare the gradient in the inefficient and efficient calculation
vertex_to_relax = random_bond[2]

# efficient calculation
gradient = NG.gradient_keating_efficient(graph_dict, vertex_to_relax)

translation_vector_eff =  NG.get_approximate_translation_vector_keating(gradient, 
        graph_dict["bond_bending_const"];
        relaxation_overshoot_factor_r = 1.5,
        relaxation_optimization_parameter_l = 1)

# inefficient calculation

# get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

# get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

# get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

# set energy, gradient and hessian for energy minimization
energy(x) = NG.energy_from_position_keating(x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr )
                                            
gradient!(gradient, x) = NG.gradient_keating!(gradient, x, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

hessian!(hessian, x) = NG.hessian_keating!(hessian, x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)
# find energy minimum
minimizer_result = Optim.optimize(
                            energy, 
                            gradient!, 
                            hessian!,
                            initial_position, 
                            Optim.Newton())

# get relaxed position and local keating energy
relaxed_position = Optim.minimizer(minimizer_result)

# calculate translation vector for relaxed vertex
translation_vector_ineff = relaxed_position .- initial_position

println("Efficient: "*string(translation_vector_eff))
println("Inefficient "*string(translation_vector_ineff))


@time graph_dict, cluster_energy_vec = NG.relax_cluster_keating!(graph_dict,
random_bond; 
nr_max_relaxation_cycles = 1,
relax_efficiently = true,
    update_total_energy = true,
    track_cluster_energy = true,
    relaxation_optimization_parameter_l=1.3)

figure = NG.plot_network(graph_dict)


random_bond = NG.get_random_bond(graph_dict)

vertex_to_relax = random_bond[2]

graph_dict, new_bond_vec = NG.switch_bond!(graph_dict, random_bond )


# efficient calculation
gradient = NG.gradient_keating_efficient(graph_dict, vertex_to_relax)

hessian = NG.hessian_keating_efficient(graph_dict, vertex_to_relax)

# calculate translation vector to approximate energy minimum
translation_vector_eff = .- LinearAlgebra.inv(hessian)*gradient

# inefficient calculation

# get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

# get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

# get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

# set energy, gradient and hessian for energy minimization
energy(x) = NG.energy_from_position_keating(x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr )
                                            
gradient!(gradient, x) = NG.gradient_keating!(gradient, x, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

hessian!(hessian, x) = NG.hessian_keating!(hessian, x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)
# find energy minimum
minimizer_result = Optim.optimize(
                            energy, 
                            gradient!, 
                            hessian!,
                            initial_position, 
                            Optim.Newton())

# get relaxed position and local keating energy
relaxed_position = Optim.minimizer(minimizer_result)

# calculate translation vector for relaxed vertex
translation_vector_ineff = relaxed_position .- initial_position

println("Efficient: "*string(translation_vector_eff))
println("Inefficient "*string(translation_vector_ineff))



# efficient calculation
neighbor_vec = collect(
        MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], central_vertex))

j=1

# get vector pointing from central vertex to neighbor j
distance_vector_j_eff = (sign(neighbor_vec[j] - central_vertex)
* graph_dict["spatial_network"][central_vertex, neighbor_vec[j]]["vector"])

bond_stretching_term_eff = ( - 3/4 * ( 
            graph_dict["spatial_network"][central_vertex, neighbor_vec[j]]["distance_squared"] - 1 
            ) ) .* distance_vector_j_eff

# inefficient calculation

# get initial position of vertex to relax 
x = graph_dict["spatial_network"][vertex_to_relax]["position"]

# get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

# get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

# set energy, gradient and hessian for energy minimization

# get vector pointing from central vertex to neighbor
distance_vector_j = neighbor_positions_mat[:,j] .- x

# get bond stretching term
bond_stretching_term_ineff = ( - 3/4 * ( LinearAlgebra.norm(distance_vector_j)^2 - 1 ) 
                                ) .* distance_vector_j

println("Efficient: "*string(distance_vector_j_eff))
println("Inefficient "*string(distance_vector_j))



switched_bond = NG.get_random_bond(graph_dict, seed = 1)

graph_dict, new_bond_vec = NG.switch_bond!(graph_dict, switched_bond )

vertexic_position_arr, cluster_energy_arr = NG.compare_relaxation_methods(graph_dict,
switched_bond,
"diamond_disordered_t10_30steps_1" )

temperature = 8

graph_dict, move_accepted, new_bond_vec = NG.monte_carlo_move!(graph_dict, 
    temperature; 
    reject_during_relaxation_cycle_threshold = 10,
        break_at_relative_cluster_energy_change = 0.0001,
        shell_nr = 3,
        relax_efficiently = true,
        thermal_fluctuations = false)

println(graph_dict["total_energy"])
println(NG.get_total_energy_keating(graph_dict))


switched_bond = NG.get_random_bond(graph_dict, seed = 9)
shell_nr = 4

initial_cluster_dict = NG.get_cluster_in_shells_dict(
                graph_dict, 
                switched_bond; 
                shell_nr = shell_nr)


# switch bond
graph_dict, new_bond_vec = NG.switch_bond!(graph_dict, switched_bond )

# get cluster after bond switch
cluster_dict = NG.get_cluster_in_shells_dict(
                graph_dict, 
                switched_bond; 
                shell_nr = shell_nr)

# relax cluster once and update cluster energy
graph_dict, cluster_dict = NG.relax_cluster_keating!(graph_dict,
cluster_dict; 
nr_max_relaxation_cycles = 25,
break_at_relative_cluster_energy_change = 0.001,
reject_during_relaxation_cycle_threshold = 10,
relax_efficiently = true,
update_total_energy = true)

# calculate new total energy and compare to actual total energy
smart_total_energy = graph_dict["total_energy"] 

actual_total_energy = NG.get_total_energy_keating(graph_dict)

println(string(smart_total_energy))
println(string(actual_total_energy))


temperature = 5
shell_nr = 3


graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network(graph_dict,
    20, 
    temperature; 
    nr_max_relaxation_cycles = 25,
        break_at_relative_cluster_energy_change = 0.001,
        reject_during_relaxation_cycle_threshold = 10,
        relax_efficiently = true,
        shell_nr = shell_nr,
    print_progress = true,
    random_evolution_seed = 3,
    thermal_fluctuations = false)

# calculate new total energy and compare to actual total energy
smart_total_energy = graph_dict["total_energy"] 

actual_total_energy = NG.get_total_energy_keating(graph_dict)

println(string(smart_total_energy))
println(string(actual_total_energy))


structure_factor_dict = NA.get_structure_factor_isotrope_by_wavenumber_vec(
        graph_dict)

hyperuniformity_parameter = NA.get_effective_hyperuniformity_parameter(structure_factor_dict)
println("hyperuniformity parameter: "*string(hyperuniformity_parameter))

local_nr_variance_dict = NA.get_local_nr_variance_by_window_radius_vec(graph_dict;
structure_factor_dict = structure_factor_dict)

Plots.plot(local_nr_variance_dict["window_radius_vec"], 
local_nr_variance_dict["local_nr_variance_vec"] ./ local_nr_variance_dict["window_radius_vec"].^3)

NG.plot_network(graph_dict)


graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence(graph_dict,
        evolution_dict;
        print_progress = true,
        print_every_nr_attempted_bond_switches = 50)

evolution_dict["temperature_vec"] = [0]
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = [30]

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence(graph_dict,
        evolution_dict;
        move_accepted_vec = move_accepted_vec,
        total_energy_vec = total_energy_vec,
        print_progress = true,
        print_every_nr_attempted_bond_switches = 50)

Plots.plot(collect(1:length(total_energy_vec)), total_energy_vec)

NG.save_graph_to_h5_and_MGformat(deepcopy(graph_dict), "example")


evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [0.1],
nr_monte_carlo_steps_per_temperature_vec = [1])

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 20)

Plots.plot(collect(1:length(total_energy_vec)), total_energy_vec)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "64_vertices_T_0.1"

save_path = raw"..\structures\random_networks\\"

NG.save_mesh_from_network(graph_dict_to_save, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict_to_save,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)
        


network_names = ["64_vertices_T_0.1", "1000_vertices_T_1_quenched", "1000_vertices_T_2_quenched", "1000_vertices_T_4_quenched"]

save_path = raw"..\structures\random_networks\\"

for name in network_names
    graph_dict_to_save = NG.load_graph_from_h5_and_MGformat(save_path*name)

    NG.save_mesh_from_network(graph_dict_to_save, name; save_path = save_path)

end


evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 ,temperature_vec = [0.25, 0],
nr_monte_carlo_steps_per_temperature_vec = [3, 50])

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
        move_accepted_vec=move_accepted_vec,
        total_energy_vec=total_energy_vec,
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

Plots.plot(collect(1:length(total_energy_vec)), total_energy_vec)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "1000_vertices_T_0.25_quenched"

save_path = raw"..\structures\random_networks\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)

        
nr_samples = 10

evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 )

evolution_dict["temperature_vec"] = zeros(nr_samples) .+ 1
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = ones(Int64, nr_samples)
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][1] = 3 

graph_dict = NG.get_periodic_network(evolution_dict)


graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200,
    save_network_after_each_step = true,
    filename = "1000_vertices_T_1",)

nr_samples = 5

evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000)

evolution_dict["temperature_vec"] = zeros(nr_samples) .+ 0.5
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = ones(Int64, nr_samples)
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][1] = 3 

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200,
    save_network_after_each_step = true,
    filename = "1000_vertices_T_0.5",)

evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000)

evolution_dict["temperature_vec"] = zeros(nr_samples) .+ 0.015625
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = ones(Int64, nr_samples)
evolution_dict["nr_monte_carlo_steps_per_temperature_vec"][1] = 3 

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200,
    save_network_after_each_step = true,
    filename = "1000_vertices_T_0.015625",)


load_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

load_name = "1000_vertices_T_1_quenched"

# load dictionary with 1000 vertices which was heated to T=1 and then quenched
graph_dict = NG.load_graph_from_h5_and_MGformat(load_path*load_name)

save_path = raw"..\structures\random_networks\\"

NG.save_mesh_from_network(graph_dict, load_name*"_thick_bonds"; save_path = load_path, bond_radius = 0.3131)

structure_factor_dict = NA.get_structure_factor_isotrope_by_wavenumber_vec(
        graph_dict)

network_names = [ "1000_vertices_T_1_quenched", "1000_vertices_T_4_quenched"]

load_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

graph_dict_1 = NG.load_graph_from_h5_and_MGformat(load_path*network_names[1])

structure_factor_dict_1 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict_1;
    sampling_distance_step_length = 0.025,
    maximal_sampling_distance = 4*graph_dict_1["supercell_edge_length"],
    save_result = false,
    save_path = raw"..\analysis_data\random_networks\1000_vertices_T_1_quenched_high_sampling_rate",
    label = nothing)


graph_dict_4 = NG.load_graph_from_h5_and_MGformat(load_path*network_names[2])

structure_factor_dict_4 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict_4;
    sampling_distance_step_length = 0.025,
    maximal_sampling_distance = 4*graph_dict_4["supercell_edge_length"],
    save_result = false,
    save_path = raw"..\analysis_data\random_networks\1000_vertices_T_4_quenched_high_sampling_rate",
    label = nothing)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [1, 0],
nr_monte_carlo_steps_per_temperature_vec = [1, 50], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "216_vertices_T_1_quenched"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 512 ,temperature_vec = [1, 0],
nr_monte_carlo_steps_per_temperature_vec = [1, 50], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "512_vertices_T_1_quenched"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)


Plots.plot(collect(1:length(total_energy_vec)) ./ (512*18), total_energy_vec, xlabel = "steps", ylabel = "total energy", label = "ylabel")


# path where structures are stored
load_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

# load graph with 512 vertices
graph_dict_5 = NG.load_graph_from_h5_and_MGformat(load_path*"512_vertices_T_1_quenched")

# determine structure factor
structure_factor_dict_5 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict_5;
    sampling_distance_step_length = 0.025,
    maximal_sampling_distance = 4*graph_dict_5["supercell_edge_length"],
    save_result = true,
    save_path = raw"..\analysis_data\random_networks\512_vertices_T_1_quenched_high_sampling_rate",
    label = "512_vertices_T_1_quenched_high_sampling_rate")


# load graph with 216 vertices
graph_dict_2 = NG.load_graph_from_h5_and_MGformat(load_path*"216_vertices_T_1_quenched")

# determine structure factor
structure_factor_dict_2 = NA.get_structure_factor_isotrope_by_wavenumber_vec(
    graph_dict_2;
    sampling_distance_step_length = 0.025,
    maximal_sampling_distance = 4*graph_dict_2["supercell_edge_length"],
    save_result = true,
    save_path = raw"..\analysis_data\random_networks\216_vertices_T_1_quenched_high_sampling_rate",
    label = "216_vertices_T_1_quenched_high_sampling_rate")


# load structure factor for 1000 vertices
dict_path_1 = raw"..\analysis_data\random_networks\1000_vertices_T_1_quenched_high_sampling_rate_structure_factor_isotrope.h5"

structure_factor_dict_1 = GU.load_h5_dict(dict_path_1)

dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.125, 0.25, 0.5, 1, 2, 4, 6, 8]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [2, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end




dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 6, 8]

for i in eachindex(temperatures)

    graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*"216_vertices_T_"*string(temperatures[i])*"_quenched")

    structure_factor_dict = NA.get_structure_factor_bartlett_isotrope_by_wavenumber_vec(
        graph_dict;
        sampling_distance_step_length = 0.05,
        maximal_sampling_distance = 4*graph_dict["supercell_edge_length"],
        save_result = true,
        save_path = raw"..\analysis_data\random_networks\216_vertices_T_"*string(temperatures[i])*"_quenched",
        label = "T = "*string(temperatures[i]))

end


evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [1],
nr_monte_carlo_steps_per_temperature_vec = [0.01], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "64_vertices_slight_disorder"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path, bond_radius = 0.3131)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [1],
nr_monte_carlo_steps_per_temperature_vec = [0.01], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

filename = "64_vertices_perfect_diamond_thick_bonds"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path, bond_radius = 0.3131)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [0.5, 0],
nr_monte_carlo_steps_per_temperature_vec = [0.1, 30], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
        total_energy_vec = total_energy_vec,
        move_accepted_vec = move_accepted_vec,
        print_progress = true,
    print_every_nr_attempted_bond_switches = 50)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "64_vertices_T_0.5_quenched_thick_bonds"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path, bond_radius = 0.3131)

NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)


evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = [1],
nr_monte_carlo_steps_per_temperature_vec = [0.03], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)

graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
        total_energy_vec = total_energy_vec,
        move_accepted_vec = move_accepted_vec,
        print_progress = true,
    print_every_nr_attempted_bond_switches = 50)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)

filename = "64_vertices_slightly_more_disorder_thick_bonds"

save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path, bond_radius = 0.3131)

NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.125, 0.25, 0.5, 1, 2, 4, 8]

evolution_dict = GU.load_h5_dict(dict_path*"216_vertices_T_"*string(temperatures[7])*"_quenched_evolution.h5")

Plots.plot(collect(1:length(evolution_dict["total_energy_vec"])) ./ (216*18), evolution_dict["total_energy_vec"], xlabel = "steps", ylabel = "total energy", label = "ylabel")




dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.1, 0.15, 0.2, 0.3, 0.4]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [2, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end


temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_cooling_gradient(temperature;
    temperature_decrease_per_monte_carlo_step = 0.1,
    quench = true )

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_cool_0.1_per_mc_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
    temperature_increase_per_monte_carlo_step = 0.1,
    quench = true )

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_cool_0.1_per_mc_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [0.01, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_heated_for_0.01_steps_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end


temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [0.05, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_heated_for_0.05_steps_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [5, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_heated_for_5.0_steps_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



temperatures = [0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5]

for temperature in temperatures

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
    nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_heated_for_10.0_steps_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end


central_vertex = 10

l_max = 12

evolution_dict = NA.get_evolution_dict(nr_vertices = 64, network_type = "diamond")
graph_dict_diamond = NG.get_periodic_network(evolution_dict)

single_vertex_q_l_diamond = NA.get_q_l_averaged_single_vertex_dict(graph_dict_diamond,
central_vertex,
l_max)

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, single_vertex_q_l_diamond[i])
    println()
end

evolution_dict = NA.get_evolution_dict(nr_vertices = 16, network_type = "bcc")
graph_dict_bcc = NG.get_periodic_network(evolution_dict)

single_vertex_q_l_bcc = NA.get_q_l_averaged_single_vertex_dict(graph_dict_bcc,
central_vertex,
l_max)

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, single_vertex_q_l_bcc[i])
    println()
end

evolution_dict = NA.get_evolution_dict(nr_vertices = 16, network_type = "fcc")
graph_dict_fcc = NG.get_periodic_network(evolution_dict)

single_vertex_q_l_fcc = NA.get_q_l_averaged_single_vertex_dict(graph_dict_fcc,
central_vertex,
l_max)

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, single_vertex_q_l_fcc[i])
    println()
end

evolution_dict = NA.get_evolution_dict(nr_vertices = 16, network_type = "primitive cubic")
graph_dict_primitive_cubic = NG.get_periodic_network(evolution_dict)

single_vertex_q_l_primitive_cubic = NA.get_q_l_averaged_single_vertex_dict(graph_dict_primitive_cubic,
central_vertex,
l_max)

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, single_vertex_q_l_primitive_cubic[i])
    println()
end


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.75, 1.0, 2.0, 4.0]

for temperature in temperatures

    for nr_heating_mc_steps in [0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 5.0, 10.0]

        evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
        nr_monte_carlo_steps_per_temperature_vec = [nr_heating_mc_steps, 50], min_ring_size = 3)

        graph_dict = NG.get_periodic_network(evolution_dict)

        graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
                evolution_dict; 
            print_progress = true,
            print_every_nr_attempted_bond_switches = 500)

        evolution_dict["total_energy_vec"] = total_energy_vec
        evolution_dict["move_accepted_vec"] = move_accepted_vec

        NG.plot_network(graph_dict)

        filename = "216_vertices_T_"*string(temperature)*"_heated_for_"*string(nr_heating_mc_steps)*"_steps_quenched"

        save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

        NG.save_graph_to_h5_and_MGformat(graph_dict,
            filename;
            evolution_dict = evolution_dict,
            save_path 
                = save_path)

    end

end


temperatures = [0.75, 1.0, 2.0, 4.0]

for temperature in temperatures

    for temperature_increase_per_monte_carlo_step in [0.025, 0.05, 0.1, 0.2]

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_increase_per_monte_carlo_step = temperature_increase_per_monte_carlo_step, 
        nr_monte_carlo_steps_per_temperature = 0.01,
        quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
        nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

        graph_dict = NG.get_periodic_network(evolution_dict)

        graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
                evolution_dict; 
            print_progress = true,
            print_every_nr_attempted_bond_switches = 500)

        evolution_dict["total_energy_vec"] = total_energy_vec
        evolution_dict["move_accepted_vec"] = move_accepted_vec

        NG.plot_network(graph_dict)

        filename = "216_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_increase_per_monte_carlo_step)*"_per_mc_quenched"

        save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

        NG.save_graph_to_h5_and_MGformat(graph_dict,
            filename;
            evolution_dict = evolution_dict,
            save_path 
                = save_path)

    end

end


temperatures = [0.75, 1.0, 2.0, 4.0]

for temperature in temperatures

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_cooling_gradient(temperature;
    temperature_decrease_per_monte_carlo_step = 0.1,
    quench = true )

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
    nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

    graph_dict = NG.get_periodic_network(evolution_dict)

    graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
            evolution_dict; 
        print_progress = true,
        print_every_nr_attempted_bond_switches = 500)

    evolution_dict["total_energy_vec"] = total_energy_vec
    evolution_dict["move_accepted_vec"] = move_accepted_vec

    NG.plot_network(graph_dict)

    filename = "216_vertices_T_"*string(temperature)*"_cool_0.1_per_mc_quenched"

    save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

    NG.save_graph_to_h5_and_MGformat(graph_dict,
        filename;
        evolution_dict = evolution_dict,
        save_path 
            = save_path)

end



dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [1.0, 2.0, 4.0, 0.75]


for nr_heating_mc_steps in [0.5, 1.0, 5.0, 10.0 ]
    for temperature in temperatures

        evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [temperature, 0],
        nr_monte_carlo_steps_per_temperature_vec = [nr_heating_mc_steps, 50], min_ring_size = 3)

        graph_dict = NG.get_periodic_network(evolution_dict)

        graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
                evolution_dict; 
            print_progress = true,
            print_every_nr_attempted_bond_switches = 500)

        evolution_dict["total_energy_vec"] = total_energy_vec
        evolution_dict["move_accepted_vec"] = move_accepted_vec

        NG.plot_network(graph_dict)

        filename = "216_vertices_T_"*string(temperature)*"_heated_for_"*string(nr_heating_mc_steps)*"_steps_quenched"

        save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

        NG.save_graph_to_h5_and_MGformat(graph_dict,
            filename;
            evolution_dict = evolution_dict,
            save_path 
                = save_path)

    end

end


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

graph_dict_low_t = NG.load_graph_from_h5_and_MGformat(dict_path*"216_vertices_T_0.1_heated_for_0.5_steps_quenched")

graph_dict_high_t = NG.load_graph_from_h5_and_MGformat(dict_path*"216_vertices_T_1.0_heated_for_0.5_steps_quenched")

l_max = 12

q_l_low_t = NA.get_q_l_total_network_mean_dict(graph_dict_low_t, l_max)

q_l_high_t = NA.get_q_l_total_network_mean_dict(graph_dict_high_t, l_max)

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, q_l_low_t[i])
    println()
end

for i in 0:l_max
    Fmt.printfmt("q_{1:d} = {2:.3f}", i, q_l_high_t[i])
    println()
end


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.1, 1.0, 2.0]


for nr_heating_mc_steps in [0.5]
    for temperature in temperatures

        evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 ,temperature_vec = [temperature, 0],
        nr_monte_carlo_steps_per_temperature_vec = [nr_heating_mc_steps, 50], min_ring_size = 3)

        graph_dict = NG.get_periodic_network(evolution_dict)

        graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
                evolution_dict; 
            print_progress = true,
            print_every_nr_attempted_bond_switches = 1000)

        evolution_dict["total_energy_vec"] = total_energy_vec
        evolution_dict["move_accepted_vec"] = move_accepted_vec

        filename = "1000_vertices_T_"*string(temperature)*"_heated_for_"*string(nr_heating_mc_steps)*"_steps_quenched"

        save_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

        NG.save_graph_to_h5_and_MGformat(graph_dict,
            filename;
            evolution_dict = evolution_dict,
            save_path 
                = save_path)

    end

end


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

temperatures = [0.1, 1.0, 2.0]

for i in eachindex(temperatures)

    graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*
    "1000_vertices_T_"
    *string(temperatures[i])*"_heated_for_0.5_steps_quenched")

    structure_factor_dict = NA.get_structure_factor_bartlett_isotrope_by_wavenumber_vec(
        graph_dict;
        sampling_distance_step_length = 0.05,
        maximal_sampling_distance = 4*graph_dict["supercell_edge_length"],
        save_result = true,
        save_path = raw"..\analysis_data\random_networks\1000_vertices_T_"
        *string(temperatures[i])*"_heated_for_0.5_steps_quenched",
        print_progress = true,
        label = "T = "*string(temperatures[i]))

end


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"

filename = "216_vertices_T_0.2_heated_for_1.0_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_MGformat(dict_path*filename)

NG.save_graph_to_h5_and_gml(graph_dict,
"my_graph")

dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\my_graph"

loaded_graph_dict = NG.load_graph_from_h5_and_gml(dict_path)

NG.plot_network(graph_dict)


directory_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\without_ring_size_limitation\\"

NG.convert_all_files_in_directory_MGformat_to_gml(directory_path)


# load some network
dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"
filename = "216_vertices_T_0.2_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)
evolution_dict = GU.load_h5_dict(dict_path*filename*"_evolution.h5")
evolution_dict["reject_during_relaxation_cycle_threshold"] = 5

# perform a bond switch
switched_chain = NG.get_random_chain(graph_dict; 
                                min_ring_size = evolution_dict["min_ring_size"], seed=15)

initial_cluster_dict = NG.get_cluster_in_shells_dict(
                                    graph_dict, 
                                    switched_chain; 
                                    shell_nr = evolution_dict["shell_nr"])

NG.switch_chain!(graph_dict,
    switched_chain )

# fully relax cluster

cluster_dict = NG.get_cluster_in_shells_dict(
                                    graph_dict, 
                                    switched_chain; 
                                    shell_nr = evolution_dict["shell_nr"])

temperature = 50
threshold_cluster_energy = initial_cluster_dict["cluster_energy"] - temperature * log(rand())

graph_dict, new_cluster_dict = NG.relax_cluster_keating!(graph_dict,
    cluster_dict, 
    evolution_dict;
    threshold_cluster_energy = threshold_cluster_energy,
    update_total_energy = false,
    print_progress = true)


    
evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 ,temperature_vec = [0.5],
    nr_monte_carlo_steps_per_temperature_vec = [1], min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)
@time graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 200)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

NG.plot_network(graph_dict)


@ProfileView.profview NG.monte_carlo_move!(graph_dict, 
evolution_dict,
temperature;
print_progress = false)

@ProfileView.profview NG.monte_carlo_move!(graph_dict, 
evolution_dict,
temperature;
print_progress = false)

@time graph_dict, move_accepted = NG.monte_carlo_move!(graph_dict, 
evolution_dict,
temperature;
print_progress = false)

@time graph_dict, move_accepted = NG.monte_carlo_move!(graph_dict, 
evolution_dict,
temperature;
print_progress = false)

@time graph_dict, move_accepted = NG.monte_carlo_move!(graph_dict, 
evolution_dict,
temperature;
print_progress = false)


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_multiple_runs/test_networks_run_1/"
save_path = "../structures/random_networks/216_vertices_multiple_runs/test_networks_run_2/"

NG.generate_graphs_from_evolution_dicts_in_directory(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 100,
    print_progress = true)


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_multiple_runs/216_vertices_run_1/"
save_path = "../structures/random_networks/216_vertices_multiple_runs/216_vertices_run_4/"

println("Starting network generation")

NG.generate_graphs_from_evolution_dicts_in_directory(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 500,
    print_progress = true)



evolution_dicts_directory_path = "../structures/random_networks/216_vertices_multiple_runs/216_vertices_run_1/"

for i in 3:5

    save_path = "../structures/random_networks/216_vertices_multiple_runs/216_vertices_run_"*string(i)*"/"

    println("Starting run "*string(i))

    NG.generate_graphs_from_evolution_dicts_in_directory(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 500,
    print_progress = true)

end
    

dict_path = raw"..\analysis_data\random_networks\\"
filename = "216_vertices_T_0.2_heated_for_0.5_steps_quenched"

structure_factor_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_array.h5")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path= dict_path*filename)


dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\216_vertices_run_4\\"
filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

NG.plot_network(graph_dict)


dict_path = raw"..\structures\random_networks\without_ring_size_limitation\\"
filename = "1000_vertices_T_0.1_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(graph_dict;
save_result = true,
save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=0.1, heated for 0.5 steps, quenched")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=0.1, heated for 0.5 steps, quenched")


filename = "1000_vertices_T_1.0_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(graph_dict;
save_result = true,
save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")


filename = "1000_vertices_perfect_diamond"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(graph_dict;
save_result = true,
save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")

dict_path = raw"..\analysis_data\random_networks\\"
filename = "1000_vertices_T_0.1_heated_for_0.5_steps_quenched"

structure_factor_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_array.h5")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
gaussian_filter_sigma_x = 2*pi/20, 
gaussian_filter_filtered_data_x_step_length = 2*pi/20,
label = "1000 vertices, T=0.1, heated for 0.5 steps, quenched")


filename = "1000_vertices_T_1.0_heated_for_0.5_steps_quenched"


structure_factor_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_array.h5")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
gaussian_filter_sigma_x = 2*pi/20, 
gaussian_filter_filtered_data_x_step_length = 2*pi/20,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")


filename = "1000_vertices_perfect_diamond"


structure_factor_dict = GU.load_h5_dict(dict_path*filename*"_structure_factor_array.h5")

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict, save_result = true, save_path = raw"..\analysis_data\random_networks\\"*filename,
gaussian_filter_sigma_x = 2*pi/20, 
gaussian_filter_filtered_data_x_step_length = 2*pi/20,
label = "1000 vertices, T=1.0, heated for 0.5 steps, quenched")



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_anneal_quench_multiple_runs\evolution_dicts\\"

randomization_temperature_vec = [2., 4., 8.]
randomization_nr_monte_carlo_steps_vec = [2., 4., 8.]
annealing_temperature_vec = [0.36, 0.5, 1.]

for randomization_temperature in randomization_temperature_vec
    for randomization_nr_monte_carlo_steps in randomization_nr_monte_carlo_steps_vec
        for annealing_temperature in annealing_temperature_vec

            temperature_vec = [randomization_temperature]
            nr_monte_carlo_steps_per_temperature_vec = [randomization_nr_monte_carlo_steps]

            for i in 1:5
                append!(temperature_vec, [annealing_temperature, 0])
                append!(nr_monte_carlo_steps_per_temperature_vec, [randomization_nr_monte_carlo_steps, 50])
            end

            evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

            filename = "216_vertices_randomization_T_"*string(randomization_temperature)*"_randomization_nr_MC_steps_"*string(randomization_nr_monte_carlo_steps)*"_annealing_T_"*string(annealing_temperature)*"_quenched_evolution.h5"

            GU.save_dict_to_h5(evolution_dict;
                        save_path=save_path*filename)

        end
    end
end


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_anneal_quench_multiple_runs/evolution_dicts/"

for i in 1:5

    save_path = "../structures/random_networks/216_vertices_anneal_quench_multiple_runs/run_"*string(i)*"/"

    println("Starting run "*string(i))

    NG.generate_graphs_from_evolution_dicts_in_directory(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 500,
    print_progress = true)

end


dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\216_vertices_run_4\\"
filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_dict = BDA.get_binary_data_from_spatial_network(graph_dict;
    bond_radius = 0.35,
    filename = filename,
    save_result=true)

save_path = raw"..\analysis_data\random_networks\\"*filename

# load dict
complete_autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction_complete.h5")

spectral_density_dict = BDA.get_spectral_density_by_wavevector_array_fft(structure_dict;
    save_complete_autocovariance_fct_direction_dict = false,
    save_result = true,
    save_path = save_path,
    complete_autocovariance_fct_direction_dict = complete_autocovariance_fct_direction_dict)


# loop through folders
for i in 1:5

    graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\run_"*string(i)*"\\"

    structure_dict_save_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\run_"*string(i)*"\\"

    # loop through all files in folder
    for filename in readdir(graph_dict_path)

        # check if file is a h5 file
        if endswith(filename, ".gml")

            # load graph dict
            graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename[1:end-4])

            # create and save voxelized data array
            structure_dict = BDA.get_binary_data_from_spatial_network(graph_dict;
            bond_radius = 0.35,
            voxel_edge_length = 0.2,
            save_path = structure_dict_save_path,
            filename = filename[1:end-4],
            save_result=true)
#
        end
    end

end


graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\216_vertices_run_4\\"
filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"

#graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename)
structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\run_4\\"
structure_dict = GU.load_h5_dict(structure_dict_path*filename*"_structure.h5")
save_path = raw"..\analysis_data\random_networks\\"*filename

# get autocovariance fct for the structure
autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction_pbc.h5")

spectral_density_dict = BDA.get_spectral_density_by_wavevector_array_fft_pbc(structure_dict;
    save_autocovariance_fct_direction_dict = false,
    save_result = true,
    save_path = save_path,
    autocovariance_fct_direction_dict = autocovariance_fct_direction_dict
    )


filename = "216_vertices_T_0.2_heat_cool_0.2_per_mc_quenched"

#graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename)
structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\run_4\\"
structure_dict = GU.load_h5_dict(structure_dict_path*filename*"_structure.h5")
save_path = raw"..\analysis_data\random_networks\\"*filename

# get autocovariance fct for the structure
autocovariance_fct_direction_dict = GU.load_h5_dict(save_path*"_autocovariance_fct_direction_pbc.h5")

# reorganize the autocovariance fct array by shifting all entries by half the length of the array
autocovariance_fct_array = cat(autocovariance_fct_direction_dict["autocovariance_fct_array"][Int(ceil(size(autocovariance_fct_direction_dict["autocovariance_fct_array"])[1]/2)):end,:,:], autocovariance_fct_direction_dict["autocovariance_fct_array"][1:Int(ceil(size(autocovariance_fct_direction_dict["autocovariance_fct_array"])[1]/2))-1,:,:], dims = 1)

autocovariance_fct_array = cat(autocovariance_fct_array[:,Int(ceil(size(autocovariance_fct_array)[2]/2)):end,:], autocovariance_fct_array[:,1:Int(ceil(size(autocovariance_fct_array)[2]/2))-1,:], dims = 2)

autocovariance_fct_array = cat(autocovariance_fct_array[:,:,Int(ceil(size(autocovariance_fct_array)[3]/2)):end], autocovariance_fct_array[:,:,1:Int(ceil(size(autocovariance_fct_array)[3]/2))-1], dims = 3)

autocovariance_fct_direction_dict["autocovariance_fct_array"] = autocovariance_fct_array

spectral_density_dict = BDA.get_spectral_density_by_wavevector_array_fft_pbc(structure_dict;
    save_autocovariance_fct_direction_dict = false,
    save_result = true,
    save_path = save_path,
    autocovariance_fct_direction_dict = autocovariance_fct_direction_dict
    )


dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\216_vertices_run_1\\"
filename = "216_vertices_T_2.0_heated_for_0.5_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(dict_path*filename)

structure_dict = NA.get_binary_data_from_spatial_network(graph_dict;
    bond_radius = 0.35,
    voxel_edge_length = 1/2, 
    filename = filename,
    save_result=false)

@ProfileView.profview autocovariance_fct_direction_dict = NA.get_autocovariance_fct_by_sampling_indices_array(structure_dict;
save_result = false,
print_progress = true)

@ProfileView.profview autocovariance_fct_direction_dict = NA.get_autocovariance_fct_by_sampling_indices_array(structure_dict;
save_result = false,
print_progress = true)



for i in 1:5

    structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_multiple_runs/run_"*string(i)*"/"
    save_path = "../analysis_data/random_networks/216_vertices_multiple_runs/run_"*string(i)*"/"

    NA.get_autocovariance_fct_direction_from_filenames_multithreading(structure_dicts_path;
    print_progress = true,
    save_path = save_path)

end


print_lock = Threads.ReentrantLock()
graph_dicts_path = "../structures/random_networks/216_vertices_multiple_runs/"

structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_multiple_runs/"

autocovariance_dicts_path = "../analysis_data/random_networks/216_vertices_multiple_runs/"

NA.get_binary_data_from_spatial_network_multithreading(graph_dicts_path,
structure_dicts_path;
print_progress = true,
print_lock = print_lock,
nr_runs = 5,
bond_radius = 0.35,
voxel_edge_length = 0.1)



print_lock = Threads.ReentrantLock()

structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_multiple_runs/"

save_path = "../analysis_data/random_networks/216_vertices_multiple_runs/"

NA.get_autocovariance_fct_direction_from_filenames_multithreading(structure_dicts_path;
print_progress = true,
save_path = save_path,
print_lock = print_lock)


evolution_dicts_directory_path = raw"..\structures\random_networks\216_vertices_multiple_runs\\"

nr_monte_carlo_steps_for_quenching_vec = NA.get_monte_carlo_steps_for_quenching_vec(evolution_dicts_directory_path)


print_lock = Threads.ReentrantLock()
evolution_dicts_directory_path = "../structures/random_networks/216_vertices_anneal_quench_multiple_runs/evolution_dicts/"

i = 1

save_path = "../structures/random_networks/216_vertices_anneal_quench_multiple_runs/run_"*string(i)*"/"
println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 500,
print_progress = true,
print_lock = print_lock)



dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

for i in 1:5

    current_dict_path = dict_path*"run_"*string(i)*"\\"

    filenames = readdir(current_dict_path)
    
    for filename in filenames
        loaded_dict = GU.load_h5_dict(current_dict_path*filename)

        new_filename = filename[1:Int((length(filename)-13)/2)]*"_structure.h5"

        GU.save_dict_to_h5(loaded_dict, current_dict_path*new_filename)

    end

end


function get_spectral_densities(structure_dict_path, analysis_data_path)

    for i in 1:5

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames
            structure_dict = GU.load_h5_dict(current_structure_dict_path*structure_dict_filename)
    
            autocovariance_fct_direction_dict = GU.load_h5_dict(current_analysis_data_path*structure_dict_filename[1:end-13]*"_autocovariance_fct_direction.h5")
    
            spectral_density_dict = NA.get_spectral_density_by_wavevector_array_fft(structure_dict;
                save_autocovariance_fct_direction_dict = false,
                save_result = true,
                save_path = current_analysis_data_path*structure_dict_filename[1:end-13],
                autocovariance_fct_direction_dict = autocovariance_fct_direction_dict
                )
    
            spectral_density_angle_averaged_dict = NA.get_spectral_density_angle_averaged(spectral_density_dict;
                gaussian_filter = true,
                gaussian_filter_sigma_x = 2*pi/25, 
                gaussian_filter_filtered_data_x_step_length = 2*pi/25,
                save_result = true,
            save_path = current_analysis_data_path*structure_dict_filename[1:end-13])
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\\"

get_spectral_densities(structure_dict_path, analysis_data_path)


randomization_temperature = 8.0
randomization_nr_monte_carlo_steps = 2.0
annealing_temperature = 1.0

i=1

temperature_vec = [randomization_temperature]
nr_monte_carlo_steps_per_temperature_vec = [randomization_nr_monte_carlo_steps]

for i in 1:5
    append!(temperature_vec, [annealing_temperature, 0])
    append!(nr_monte_carlo_steps_per_temperature_vec, [randomization_nr_monte_carlo_steps, 50])
end

evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

graph_dict = NG.get_periodic_network(evolution_dict)
graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 500)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec


filename = "216_vertices_randomization_T_"*string(randomization_temperature)*"_randomization_nr_MC_steps_"*string(randomization_nr_monte_carlo_steps)*"_annealing_T_"*string(annealing_temperature)*"_quenched_evolution.h5"

save_path = "../structures/random_networks/216_vertices_anneal_quench_multiple_runs/run_"*string(i)*"/"

NG.save_graph_to_h5_and_MGformat(graph_dict,
            filename;
            evolution_dict = evolution_dict,
            save_path 
                = save_path)



filename_disorder = "216_vertices_T_2.0_heated_for_0.5_steps_quenched"
filename_order = "216_vertices_T_0.4_heated_for_0.25_steps_quenched"

filename_other = "216_vertices_T_0.25_heated_for_0.05_steps_quenched"

data_path_order = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"*filename_order
data_path_disorder = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"*filename_disorder

data_path_other = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"*filename_other

spectral_density_angle_averaged_dict = GU.load_h5_dict(data_path_other*"_spectral_density_angle_averaged.h5")

function two_gaussians(x, p)
    a2, b2, c2, a3, b3, c3 = p
    return a2 * exp.(-b2 .* (x .- c2).^2) .+ a3 * exp.(-b3 .* (x .- c3).^2)
end

function exp_gaussian(x, p)
    a1, b1, a2, b2, c2 = p
    return a1 * exp.(-b1 .* x) .+ a2 * exp.(-b2 .* (x .- c2).^2)
end

function exp_decay(x, p)
    a1, b1 = p
    return a1 * exp.(-b1 .* x)
end


# define function to estimate the fit parameters from data
function estimate_fit_parameters(x_vec, y_vec)

    # get first parameters by assuming that the exponential decay is the dominant feature
    # and linear close to x=0
    a1 = ( x_vec[2]*y_vec[1] - x_vec[1]*y_vec[2] )/(x_vec[2] - x_vec[1])
    b1 = 3*(y_vec[1] - y_vec[2])/( x_vec[2]*y_vec[1] - x_vec[1]*y_vec[2] )

    # locate peaks of spectral density
    pks, vals = Peaks.findmaxima( y_vec )

    # get parameters for two gaussians
    a2 = vals[1]
    b2 = 2.
    c2 = x_vec[pks[1]]
    a3 = vals[2]
    b3 = 2.
    c3 = x_vec[pks[2]]

    # if the exponential decay is not the dominant feature, get parameters for gaussians
    if a1 < 0
        return [a2, b2, c2, a3, b3, c3]
    else
        # if there is no dominant peak, use only exponential decay
        if a2 < 0.5
            return [a1, b1]
        else
            return [a1, b1, a2, b2, c2 ]
        end
    end

end

# get wavenumber vector and structure factor vector
x_vec = spectral_density_angle_averaged_dict["wavenumber_vec"]
y_vec = Measurements.value.(spectral_density_angle_averaged_dict["spectral_density_vec"])

# estimate fit parameters
p0 = estimate_fit_parameters(x_vec, y_vec)

# get function to fit
if length(p0) == 6
    fit_function = two_gaussians
elseif length(p0) == 5
    fit_function = exp_decay_gaussian
else
    fit_function = exp_decay
end

# fit function to data
fit_result = LsqFit.curve_fit(fit_function, x_vec, y_vec, p0)

# get fit parameters and their uncertainties
fit_param = Measurements.measurement.(fit_result.param, 
    sqrt.(LinearAlgebra.diag(LsqFit.estimate_covar(fit_result))) )

Plots.plot(x_vec, y_vec, label="data")
Plots.plot!(x_vec, fit_function(x_vec, p0), label="initial guess", ls=:dot)
Plots.plot!(x_vec, fit_function(x_vec, Measurements.value.(fit_param)), label="fit", ls=:dash)
Plots.ylims!(0,100)


graph_dict_path = raw"..\structure_analysis\structures\random_networks\216_vertices_anneal_quench_multiple_runs\run_1\\"


all_filenames = readdir(graph_dict_path)

filenames_filtered = filter(filename -> endswith(filename, ".gml"), all_filenames)

filenames = [filename[1:end-4] for filename in filenames_filtered]

for filename in filenames
    println(filename)

    # Use a regular expression to match all numbers (both integers and floating point numbers)
    regex = r"\d+\.?\d*"

    # Find all matches in the string
    matches = eachmatch(regex, filename)

    # Convert the matches to numbers (Float64)
    numbers = [parse(Float64, match.match) for match in matches]

    graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename)

    new_filename = "216_vertices_rand_T_"*string(numbers[2])*"_rand_nr_MC_steps_"*string(numbers[3])*"_anneal_T_"*string(numbers[4])*"_quenched"

    NG.save_graph_to_h5_and_gml(graph_dict,
    new_filename;
    evolution_dict = GU.load_h5_dict(graph_dict_path*filename*"_evolution.h5"),
    save_path  = graph_dict_path)
end


graph_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_anneal_quench_multiple_runs\run_1\\"

all_filenames = readdir(graph_dict_path)
filenames = filter(filename -> endswith(filename, "_evolution.h5"), all_filenames)
final_energy_vec = Float64[]

for filename in filenames
    println(filename)

    evolution_dict = GU.load_h5_dict(graph_dict_path*filename)

    push!(final_energy_vec, evolution_dict["total_energy_vec"][end])
end

filenames_sorted = filenames[sortperm(final_energy_vec)   ]
sort!(final_energy_vec)

graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filenames_sorted[26][1:end-13])
NG.plot_spatial_network(graph_dict)



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_anneal_quench_multiple_runs\evolution_dicts_2\\"

randomization_temperature_vec = [2., 4., 8.]
randomization_nr_monte_carlo_steps_vec = [2., 4., 8.]
annealing_temperature_vec = [0.36, 0.5, 1.]

for randomization_temperature in randomization_temperature_vec
    for randomization_nr_monte_carlo_steps in randomization_nr_monte_carlo_steps_vec
        for annealing_temperature in annealing_temperature_vec

            temperature_vec = [randomization_temperature]
            nr_monte_carlo_steps_per_temperature_vec = [randomization_nr_monte_carlo_steps]

            for i in 1:5
                append!(temperature_vec, [annealing_temperature, 0])
                append!(nr_monte_carlo_steps_per_temperature_vec, [6, 50])
            end

            evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

            filename = "216_vertices_rand_T_"*string(randomization_temperature)*"_rand_nr_MC_steps_"*string(randomization_nr_monte_carlo_steps)*"_anneal_T_"*string(annealing_temperature)*"_quenched_evolution.h5"

            GU.save_dict_to_h5(evolution_dict, save_path*filename)

        end
    end
end



print_lock = Threads.ReentrantLock()
evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts/"

i = 2

save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
print_lock = print_lock)


graph_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\run_2\\"

filename = "216_vertices_rand_T_8.0_rand_nr_MC_steps_4.0_anneal_T_0.36_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename)
NG.plot_spatial_network(graph_dict)

evolution_dict = GU.load_h5_dict(graph_dict_path*filename*"_evolution.h5")

Plots.plot(collect(1:length(evolution_dict["total_energy_vec"]))./(216*18), evolution_dict["total_energy_vec"]./216, xlabel="MC steps", ylabel="Energy", label="energy per vertex")


temperature_vec = [2,0,2,0]
nr_monte_carlo_steps_per_temperature_vec = [0.01, 50, 0.01, 50]

evolution_dict = NA.get_evolution_dict(;nr_vertices = 64 ,temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

graph_dict = NG.get_periodic_network(evolution_dict)
graph_dict, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(graph_dict,
        evolution_dict; 
    print_progress = true,
    print_every_nr_attempted_bond_switches = 50)

evolution_dict["total_energy_vec"] = total_energy_vec
evolution_dict["move_accepted_vec"] = move_accepted_vec

Plots.plot(collect(1:length(evolution_dict["total_energy_vec"]))./(64*18), evolution_dict["total_energy_vec"]./64, xlabel="MC steps", ylabel="Energy", label="energy per vertex")



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_3\\"

randomization_temperature = 1.45

randomization_nr_monte_carlo_steps_vec = [4., 6., 8., 10.]
cooling_nr_monte_carlo_steps_vec = [0.5, 1., 2., 4.]

temperature_vec = vcat(collect(1.45:-0.1:0.35), [0])

for randomization_nr_monte_carlo_steps in randomization_nr_monte_carlo_steps_vec
    for cooling_nr_monte_carlo_steps in cooling_nr_monte_carlo_steps_vec
        
        nr_monte_carlo_steps_per_temperature_vec = vcat([randomization_nr_monte_carlo_steps], cooling_nr_monte_carlo_steps .* ones(11), [50] )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
        nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

        filename = "216_vertices_rand_nr_MC_steps_"*string(randomization_nr_monte_carlo_steps)*"_cool_nr_MC_steps_"*string(cooling_nr_monte_carlo_steps)*"_quenched_evolution.h5"

        GU.save_dict_to_h5(evolution_dict, save_path*filename)
    end
end

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_3\\"

randomization_temperature_vec = [1.45]
randomization_nr_monte_carlo_steps_vec = [4., 6., 8., 10.]
annealing_temperature_vec = [0.36]

for randomization_temperature in randomization_temperature_vec
    for randomization_nr_monte_carlo_steps in randomization_nr_monte_carlo_steps_vec
        for annealing_temperature in annealing_temperature_vec

            temperature_vec = [randomization_temperature]
            nr_monte_carlo_steps_per_temperature_vec = [randomization_nr_monte_carlo_steps]

            for i in 1:5
                append!(temperature_vec, [annealing_temperature, 0])
                append!(nr_monte_carlo_steps_per_temperature_vec, [6, 50])
            end

            evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

            filename = "216_vertices_rand_T_"*string(randomization_temperature)*"_rand_nr_MC_steps_"*string(randomization_nr_monte_carlo_steps)*"_anneal_T_"*string(annealing_temperature)*"_quenched_evolution.h5"

            GU.save_dict_to_h5(evolution_dict, save_path*filename)

        end
    end
end


print_lock = Threads.ReentrantLock()

i = 3

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_3/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"



println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
print_lock = print_lock)


save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_4\\"

evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [1.45],
nr_monte_carlo_steps_per_temperature_vec = [10], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(10.0)*"_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [1.45, 0],
nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(10.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = collect(1.45:-0.05:0)
nr_monte_carlo_steps_per_temperature_vec = vcat([10], 0.5 .* ones(length(temperature_vec)-2), [50])

evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(10.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = collect(1.45:-0.05:0)
nr_monte_carlo_steps_per_temperature_vec = vcat([10], 1.0 .* ones(length(temperature_vec)-2), [50])

evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(10.0)*"_cool_nr_MC_steps_"*string(1.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


print_lock = Threads.ReentrantLock()

i = 4

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_4/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"


println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
print_lock = print_lock)



graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_multiple_runs\run_4\\"

evolution_dicts_directory_path = raw"..\structures\random_networks\216_vertices_multiple_runs\evolution_dicts_check\\"

filenames = ["216_vertices_T_0.1_cool_0.1_per_mc_quenched", "216_vertices_T_0.5_heated_for_1.0_steps_quenched", "216_vertices_T_0.25_heated_for_5.0_steps_quenched"]

for filename in filenames

    # load evolution dict
    evolution_dict = GU.load_h5_dict(graph_path * filename * "_evolution.h5")

    # add missing keys to the evolution dict
    evolution_dict["relax_globally_after_threshold_cycle"] = true
    evolution_dict["mean_nr_monte_carlo_steps_for_quenching"] = 13.7

    # save the evolution dict to new folder
    GU.save_dict_to_h5(evolution_dict, evolution_dicts_directory_path * filename * "_evolution.h5")

end


print_lock = Threads.ReentrantLock()


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_multiple_runs/evolution_dicts_check/"
save_path = "../structures/random_networks/216_vertices_multiple_runs/run_check/"


NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
print_lock = print_lock)


graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_multiple_runs\run_2\\"

filenames = ["216_vertices_T_0.1_heated_for_0.5_steps_quenched",
"216_vertices_T_0.2_heated_for_0.5_steps_quenched",
"216_vertices_T_0.4_heated_for_0.5_steps_quenched",
"216_vertices_T_0.5_heated_for_0.5_steps_quenched",
"216_vertices_T_1.0_heated_for_0.5_steps_quenched",
"216_vertices_T_2.0_heated_for_0.5_steps_quenched",
"216_vertices_T_0.1_heated_for_1.0_steps_quenched",
"216_vertices_T_0.1_heated_for_10.0_steps_quenched",
"216_vertices_T_0.4_heated_for_0.1_steps_quenched",
"216_vertices_T_0.4_heated_for_0.25_steps_quenched",
"216_vertices_T_0.5_heated_for_10.0_steps_quenched",
]


evolution_dicts_directory_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\1728_vertices\evolution_dicts\\"

for filename in filenames

    # load evolution dict
    evolution_dict = GU.load_h5_dict(graph_path * filename * "_evolution.h5")

    # add missing keys to the evolution dict
    evolution_dict["relax_globally_after_threshold_cycle"] = true
    evolution_dict["mean_nr_monte_carlo_steps_for_quenching"] = 13.7
    evolution_dict["nr_vertices"] = 1728

    # save the evolution dict to new folder
    GU.save_dict_to_h5(evolution_dict, evolution_dicts_directory_path * "1728" *filename[4:end] * "_evolution.h5")

end


evolution_dicts_directory_path = raw"../structures/random_networks/1728_vertices/evolution_dicts/"

print_lock = Threads.ReentrantLock()

save_path = "../structures/random_networks/1728_vertices/run_cubic/"


NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = true,
print_lock = print_lock)


function get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)

    for i in 1:5

        current_graph_dict_path = graph_dict_path*"run_"*string(i)*"\\"

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames

            structure_dict = GU.load_h5_dict(current_structure_dict_path*structure_dict_filename)
    
            autocovariance_fct_direction_dict = GU.load_h5_dict(current_analysis_data_path*structure_dict_filename[1:end-13]*"_autocovariance_fct_direction.h5")
    
            volume_fract_variance_dict = NA.get_volume_fract_variance(autocovariance_fct_direction_dict;
            save_result = true,
            save_path = current_analysis_data_path*structure_dict_filename[1:end-13])
            
            graph_dict = NG.load_graph_from_h5_and_gml(current_graph_dict_path*structure_dict_filename[1:end-13])

            structure_factor_dict = NA.get_structure_factor_by_wavevector_array(graph_dict;
                save_result = true,
                save_path = current_analysis_data_path*structure_dict_filename[1:end-13])

            structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(structure_factor_dict;
                gaussian_filter = true,
                gaussian_filter_sigma_x = 2*pi/25, 
                gaussian_filter_filtered_data_x_step_length = 2*pi/25,
                save_result = true,
                save_path = current_analysis_data_path*structure_dict_filename[1:end-13])

            #local_nr_variance_dict = NA.get_local_nr_variance_by_window_radius_vec(
            #    graph_dict;
            #    structure_factor_dict = structure_factor_angle_averaged_dict,
            #    window_radius_step_length = 0.2,
            #    save_result = true,
            #    save_path = current_analysis_data_path*structure_dict_filename[1:end-13])
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\\"

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\\"

get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)


graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_multiple_runs\run_2\\"

data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"

filename = "216_vertices_T_0.2_heated_for_0.05_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(graph_path*filename)

structure_factor_angle_averaged_dict = GU.load_h5_dict(data_path*filename*"_structure_factor_angle_averaged.h5")

Plots.plot(structure_factor_angle_averaged_dict["wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict["structure_factor_vec"]) , xlims = (0,10))


graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_multiple_runs\run_2\\"

data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"

filename = "216_vertices_T_0.1_heated_for_5.0_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(graph_path*filename)

spectral_density_angle_averaged_dict = GU.load_h5_dict(data_path*filename*"_spectral_density_angle_averaged.h5")

anisotropy_metric_from_spectral_density = NA.get_anisotropy_metric_from_spectral_density(spectral_density_angle_averaged_dict)


function get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)

    for i in 1:5

        current_graph_dict_path = graph_dict_path*"run_"*string(i)*"\\"

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames
            
            graph_dict = NG.load_graph_from_h5_and_gml(current_graph_dict_path*structure_dict_filename[1:end-13])

            correlation_functions_dict = NA.get_correlation_functions(graph_dict;
                distance_histogram_bin_width = 0.02,
                save_result = true,
                save_path = current_analysis_data_path*structure_dict_filename[1:end-13])
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\\"

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\\"

get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)



graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\run_2\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"

filename = "216_vertices_T_0.2_heated_for_0.01_steps_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*filename)

small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(filename,
    graph_path,
    analysis_data_path;
    save_result = false,
    )


function get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)

    for i in 1:5

        current_graph_dict_path = graph_dict_path*"run_"*string(i)*"\\"

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames

            small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(structure_dict_filename[1:end-13],
            current_graph_dict_path,
            current_analysis_data_path,
            save_result = true)
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

graph_dict_path = raw"..\structures\random_networks\216_vertices_multiple_runs\\"

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_multiple_runs\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\\"

get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)



analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_multiple_runs\run_2\\"

all_filenames = readdir(analysis_data_path)

# get filenames of small scale order anisotropy_metric_from_spectral_density
order_metrics_filenames = [filename for filename in all_filenames if occursin("small_scale_order_metrics", filename)]

# initialize vectors for all order metrics
total_keating_energy_vec = Vector{Float64}(undef, length(order_metrics_filenames))
bond_length_std_vec = Vector{Float64}(undef, length(order_metrics_filenames))
bond_angle_std_vec = Vector{Float64}(undef, length(order_metrics_filenames))
dihedral_angle_std_vec = Vector{Float64}(undef, length(order_metrics_filenames))
q_l_vec_vec = Vector{Vector{Measurements.Measurement{Float64}}}(undef, length(order_metrics_filenames))
cluster_metric_vec = Vector{Float64}(undef, length(order_metrics_filenames))
anisotropy_metric_from_structure_factor_vec = Vector{Float64}(undef, length(order_metrics_filenames))
anisotropy_metric_from_spectral_density_vec = Vector{Float64}(undef, length(order_metrics_filenames))


# loop through order metric filenames
for i in eachindex(order_metrics_filenames)

    # load order metrics
    order_metrics_dict = GU.load_h5_dict(analysis_data_path*order_metrics_filenames[i])

    # get all order metrics and save them to the corresponding vectors
    total_keating_energy_vec[i] = order_metrics_dict["total_keating_energy"]
    bond_length_std_vec[i] = order_metrics_dict["bond_length_std"]
    bond_angle_std_vec[i] = order_metrics_dict["bond_angle_std"]
    dihedral_angle_std_vec[i] = order_metrics_dict["dihedral_angle_std"]
    q_l_vec_vec[i] = order_metrics_dict["q_l_vec"]
    cluster_metric_vec[i] = order_metrics_dict["cluster_metric"]
    anisotropy_metric_from_structure_factor_vec[i] = order_metrics_dict["anisotropy_metric_from_structure_factor"]
    anisotropy_metric_from_spectral_density_vec[i] = order_metrics_dict["anisotropy_metric_from_spectral_density"]

end

# sort all vectors with respect to the keating energy
sort!(total_keating_energy_vec)

order_metrics_filenames = order_metrics_filenames[sortperm(total_keating_energy_vec)]

bond_length_std_vec = bond_length_std_vec[sortperm(total_keating_energy_vec)]
bond_angle_std_vec = bond_angle_std_vec[sortperm(total_keating_energy_vec)]
dihedral_angle_std_vec = dihedral_angle_std_vec[sortperm(total_keating_energy_vec)]
q_l_vec_vec = q_l_vec_vec[sortperm(total_keating_energy_vec)]
cluster_metric_vec = cluster_metric_vec[sortperm(total_keating_energy_vec)]
anisotropy_metric_from_structure_factor_vec = anisotropy_metric_from_structure_factor_vec[sortperm(total_keating_energy_vec)]
anisotropy_metric_from_spectral_density_vec = anisotropy_metric_from_spectral_density_vec[sortperm(total_keating_energy_vec)]

Plots.scatter(bond_angle_std_vec[1:end-2], bond_length_std_vec[1:end-2])
Plots.scatter(anisotropy_metric_from_structure_factor_vec, anisotropy_metric_from_spectral_density_vec)
Plots.scatter(total_keating_energy_vec, anisotropy_metric_from_spectral_density_vec)

Plots.scatter(total_keating_energy_vec[1:end-2], cluster_metric_vec[1:end-2])


filename = order_metrics_filenames[80][1:end-29]
graph_dict = NG.load_graph_from_h5_and_gml(graph_path*filename)

evolution_dict = GU.load_h5_dict(graph_path*filename*"_evolution.h5")
evolution_dict["relax_globally_after_threshold_cycle"] = true
evolution_dict["reject_during_relaxation_cycle_threshold"]  = 5
random_chain = NG.get_random_chain(graph_dict)

total_energy = NG.get_total_energy_keating(graph_dict)

original_graph_dict = deepcopy(graph_dict)

graph_dict = NG.relax_network_keating!(graph_dict,
random_chain,
evolution_dict;
threshold_total_energy = Inf,
update_total_energy = true,
print_progress = true)

println(graph_dict["total_energy"])


function relax_all_networks_globally(graph_path)

    for i in 1:5
        # get current path
        current_path = graph_path * "run_" * string(i) * "\\"

        # get all files in directory
        filenames = readdir(current_path)

        # get filenames of all evolution dicts
        filenames_evolution_dicts = filter(filename -> endswith(filename, "_evolution.h5"), filenames)

        for filename in filenames_evolution_dicts

            println(filename)

            # load evolution dict
            evolution_dict = GU.load_h5_dict(current_path * filename)

            # load graph
            graph_dict = NG.load_graph_from_h5_and_gml(current_path * filename[1:end-13])

            evolution_dict["relax_globally_after_threshold_cycle"] = true
            evolution_dict["reject_during_relaxation_cycle_threshold"]  = 5
            evolution_dict["mean_nr_monte_carlo_steps_for_quenching"]  = 13.7

            # relax network
            graph_dict = NG.relax_network_keating!(graph_dict,
            NG.get_random_chain(graph_dict),
            evolution_dict;
            threshold_total_energy = Inf,
            update_total_energy = true,
            print_progress = false)

            # save graph
            NG.save_graph_to_h5_and_gml(graph_dict, filename[1:end-13]; evolution_dict = evolution_dict,
            save_path = current_path)
        end

    end

    return
end

graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_globally_relaxed\\"

relax_all_networks_globally(graph_path)



graph_dicts_path = "../structures/random_networks/216_vertices_globally_relaxed/"
structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_globally_relaxed/"
analysis_data_path = "../analysis_data/random_networks/216_vertices_globally_relaxed/"

print_lock = Threads.ReentrantLock()

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
structure_dicts_path,
analysis_data_path,
print_progress = true,
runs_vec = [2],
print_lock = print_lock)


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_2\\"

order_metrics_dict = NA.get_small_length_scale_order_metrics_all_files(analysis_data_path;
    l_max_steinhardt_q_l = 12,
    save_result = true,)


for i in [1,3,4,5]
    analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_"*string(i)*"\\"

    order_metrics_dict = NA.get_small_length_scale_order_metrics_all_files(analysis_data_path;
        l_max_steinhardt_q_l = 12,
        save_result = true,)

    println("run $(i) done")

end



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_5\\"

temperature_vec = vcat(0.07 .* ones(20), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(20), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(80.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.07 .* ones(20), collect(0.07:-0.01:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(20), ones(length(collect(0.07:-0.01:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [0.07, 0],
nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(80.0)*"_cool_nr_MC_steps_"*string(1.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.07 .* ones(20), collect(0.07:-0.01:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(20), 0.5 .* ones(length(collect(0.07:-0.01:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [0.07, 0],
nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(80.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


temperature_vec = vcat(0.07 .* ones(10), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(10), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(40.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.07 .* ones(10), collect(0.07:-0.01:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(10), ones(length(collect(0.07:-0.01:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [0.07, 0],
nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(40.0)*"_cool_nr_MC_steps_"*string(1.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.07 .* ones(10), collect(0.07:-0.01:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(10), 0.5 .* ones(length(collect(0.07:-0.01:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 ,temperature_vec = [0.07, 0],
nr_monte_carlo_steps_per_temperature_vec = [10, 50], min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(40.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


print_lock = Threads.ReentrantLock()

i = 5

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_"*string(i)*"/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = true,
print_lock = print_lock)


graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_multiple_runs\run_2\\"

filenames = ["216_vertices_T_0.125_heat_cool_0.1_per_mc_quenched",
"216_vertices_T_0.2_heat_cool_0.1_per_mc_quenched",
"216_vertices_T_0.3_heat_cool_0.1_per_mc_quenched",
"216_vertices_T_0.1_cool_0.1_per_mc_quenched",
"216_vertices_T_0.15_cool_0.1_per_mc_quenched",
"216_vertices_T_0.2_cool_0.1_per_mc_quenched",
]

evolution_dicts_directory_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\1728_vertices\evolution_dicts_2\\"

for filename in filenames

    # load evolution dict
    evolution_dict = GU.load_h5_dict(graph_path * filename * "_evolution.h5")

    # add missing keys to the evolution dict
    evolution_dict["relax_globally_after_threshold_cycle"] = true
    evolution_dict["mean_nr_monte_carlo_steps_for_quenching"] = 13.7
    evolution_dict["reject_during_relaxation_cycle_threshold"] = 5
    evolution_dict["nr_vertices"] = 1728

    # save the evolution dict to new folder
    GU.save_dict_to_h5(evolution_dict, evolution_dicts_directory_path * "1728" *filename[4:end] * "_evolution.h5")

end


save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_6\\"


temperature_vec = vcat(0.065 .* ones(10), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(10), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_T_"*string(0.065)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.062 .* ones(10), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(10), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_T_"*string(0.062)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


# 2 threads
print_lock = Threads.ReentrantLock()

i = 6

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_"*string(i)*"/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = true,
print_lock = print_lock)


# 6 threads
evolution_dicts_directory_path = raw"../structures/random_networks/1728_vertices/evolution_dicts_2/"

print_lock = Threads.ReentrantLock()

save_path = "../structures/random_networks/1728_vertices/run_2/"


NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = true,
print_lock = print_lock)


analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_globally_relaxed\run_2\\"

all_filenames = readdir(analysis_data_path)

# get filenames of small scale order anisotropy_metric_from_spectral_density
order_metrics_filenames = [filename for filename in all_filenames if occursin("small_scale_order_metrics", filename)]

# initialize vectors for all order metrics
anisotropy_metric_from_structure_factor_vec = Vector{Float64}(undef, length(order_metrics_filenames))
anisotropy_metric_from_spectral_density_vec = Vector{Float64}(undef, length(order_metrics_filenames))


# loop through order metric filenames
for i in eachindex(order_metrics_filenames)

    # load order metrics
    order_metrics_dict = GU.load_h5_dict(analysis_data_path*order_metrics_filenames[i])

    # get all order metrics and save them to the corresponding vectors
    total_keating_energy_vec[i] = order_metrics_dict["total_keating_energy"]
    anisotropy_metric_from_structure_factor_vec[i] = order_metrics_dict["anisotropy_metric_from_structure_factor"]
    anisotropy_metric_from_spectral_density_vec[i] = order_metrics_dict["anisotropy_metric_from_spectral_density"]

end

# sort all vectors with respect to the keating energy
sort!(total_keating_energy_vec)

anisotropy_metric_from_structure_factor_vec = anisotropy_metric_from_structure_factor_vec[sortperm(total_keating_energy_vec)]
anisotropy_metric_from_spectral_density_vec = anisotropy_metric_from_spectral_density_vec[sortperm(total_keating_energy_vec)]

Plots.scatter(bond_angle_std_vec[1:end-2], bond_length_std_vec[1:end-2])
Plots.scatter(anisotropy_metric_from_structure_factor_vec, anisotropy_metric_from_spectral_density_vec)
Plots.scatter(total_keating_energy_vec, anisotropy_metric_from_spectral_density_vec)

Plots.scatter(total_keating_energy_vec[1:end-2], cluster_metric_vec[1:end-2])



filename_diamond = "216_vertices_T_0.1_heated_for_0.01_steps_quenched"

graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_globally_relaxed\run_2\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_2\\"

graph_dict = NG.load_graph_from_h5_and_gml(graph_path*filename_diamond) 
NG.plot_spatial_network(graph_dict)

structure_factor_angle_averaged_dict_diamond = GU.load_h5_dict(analysis_data_path*filename_diamond*"_structure_factor_angle_averaged.h5")

spectral_density_angle_averaged_dict_diamond = GU.load_h5_dict(analysis_data_path*filename_diamond*"_spectral_density_angle_averaged.h5")

wavenumbers_to_check_vec = 2*pi*collect(0.2:0.1:1.0)

# get the wavenumbers that lie the clostest to the wavenumbers_to_check_vec
index_vec = [argmin(abs.(spectral_density_angle_averaged_dict_diamond["wavenumber_vec"] .- wavenumbers_to_check_vec[i])) for i in eachindex(wavenumbers_to_check_vec)]

summed_spectral_density_to_check = sum( spectral_density_angle_averaged_dict_diamond["spectral_density_vec"][index_vec] )

std_value_ratio = (Measurements.uncertainty(summed_spectral_density_to_check) / Measurements.value(summed_spectral_density_to_check))



function get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)

    for i in 1:5

        current_graph_dict_path = graph_dict_path*"run_"*string(i)*"\\"

        current_structure_dict_path = structure_dict_path*"run_"*string(i)*"\\"
    
        current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"
    
        structure_dict_filenames = readdir(current_structure_dict_path)
        
        for structure_dict_filename in structure_dict_filenames
            small_scale_order_metrics_dict = GU.load_h5_dict(current_analysis_data_path*structure_dict_filename[1:end-13]*"_small_scale_order_metrics.h5")

            structure_factor_angle_averaged_dict = GU.load_h5_dict(current_analysis_data_path*structure_dict_filename[1:end-13]*"_structure_factor_angle_averaged.h5")

            anisotropy_metric_from_structure_factor = NA.get_anisotropy_metric_from_structure_factor(
    structure_factor_angle_averaged_dict)

            spectral_density_angle_averaged_dict = GU.load_h5_dict(current_analysis_data_path*structure_dict_filename[1:end-13]*"_spectral_density_angle_averaged.h5")

            anisotropy_metric_from_spectral_density = NA.get_anisotropy_metric_from_spectral_density(
                spectral_density_angle_averaged_dict)

            small_scale_order_metrics_dict["anisotropy_metric_from_structure_factor"] = anisotropy_metric_from_structure_factor
            small_scale_order_metrics_dict["anisotropy_metric_from_spectral_density"] = anisotropy_metric_from_spectral_density

            GU.save_dict_to_h5(small_scale_order_metrics_dict, current_analysis_data_path*structure_dict_filename[1:end-13]*"_small_scale_order_metrics.h5")
    
            println(structure_dict_filename*" done")
    
        end
    
    end
end

graph_dict_path = raw"..\structures\random_networks\216_vertices_globally_relaxed\\"

structure_dict_path = raw"..\structures\random_networks\binary_structures\216_vertices_globally_relaxed\\"

analysis_data_path = raw"..\analysis_data\random_networks\216_vertices_globally_relaxed\\"

get_different_functions(graph_dict_path, structure_dict_path, analysis_data_path)


function my_function()

    for i in 1:5
        analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_"*string(i)*"\\"

        order_metrics_dict = NA.get_small_length_scale_order_metrics_all_files(analysis_data_path;
            l_max_steinhardt_q_l = 12,
            save_result = true,)

        println("run $(i) done")

    end

    return
end

my_function()


# 6 threads
evolution_dicts_directory_path = raw"../structures/random_networks/1728_vertices/evolution_dicts_2/"

print_lock = Threads.ReentrantLock()

save_path = "../structures/random_networks/1728_vertices/run_2/"


NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = true,
further_evolve_previous_networks = true,
print_lock = print_lock)


analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\\"

order_metrics_dict = Dict()

# loop through folders and append all order metrics to the order_metrics_dict
for i in 1:5
    
    current_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_globally_relaxed\run_"*string(i)*"\\"

    current_order_metrics_dict = GU.load_h5_dict(current_analysis_data_path*"all_order_metrics.h5")

    for (key, value) in current_order_metrics_dict
        if haskey(order_metrics_dict, key)
            order_metrics_dict[key] = vcat(order_metrics_dict[key], value)
        else
            order_metrics_dict[key] = value
        end
    end
    
end

order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec", "cluster_metric_vec"]

# sort all vectors in order of the total keating energy
for order_metric_name in order_metrics_names
    order_metrics_dict[order_metric_name] = order_metrics_dict[order_metric_name][sortperm(order_metrics_dict["total_keating_energy_vec"])]
end
order_metrics_dict["filenames_vec"] = order_metrics_dict["filenames_vec"][sortperm(order_metrics_dict["total_keating_energy_vec"])]
sort!(order_metrics_dict["total_keating_energy_vec"])

mask_vec = [contains.(order_metrics_dict["filenames_vec"], "heated_for_0.1_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_0.25_steps"),
    contains.(order_metrics_dict["filenames_vec"], "heated_for_0.5_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_1.0_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_5.0_steps"),
contains.(order_metrics_dict["filenames_vec"], "heated_for_10.0_steps"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.025"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.05"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.1"),
contains.(order_metrics_dict["filenames_vec"], "heat_cool_0.2"),
( contains.(order_metrics_dict["filenames_vec"], "cool_0.1")
    .& .!(contains.(order_metrics_dict["filenames_vec"], "heat")) ),
]

filename_vec = ["heated_for_0.1_steps", "heated_for_0.25_steps", "heated_for_0.5_steps", "heated_for_1.0_steps", "heated_for_5.0_steps", "heated_for_10.0_steps", "heat_cool_0.025", "heat_cool_0.05", "heat_cool_0.1", "heat_cool_0.2", "cool_0.1"]

title_vec = ["heated for 0.1 steps", "heated for 0.25 steps", "heated for 0.5 steps", "heated for 1.0 steps", "heated for 5.0 steps", "heated for 10.0 steps", "heat and cool 0.025/MC step", "heat and cool 0.05/MC step", "heat and cool 0.1/MC step", "heat and cool 0.2/MC step", "cool 0.1/MC step"]

i=9

mask = mask_vec[i]
filtered_filenames_vec = order_metrics_dict["filenames_vec"][mask]
filtered_total_keating_energy_vec = order_metrics_dict["total_keating_energy_vec"][mask]
filtered_bond_length_std_vec = order_metrics_dict["bond_length_std_vec"][mask]
filtered_bond_angle_std_vec = order_metrics_dict["bond_angle_std_vec"][mask]
filtered_dihedral_angle_std_vec = order_metrics_dict["dihedral_angle_std_vec"][mask]
filtered_anisotropy_metric_from_structure_factor_vec = order_metrics_dict["anisotropy_metric_from_structure_factor_vec"][mask]
filtered_anisotropy_metric_from_spectral_density_vec = order_metrics_dict["anisotropy_metric_from_spectral_density_vec"][mask]
filtered_cluster_metric_vec = order_metrics_dict["cluster_metric_vec"][mask]

pattern = r"T_([0-9\.]+)"
extracted_numbers = [match(pattern, s).captures[1] for s in filtered_filenames_vec]
temperatures = parse.(Float64, extracted_numbers)
min_temp = minimum(temperatures)
max_temp = maximum(temperatures)
normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

graph_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_globally_relaxed\run_2\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_globally_relaxed\run_2_thick_bonds\\"

run_2_index_vec = []

for temperature in [0.1, 0.125, 0.15,  0.2, 0.25, 0.3, 0.4, 0.5]

    current_filename = "216_vertices_T_"*string(temperature)*"_heat_cool_0.1_per_mc_quenched"

    graph_dict = NG.load_graph_from_h5_and_gml(graph_dict_path*current_filename)

    run_2_index = argmin(abs.(graph_dict["total_energy"] .- filtered_total_keating_energy_vec))
    push!(run_2_index_vec, run_2_index)

    NG.save_mesh_from_spatial_network(graph_dict, current_filename;
    bond_radius = 0.3131,
    save_path=save_path)
end

println("run_2")
println(run_2_index_vec)
filtered_total_keating_energy_vec[run_2_index_vec] ./216



graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\run_6\\"

filename = "216_vertices_rand_T_0.062_quenched_3"

graph_dict = NG.load_graph_from_h5_and_gml(graph_path * filename)
evolution_dict = GU.load_h5_dict(graph_path * filename * "_evolution.h5")
NG.plot_spatial_network(graph_dict)


save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\anneal_quench\evolution_dicts_7\\"


temperature_vec = vcat(0.062 .* ones(2), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(2), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(8.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.062 .* ones(3), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(3), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(12.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.062 .* ones(4), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(4), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(16.0)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat(0.062 .* ones(2), collect(0.06:-0.02:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(2), 0.5 .* ones(length(collect(0.06:-0.02:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(8.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)

temperature_vec = vcat( 0.062 .* ones(3), collect(0.06:-0.02:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(3), 0.5 .* ones(length(collect(0.06:-0.02:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(12.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


temperature_vec = vcat(0.062 .* ones(4), collect(0.06:-0.02:0), [0])
nr_monte_carlo_steps_per_temperature_vec = vcat(4 .* ones(4), 0.5 .* ones(length(collect(0.06:-0.02:0))), [50])
evolution_dict = NA.get_evolution_dict(;nr_vertices = 216 , temperature_vec = temperature_vec,
nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 4)

filename = "216_vertices_rand_nr_MC_steps_"*string(16.0)*"_cool_nr_MC_steps_"*string(0.5)*"_quenched_evolution.h5"

GU.save_dict_to_h5(evolution_dict, save_path*filename)


i = 7

print_lock = Threads.ReentrantLock()

evolution_dicts_directory_path = "../structures/random_networks/anneal_quench/evolution_dicts_"*string(i)*"/"
save_path = "../structures/random_networks/anneal_quench/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = false,
print_lock = print_lock)



dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.285\run_1\\"

save_path_1 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\evolution_dicts\\"

save_path_2 = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.36\evolution_dicts\\"

# load all files in dict path that end with "evolution.h5"
# and save them to save_path

# get all files in dict_path
files = readdir(dict_path)

# filter files that end with "evolution.h5"
files = filter(x -> occursin("evolution.h5", x), files)

# load all files and save them to save_path
for file in files

    # load the file
    evolution_dict = GU.load_h5_dict(dict_path*file)

    evolution_dict["bond_bending_const"] = 0.21
    
    # save the file
    GU.save_dict_to_h5(evolution_dict, save_path_1*file)

    evolution_dict["bond_bending_const"] = 0.36

    # save the file
    GU.save_dict_to_h5(evolution_dict, save_path_2*file)
end



i = 5

print_lock = Threads.ReentrantLock()

evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.21/evolution_dicts/"

save_path = "../structures/random_networks/216_vertices_bond_bending_0.21/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = false,
print_lock = print_lock)



graph_dict_path = "../structures/random_networks/1728_vertices/run_2/"

structure_dict_path = "../structures/random_networks/binary_structures/1728_vertices/"

analysis_data_path = "../analysis_data/random_networks/1728_vertices/run_2/"

filename = "1728_vertices_T_0.2_heat_cool_0.1_per_mc_quenched"


NA.get_all_dicts_from_graph_single_file(filename,
    graph_dict_path,
    structure_dict_path,
    analysis_data_path;
    print_progress = true)



graph_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\1728_vertices\run_2\\"
analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1728_vertices\run_2\\"

filename = "1728_vertices_T_0.2_heat_cool_0.1_per_mc_quenched"

small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(filename,
    graph_path,
    analysis_data_path;
    save_result = true,
    )

filename = "1728_vertices_T_0.125_heat_cool_0.1_per_mc_quenched"

small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(filename,
        graph_path,
        analysis_data_path;
        save_result = true,
        )



save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\1000_vertices_bond_bending_0.285\evolution_dicts\\"

temperatures = [0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24]

temperature_increase_per_monte_carlo_step = 0.1

for temperature in temperatures

    temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
        temperature_increase_per_monte_carlo_step = temperature_increase_per_monte_carlo_step, 
        nr_monte_carlo_steps_per_temperature = 0.01,
        quench = true )

    evolution_dict = NA.get_evolution_dict(;nr_vertices = 1000 ,temperature_vec = temperature_vec,
        nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3)

    filename = "1000_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_increase_per_monte_carlo_step)*"_per_mc_quenched"

    GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

end


i = 1

print_lock = Threads.ReentrantLock()

evolution_dicts_directory_path = "../structures/random_networks/1000_vertices_bond_bending_0.285/evolution_dicts/"

save_path = "../structures/random_networks/1000_vertices_bond_bending_0.285/run_"*string(i)*"/"

println("Starting run "*string(i))

NG.generate_graphs_from_evolution_dicts_in_directory(
evolution_dicts_directory_path,
save_path;
print_every_nr_attempted_bond_switches = 200,
print_progress = true,
save_network_after_each_temperature = false,
print_lock = print_lock)



temperatures = [0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24]

temperature_increase_per_monte_carlo_step = 0.1

nr_vertices = 512

bond_bending_const_vec = [0.21, 0.28, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
            temperature_increase_per_monte_carlo_step = temperature_increase_per_monte_carlo_step, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_increase_per_monte_carlo_step)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end


nr_vertices = 1000

bond_bending_const_vec = [0.21, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
            temperature_increase_per_monte_carlo_step = temperature_increase_per_monte_carlo_step, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_increase_per_monte_carlo_step)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end



nr_vertices = 216

bond_bending_const_vec = [0.21, 0.285, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"_heat_cool\\evolution_dicts\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
            temperature_increase_per_monte_carlo_step = temperature_increase_per_monte_carlo_step, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_increase_per_monte_carlo_step)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.36_heat_cool/evolution_dicts/"
save_path = "../structures/random_networks/216_vertices_bond_bending_0.36_heat_cool/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


evolution_dicts_directory_path = "../structures/random_networks/512_vertices_bond_bending_0.285/evolution_dicts/"
save_path = "../structures/random_networks/512_vertices_bond_bending_0.285/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:2),
    print_lock = Threads.ReentrantLock())



temperatures = [0.11, 0.13, 0.15, 0.17]

temperature_increase_per_monte_carlo_step = 0.1

nr_vertices = 512

bond_bending_const_vec = [0.21, 0.285, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts_2\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
            temperature_increase_per_monte_carlo_step = temperature_increase_per_monte_carlo_step, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_increase_per_monte_carlo_step)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end


nr_vertices = 1000

bond_bending_const_vec = [0.21, 0.285, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\evolution_dicts_2\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
            temperature_increase_per_monte_carlo_step = temperature_increase_per_monte_carlo_step, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_increase_per_monte_carlo_step)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end



nr_vertices = 216

bond_bending_const_vec = [0.21, 0.285, 0.36]

for bond_bending in bond_bending_const_vec

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"* string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"_heat_cool\\evolution_dicts_2\\"

    for temperature in temperatures

        temperature_vec, nr_monte_carlo_steps_per_temperature_vec = NA.get_temperature_sequence_heating_cooling_gradient(temperature;
            temperature_increase_per_monte_carlo_step = temperature_increase_per_monte_carlo_step, 
            nr_monte_carlo_steps_per_temperature = 0.01,
            quench = true )

        evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices ,temperature_vec = temperature_vec,
            nr_monte_carlo_steps_per_temperature_vec = nr_monte_carlo_steps_per_temperature_vec, min_ring_size = 3,
            bond_bending_const = bond_bending)

        filename = string(nr_vertices)*"_vertices_T_"*string(temperature)*"_heat_cool_"*string(temperature_increase_per_monte_carlo_step)*"_per_mc_quenched"

        GU.save_dict_to_h5(evolution_dict, save_path*filename*"_evolution.h5")

    end
end



evolution_dicts_directory_path = "../structures/random_networks/512_vertices_bond_bending_0.285/evolution_dicts_2/"
save_path = "../structures/random_networks/512_vertices_bond_bending_0.285/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())
    

evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.285_heat_cool/evolution_dicts_2/"
save_path = "../structures/random_networks/216_vertices_bond_bending_0.285_heat_cool/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


evolution_dicts_directory_path = "../structures/random_networks/1000_vertices_bond_bending_0.21/evolution_dicts/"
save_path = "../structures/random_networks/1000_vertices_bond_bending_0.21/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.21_heat_cool/evolution_dicts/"
save_path = "../structures/random_networks/216_vertices_bond_bending_0.21_heat_cool/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


graph_dict_path = "../structures/random_networks/1000_vertices_bond_bending_0.285/run_1/"

structure_dict_path = "../structures/random_networks/binary_structures/1000_vertices_bond_bending_0.285/run_1/"

analysis_data_path = "../analysis_data/random_networks/1000_vertices_bond_bending_0.285/run_1/"

# loop through all files in folder
for filename_with_format in readdir(graph_dict_path)

    # check if file is a h5 file
    if endswith(filename_with_format, ".gml")

        filename = filename_with_format[1:end-4]

        NA.get_all_dicts_from_graph_single_file(filename,
            graph_dict_path,
            structure_dict_path,
            analysis_data_path;
            print_progress = true)

    end
end


network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285\run_2\216_vertices_T_0.15_heat_cool_0.1_per_mc_quenched_structure.h5"

structure_dict_network = GU.load_h5_dict(network_path)

# load pachy weevil data

pachy_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\pachy_blue_structure.h5"

structure_dict_pachy = GU.load_h5_dict(pachy_path)

# adjust pachy weevil data to use same keys as random network

structure_dict_pachy["data_binary"] = structure_dict_pachy["data_binary"][:,1:308,1:308]
structure_dict_pachy["size_data"] = (308,308,308)


structure_dict_pachy["volume_fract_tot"] = sum(structure_dict_pachy["data_binary"])/structure_dict_pachy["size_data"][1]^3

# voxel size = 10 nm
# mean bond length =~ 160 nm
# voxel edge length in units of bond length = 0.0625

structure_dict_pachy["voxel_edge_length"] = 0.0625

structure_dict_pachy["mean_edge_length_data"] = structure_dict_pachy["size_data"][1] * structure_dict_pachy["voxel_edge_length"]

# save pachy weevil data

GU.save_dict_to_h5(structure_dict_pachy, raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\pachy_blue_structure_adjusted.h5")


pachy_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\pachy_red_structure.h5"

structure_dict_pachy = GU.load_h5_dict(pachy_path)


structure_dict_pachy["data_binary"] = structure_dict_pachy["data_binary"][:,1:303,1:303]
structure_dict_pachy["size_data"] = (303,303,303)


structure_dict_pachy["volume_fract_tot"] = sum(structure_dict_pachy["data_binary"])/structure_dict_pachy["size_data"][1]^3

# voxel size = 9 nm
# mean bond length =~ 185 nm
# voxel edge length in units of bond length = 0.04865

structure_dict_pachy["voxel_edge_length"] = 0.04865

structure_dict_pachy["mean_edge_length_data"] = structure_dict_pachy["size_data"][1] * structure_dict_pachy["voxel_edge_length"]

# save pachy weevil data

GU.save_dict_to_h5(structure_dict_pachy, raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\pachy_red_structure_adjusted.h5")



structure_dict_path = "../structures/biological/"

analysis_data_path = "../analysis_data/biological/"

filename = "pachy_blue"

NA.get_all_dicts_from_voxelized_structure(filename,
    structure_dict_path,
    analysis_data_path;
    print_progress = true,
    print_lock = Threads.ReentrantLock())

filename = "pachy_red"

NA.get_all_dicts_from_voxelized_structure(filename,
    structure_dict_path,
    analysis_data_path;
    print_progress = true,
    print_lock = Threads.ReentrantLock())


network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285\run_2\216_vertices_T_0.15_heat_cool_0.1_per_mc_quenched_structure.h5"

structure_dict_network = GU.load_h5_dict(network_path)



evolution_dicts_directory_path = "../structures/random_networks/512_vertices_bond_bending_0.21/evolution_dicts/"
save_path = "../structures/random_networks/512_vertices_bond_bending_0.21/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:2),
    print_lock = Threads.ReentrantLock())


network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285\run_2\216_vertices_T_0.4_heat_cool_0.1_per_mc_quenched_structure.h5"

structure_dict_network = GU.load_h5_dict(network_path)

pore_size_distribution_dict = NA.get_pore_size_distribution(structure_dict_network)

pore_size_distribution_second_moment = NA.get_pore_size_distribution_second_moment(pore_size_distribution_dict)

function get_pore_size_distributions()
    file_count = 0

    structure_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285\run_"

    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.285\run_"

    for i in 1:5

        current_structure_path = structure_path * string(i) * "\\"
        current_save_path = save_path * string(i) * "\\"

        # read all files in the current structure path
        structure_files = readdir(current_structure_path)

        for file in structure_files
            # load structure
            structure_dict = GU.load_h5_dict(current_structure_path * file)

            filename = file[1:end-13]

            # get pore size distribution
            pore_size_distribution_dict = NA.get_pore_size_distribution(structure_dict;
            save_result = true,
            save_path = current_save_path * filename)

            # get pore size distribution second moment
            pore_size_distribution_second_moment = NA.get_pore_size_distribution_second_moment(pore_size_distribution_dict)

            # load small scale order metrics dict 
            small_scale_order_metrics_dict = GU.load_h5_dict(current_save_path * filename * "_small_scale_order_metrics.h5")

            # save small scale order metrics dict with pore size distribution second moment
            small_scale_order_metrics_dict["pore_size_distribution_second_moment"] = pore_size_distribution_second_moment

            GU.save_dict_to_h5(small_scale_order_metrics_dict,
            current_save_path * filename*"_small_scale_order_metrics.h5")

            file_count += 1
            println("File ", file_count, " done.")
            
        end
        
    end

    return

end

get_pore_size_distributions()



# load pachy weevil data

pachy_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\pachy_red_autocovariance_fct_direction_complete_unadjusted.h5"

autocovariance_dict_pachy = GU.load_h5_dict(pachy_path)

# cut all arrays to cubic shape

previous_size = size(autocovariance_dict_pachy["autocovariance_fct_array"])[1:3]
central_indices = Int.( (previous_size .+ 1) ./ 2)
half_cubic_edge_length = Int( ( minimum(previous_size[1:3] ) - 1) /2)


autocovariance_dict_pachy["autocovariance_fct_array"] = autocovariance_dict_pachy["autocovariance_fct_array"][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length, central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length, central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length, :]

autocovariance_dict_pachy["autocovariance_fct_array_uncertainty"] = Measurements.uncertainty.(autocovariance_dict_pachy["autocovariance_fct_array"])
autocovariance_dict_pachy["autocovariance_fct_array"] = Measurements.value.(autocovariance_dict_pachy["autocovariance_fct_array"])

autocovariance_dict_pachy["sampling_index_array"] = autocovariance_dict_pachy["sampling_vec_array"][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length, central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length, central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length, :]

autocovariance_dict_pachy["sampling_index_vec_vec"] = [autocovariance_dict_pachy["sampling_distance_vec_vec"][1][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length], autocovariance_dict_pachy["sampling_distance_vec_vec"][2][central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length], autocovariance_dict_pachy["sampling_distance_vec_vec"][3][central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length]]

# voxel size = 9 nm
# mean bond length =~ 185 nm
# voxel edge length in units of bond length = 0.04865
autocovariance_dict_pachy["voxel_edge_length"] = 0.04865
autocovariance_dict_pachy["sampling_distance_vec_vec"] = autocovariance_dict_pachy["sampling_index_vec_vec"] .* autocovariance_dict_pachy["voxel_edge_length"]
autocovariance_dict_pachy["sampling_distance_array"] = autocovariance_dict_pachy["sampling_index_array"] .* autocovariance_dict_pachy["voxel_edge_length"]

# save pachy weevil data

GU.save_dict_to_h5(autocovariance_dict_pachy, raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\pachy_red_autocovariance_fct_direction.h5")



# load pachy weevil data

pachy_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\pachy_blue_autocovariance_fct_direction_complete_unadjusted.h5"

autocovariance_dict_pachy = GU.load_h5_dict(pachy_path)

# cut all arrays to cubic shape

previous_size = size(autocovariance_dict_pachy["autocovariance_fct_array"])[1:3]
central_indices = Int.( (previous_size .+ 1) ./ 2)
half_cubic_edge_length = Int( ( minimum(previous_size[1:3] ) - 1) /2)


autocovariance_dict_pachy["autocovariance_fct_array"] = autocovariance_dict_pachy["autocovariance_fct_array"][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length, central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length, central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length, :]

autocovariance_dict_pachy["autocovariance_fct_array_uncertainty"] = Measurements.uncertainty.(autocovariance_dict_pachy["autocovariance_fct_array"])
autocovariance_dict_pachy["autocovariance_fct_array"] = Measurements.value.(autocovariance_dict_pachy["autocovariance_fct_array"])

autocovariance_dict_pachy["sampling_index_array"] = autocovariance_dict_pachy["sampling_vec_array"][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length, central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length, central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length, :]

autocovariance_dict_pachy["sampling_index_vec_vec"] = [autocovariance_dict_pachy["sampling_distance_vec_vec"][1][central_indices[1]-half_cubic_edge_length:central_indices[1]+half_cubic_edge_length], autocovariance_dict_pachy["sampling_distance_vec_vec"][2][central_indices[2]-half_cubic_edge_length:central_indices[2]+half_cubic_edge_length], autocovariance_dict_pachy["sampling_distance_vec_vec"][3][central_indices[3]-half_cubic_edge_length:central_indices[3]+half_cubic_edge_length]]

# voxel size = 10 nm
# mean bond length =~ 160 nm
# voxel edge length in units of bond length = 0.0625

autocovariance_dict_pachy["voxel_edge_length"] = 0.0625
autocovariance_dict_pachy["sampling_distance_vec_vec"] = autocovariance_dict_pachy["sampling_index_vec_vec"] .* autocovariance_dict_pachy["voxel_edge_length"]
autocovariance_dict_pachy["sampling_distance_array"] = autocovariance_dict_pachy["sampling_index_array"] .* autocovariance_dict_pachy["voxel_edge_length"]

# save pachy weevil data

GU.save_dict_to_h5(autocovariance_dict_pachy, raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\pachy_blue_autocovariance_fct_direction.h5")


function get_spectral_density_from_autocovariance_fct(filename)

    # set path to analysis data
    analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\\"

    # load structure dict
    pachy_structure_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\biological\\"*filename*"_structure.h5"

    # load autocovariance function dict
    pachy_auto_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\biological\\"*filename*"_autocovariance_fct_direction.h5"

    # load structure dict
    structure_dict = GU.load_h5_dict(pachy_structure_path)

    # load autocovariance function dict
    autocovariance_fct_direction_dict = GU.load_h5_dict(pachy_auto_path)

    # get spectral density by wavevector array
    spectral_density_dict = NA.get_spectral_density_by_wavevector_array_fft(structure_dict;
    save_autocovariance_fct_direction_dict = false,
    save_result = true,
    save_path = analysis_data_path*filename,
    autocovariance_fct_direction_dict = autocovariance_fct_direction_dict)

    # get angle averaged spectral density
    spectral_density_angle_averaged_dict = NA.get_spectral_density_angle_averaged(spectral_density_dict;
    gaussian_filter = true,
    gaussian_filter_sigma_x = 2*pi/25, 
    gaussian_filter_filtered_data_x_step_length = 2*pi/25,
    save_result = true,
    save_path = analysis_data_path*filename)

    return

end

get_spectral_density_from_autocovariance_fct("pachy_blue")

get_spectral_density_from_autocovariance_fct("pachy_red")



evolution_dicts_directory_path = "../structures/random_networks/512_vertices_bond_bending_0.21/evolution_dicts/"
save_path = "../structures/random_networks/512_vertices_bond_bending_0.21/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(3:5),
    print_lock = Threads.ReentrantLock())


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"

filename = "512_vertices_perfect_diamond"

# load the network
evolution_dict = NA.get_evolution_dict(;nr_vertices = 512)

graph_dict = NG.get_periodic_network(evolution_dict)

NG.save_graph_to_h5_and_gml(graph_dict, filename; 
            save_path = path)

graph_dict = NG.load_graph_from_h5_and_gml(path * filename)

NG.plot_spatial_network(graph_dict)

NG.save_mesh_from_spatial_network(graph_dict, filename; save_path = path, bond_radius = 0.3131)


graph_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"

structure_dict_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\diamonds\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

filename = "216_vertices_perfect_diamond"

NA.get_all_dicts_from_graph_single_file(filename,
    graph_dict_path,
    structure_dict_path,
    analysis_data_path;
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    print_lock = Threads.ReentrantLock())



graph_dicts_path = "../structures/random_networks/1000_vertices_bond_bending_0.21/"
structure_dicts_path = "../structures/random_networks/binary_structures/1000_vertices_bond_bending_0.21/"
analysis_data_path = "../analysis_data/random_networks/1000_vertices_bond_bending_0.21/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
    structure_dicts_path,
    analysis_data_path::String;
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


graph_dicts_path = "../structures/random_networks/512_vertices_bond_bending_0.21/"
structure_dicts_path = "../structures/random_networks/binary_structures/512_vertices_bond_bending_0.21/"
analysis_data_path = "../analysis_data/random_networks/512_vertices_bond_bending_0.21/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
        structure_dicts_path,
        analysis_data_path::String;
        structure_factor_diamond_std_value_ratio = 1,
        spectral_density_diamond_std_value_ratio = 1,
        pore_size_distribution_nr_sampled_voxels = 20000,
        print_progress = true,
        runs_vec = collect(1:5),
        print_lock = Threads.ReentrantLock())


graph_dicts_path = "../structures/random_networks/216_vertices_bond_bending_0.21_heat_cool/"
structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_bond_bending_0.21_heat_cool/"
analysis_data_path = "../analysis_data/random_networks/216_vertices_bond_bending_0.21_heat_cool/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
                structure_dicts_path,
                analysis_data_path::String;
                structure_factor_diamond_std_value_ratio = 1,
                spectral_density_diamond_std_value_ratio = 1,
                pore_size_distribution_nr_sampled_voxels = 20000,
                print_progress = true,
                runs_vec = collect(1:5),
                print_lock = Threads.ReentrantLock())


graph_dicts_path = "../structures/random_networks/1000_vertices_bond_bending_0.285/"
structure_dicts_path = "../structures/random_networks/binary_structures/1000_vertices_bond_bending_0.285/"
analysis_data_path = "../analysis_data/random_networks/1000_vertices_bond_bending_0.285/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
    structure_dicts_path,
    analysis_data_path::String;
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


graph_dicts_path = "../structures/random_networks/512_vertices_bond_bending_0.285/"
structure_dicts_path = "../structures/random_networks/binary_structures/512_vertices_bond_bending_0.285/"
analysis_data_path = "../analysis_data/random_networks/512_vertices_bond_bending_0.285/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
        structure_dicts_path,
        analysis_data_path::String;
        structure_factor_diamond_std_value_ratio = 1,
        spectral_density_diamond_std_value_ratio = 1,
        pore_size_distribution_nr_sampled_voxels = 20000,
        print_progress = true,
        runs_vec = collect(1:5),
        print_lock = Threads.ReentrantLock())


graph_dicts_path = "../structures/random_networks/216_vertices_bond_bending_0.285_heat_cool/"
structure_dicts_path = "../structures/random_networks/binary_structures/216_vertices_bond_bending_0.285_heat_cool/"
analysis_data_path = "../analysis_data/random_networks/216_vertices_bond_bending_0.285_heat_cool/"

NA.get_all_dicts_from_graphs_multithreading(graph_dicts_path,
                structure_dicts_path,
                analysis_data_path::String;
                structure_factor_diamond_std_value_ratio = 1,
                spectral_density_diamond_std_value_ratio = 1,
                pore_size_distribution_nr_sampled_voxels = 20000,
                print_progress = true,
                runs_vec = collect(1:5),
                print_lock = Threads.ReentrantLock())


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.21\run_1\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\geometrical_models\216_vertices_bond_bending_0.21\run_1\\"

filename = "216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(path * filename)

#NG.plot_spatial_network(graph_dict)

NG.save_mesh_from_spatial_network(graph_dict, filename; save_path = save_path, bond_radius = 0.35,
vector_out_of_supercell_length = 1)


paths = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\216_vertices_bond_bending_0.285_heat_cool\run_1\216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_structure.h5",
    raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\512_vertices_bond_bending_0.285\run_1\512_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_structure.h5",
    raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\binary_structures\1000_vertices_bond_bending_0.285\run_1\1000_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_structure.h5"]


for path in paths

    structure_dict = GU.load_h5_dict(path)

    NA.get_digital_sphere_mask_dict(structure_dict["size_data"];
        save_result = true,
        save_path = raw"..\analysis_data\random_networks\digital_sphere_masks\\")

    println("done")

end



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

files = ["216_vertices_perfect_diamond_small_scale_order_metrics.h5", "512_vertices_perfect_diamond_small_scale_order_metrics.h5", "1000_vertices_perfect_diamond_small_scale_order_metrics.h5"]

anisotropy_metric_from_structure_factor_vec = Vector{Float64}(undef, 3)
anisotropy_metric_from_spectral_density_vec = Vector{Float64}(undef, 3)

for i in eachindex(files)

    current_dict = GU.load_h5_dict(path * files[i])

    anisotropy_metric_from_structure_factor_vec[i] = current_dict["anisotropy_metric_from_structure_factor"]
    anisotropy_metric_from_spectral_density_vec[i] = current_dict["anisotropy_metric_from_spectral_density"]

    println("Anisotropy metric from structure factor for ", files[i], ": ", anisotropy_metric_from_structure_factor_vec[i])
    println("Anisotropy metric from spectral density for ", files[i], ": ", anisotropy_metric_from_spectral_density_vec[i])
end


analysis_data_paths = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.21_heat_cool\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.285_heat_cool\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.285\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.285\\"]

for analysis_data_path in analysis_data_paths
    for i in 1:5
        NA.get_small_length_scale_order_metrics_all_files(analysis_data_path*"run_$i\\";
        l_max_steinhardt_q_l = 12,
        save_result = true,)
    end
    println("done with $analysis_data_path")
end


evolution_dicts_directory_path = "../structures/random_networks/216_vertices_bond_bending_0.36/evolution_dicts/"

save_path = "../structures/random_networks/216_vertices_bond_bending_0.36/"

NG.generate_graphs_from_evolution_dicts_in_directory_multiple_runs(
    evolution_dicts_directory_path,
    save_path;
    print_every_nr_attempted_bond_switches = 200,
    print_progress = true,
    save_network_after_each_temperature = false,
    further_evolve_previous_networks = false,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\geometrical_models\diamonds\\"

nr_vertices_vec = [216, 512, 1000]

bond_radius_vec = [0.26, 0.35]

for nr_vertices_vec in nr_vertices_vec
    filename = string(nr_vertices_vec, "_vertices_perfect_diamond")
    graph_dict = NG.load_graph_from_h5_and_gml(path*filename)

    for bond_radius in bond_radius_vec

        save_filename = string(filename, "_bond_radius_", bond_radius)
        
        NG.save_mesh_from_spatial_network(graph_dict, save_filename; save_path = save_path, bond_radius = bond_radius,
        vector_out_of_supercell_length = 1)
    end
end


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\1000_vertices_bond_bending_0.285\run_2\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\geometrical_models\1000_vertices_bond_bending_0.285\run_2\\"

temperatures = vcat(collect(0.1:0.01:0.18), collect(0.2:0.02:0.24))

for temperature in temperatures
    filename = string("1000_vertices_T_", temperature, "_heat_cool_0.1_per_mc_quenched")

    graph_dict = NG.load_graph_from_h5_and_gml(path*filename)

    save_filename = string(filename, "_br_0.26")
        
    NG.save_mesh_from_spatial_network(graph_dict, save_filename; save_path = save_path, bond_radius = 0.26,
        vector_out_of_supercell_length = 1)
end



path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.285_heat_cool\run_2\\"

save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"


filename = "216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched"

graph_dict = NG.load_graph_from_h5_and_gml(path*filename)

NG.save_spatial_network_to_gml(graph_dict["spatial_network"], filename, save_path=save_path)


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"

load_path = raw"C:\\Users\\HemmannF\\OneDrive - Université de Fribourg\\structure_analysis\\structures\\random_networks\\216_vertices_bond_bending_0.285\\run_1\\216_vertices_T_0.25_heated_for_0.1_steps_quenched.gml"
path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"

spatial_network = NG.load_spatial_network_from_gml(load_path)

NG.save_spatial_network_to_gml(spatial_network,
    "216_vertices_T_0.25_heated_for_0.1_steps_quenched";
    save_path=path)

path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"
function myfunc()
    count = 0
    filepath_list = []
    # Iterate through all directories and subdirectories
    for (root, dirs, files) in walkdir(path)
        for file in files
            if endswith(file, ".gml")

                joined_path = joinpath(root, file)
                spatial_network = NG.load_spatial_network_from_gml(joined_path)
                NG.save_spatial_network_to_gml(spatial_network,
                file[1:end-4];
                    save_path=root*"\\")

                #println(joinpath(root, file))  # Print the full path to the .gml file
                push!(filepath_list, joinpath(root, file))

                #print every 50th file
                count += 1
                if count % 50 == 0
                    println(count, " ", joined_path)
                end
            end
        end
    end

    return
end

myfunc()


path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.285\run_2\\"

filename = "216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched"

# load the network
graph_dict = NG.load_graph_from_h5_and_gml(path*filename)

my_graph_dict = MetaGraphsNext.MetaGraph(
    Graphs.Graph();  
    label_type=Int64, 
    vertex_data_type=Dict{String, Any},  
    edge_data_type=Dict{String, Any}, 
    graph_data=Dict{String, Any}("coordination_nr"   => 4)
)



function my_func()

    nr_vertices_vec = [216, 512, 1000]
    bond_bending_vec = [0.21, 0.285, 0.36]

    # loop through each number of vertices
    for nr_vertices in nr_vertices_vec

        # loop through each bond bending
        for bond_bending in bond_bending_vec

            path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"*string(nr_vertices, "_vertices_bond_bending_", bond_bending, "\\")

            for run in 1:5
                current_path = path*"run_"*string(run)*"\\"

                # get all files in the current path
                files = readdir(current_path)

                # get vector of all gml files
                gml_files = [file for file in files if endswith(file, ".gml")]

                # loop through each file that is a gml file
                for gml_file in gml_files
                    filename  = gml_file[1:end-4]

                    spatial_network = NG.load_spatial_network_from_h5_and_gml(current_path*filename)

                    # save spatial network
                    NG.save_spatial_network_to_gml(
                        spatial_network,
                        filename;
                        save_path = current_path)

                    # delete h5 file
                    rm(current_path*filename*".h5")
                end
            end
        end
    end
end

my_func()



nr_vertices = 512

spatial_network_path = "../structures/random_networks/"*string(nr_vertices)*"_vertices_bond_bending_0.36/"
structure_dict_path = "../structures/random_networks/binary_structures/"*string(nr_vertices)*"_vertices_bond_bending_0.36/"
analysis_data_path = "../analysis_data/random_networks/"*string(nr_vertices)*"_vertices_bond_bending_0.36/"

NA.get_all_dicts_from_networks_multithreading(
    spatial_network_path,
    structure_dict_path,
    analysis_data_path;
    bond_radius = 0.35,
    voxel_edge_length = 0.1,
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())



spatial_network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\216_vertices_bond_bending_0.36\run_1\\"

filename = "216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched"

structure_dict_path = raw"C:\Users\HemmannF\Downloads\\"

analysis_data_path = structure_dict_path

NA.get_all_dicts_from_network_single_file(
    filename,
    spatial_network_path,
    structure_dict_path,
    analysis_data_path;
    pore_size_distribution_nr_sampled_voxels = 1000,
    print_progress = true)



analysis_data_paths = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.36\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.36\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.36\\"]

for analysis_data_path in analysis_data_paths
    for i in 1:5
        NA.get_small_length_scale_order_metrics_all_files(analysis_data_path*"run_$i\\";
        l_max_steinhardt_q_l = 12,
        save_result = true,)
    end
    println("done with $analysis_data_path")
end

analysis_data_paths = [raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.21\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\216_vertices_bond_bending_0.285\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\512_vertices_bond_bending_0.285\\",
raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\1000_vertices_bond_bending_0.285\\"]

for analysis_data_path in analysis_data_paths
    for i in 1:5
        NA.get_small_length_scale_order_metrics_all_files(analysis_data_path*"run_$i\\";
        l_max_steinhardt_q_l = 12,
        save_result = true,)
    end
    println("done with $analysis_data_path")
end


function my_func()

    nr_vertices_vec = [512, 1000]
    bond_bending_vec = [0.21, 0.285, 0.36]

    # loop through each number of vertices
    for nr_vertices in nr_vertices_vec

        # loop through each bond bending
        for bond_bending in bond_bending_vec

            path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"*string(nr_vertices, "_vertices_bond_bending_", bond_bending, "\\")

            analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\\"*string(nr_vertices, "_vertices_bond_bending_", bond_bending, "\\")

            for run in 1:5

                println("Starting run ", run, " with ", nr_vertices, " vertices and bond bending ", bond_bending)

                current_path = path*"run_"*string(run)*"\\"
                current_analysis_data_path = analysis_data_path*"run_"*string(run)*"\\"

                # get all files in the current path
                files = readdir(current_path)

                # get vector of all gml files
                gml_files = [file for file in files if endswith(file, ".gml")]

                # loop through each file that is a gml file
                for gml_file in gml_files
                    filename  = gml_file[1:end-4]

                    spatial_network = NG.load_spatial_network_from_gml(current_path*filename*".gml")

                    # get correlation functions
                    correlation_functions_dict = NA.get_correlation_functions(
                        spatial_network;
                        distance_histogram_bin_width = 0.02,
                        save_result = true,
                        save_path = current_analysis_data_path*filename,
                        label = nothing)

                    # get all order metrics that contain information about small length scales
                    small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(
                        filename,
                        current_path,
                        current_analysis_data_path;
                        l_max_steinhardt_q_l = 12,
                        structure_factor_diamond_std_value_ratio 
                            = 1,
                        spectral_density_diamond_std_value_ratio 
                            = 1,
                        save_result = true,
                        )
                end
            end
        end
    end
end

my_func()


current_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"
current_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

# get all files in the current path
files = readdir(current_path)

# get vector of all gml files
gml_files = [file for file in files if endswith(file, ".gml")]

# loop through each file that is a gml file
for gml_file in gml_files
    filename  = gml_file[1:end-4]
    spatial_network = NG.load_spatial_network_from_h5_and_gml(current_path*filename)

    NG.save_spatial_network_to_gml(
        spatial_network,
        filename;
        save_path = current_path)

    # delete h5 file
    rm(current_path*filename*".h5")
end


function my_func()
    
    current_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"
    current_analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\diamonds\\"

    # get all files in the current path
    files = readdir(current_path)
    # get vector of all gml files
    gml_files = [file for file in files if endswith(file, ".gml")]
    
    # loop through each file that is a gml file
    for gml_file in gml_files
        filename  = gml_file[1:end-4]
        spatial_network = NG.load_spatial_network_from_gml(current_path*filename*".gml")

        # get correlation functions
        correlation_functions_dict = NA.get_correlation_functions(
            spatial_network;
            distance_histogram_bin_width = 0.02,
            save_result = true,
            save_path = current_analysis_data_path*filename,
            label = nothing)
        # get all order metrics that contain information about small length scales
        small_scale_order_metrics_dict = NA.get_small_length_scale_order_metrics(
            filename,
            current_path,
            current_analysis_data_path;
            l_max_steinhardt_q_l = 12,
            structure_factor_diamond_std_value_ratio 
                = 1,
            spectral_density_diamond_std_value_ratio 
                = 1,
            save_result = true,
            )
    end
end

my_func()


nr_vertices = 512

spatial_network_path = "../structures/random_networks/"*string(nr_vertices)*"_vertices_bond_bending_0.21/"
structure_dict_path = "../structures/random_networks/binary_structures/"*string(nr_vertices)*"_vertices_bond_bending_0.21/"
analysis_data_path = "../analysis_data/random_networks/"*string(nr_vertices)*"_vertices_bond_bending_0.21/"

NA.get_all_dicts_from_networks_multithreading(
    spatial_network_path,
    structure_dict_path,
    analysis_data_path;
    bond_radius = 0.35,
    voxel_edge_length = 0.1,
    structure_factor_diamond_std_value_ratio = 1,
    spectral_density_diamond_std_value_ratio = 1,
    pore_size_distribution_nr_sampled_voxels = 20000,
    print_progress = true,
    runs_vec = collect(1:5),
    print_lock = Threads.ReentrantLock())


function my_func()
    
    path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"

    nr_vertices_vec = [216, 512, 1000]
    bond_bending_vec = [0.21, 0.285, 0.36]

    for nr_vertices in nr_vertices_vec
        for bond_bending in bond_bending_vec
            for run in 1:5
                current_path = path*string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\run_"*string(run)*"\\"

                save_path = path*"networks_for_simulation\\"*string(nr_vertices)*"_vertices_bond_bending_"*string(bond_bending)*"\\run_"*string(run)*"\\"

                # get all files in the current path
                files = readdir(current_path)
                # get vector of all gml files
                gml_files = [file for file in files if endswith(file, ".gml")]

                # loop through each file that is a gml file
                for gml_file in gml_files
                    filename  = gml_file[1:end-4]
                    spatial_network = NG.load_spatial_network_from_gml(current_path*filename*".gml")

                    spatial_network_for_simulation = NG.get_spatial_network_for_simulation(
                        spatial_network;
                        vector_out_of_supercell_length = 1,
                        duplicate_bonds_close_to_supercell_edge = true,
                        bond_radius = 0.35,
                        save_result = true,
                        filename = filename,
                        save_path = save_path)
                end
            end

            println("finished run for ", nr_vertices, " vertices and bond bending ", bond_bending)
        end
    end

    current_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\diamonds\\"
    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\networks_for_simulation\diamonds\\"

    # get all files in the current path
    files = readdir(current_path)
    # get vector of all gml files
    gml_files = [file for file in files if endswith(file, ".gml")]

    # loop through each file that is a gml file
    for gml_file in gml_files
        filename  = gml_file[1:end-4]
        spatial_network = NG.load_spatial_network_from_gml(current_path*filename*".gml")

        spatial_network_for_simulation = NG.get_spatial_network_for_simulation(
            spatial_network;
            vector_out_of_supercell_length = 1,
            duplicate_bonds_close_to_supercell_edge = true,
            bond_radius = 0.35,
            save_result = true,
            filename = filename,
            save_path = save_path)
    end

end

my_func()