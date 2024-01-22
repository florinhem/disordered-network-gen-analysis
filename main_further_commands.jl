
"""
These are the calculations for the pachy weevil from the 10.1002advs.202202145 paper
"""

#set raw data path
data_path_raw = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy_10.1002advs.202202145\3Dvolumes\SD_ff04_a51_box5_vox1.tif"


#load data, correct voxel size anisotropy and save data
data_binary = BDA.get_binary_data_from_colorscale(data_path_raw; 
                                        voxel_size=(1,1,1), 
                                        save_data=true, 
                                        data_path_corrected=data_path_corrected)


#set path to voxel size corrected data
data_path_corrected = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy_10.1002advs.202202145\pachy_blue.h5"
      
#compare hyperuniformity criterion for red and blue patches
data_path_corrected_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy_10.1002advs.202202145\pachy_blue.h5"]

#set labels for plotting 
label_vec = ["blue patch"] #, "red patch", "simple diamond"

#create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


#loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    #load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )

    #determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 10,
                                                                                                    window_positioning="scanned",
                                                                                                    window_shape="spherical")

    #create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


#plot window volume times local volume fraction uncertainty as a function of window edge length
#to determine whether structure is hyperuniform
BDA.plot_volume_fraction_variance_times_window_volume(nr_dimensions_data,
                            plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of scanned positions",
                            save_filename="pachy_scanned_spherical_windows_blue_red_sd")





#compare hyperuniformity criterion for red and blue patches
data_path_corrected_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy_10.1002advs.202202145\pachy_blue.h5",
                        raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy_10.1002advs.202202145\pachy_red.h5",
                        raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy_10.1002advs.202202145\simple_diamond.h5"]

#set labels for plotting 
label_vec = ["blue patch", "red patch", "simple diamond"] #, "red patch", "simple diamond"

#in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
#they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

#create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


#loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    #load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )

    #determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 100,
                                                                                                    window_positioning="random",
                                                                                                    window_shape="spherical")

    #create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "voxel_edge_length" => voxel_edge_length_vec[i],
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)

    
    #save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="pachy_volume_fraction_variance_random_spherical_"*label_vec[i])
                                                                                            
end


#plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="pachy_volume_fraction_variance_random_spherical_windows_blue_red_sd")




#compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy_volume_fraction_variance_random_spherical_blue patch.h5",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy_volume_fraction_variance_random_spherical_red patch.h5",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\volume_fraction_variance_random_spherical_D_surface.h5"]

#set labels for plotting 
label_vec = ["blue patch", "red patch", "perfect diamond"] #, "red patch", "simple diamond"

#in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
#they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

#create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []

#loop through data that will be analyzed
for i in eachindex(data_path_vec)

    #load dictionary that contains the following keys:
    #"window_edge_length_vec"
    #"local_volume_fract_variance_vec"
    #"voxel_edge_length"
    #"label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    #adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


#plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="pachy_volume_fraction_variance_random_spherical_windows_blue_red_sd")


#perform convergence analysis to determine the nr of measurements per distance 
#when calculating the autocovariance function
BDA.convergence_analysis_autocovariance_fct_nr_measurements_per_distance(size_data,
                            volume_fract_tot,
                            data_binary;
                            save_plot = true )




#loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    #load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )


    #get autocovariance function as a function of sampling distance
    sampling_distance_vec, autocovariance_fct_vec = BDA.get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                                                                    size_data, 
                                                                                    volume_fract_tot,
                                                                                    data_binary;
                                                                                    nr_measurements_per_distance = 5000)

    #create dictionary for current plot
    plot_dict = Dict("sampling_distance_vec" => sampling_distance_vec,
                    "autocovariance_fct_vec" => autocovariance_fct_vec,
                    "voxel_edge_length" => voxel_edge_length_vec[i],
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)

    
    #save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="pachy_autocovariance_fct_"*label_vec[i])
                                                                                            
end


#plot real part, imaginary part and absolute value of spectral 
#density as a function of the wavenumber
BDA.plot_autocovariance_fct(plot_dict_vec;
                        title="Autocovariance function",
                        save_plot = true,
                        save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\\",
                        save_filename="pachy_autocovariance_fct_blue_red_sd")


#compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy_autocovariance_fct_blue patch.h5",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy_autocovariance_fct_red patch.h5",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\autocovariance_fct_D_surface.h5"]

#set labels for plotting 
label_vec = ["blue patch", "red patch", "perfect diamond"] #, "red patch", "simple diamond"

#in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
#they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

#create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []

#loop through data that will be analyzed
for i in eachindex(data_path_vec)

    #load dictionary that contains the following keys:
    #"sampling_distance_vec"
    #"autocovariance_fct_vec"
    #"voxel_edge_length"
    #"label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    #adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


#plot local volume fraction variance as a function of window edge length
BDA.plot_autocovariance_fct(plot_dict_vec;
                            title="Autocovariance function",
                            save_plot = true,
                            save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\\",
                            save_filename="pachy_autocovariance_fct_blue_red_sd_zoom")



#loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    #load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )


    #get spectral density as a function of the wavenumber
    wavenumber_vec, spectral_density_vec = BDA.get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data,
                                                                                        size_data,
                                                                                        volume_fract_tot,
                                                                                        data_binary;
                                                                                        nr_measurements_per_distance = 5000,
                                                                                        nr_wavenumbers=50) 

    #create dictionary for current plot
    plot_dict = Dict("wavenumber_vec" => wavenumber_vec,
                    "spectral_density_vec" => spectral_density_vec,
                    "voxel_edge_length" => voxel_edge_length_vec[i],
                    "label" => label_vec[i] )

    push!(plot_dict_vec, plot_dict)
    
    #save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="pachy_spectral_density_fct_"*label_vec[i])
                                                                                            
end


#plot real part, imaginary part and absolute value of spectral 
#density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\\",
                        save_filename="pachy_spectral_density_blue_red_sd")


#compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy_spectral_density_fct_blue patch.h5",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy_spectral_density_fct_red patch.h5",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\spectral_density_D_surface.h5"]

#set labels for plotting 
label_vec = ["blue patch", "red patch", "perfect diamond"] #, "red patch", "simple diamond"

#in order to properly scale the x axis, save voxel edge lengths of the anisotropy corrected voxels
#they are: blue 10nm, red 9nm, simple diamond 8.5nm (estimately)
voxel_edge_length_vec = [10, 9, 8.5]

#create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []

#loop through data that will be analyzed
for i in eachindex(data_path_vec)

    #load dictionary that contains the following keys:
    #"wavenumber_vec"
    #"spectral_density_vec"
    #"voxel_edge_length"
    #"label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    #adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


#plot real part, imaginary part and absolute value of spectral 
#density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\\",
                        save_filename="pachy_spectral_density_blue_red_sd")


#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_blue.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "blue", 
                                save_path)


#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_red.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                9, 
                                "red", 
                                save_path)

                    

#set paths where statistical data is stored
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\D_surface"]

#set path where plot will be stored
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\pachy\pachy_blue_red_d"

#plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            voxel_edge_length_vec=[10, 9, 8.5],
            label_vec=["P. c. mirabilis blue", "P. c. mirabilis red", "perfect diamond"]
            )


#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_blue.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "blue", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)


#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_red.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                9, 
                                "red", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)



#set data path
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red"]

#set path to save plot
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\pachy\\"

BDA.plot_statistical_measures(data_path_vec,
            save_path;
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = true
            )



data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_red.h5"


#load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

#set data path
dict_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

#load dict
data_dict = BDA.load_h5_dict(dict_path)

#calculate spectral density
sampled_wavenumbers_vec_vec, sampled_wavevectors_array, spectral_density_array = BDA.get_spectral_density(size_data, 
                                                volume_fract_tot,
                                                data_binary;
                                                nr_measurements_per_direction = 1000,
                                                sampling_vec_array = data_dict["sampling_vec_array"],
                                                autocovariance_fct_array = data_dict["autocovariance_fct_array"])

#path where spectral density is saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red_spectral_density_direction.h5"

#create dict to save
saving_dict = Dict("sampled_wavevectors_array" => sampled_wavevectors_array,
                    "sampled_wavenumbers_vec_vec" => sampled_wavenumbers_vec_vec,
                    "spectral_density_array" => spectral_density_array,
                    "voxel_edge_length" => data_dict["voxel_edge_length"],
                    "label" => data_dict["label"])

BDA.save_dict_to_h5(saving_dict; save_path)



#path where spectral density data is saved
dict_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red_spectral_density_direction.h5"


#path where plot is saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\pachy\pachy_red_"

spectral_density_dict = BDA.load_h5_dict(dict_path)

BDA.plot_spectral_density_heatmap(spectral_density_dict,
    save_path;
    save_plot = false,
    clims = (0,200),
    wavevector_component_to_fix = 3,
    wavevector_value_fixed = 0)



#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_red.h5"

#load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

#set data path
dict_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

#load dict
data_dict = BDA.load_h5_dict(dict_path)

#set matrix of measured direction vectors (in this case the identity matrix)
direction_vec_mat = [1 0 0; 0 1 0; 0 0 1]

#set vector of labels
label_vec = ["pachy red [1,0,0]", "pachy red [0,1,0]", "pachy red [0,0,1]" ]

for i in 1:3

    direction_vec = direction_vec_mat[:,i]

    #determine spectral density
    sampled_wavenumbers_vec, spectral_density_vec = BDA.get_spectral_density_along_direction(size_data, 
                volume_fract_tot,
                data_binary,
                direction_vec;
                sampling_vec_array = data_dict["sampling_vec_array"],
                autocovariance_fct_array = data_dict["autocovariance_fct_array"])
    
    #path where spectral density is saved
    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\\"* label_vec[i] *"_spectral_density_direction.h5"

    #create dict to save
    saving_dict = Dict("wavenumber_vec" => sampled_wavenumbers_vec,
                        "spectral_density_vec" => spectral_density_vec,
                        "voxel_edge_length" => data_dict["voxel_edge_length"],
                        "label" => label_vec[i])

    BDA.save_dict_to_h5(saving_dict; save_path)
end


#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_red.h5"

#load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

#set data path
dict_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

#load dict
data_dict = BDA.load_h5_dict(dict_path)

#set matrix of measured direction vectors (in this case the identity matrix)
direction_vec_mat = [1/sqrt(2) 1/sqrt(3) 1/sqrt(6); 
                    -1/sqrt(2) 1/sqrt(3) 1/sqrt(6); 
                    0 1/sqrt(3) -2/sqrt(6)]

#set vector of labels
label_vec = ["pachy red 1/sqrt(2)*[1,-1,0]", "pachy red 1/sqrt(3)*[1,1,1]", "pachy red 1/sqrt(6)*[1,1,-2]" ]


naming_vec = ["pachy red [1,-1,0]", "pachy red [1,1,1]", "pachy red [1,1,-2]" ]

for i in 1:3

    direction_vec = direction_vec_mat[:,i]

    #determine spectral density
    sampled_wavenumbers_vec, spectral_density_vec = BDA.get_spectral_density_along_direction(size_data, 
                volume_fract_tot,
                data_binary,
                direction_vec;
                sampling_vec_array = data_dict["sampling_vec_array"],
                autocovariance_fct_array = data_dict["autocovariance_fct_array"])
    
    #path where spectral density is saved
    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\\"* naming_vec[i] *"_spectral_density_direction.h5"

    #create dict to save
    saving_dict = Dict("wavenumber_vec" => sampled_wavenumbers_vec,
                        "spectral_density_vec" => spectral_density_vec,
                        "voxel_edge_length" => data_dict["voxel_edge_length"],
                        "label" => label_vec[i])

    BDA.save_dict_to_h5(saving_dict; save_path)
end


naming_vec = ["pachy red [1,-1,0]", "pachy red [1,1,1]", "pachy red [1,1,-2]" ]


#initialize plot dict vec 
plot_dict_vec = []

for i in 1:3
    
    #path where spectral density is saved
    load_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\\"* naming_vec[i] *"_spectral_density_direction.h5"
    
    #load dict
    data_dict = BDA.load_h5_dict(load_path)

    #add current dict to plot dict vector
    push!(plot_dict_vec, data_dict)
end

#path where plot will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\pachy\pachy_red_direction_rotated_axes"

#plot the spectral densities
BDA.plot_spectral_density(plot_dict_vec,
                                save_path,
                                save_plot = true,
                                xlims=[0,0.1])


                    
#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_blue.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "P. c. mirabilis blue", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)


#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_red.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                9, 
                                "P. c. mirabilis red", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)


#set paths where statistical data is stored
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue"]

#set path where plot will be stored
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\pachy\\"

#plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            spectral_density_xlims = [0,0.1],
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = false,
            plot_spectral_density_direction_bool = true
            )


data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_red.h5"

#load data and get all its essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

#set data path
dict_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"

#load dict
data_dict = BDA.load_h5_dict(dict_path)


#path where spectral density is saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red_spectral_density_direction.h5"

#calculate spectral density
wavenumber_vec_vec, wavevector_array, spectral_density_array = BDA.get_spectral_density_by_wavevector_array(size_data, 
                                                volume_fract_tot,
                                                data_binary;
                                                nr_measurements_per_direction = 1000,
                                                sampling_distance_vec_vec = data_dict["sampling_distance_vec_vec"],
                                                sampling_vec_array = data_dict["sampling_vec_array"],
                                                autocovariance_fct_array = data_dict["autocovariance_fct_array"],
                                                save_result = false,
                save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red",
                voxel_edge_length = 9,
                label = "P. c. mirabilis red")



#set paths for structure dict, autocovariance fct dict and where spectral_density along direction is saved 
structure_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_blue_structure.h5"
autocovariance_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue_autocovariance_fct.h5"
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue"
plot_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\pachy\pachy_blue"

#load structure dict
structure_dict = BDA.load_h5_dict(structure_path)

#load autocovariance fct direction dict
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



#set paths for structure dict, autocovariance fct dict and where spectral_density along direction is saved 
structure_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_red_structure.h5"
autocovariance_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red"
plot_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\pachy\pachy_red"

#load structure dict
structure_dict = BDA.load_h5_dict(structure_path)

#load autocovariance fct direction dict
autocovariance_fct_direction_dict = BDA.load_h5_dict(autocovariance_path)

#calculate and save complete autocovariance fct 
complete_autocovariance_dict = BDA.get_complete_autocovariance_fct_by_sampling_vec_array(
    autocovariance_fct_direction_dict;
    save_result = true,
    save_path = save_path)

println("complete dict calculated")

spectral_density_dict = BDA.get_spectral_density_array_by_fft(complete_autocovariance_dict;
        save_result = true,
        save_path = save_path)

println("spectral density dict calculated")

#plot spectral density
BDA.plot_spectral_density_heatmap(spectral_density_dict,
    plot_path;
    title="Spectral density",
    save_plot = true,
    clims = nothing,
    wavevector_component_to_fix = 3,
    wavevector_value_fixed = 0)


data_path_raw = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\3Dvolumes\Blue_SI.tif"


structure_dict = BDA.get_structure_dict_from_colorscale(data_path_raw; 
    voxel_size=(10,12,10), 
    label = "P. c. mirabilis blue",
    save_result=true, 
    save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_blue")


"""
This is where the calculations for nodal surfaces begin
"""


#set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

#set properties of generated data
unit_cell_length = 500
nr_unit_cells = 10


#loop through surfaces
for label in label_vec

    #generate 3d binary data for current nodal surface
    data_binary = BDA.get_binary_data_from_nodal_eqn(unit_cell_length, 
                                                nr_unit_cells,
                                                label)


    #save current nodal surface to h5 file
    BDA.save_nodal_surface_data(data_binary,
                            unit_cell_length, 
                            nr_unit_cells,
                            label)
                                                                                            
end



#compare hyperuniformity criterion for red and blue patches
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\\"

#set surfaces that are analyzed
label_vec = ["D", "G", "P", "I-WP"] 

#create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


#loop through data that will be analyzed
for label in label_vec

    #determine data path of binary data for current nodal structure
    current_path = data_path*label*"_surface.h5"

    #load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(current_path )

    #determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 100,
                                                                                                    window_positioning="random",
                                                                                                    window_shape="spherical")

    #create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "voxel_edge_length" => 10,
                    "label" => label )

    push!(plot_dict_vec, plot_dict)

    
    #save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="volume_fraction_variance_random_spherical_"*label*"_surface")
                                                                                            
end


#plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="volume_fraction_variance_random_spherical_nodal_surfaces")



#loop through data that will be analyzed
for label in label_vec

    #determine data path of binary data for current nodal structure
    current_path = data_path*label*"_surface.h5"

    #load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(current_path )

    #get autocovariance function as a function of sampling distance
    sampling_distance_vec, autocovariance_fct_vec = BDA.get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                                                                    size_data, 
                                                                                    volume_fract_tot,
                                                                                    data_binary;
                                                                                    nr_measurements_per_distance = 20000)

    #create dictionary for current plot
    plot_dict = Dict("sampling_distance_vec" => sampling_distance_vec,
                    "autocovariance_fct_vec" => autocovariance_fct_vec,
                    "voxel_edge_length" => 10,
                    "label" => label )

    push!(plot_dict_vec, plot_dict)

    
    #save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="autocovariance_fct_"*label*"_surface")
                                                                                            
end


#plot local volume fraction variance as a function of window edge length
BDA.plot_autocovariance_fct(plot_dict_vec;
                            title="Autocovariance function",
                            save_plot = true,
                            save_filename="autocovariance_fct_nodal_surfaces")


                            
#create vector for plot_dicts
plot_dict_vec = []

#loop through data that will be analyzed
for label in label_vec

    #determine data path of binary data for current nodal structure
    current_path = data_path*label*"_surface.h5"

    #load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(current_path )

    #load dictionary that contains the following keys:
    #"sampling_distance_vec"
    #"autocovariance_fct_vec"
    #"voxel_edge_length"
    #"label"
    dict = BDA.load_h5_dict("autocovariance_fct_"*label*"_surface")

    #calculate spectral density for loaded autocovariance function
    wavenumber_vec, spectral_density_vec = BDA.get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data,
                size_data, 
                volume_fract_tot,
                data_binary;
                nr_measurements_per_distance = 20000,
                nr_wavenumbers=200,
                sampling_distance_vec = dict["sampling_distance_vec"],
                autocovariance_fct_vec = dict["autocovariance_fct_vec"])

    #create dictionary for current plot
    plot_dict = Dict("wavenumber_vec" => wavenumber_vec,
                    "spectral_density_vec" => spectral_density_vec,
                    "voxel_edge_length" => 10,
                    "label" => label )

    push!(plot_dict_vec, plot_dict)

    
    #save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
                        save_filename="spectral_density_"*label*"_surface")
                                                                                            
end

#plot real part, imaginary part and absolute value of spectral 
#density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_filename="spectral_density_nodal_surfaces")


#set paths where statistical data is stored
data_path_vec = (raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\\"
                    .* ["D", "I-WP", "P", "G"] .* "_surface" )

#set path where plot will be stored
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\nodal_surfaces\nodal_surfaces"

#plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path)



#get data essentials of stervi data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_green.h5"
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )
nr_sampling_distances = BDA.get_nr_sampling_distances(mean_edge_length_data)


#now analyze I-WP with the nr of sampling distances of the stervi weevil
#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\I-WP_surface.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\I-WP_surface_fewer_sampling_distances"
    

#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                10, 
                                "I-WP", 
                                save_path;
                                nr_sampling_distances = nr_sampling_distances,
                                save_local_volume_fraction_variance=false)



label_vec = ["D", "I-WP", "P", "G"]

for label in label_vec

    #path of original data
    data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\\"*label*"_surface.h5"

    #path where analysis data will be saved
    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\\"*label*"_surface"

    #calculate and save all statistical measures
    BDA.save_statistical_measures(data_path, 
                                    10, 
                                    label, 
                                    save_path;
                                    nr_sampling_distances = 500,
                                    save_local_volume_fraction_variance=false)

end


#set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

#set properties of generated data
unit_cell_length = 500
nr_unit_cells = 1


#loop through surfaces
for label in label_vec

    #generate 3d binary data for current nodal surface
    data_binary = BDA.get_binary_data_from_nodal_eqn(unit_cell_length, 
                                                nr_unit_cells,
                                                label)


    #save current nodal surface to h5 file
    BDA.save_nodal_surface_data(data_binary,
                            unit_cell_length, 
                            nr_unit_cells,
                            "single_unit_cell_"*label)
                                                                                            
end


#set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

#set data path
data_path_vec = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\single_unit_cell_" .* label_vec .* "_surface.h5" 


#set path where autocovariance dict will be stored
save_path_vec = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\single_unit_cell_" .* label_vec .* "_surface"


for i in eachindex(data_path_vec)

    #get data and essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_vec[i])

    #get autocovariance function array
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



#set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 

#set path where autocovariance dict of single unit cell is stored
suc_path_vec = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\single_unit_cell_" .* label_vec .* "_surface_autocovariance_fct_direction.h5"

save_path_vec = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\\" .* label_vec .* "_surface"


for i in eachindex(suc_path_vec)

    #load dict
    suc_dict = BDA.load_h5_dict(suc_path_vec[i])

    #get autocovariance function array
    sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.extrapolate_periodic_data_autocovariance_fct_by_sampling_vec_array(
                suc_dict;
                save_result = true,
                save_path = save_path_vec[i])

    println(label_vec[i]*" done")

end



#set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"] 


#path where analysis data will be saved
data_path_vec = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\\" .* label_vec .* "_surface.h5"

#path where analysis data will be saved
save_path_vec = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\\" .* label_vec .* "_surface"

for i in eachindex(data_path_vec)
    
    #calculate and save all statistical measures
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


#set paths where statistical data is stored
plot_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\nodal_surfaces\\"

#plot all statistical measures
BDA.plot_statistical_measures(save_path_vec,
                    plot_path;
                    spectral_density_xlims = [0,0.1],
                    plot_autocovariance_fct_bool = false,
                    plot_spectral_density_bool = false,
                    plot_local_volume_fraction_variance_bool = false,
                    plot_autocovariance_fct_direction_bool = false,
                    plot_spectral_density_direction_bool = true
                    )



#compare hyperuniformity criterion for red and blue patches
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\\"

#set surfaces that are analyzed
label_vec = ["D", "G", "P", "I-WP"] 

structure_dict_path_vec = data_path .* label_vec .* "_surface"

voxel_edge_length_vec = [10, 10, 10, 10]

for i in eachindex(structure_dict_path_vec)

    #load dictionary
    structure_dict = BDA.load_h5_dict(structure_dict_path_vec[i]* ".h5")

    #get essential information of data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(structure_dict["data"])

    #add this info to dictionary and save it
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

#compare hyperuniformity criterion for red and blue patches
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\\single_unit_cell_"

structure_dict_path_vec = data_path .* label_vec .* "_surface"

for i in eachindex(structure_dict_path_vec)

    #load dictionary
    structure_dict = BDA.load_h5_dict(structure_dict_path_vec[i]* ".h5")

    #get essential information of data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(structure_dict["data"])

    #add this info to dictionary and save it
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


#set which surfaces will be generated
label_vec = ["D", "G", "P", "I-WP"]  #

#set properties of generated data
unit_cell_length = 500
nr_unit_cells = 10

for label in label_vec

    #generate 3d binary data for current nodal surface
    structure_dict = BDA.get_binary_data_from_nodal_eqn(unit_cell_length, 
                                                nr_unit_cells,
                                                label;
    save_result=true, 
    save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\\"*label*"_surface_structure.h5")

    println(label*" done")
    
end


"""
These are the calculations for the stervi beetle from the 10.1002/adfm.202302720 paper.
The data was sent by Viola and is not directly taken from Zenodo
"""


#set raw data path
data_path_raw_prefix = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\2d_images_green\slice_"
data_path_raw_suffix = "_max11.tif"

data_path_corrected = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_green.h5"

#load data, correct voxel size anisotropy and save data
data_binary = BDA.get_binary_data_from_colorscale_stack(data_path_raw_prefix,
                                        data_path_raw_suffix,
                                        301; 
                                        save_data=true, 
                                        data_path_corrected=data_path_corrected)


                                        
#set raw data path
data_path_raw_prefix = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\2d_images_blue\slice_"
data_path_raw_suffix = "_max11.tif"

data_path_corrected = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_blue.h5"

#load data, correct voxel size anisotropy and save data
data_binary = BDA.get_binary_data_from_colorscale_stack(data_path_raw_prefix,
                                        data_path_raw_suffix,
                                        250; 
                                        save_data=true, 
                                        data_path_corrected=data_path_corrected)


#loop through data that will be analyzed
for i in eachindex(data_path_corrected_vec)

    #load data and get all its essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_corrected_vec[i] )

    
    #get autocovariance function as a function of sampling distance
    sampling_distance_vec, autocovariance_fct_vec = BDA.get_autocovariance_fct_isotrope_by_sampling_distance_vec(mean_edge_length_data,
                                                                                    size_data, 
                                                                                    volume_fract_tot,
                                                                                    data_binary;
                                                                                    nr_measurements_per_distance = 10000)

    #create dictionary for current plot
    plot_dict = Dict("sampling_distance_vec" => sampling_distance_vec,
                    "autocovariance_fct_vec" => autocovariance_fct_vec,
                    "voxel_edge_length" => 11,
                    "label" => "green" )

    
    #save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="stervi_autocovariance_fct_green")


    #calculate spectral density for loaded autocovariance function
    wavenumber_vec, spectral_density_vec = BDA.get_spectral_density_isotrope_by_wavenumber_vec(mean_edge_length_data,
                size_data, 
                volume_fract_tot,
                data_binary;
                nr_measurements_per_distance = 10000,
                nr_wavenumbers=200,
                sampling_distance_vec = sampling_distance_vec,
                autocovariance_fct_vec = autocovariance_fct_vec)

    #create dictionary for current plot
    plot_dict = Dict("wavenumber_vec" => wavenumber_vec,
                    "spectral_density_vec" => spectral_density_vec,
                    "voxel_edge_length" => 11,
                    "label" => "green" )

    
    #save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
                        save_filename="stervi_spectral_density_green")


    #determine local volume fraction variance vector by using measuring windows
    window_edge_length_vec, local_volume_fract_variance_vec = BDA.get_local_volume_fract_variance_by_window_vec(nr_dimensions_data, 
                                                                                                    mean_edge_length_data,
                                                                                                    size_data, 
                                                                                                    volume_fract_tot,
                                                                                                    data_binary;
                                                                                                    nr_window_sizes = 100,
                                                                                                    window_positioning="random",
                                                                                                    window_shape="spherical")

    #create dictionary for current plot
    plot_dict = Dict("window_edge_length_vec" => window_edge_length_vec,
                    "local_volume_fract_variance_vec" => local_volume_fract_variance_vec,
                    "voxel_edge_length" => 11,
                    "label" => "green" )

    
    #save the plot_dict to a H5 file
    BDA.save_dict_to_h5(copy(plot_dict);
    save_filename="stervi_volume_fraction_variance_random_spherical_green")

                                                                                            
end


#compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi_volume_fraction_variance_random_spherical_green.h5",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\volume_fraction_variance_random_spherical_I-WP_surface.h5"]

#set surfaces that are analyzed
label_vec = ["S. virescens green", "perfect I-WP"]

#set edge lengths
voxel_edge_length_vec = [11, 6]

#create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


#loop through data that will be analyzed
for i in eachindex(data_path_vec)
    
    #load dictionary that contains the following keys:
    #"window_edge_length_vec"
    #"local_volume_fract_variance_vec"
    #"voxel_edge_length"
    #"label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    #adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end

#plot local volume fraction variance as a function of window edge length
BDA.plot_volume_fraction_variance(plot_dict_vec,
                            save_plot=true,
                            title="Spherical measuring windows of random positions",
                            save_filename="stervi_volume_fraction_variance_random_spherical_windows_green_i_wp")



#compare hyperuniformity criterion for red and blue patches
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi_spectral_density_green.h5",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\spectral_density_I-WP_surface.h5"]

#set surfaces that are analyzed
label_vec = ["S. virescens green", "perfect I-WP"]

#set edge lengths
voxel_edge_length_vec = [11, 6]

#create empty vector where plot dictionaries will be stored in            
plot_dict_vec = []


#loop through data that will be analyzed
for i in eachindex(data_path_vec)
    
    #load dictionary that contains the following keys:
    #"wavenumber_vec"
    #"spectral_density_vec"
    #"voxel_edge_length"
    #"label"
    plot_dict = BDA.load_h5_dict(data_path_vec[i])

    #adjust label and voxel edge length
    plot_dict["label"] = label_vec[i]
    plot_dict["voxel_edge_length"] = voxel_edge_length_vec[i]

    push!(plot_dict_vec, plot_dict)
                                                                                            
end


#plot real part, imaginary part and absolute value of spectral 
#density as a function of the wavenumber
BDA.plot_spectral_density(plot_dict_vec;
                        title="Spectral density",
                        save_plot = true,
                        save_filename="stervi_spectral_density_green_i_wp")


#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_green.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_green"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "green", 
                                save_path)

#set paths where statistical data is stored
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_blue",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\I-WP_surface",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_green"]

#set path where plot will be stored
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\stervi\stervi_blue_green_i_wp"

#plot all statistical measures
BDA.plot_statistical_measures(data_path_vec,
            save_path;
            voxel_edge_length_vec=[11,6,11],
            label_vec=["S. virescens blue", "perfect I-WP", "S. virescens green"]
            )



#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_green.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_green"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "green", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)


#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_blue.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_blue"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "blue", 
                                save_path;
                                nr_sampling_distances = 500,
                                save_local_volume_fraction_variance=false)



#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_blue.h5"

save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"

#get essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path )

#get autocovariance function array
sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(size_data,
                       volume_fract_tot,
                       data_binary,
                       nr_measurements_per_direction=10)

#create dict to save
saving_dict = Dict("sampling_vec_array" => sampling_vec_array,
                    "autocovariance_fct_array" => autocovariance_fct_array,
                    "voxel_edge_length" => 11,
                    "label" => "stervi blue")

save_dict_to_h5(saving_dict; save_path)


#path where autocovariance fct data is saved
dict_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"

#load dict
data_dict = BDA.load_h5_dict(dict_path)

#get sampling vec array and autocovariance fct array from dict
sampling_vec_array = data_dict["sampling_vec_array"]
autocovariance_fct_array = data_dict["autocovariance_fct"]

#calculate spectral density
sampled_wavevectors_array, spectral_density_array = get_spectral_density(size_data, 
                                                volume_fract_tot,
                                                data_binary;
                                                sampling_vec_array = sampling_vec_array,
                                                autocovariance_fct_array)

#path where spectral density is saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_spectral_density_direction_blue_small_sampling.h5"

#create dict to save
saving_dict = Dict("sampled_wavevectors_array" => sampled_wavevectors_array,
                    "spectral_density_array" => spectral_density_array,
                    "voxel_edge_length" => 11,
                    "label" => "stervi blue")

save_dict_to_h5(saving_dict; save_path)



#set data path
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_blue.h5"

#get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)


#path of dict
load_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"

#load dict
data_dict = BDA.load_h5_dict(load_path)

autocovariance_fct_array = data_dict["autocovariance_fct_array"]

#calculate spectral density
sampled_wavenumbers_vec_vec, sampled_wavevectors_array, spectral_density_array = BDA.get_spectral_density(size_data, 
volume_fract_tot,
data_binary;
nr_measurements_per_direction = 50,
autocovariance_fct_array = autocovariance_fct_array)


#create dictionary for current plot
plot_dict = Dict("sampled_wavenumbers_vec_vec" => sampled_wavenumbers_vec_vec,
                "sampled_wavevectors_array" => sampled_wavevectors_array,
                "spectral_density_array" => spectral_density_array,
                "voxel_edge_length" => 11,
                "label" => data_dict["label"])


#path of dict
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_spectral_density_direction_blue_small_sampling.h5"

#save the plot_dict to a H5 file
BDA.save_dict_to_h5(copy(plot_dict); save_path=save_path)



#plot heat map of spectral density in x-y-plane for k_z=0
BDA.plot_spectral_density_heatmap(plot_dict,
                                save_path;
                                save_plot = false)

#wait until key is pressed
readline()




autocovariance_fct_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_autocovariance_fct_direction_blue_small_sampling.h5"
autocovariance_fct_dict = BDA.load_h5_dict(autocovariance_fct_path)

spectral_density_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_spectral_density_direction_blue_small_sampling.h5"
spectral_density_dict = BDA.load_h5_dict(spectral_density_path)


#set data path
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_blue.h5"

#get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)

sampling_distance_vec_vec = BDA.get_sampling_distance_vec_vec(size_data)

autocovariance_fct_dict["sampling_distance_vec_1"] = sampling_distance_vec_vec[1]
autocovariance_fct_dict["sampling_distance_vec_2"] = sampling_distance_vec_vec[2]
autocovariance_fct_dict["sampling_distance_vec_3"] = sampling_distance_vec_vec[3]

#save the plot_dict to a H5 file
BDA.save_dict_to_h5(copy(autocovariance_fct_dict); save_path=autocovariance_fct_path)



#set data path
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_blue.h5"

#set path where autocovariance dict will be stored
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_autocovariance_fct_direction_blue.h5"

#get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)

#get autocovariance function array
sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(
                    size_data,
                    volume_fract_tot,
                    data_binary,
                    nr_measurements_per_direction=100)

#create dict to save
saving_dict = Dict("sampling_vec_array" => sampling_vec_array,
                    "sampling_distance_vec_vec" => sampling_distance_vec_vec,
                    "autocovariance_fct_array" => autocovariance_fct_array,
                    "voxel_edge_length" => 11,
                    "label" => "stervi blue")

BDA.save_dict_to_h5(saving_dict; save_path)


#set data path
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_blue.h5"

#set path where autocovariance dict will be stored
dict_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_autocovariance_fct_direction_blue.h5"

#get data and essential information
data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path)

#load autocovariance array
loaded_dict = BDA.load_h5_dict(dict_path)

plot_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\stervi\stervi_direction_blue"


BDA.plot_autocovariance_fct_heatmap(loaded_dict,
    plot_path;
    save_plot = false,
    sampling_vector_component_to_fix = 3,
    sampling_vector_value_fixed = 0 )


    
#set data path
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_blue.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_green.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_blue.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_red.h5"
]

#set path where autocovariance dict will be stored
save_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_blue_autocovariance_fct_direction.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_green_autocovariance_fct_direction.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue_autocovariance_fct_direction.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red_autocovariance_fct_direction.h5"
]

voxel_edge_length_vec = [11,11,10,9]

label_vec = ["stervi blue", "stervi green", "pachy blue","pachy red"]

for i in eachindex(data_path_vec)

    #get data and essential information
    data_binary, volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = BDA.get_data_essentials(data_path_vec[i])

    #get autocovariance function array
    sampling_distance_vec_vec, sampling_vec_array, autocovariance_fct_array = BDA.get_autocovariance_fct_by_sampling_vec_array(
                        size_data,
                        volume_fract_tot,
                        data_binary,
                        nr_measurements_per_direction=1000)

    #create dict to save
    saving_dict = Dict("sampling_vec_array" => sampling_vec_array,
                        "sampling_distance_vec_vec" => sampling_distance_vec_vec,
                        "autocovariance_fct_array" => autocovariance_fct_array,
                        "voxel_edge_length" => voxel_edge_length_vec[i],
                        "label" => label_vec[i])

    BDA.save_dict_to_h5(saving_dict; save_path=save_path_vec[i])

    println(label_vec[i]*" done")

end


#set data path
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_blue",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_green"]

#set path to save plot
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\stervi\\"

BDA.plot_statistical_measures(data_path_vec,
            save_path;
            plot_autocovariance_fct_bool = false,
            plot_spectral_density_bool = false,
            plot_local_volume_fraction_variance_bool = false,
            plot_autocovariance_fct_direction_bool = true
            )

            
#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_blue.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_blue"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "S. virescens blue", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)

#path of original data
data_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_green.h5"

#path where analysis data will be saved
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_green"
    
#calculate and save all statistical measures
BDA.save_statistical_measures(data_path, 
                                11, 
                                "S. virescens green", 
                                save_path;
                                save_autocovariance_fct = false,
                                save_spectral_density = false,
                                save_local_volume_fraction_variance = false,
                                save_autocovariance_fct_direction = false,
                                save_spectral_density_along_directions = true)


#set paths where statistical data is stored
data_path_vec = [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_green",
                raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_blue"]

#set path where plot will be stored
save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\stervi\\"

#plot all statistical measures
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


data_path_vec = [ raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_red_structure.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_blue_structure.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_green_structure.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_blue_structure.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\D_surface_structure.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\G_surface_structure.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\P_surface_structure.h5",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\nodal_surfaces\I-WP_surface_structure.h5"
]

save_path_vec = [ raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_green",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_blue",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\D_surface",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\G_surface",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\P_surface",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\I-WP_surface"
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


data_path_vec_vec =  [ [raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_blue",
    raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\pachy\pachy_red",
    raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\D_surface"],
[
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_blue",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\I-WP_surface",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\stervi\stervi_green"],
[ raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\P_surface",
    raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\I-WP_surface",
    raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\D_surface",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\nodal_surfaces\G_surface"
]
] 

save_path_vec = [ raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\pachy\pachy_blue_red_d",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\stervi\stervi_blue_green_i_wp",
raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\nodal_surfaces\nodal_surfaces"
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


#set raw data path
data_path_raw_prefix = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\2d_images_green\slice_"
data_path_raw_suffix = "_max11.tif"

data_path_corrected = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_green.h5"

structure_dict = BDA.get_structure_dict_from_colorscale_stack(data_path_raw_prefix,
    data_path_raw_suffix,
    301; 
    voxel_size=(11,11,11), 
    label = "S. virescens green",
    save_result=true, 
    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\stervi\stervi_green")



"""
This is where commands for network generation begin
"""


#import my module that contains all functions for the analysis of binary structure data
import .NetworkGeneration as NG
import MetaGraphsNext as MGN


network_dict = NG.get_periodic_network( ; nr_vertices = 150 , 
                            nr_dimensions = 3, 
                            network_type = "diamond")

network_dict["bond_bending_const"] = 0.285

#get list of bonds (edges)
edges_vec = collect(MGN.edge_labels(network_dict["network_graph"]))

#break a random bond
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

#get and print neighbors
vertex_neighbors = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], vertex))

println(vertex_neighbors)

#get position of some vertex
println(graph_dict["spatial_network"][vertex]["position"] )
println(graph_dict["spatial_network"][vertex]["local_energy"] )

#move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex, 
[1,-0.5,0.3] )

println(graph_dict["spatial_network"][vertex]["position"] )
println(graph_dict["spatial_network"][vertex]["local_energy"] )


#relax vertex
graph_dict = NG.relax_single_vertex_keating!(graph_dict, vertex)

println(graph_dict["spatial_network"][vertex]["position"] )
println(graph_dict["spatial_network"][vertex]["local_energy"] )

graph_dict = NG.get_periodic_network( ; nr_vertices = 1500 ,
                            network_type = "diamond")

                            
vertex_vec = [5,12]

#get and print neighbors
neighbors_in_shells_dict, all_vertices_vec = NG.get_neighbors_in_shells_dict(graph_dict, 
                                    vertex_vec; 
                                    shell_nr = 4)


vertex_vec = [5,12]


#get position of some vertex
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

#move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[1], 
[1,-0.5,0.3] )

#move vertex
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


#get position of some vertex
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["local_energy"] )

#move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[1], 
[1,-0.5,0.3] )

#move vertex
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


#get position of some vertex
println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex_vec[1]]["position"] )

#move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_vec[1], 
[1,-0.5,0.3] )

#move vertex
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

#move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex_to_relax, 
[0.2,0.5,-0.9] )

#get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

#get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

#get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

#set energy, gradient and hessian for energy minimization
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

#move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex, 
[1,-0.5,0.3]; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )


#relax vertex
graph_dict = NG.relax_single_vertex_keating!(graph_dict, vertex; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )

println("_______")

vertex = 3

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )

#move vertex
graph_dict = NG.move_vertex!(graph_dict, 
vertex, 
[1,-0.5,0.3]; update_total_energy=true)

println(graph_dict["total_energy"] )
println(graph_dict["spatial_network"][vertex]["position"] )


#relax vertex
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

    #get central vertices neighbors 
    neighbor_vec = collect(MetaGraphsNext.neighbor_labels(
                                graph_dict["spatial_network"], central_vertex))
#create array to store next to nearest neighbors coordinates
#The first array index labels the number of the direct neighbor
next_neighbor_positions_arr = Array{Float64}(undef, 
                                            graph_dict["coordination_nr"],
                                            graph_dict["nr_dimensions"],
                                            graph_dict["coordination_nr"]-1)

#loop through central vertices neighbors
for i in 1:graph_dict["coordination_nr"]
    current_next_neighbor = 1
    #loop through the current neighbor's neighbors
    for next_neighbor in MetaGraphsNext.neighbor_labels(
                                    graph_dict["spatial_network"], neighbor_vec[i])
        if next_neighbor !== central_vertex
            #get next neighbor's virtual coordinates which might be outside of the 
            #supercell if periodic boundary conditions play a role
            next_neighbor_positions_arr[i,:,current_next_neighbor] = NG.get_virtual_position(
                        central_vertex_position,
                        graph_dict["spatial_network"][next_neighbor]["position"],
                        graph_dict["supercell_edge_length"] )
            current_next_neighbor += 1
        end
    end
end



#get cluster after bond switch
cluster_dict = NG.get_cluster_in_shells_dict(
                                graph_dict, 
                                switched_bond; 
                                shell_nr = shell_nr)

#relax cluster around switched bond
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

#figure = NG.plot_network(graph_dict)

switched_bond = (5,12) #NG.get_random_bond(graph_dict)

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


#compare the gradient in the inefficient and efficient calculation
vertex_to_relax = random_bond[2]

#efficient calculation
gradient_eff = NG.gradient_keating_efficient(graph_dict, vertex_to_relax)

#inefficient calculation

#get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

#get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

#get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

#set energy, gradient and hessian for energy minimization

gradient_ineff = zeros(3)
                                        
NG.gradient_keating!(gradient_ineff, initial_position, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

println("Efficient: "*string(gradient_eff))
println("Inefficient "*string(gradient_ineff))


random_bond = NG.get_random_bond(graph_dict)

graph_dict, new_bond_vec = NG.switch_bond!(graph_dict, random_bond )


#compare the gradient in the inefficient and efficient calculation
vertex_to_relax = random_bond[2]

#efficient calculation
gradient = NG.gradient_keating_efficient(graph_dict, vertex_to_relax)

translation_vector_eff =  NG.get_approximate_translation_vector_keating(gradient, 
        graph_dict["bond_bending_const"];
        relaxation_overshoot_factor_r = 1.5,
        relaxation_optimization_parameter_l = 1)

#inefficient calculation

#get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

#get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

#get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

#set energy, gradient and hessian for energy minimization
energy(x) = NG.energy_from_position_keating(x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr )
                                            
gradient!(gradient, x) = NG.gradient_keating!(gradient, x, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

hessian!(hessian, x) = NG.hessian_keating!(hessian, x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)
#find energy minimum
minimizer_result = Optim.optimize(
                            energy, 
                            gradient!, 
                            hessian!,
                            initial_position, 
                            Optim.Newton())

#get relaxed position and local keating energy
relaxed_position = Optim.minimizer(minimizer_result)

#calculate translation vector for relaxed vertex
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


#efficient calculation
gradient = NG.gradient_keating_efficient(graph_dict, vertex_to_relax)

hessian = NG.hessian_keating_efficient(graph_dict, vertex_to_relax)

#calculate translation vector to approximate energy minimum
translation_vector_eff = .- LinearAlgebra.inv(hessian)*gradient

#inefficient calculation

#get initial position of vertex to relax 
initial_position = graph_dict["spatial_network"][vertex_to_relax]["position"]

#get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

#get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

#set energy, gradient and hessian for energy minimization
energy(x) = NG.energy_from_position_keating(x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr )
                                            
gradient!(gradient, x) = NG.gradient_keating!(gradient, x, graph_dict, 
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)

hessian!(hessian, x) = NG.hessian_keating!(hessian, x, graph_dict,
                                            neighbor_positions_mat,
                                            next_neighbor_positions_arr)
#find energy minimum
minimizer_result = Optim.optimize(
                            energy, 
                            gradient!, 
                            hessian!,
                            initial_position, 
                            Optim.Newton())

#get relaxed position and local keating energy
relaxed_position = Optim.minimizer(minimizer_result)

#calculate translation vector for relaxed vertex
translation_vector_ineff = relaxed_position .- initial_position

println("Efficient: "*string(translation_vector_eff))
println("Inefficient "*string(translation_vector_ineff))



#efficient calculation
neighbor_vec = collect(
        MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], central_vertex))

j=1

#get vector pointing from central vertex to neighbor j
distance_vector_j_eff = (sign(neighbor_vec[j] - central_vertex)
* graph_dict["spatial_network"][central_vertex, neighbor_vec[j]]["vector"])

bond_stretching_term_eff = ( - 3/4 * ( 
            graph_dict["spatial_network"][central_vertex, neighbor_vec[j]]["distance_squared"] - 1 
            ) ) .* distance_vector_j_eff

#inefficient calculation

#get initial position of vertex to relax 
x = graph_dict["spatial_network"][vertex_to_relax]["position"]

#get matrix of the vertex's neighbors' positions 
neighbor_positions_mat = NG.get_neighbor_positions_mat(graph_dict, vertex_to_relax)

#get next to nearest neighbors' positions
next_neighbor_positions_arr = NG.get_next_neighbor_positions_arr(graph_dict, vertex_to_relax)

#set energy, gradient and hessian for energy minimization

#get vector pointing from central vertex to neighbor
distance_vector_j = neighbor_positions_mat[:,j] .- x

#get bond stretching term
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


#switch bond
graph_dict, new_bond_vec = NG.switch_bond!(graph_dict, switched_bond )

#get cluster after bond switch
cluster_dict = NG.get_cluster_in_shells_dict(
                graph_dict, 
                switched_bond; 
                shell_nr = shell_nr)

#relax cluster once and update cluster energy
graph_dict, cluster_dict = NG.relax_cluster_keating!(graph_dict,
cluster_dict; 
nr_max_relaxation_cycles = 25,
break_at_relative_cluster_energy_change = 0.001,
reject_during_relaxation_cycle_threshold = 10,
relax_efficiently = true,
update_total_energy = true)

#calculate new total energy and compare to actual total energy
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

#calculate new total energy and compare to actual total energy
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

save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\\"

NG.save_mesh_from_network(graph_dict_to_save, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict_to_save,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)
        


network_names = ["64_vertices_T_0.1", "1000_vertices_T_1_quenched", "1000_vertices_T_2_quenched", "1000_vertices_T_4_quenched"]

save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\\"

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

save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\\"

NG.save_mesh_from_network(graph_dict, filename; save_path = save_path)

NG.save_graph_to_h5_and_MGformat(graph_dict,
    filename;
    evolution_dict = evolution_dict,
    save_path 
        = save_path)