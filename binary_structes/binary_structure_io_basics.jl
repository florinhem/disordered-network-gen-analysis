"""
These functions are required to load data, convert and correct it and to extract
some basic measures.

This code requires the following Packages:
import FileIO   #for file loading and saving
import Images #for image conversion, in our case from color to grayscale
import Statistics   #for statistical operations like mean()
"""


"""
function to add a singleton dimension to an array, taken from
https://stackoverflow.com/questions/42312319/how-do-i-add-a-dimension-to-an-array-opposite-of-squeeze
"""
function extend_dims(array, dim)

    #get current shape
    shape = [size(array)...]

    #set new shape by adding singleton dimension
    insert!(shape,dim,1)

    #reshape array
    return reshape(array, shape...)
end


"""
correct voxel size in data array in case of anisotropic voxel size
"""
function correct_voxel_size!(data_binary::Array{Bool}, voxel_size::Tuple) 

    #determine smallest voxel size
    min_voxel_size = minimum(voxel_size)

    #get sizes of voxel corrected (stretched) data array
    compensated_data_size = ( Int( round( size(data_binary)[1]*voxel_size[1]/min_voxel_size ) ),
                            Int( round( size(data_binary)[2]*voxel_size[2]/min_voxel_size ) ),
                            Int( round( size(data_binary)[3]*voxel_size[3]/min_voxel_size ) ) )

    #loop through all dimensions to insert sclices along those directions that
    #need to be stretched
    for i in 1:3

        #check if array needs to be stretched along this direction
        if min_voxel_size !== voxel_size[i]
            println("dimension to be corrected: "*string(i))

            #find the nr of slices that need to be inserted
            nr_insertion_slices = compensated_data_size[i] - size(data_binary)[i]

            #determine indices of the compensated array that will be duplicates along given direction
            insertion_indices = Int.( round.( (compensated_data_size[i] / nr_insertion_slices) 
                                            .* (collect(1:nr_insertion_slices) .- 1/2 ) 
                                ))

            #loop through data array along current direction and insert slices at insertion indices
            for j in insertion_indices
                println("slice to be inserted: "*string(j))

                #println("shape of insertion slice: "*string(size(insertion_slice)))

                #insert slice at index j along current direction
                if i == 1
                    #reshape slice that will be inserted by adding a singleton dimension
                    insertion_slice = extend_dims(data_binary[j,:,:], i)

                    #insert slice
                    data_binary = cat(data_binary[1:j,:,:], insertion_slice, data_binary[j+1:end,:,:], dims=i )

                elseif i == 2
                    #reshape slice that will be inserted by adding a singleton dimension
                    insertion_slice = extend_dims(data_binary[:,j,:], i)

                    #insert slice
                    data_binary = cat(data_binary[:,1:j,:], insertion_slice, data_binary[:,j+1:end,:], dims=i )
                else
                    #reshape slice that will be inserted by adding a singleton dimension
                    insertion_slice = extend_dims(data_binary[:,:,j], i)

                    #insert slice
                    data_binary = cat(data_binary[:,:,1:j], insertion_slice, data_binary[:,:,j+1:end], dims=i )
                end

                #println("shape after insertion: "*string(size(data_binary)))

            end

        end

    end

    return data_binary 

end


"""
determine volume fraction within given data array
"""
function get_volume_fraction(data_binary::Array{Bool})
    
    volume_fraction = Statistics.mean(data_binary)

    return volume_fraction

end


"""
load data and get all its essential properties
"""
function get_data_essentials(data_binary::Array{Bool})

    #get total volume fraction
    volume_fract_tot = get_volume_fraction(data_binary)

    #get size and mean edge length of data array. Since data array is expected to be roughly 
    #cubic, mean should yield a sensible window length scale
    size_data = size(data_binary)
    mean_edge_length_data = Int(round( Statistics.mean(size_data) )) 

    #get the number of dimensions of data array
    nr_dimensions_data = ndims(data_binary)

    #warn if not 3d, since some functions only work for 3d data
    if nr_dimensions_data !== 3
        @warn "Data is not 3D. Most functions thus won't work!"
    end

    return [volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data]

end


"""
Convert a colorscale array to binary data
"""
function convert_colorscale_to_binary(data_colorscale::Array)

    #convert colorscale data to grayscale, then to float and then to binary data
    data_gray = Images.Gray.(data_colorscale)
    data_float = Float64.( data_gray ) 
    data_binary = Array(Bool.(round.( (1 / maximum( data_float ) ) .* data_float )))

    return data_binary

end


"""
load structure data, bring to binary form, correct asymmetric voxels and
if desired save to dictionary
"""
function get_structure_dict_from_colorscale(data_path_raw::String; 
    voxel_size::Tuple=(1,1,1), 
    label = "some structure",
    save_result::Bool=false, 
    save_path::String=raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\binary_data")

    #load colorscale structure data
    data_colorscale = FileIO.load(data_path_raw)

    #convert colorscale data to grayscale, then to float and then to binary data
    data_binary = convert_colorscale_to_binary(data_colorscale)
    
    #in case of anisotropic voxel size, correct anisotropy by stretching the array
    if voxel_size !== (1,1,1)
        data_binary_corrected_voxel_size = correct_voxel_size!(data_binary, voxel_size)

    else
        data_binary_corrected_voxel_size = data_binary
    end

    #get essential information about the structure data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = get_data_essentials(
                                                                data_binary_corrected_voxel_size)

    #save everything in dictionary
    structure_dict = Dict("data_binary" => data_binary_corrected_voxel_size, 
                            "volume_fract_tot" => volume_fract_tot, 
                            "size_data" => size_data, 
                            "mean_edge_length_data" => mean_edge_length_data, 
                            "nr_dimensions_data" => nr_dimensions_data,
                            "voxel_edge_length" => minimum(voxel_size) ,
                            "label" => label )

    #if desired, save corrected data
    if save_result
        save_dict_to_h5(structure_dict; save_path = save_path*"_structure.h5")

    end

    return structure_dict

end


"""
load 3d data from a stack of 2d images and bring to binary form
"""
function get_structure_dict_from_colorscale_stack(data_path_raw_prefix::String,
    data_path_raw_suffix::String,
    nr_images::Int64; 
    voxel_size::Tuple=(1,1,1), 
    label = "some structure",
    save_result::Bool=false, 
    save_path::String=raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\\binary_data")

    #load first image to get array dimensions
    image_1_colorscale = FileIO.load(data_path_raw_prefix*string(0)*data_path_raw_suffix)

    #initialize empty array where data will be stored in
    data_binary = Array{Bool}(undef, 
                                size(image_1_colorscale)..., 
                                nr_images)

    #loop through all images and save them to array
    for i in 1:nr_images

        #load colorscale structure data
        data_colorscale = FileIO.load(data_path_raw_prefix
                                    *string(i-1)
                                    *data_path_raw_suffix)

        #convert colorscale data to grayscale, then to float and then to binary data
        data_binary[:,:,i] = convert_colorscale_to_binary(data_colorscale)

        println("slice "*string(i)*" done")

    end

    #in case of anisotropic voxel size, correct anisotropy by stretching the array
    if voxel_size !== (1,1,1)
        data_binary_corrected_voxel_size = correct_voxel_size!(data_binary, voxel_size)

    else
        data_binary_corrected_voxel_size = data_binary
    end

    #get essential information about the structure data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = get_data_essentials(
                                                                data_binary_corrected_voxel_size)

    #save everything in dictionary
    structure_dict = Dict("data_binary" => data_binary_corrected_voxel_size, 
                            "volume_fract_tot" => volume_fract_tot, 
                            "size_data" => size_data, 
                            "mean_edge_length_data" => mean_edge_length_data, 
                            "nr_dimensions_data" => nr_dimensions_data,
                            "voxel_edge_length" => minimum(voxel_size) ,
                            "label" => label )

    #if desired, save corrected data
    if save_result
        save_dict_to_h5(structure_dict; save_path = save_path*"_structure.h5")

    end

    return structure_dict

end



"""
decompose measurement types in dictionary and store them as new keys
"""
function decompose_measurements_in_dict(dict::Dict)

    #loop through each key of the dict
    for (key, value) in dict

        #if key if of type Measurements.Measurement, save values and uncertainties
        #to two seperate keys and delete previous key
        if typeof(value) in [Vector{Measurements.Measurement},
                            Vector{Measurements.Measurement{Float64}}, 
                            Vector{Complex{Measurements.Measurement{Float64}}},
                            Array{Measurements.Measurement{Float64}, 3},
                            Array{Measurements.Measurement{Float64}, 3}, 
                            Array{Complex{Measurements.Measurement{Float64}}, 3}]

            dict[key*"_values"] = Measurements.value.( value )
            dict[key*"_uncertainties"] = Measurements.uncertainty.( value )
            delete!(dict, key)

        end
    end

    return dict

end


"""
decompose vectors of vectors into seperate vectors
"""
function decompose_vec_vecs_in_dict(dict::Dict)

    #loop through each key of the dict
    for (key, value) in dict

        #if key ends with _vec_vec
        if endswith(key, "_vec_vec")

            #cut of the last "vec" from key
            core_key = key[1:end-3]

            #loop through vector of vectors and save individual vectors
            for i in eachindex(value)
                dict[core_key*string(i)] = value[i]
            end

            #delete vector of vectors
            delete!(dict, key)

        end
    end

    return dict
end


"""
turn tuples into vectors
"""
function tuples_to_vectors_in_dict(dict::Dict)

    #loop through each key of the dict
    for (key, value) in dict

        #check if value is a tuple
        if typeof(value) in [Tuple{Int64, Int64}, Tuple{Int64, Int64, Int64}]

            #save vector to dict
            dict[key*"_tuple"] = collect(value)

            #delete vector of vectors
            delete!(dict, key)

        end
    end

    return dict
end


"""
save dict to H5 file.
The dicts can not be stored right away, because they contain some variables of type
Measurements.Measurement. These need to be decomposed into value and uncertainty first.
Also vectors of vectors need to be decomposed and the vectors have to be stored individually
and tuples need to be converted into vectors
"""
function save_dict_to_h5(dict::Dict;
                        save_path::String)

    #decompose those keys of the dict that are of measurement type
    decomposed_measurements_dict = decompose_measurements_in_dict(dict)

    #decompose vectors of vectors into seperate vectors
    decomposed_vec_vec_dict = decompose_vec_vecs_in_dict(decomposed_measurements_dict)

    #turn tuples into vectors
    saving_dict = tuples_to_vectors_in_dict(decomposed_vec_vec_dict)

    #save dict
    FileIO.save(save_path, saving_dict)

    return

end


"""
restore measurement types in dictionary
"""
function restore_measurement_types(dict::Dict)

    #loop through each key of the dict
    for (key, value) in dict

        #if key ends on "_values", then merge it with corresponding uncertainties
        #into measurement type
        if endswith(key, "_values")

            #get original key
            original_key = first(key, findfirst("_values", key)[1]-1)

            #get key of uncertainties
            uncertainty_key = original_key*"_uncertainties"

            #create vector of type Measurements.Measurement after checking whether
            #values are complex or real
            if typeof(value[1]) == ComplexF64
                measurement_vector = Complex.( Measurements.measurement.( real.( value ), real.( dict[uncertainty_key] ) ),
                                                Measurements.measurement.( imag.( value ), imag.( dict[uncertainty_key] ) )  )

            else
                measurement_vector = Measurements.measurement.( dict[key], dict[uncertainty_key] )

            end

            #save measurement vector to dict and delete value and uncertainty keys
            dict[original_key] = measurement_vector
            delete!(dict, key)
            delete!(dict, uncertainty_key)

        end
    end

    return dict

end


"""
restore vectors of vectors that were decomposed to save the dict
"""
function restore_vec_vecs_in_dict(dict::Dict)

    #loop through each key of the dict
    for (key, value) in dict

        #if key ends on "_1", then merge it with the other vectors
        #ending on _2 and _3
        if endswith(key, "_1")

            #get core key
            core_key = key[1:end-1]

            vec_vec = [value, dict[core_key*"2"], dict[core_key*"3"]]

            #save vec vec to dict
            dict[core_key*"vec"] = vec_vec

            #delete single vectors
            for i in 1:3
                delete!(dict, core_key*string(i))

            end

        end
    end

    return dict
end


"""
restore tuples in a dictionary
"""
function restore_tuples_in_dict(dict::Dict)
    
    #loop through each key of the dict
    for (key, value) in dict

        #if key ends on "_tuple", then convert it to a vector
        if endswith(key, "_tuple")

            #get core key
            dict[key[1:end-6]] = Tuple(value )
            
            delete!(dict, key)

        end
    end

    return dict

end


"""
load dict from h5 file.
When loaded, variables of type Measurements.Measurement whcih were decomposed are
restored again. Also vectors of vectors had to be decomposed for saving and are restored
"""
function load_h5_dict(dict_path)

    #save dict
    loaded_dict = FileIO.load(dict_path)

    #restore vectors of type Measurements.Measurement
    measurements_restored_dict = restore_measurement_types(loaded_dict)

    #restore vectors of vectors
    vec_vecs_restored_dict = restore_vec_vecs_in_dict(measurements_restored_dict)

    #restore tuples
    restored_dict = restore_tuples_in_dict(vec_vecs_restored_dict)

    return restored_dict

end


"""
modify the passed keys of a dictionary if they are not nothing
"""
function modify_keys_in_dict(dict::Dict, voxel_edge_length, label)

    if voxel_edge_length !== nothing
        dict["voxel_edge_length"] = voxel_edge_length
    end
    if label !== nothing
        dict["label"] = label
    end

    return dict

end



"""
function to calculate and save the following measures for a given 3d data set
- local volume fraction variance
- autocovariance function as a function of sampling distance (assuming an isotropic medium)
- spectral density as a function of sampling distance (assuming an isotropic medium)
- autocovariance function as a function of sampling vector (not assuming an isotropic medium)
- spectral density along 6 different directions (not assuming an isotropic medium)
"""
function save_statistical_measures(data_path::String,
                                    save_path::String;
                                    voxel_edge_length = nothing, 
                                    label = nothing,
                                    nr_sampling_distances = nothing,
                                    nr_measurements_per_distance = 10000,
                                    nr_window_sizes = 100,
                                    nr_measurements_per_direction = 1000,
                                    save_autocovariance_fct = true,
                                    save_spectral_density = true,
                                    save_local_volume_fraction_variance = true,
                                    save_autocovariance_fct_direction = true,
                                    save_spectral_density_along_directions = true,
                                    save_complete_autocovariance_fct_direction = true,
                                    save_spectral_density_array = true)

    #load structure dictionary which contains all essential information about the structure
    structure_dict = load_h5_dict(data_path)

    #set voxel edge length and label from structure dict if not specified in the arguments
    if voxel_edge_length === nothing
        voxel_edge_length = structure_dict["voxel_edge_length"]
    end
    if label === nothing
        label = structure_dict["label"]
    end

    #if not specified, determine nr of sampling distances for autocovariance function
    #and spectral density
    if nr_sampling_distances === nothing
        nr_sampling_distances = get_nr_sampling_distances(structure_dict["mean_edge_length_data"] )
    end

    #save autocovariance function if desired
    if save_autocovariance_fct

        #get autocovariance function as a function of sampling distance
        autocovariance_fct_isotrope_dict = get_autocovariance_fct_isotrope_by_sampling_distance_vec(
                                                structure_dict;
                                                nr_sampling_distances = nr_sampling_distances,
                                                nr_measurements_per_distance = nr_measurements_per_distance,
                                                save_result = true,
                                                save_path = save_path,
                                                voxel_edge_length = voxel_edge_length,
                                                label = label)

    end


    #save spectral density if desired
    if save_spectral_density

        #check if autocovariance fct was calculated within this function
        if !save_autocovariance_fct

            #if autocovariance function was not calculated within this function, check if it was
            #calculated before and if so, load the corresponding dictionary
            if isfile(save_path*"_autocovariance_fct.h5")

                #load autocovariance function dict
                autocovariance_fct_isotrope_dict = load_h5_dict(save_path*"_autocovariance_fct.h5")

            else
                autocovariance_fct_isotrope_dict = get_autocovariance_fct_isotrope_by_sampling_distance_vec(
                                                structure_dict;
                                                nr_sampling_distances = nr_sampling_distances,
                                                nr_measurements_per_distance = nr_measurements_per_distance,
                                                save_result = false,
                                                save_path = save_path,
                                                voxel_edge_length = voxel_edge_length,
                                                label = label)
            end
        end
    
        #calculate spectral_density
        spectral_density_isotrope_dict = get_spectral_density_isotrope_by_wavenumber_vec(
                structure_dict;
                nr_sampling_distances = length(autocovariance_fct_isotrope_dict["sampling_distance_vec"]),
                nr_measurements_per_distance = autocovariance_fct_isotrope_dict["nr_measurements_per_distance"],
                sampling_distance_vec = autocovariance_fct_isotrope_dict["sampling_distance_vec"],
                autocovariance_fct_dict = autocovariance_fct_isotrope_dict,
                save_result = true,
                save_path = save_path,
                voxel_edge_length = voxel_edge_length,
                label = label)
                    
    end


    #save local volume fraction if desired
    if save_local_volume_fraction_variance

        #determine local volume fraction variance vector by using measuring windows
        local_volume_fract_variance_dict = get_local_volume_fract_variance_by_window_vec(
                                                            structure_dict;
                                                            nr_window_sizes = nr_window_sizes,
                                                            window_positioning="random",
                                                            window_shape="spherical",
                                                            save_result = true,
                                                            save_path = save_path,
                                                            voxel_edge_length = voxel_edge_length,
                                                            label = label)
    
    end


    #save autocovariance function as a function of sampling direction
    if save_autocovariance_fct_direction

        #get autocovariance function array
        autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                            structure_dict;
                            nr_measurements_per_direction = nr_measurements_per_direction,
                            save_result = true,
                            save_path = save_path,
                            voxel_edge_length = voxel_edge_length,
                            label = label)
    
    end


    #save spectral density along certain directions if desired
    if save_spectral_density_along_directions

        #set vectors along which spectral density will be measured
        direction_vec_vec = [[1,0,0], 
                            [0,1,0],
                            [0,0,1],
                            (1/sqrt(2)) .* [1,-1,0],
                            (1/sqrt(3)) .* [1,1,1],
                            (1/sqrt(6)) .* [1,1,-2]]

        #set vector with names to save the files
        naming_vec = string.( [[1,0,0], 
                                [0,1,0],
                                [0,0,1],
                                [1,-1,0],
                                [1,1,1],
                                [1,1,-2]] )

        #check if autocovariance function was previously calculated in this function
        if !save_autocovariance_fct_direction

            #check if data can be loaded from file
            if isfile( save_path*"_autocovariance_fct_direction.h5" )
            
                #load autocovariance function per direction dict
                autocovariance_fct_direction_dict = load_h5_dict(save_path*"_autocovariance_fct_direction.h5")
            
            else
                #get autocovariance function array
                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                    structure_dict;
                    nr_measurements_per_direction = nr_measurements_per_direction,
                    save_result = false,
                    save_path = save_path,
                    voxel_edge_length = voxel_edge_length,
                    label = label)

            end
        end

        for i in eachindex(direction_vec_vec)

            #determine spectral density
            spectral_density_along_direction_dict = get_spectral_density_along_direction_by_wavenumber_vec(
                    structure_dict,
                    direction_vec_vec[i];
                    nr_measurements_per_direction = nr_measurements_per_direction,
                    sampling_distance_vec_vec = autocovariance_fct_direction_dict["sampling_distance_vec_vec"],
                    autocovariance_fct_dict = autocovariance_fct_direction_dict,
                    save_result = true,
                    save_path = save_path*"_"*naming_vec[i],
                    voxel_edge_length = voxel_edge_length,
                    label = label*" "*naming_vec[i])


        end
    end


    #save complete autocovariance function as a function of sampling direction for all spacial directions,
    #not only the half space considered previously
    if save_complete_autocovariance_fct_direction

        #check if autocovariance function was previously calculated in this function
        if !save_autocovariance_fct_direction

            #check if data can be loaded from file
            if isfile( save_path*"_autocovariance_fct_direction.h5" )
            
                #load autocovariance function per direction dict
                autocovariance_fct_direction_dict = load_h5_dict(save_path*"_autocovariance_fct_direction.h5")
            
            else
                #get autocovariance function array
                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                    structure_dict;
                    nr_measurements_per_direction = nr_measurements_per_direction,
                    save_result = false,
                    save_path = save_path,
                    voxel_edge_length = voxel_edge_length,
                    label = label)

            end
        end
    
        #calculate complete autocovariance function
        complete_autocovariance_fct_direction_dict = get_complete_autocovariance_fct_by_sampling_vec_array(
            autocovariance_fct_direction_dict;
            save_result = true,
            save_path = save_path)
                    
    end


    #save complete autocovariance function as a function of sampling direction for all spacial directions,
    #not only the half space considered previously
    if save_spectral_density_array

        #check if complete autocovariance function was previously calculated in this function
        if save_complete_autocovariance_fct_direction

            #check if data can be loaded from file
            if isfile( save_path*"_autocovariance_fct_direction_complete.h5" )
            
                #load complete autocovariance function per direction dict
                complete_autocovariance_fct_direction_dict = load_h5_dict(save_path*"_autocovariance_fct_direction_complete.h5")
            
            else

                #check if autocovariance function was previously calculated in this function
                if !save_autocovariance_fct_direction
                
                    #check if data can be loaded from file
                    if isfile( save_path*"_autocovariance_fct_direction.h5" )
                    
                        #load autocovariance function per direction dict
                        autocovariance_fct_direction_dict = load_h5_dict(save_path*"_autocovariance_fct_direction.h5")
                    
                    else
                        #get autocovariance function array
                        autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_vec_array(
                            structure_dict;
                            nr_measurements_per_direction = nr_measurements_per_direction,
                            save_result = false,
                            save_path = save_path,
                            voxel_edge_length = voxel_edge_length,
                            label = label)
                    
                    end
                end

                #get complete autocovariance function array
                complete_autocovariance_fct_direction_dict = get_complete_autocovariance_fct_by_sampling_vec_array(
                    autocovariance_fct_direction_dict;
                    save_result = false,
                    save_path = save_path)

            end
        end
    
        #calculate spectral_density
        spectral_density_array_dict = get_spectral_density_array_by_fft(complete_autocovariance_fct_direction_dict;
                                                                    save_result = true,
                                                                    save_path = save_path,
                                                                    voxel_edge_length = voxel_edge_length,
                                                                    label = label)
                    
    end

    return
    
end