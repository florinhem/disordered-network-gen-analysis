"""
Functions for general IO utilities
"""


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
function load_h5_dict(dict_path::String)

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