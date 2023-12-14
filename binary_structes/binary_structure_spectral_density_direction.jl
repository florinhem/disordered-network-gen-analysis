"""
Functions to calculate the autocovariance function and the spectral
density for 3d data as a function of distance vector and momentum vector

This code requires the following Packages:
import Statistics   #for statistical operations like mean()
import Measurements    #for handling data with uncertainty and error propagation
import FFTW     #to calculate the Fast Fourier Transform
"""


"""
Get vector of vectors containing the vector components at which 
the autocovariance function will be calculated
"""
function get_sampling_distance_vec_vec(size_data::Tuple)

    #determine the maximal sampling distances along the three axes
    max_sampling_distances = Int.( ceil.( (size_data ) ./ 2 ))

    #get sampling distance vec vec
    #Along one axis (z direction is chosen here) only positive directions are considered,
    #because negative ones would yield redundant information
    sampling_distance_vec_vec = [
                collect(-max_sampling_distances[1]:max_sampling_distances[1]),
                collect(-max_sampling_distances[2]:max_sampling_distances[2]),
                collect(0:max_sampling_distances[3])]

    return sampling_distance_vec_vec

end



"""
Get 2 random coordinates that lie inside the structure at a given direction to each
other
"""
function get_random_coordinates_at_direction(sampling_vec::Vector,
                                size_data::Tuple)

    #get first random coordinate x1 such that second random coordinate at
    #given direction is still inside array
    x1_vec = Vector{Int64}(undef, 3)

    for i in 1:3
        if sampling_vec[i] < 0
            x1_vec[i] = rand(1-sampling_vec[i]:size_data[i])

        else
            x1_vec[i] = rand(1:size_data[i]-sampling_vec[i])
        end

    end

    #convert vectors to tuples
    x1 = Tuple(x1_vec)

    #get second coordinate at given direction
    x2 = Tuple(x1_vec .+ sampling_vec)

    return [x1, x2]

end



"""
Get array vectors out of three different vectors containing the
x, y and z components
"""
function get_vector_array(component_vec_vec::Vector)

    #initialize array where vectors will be stored
    vector_array = Array{typeof(component_vec_vec[1][1])}(undef, 
                                    length.(component_vec_vec)..., 3 )

    #save vectors to array
    for i in eachindex(component_vec_vec[1])
        for j in eachindex(component_vec_vec[2])
            for k in eachindex(component_vec_vec[3])

                #save current vector
                vector_array[i,j,k, :] = [component_vec_vec[1][i],
                                        component_vec_vec[2][j],
                                        component_vec_vec[3][k]]

            end
        end
    end

    return vector_array
end



"""
Get the autocovariance function for 3d media
"""
function get_autocovariance_fct(sampling_vec::Vector{Int64},
                                structure_dict::Dict;
                                nr_measurements_per_direction::Int64 = 1000)

    #initialize vector from which the two point prob. fct. will be calculated later
    two_point_prob_fct_summand_vec = Vector{Float64}(undef, nr_measurements_per_direction)

    #generate given number of random vectors to sample 2 point prob. function
    for i in 1:nr_measurements_per_direction

        #get two random coordinates at given sampling distance to each other
        x1, x2 = get_random_coordinates_at_direction(sampling_vec, structure_dict["size_data"])

        #calculate the contribution to the two point prob. fct. from these coodinates
        two_point_prob_fct_summand_vec[i] = structure_dict["data_binary"][x1...] * structure_dict["data_binary"][x2...]

    end

    #calculate ensemble averaged 2 point prob. function
    #the error is calculated according to the equation given at
    #https://en.wikipedia.org/wiki/Checking_whether_a_coin_is_fair
    #below the sentence
    #"In statistics, the estimate of a proportion of a sample "
    #(denoted by p) has a standard error given by
    two_point_prob_fct = Measurements.measurement( 
                                    Statistics.mean( two_point_prob_fct_summand_vec ),
                                    sqrt( structure_dict["volume_fract_tot"]^2 * (1 - structure_dict["volume_fract_tot"]^2) 
                                            / nr_measurements_per_direction  ) )

    #determine autocovariance function
    autocovariance_fct = two_point_prob_fct - structure_dict["volume_fract_tot"]^2

    return autocovariance_fct
    
end



"""
determine autocovariance function as a function of sampling direction.
The argument structure_dict is a dict with the following keys:
- data_binary
- volume_fract_tot
- size_data
"""
function get_autocovariance_fct_by_sampling_vec_array(structure_dict::Dict;
                sampling_distance_vec_vec = get_sampling_distance_vec_vec(structure_dict["size_data"]),
                sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
                nr_measurements_per_direction::Int64 = 1000,
                save_result = false,
                save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\sample_name",
                voxel_edge_length = nothing,
                label = nothing)

    #get size of autocovariance fct array
    autocovariance_fct_array_size = size(sampling_vec_array)[1:3]

    #create vector where for each sampling distance the autocovariance function and its
    #uncertainty will be stored
    autocovariance_fct_array = Array{Measurements.Measurement}(undef, autocovariance_fct_array_size...)

    #for each sampling distance get autocovariance function and its uncertainty
    for i in 1:autocovariance_fct_array_size[1]
        for j in 1:autocovariance_fct_array_size[2]
            for k in 1:autocovariance_fct_array_size[3]

                autocovariance_fct_array[i,j,k] = get_autocovariance_fct(sampling_vec_array[i,j,k,:],
                                        structure_dict;
                                        nr_measurements_per_direction = nr_measurements_per_direction)
            end
        end

        println("layer "*string(i)*" along x direction done" )

    end
^
    #create dict
    autocovariance_fct_dict = Dict("sampling_vec_array" => sampling_vec_array,
                        "sampling_distance_vec_vec" => sampling_distance_vec_vec,
                        "autocovariance_fct_array" => autocovariance_fct_array,
                        "nr_measurements_per_direction" => nr_measurements_per_direction,
                        "voxel_edge_length" => structure_dict["voxel_edge_length"] ,
                        "label" => structure_dict["label"])
                  
    #if desired, adjust voxel edge length and label
    autocovariance_fct_dict = modify_keys_in_dict(autocovariance_fct_dict, voxel_edge_length, label)

    #save results if desired
    if save_result

        save_dict_to_h5(copy(autocovariance_fct_dict);
                        save_path=save_path*"_autocovariance_fct_direction.h5")

    end

    return autocovariance_fct_dict
    
end


"""
Mirror the autocovariance function as a function of direction on the x-y-plane.
When autocovariance function is calculated as a function of sampling vectors, only
half of the 3d space with with z>=0 is considered. Due to the mirror symmetry of the
autocovariance function, this contains the full information. I think for the Fast
Fourier Transform the full data is needed however
"""
function get_complete_autocovariance_fct_by_sampling_vec_array(
        autocovariance_fct_direction_dict::Dict;
        save_result = false,
        save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\sample_name")

    #correct the sampling distance vec along the third dimension, where due to the mirror
    #symmetry of the autocovariance fct only positive z values where considered
    complete_sampling_distance_vec_vec = [autocovariance_fct_direction_dict["sampling_distance_vec_vec"][1], 
                                            autocovariance_fct_direction_dict["sampling_distance_vec_vec"][2],
                                            vcat(
                    .- reverse(autocovariance_fct_direction_dict["sampling_distance_vec_vec"][3][2:end]),
                    autocovariance_fct_direction_dict["sampling_distance_vec_vec"][3]
                                                ) ]

    #create array of sampling vectors out of sampling distances
    complete_sampling_vec_array = get_vector_array(complete_sampling_distance_vec_vec)

    #create an autocovariance function array that is mirrored in the first two dimensions
    autocovariance_fct_array_mirrored = reverse( reverse( 
                            reverse(autocovariance_fct_direction_dict["autocovariance_fct_array"], dims=1), 
                                                        dims=2), dims=3) 

    #get mirrored autocovariance function
    complete_autocovariance_fct_array = cat(autocovariance_fct_array_mirrored[:,:,1:end-1],
                                            autocovariance_fct_direction_dict["autocovariance_fct_array"],
                                            dims = 3)

    
    #create dict to save
    complete_autocovariance_fct_dict = Dict("sampling_vec_array" => complete_sampling_vec_array,
                    "sampling_distance_vec_vec" => complete_sampling_distance_vec_vec,
                    "autocovariance_fct_array" => complete_autocovariance_fct_array,
                    "nr_measurements_per_direction" 
                            => autocovariance_fct_direction_dict["nr_measurements_per_direction"],
                    "voxel_edge_length" => autocovariance_fct_direction_dict["voxel_edge_length"],
                    "label" => autocovariance_fct_direction_dict["label"])

    #save results if desired
    if save_result
        save_dict_to_h5(copy(complete_autocovariance_fct_dict);
                        save_path=save_path*"_autocovariance_fct_direction_complete.h5")

    end

    return complete_autocovariance_fct_dict
end


"""
if the autocovariance function was measured as a function of all possible directions for
a single unit cell of a perfectly periodic medium (e.g. one created from a nodal equation),
extrapolate the autocovariance function of many unit cells (50x50x50 as a standard).
The dictionary that is passed as an argument, should contain the following keys:
- sampling_vec_array
- sampling_distance_vec_vec
- autocovariance_fct_array
- voxel_edge_length
- label
"""
function extrapolate_periodic_data_autocovariance_fct_by_sampling_vec_array(
            autocovariance_fct_direction_dict::Dict;
            size_data_single_unit_cell::Tuple = (50,50,50),
            nr_unit_cells=10,
            save_result = false,
            save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\sample_name")

    #get vector of sampling distances along the three coordinate axes for extrapolated data
    extrapolated_sampling_distance_vec_vec = get_sampling_distance_vec_vec(size_data_single_unit_cell 
                                                                            .* nr_unit_cells)
    
    #get the corresponding array of vectors, where along the fourth dimension the three vector
    #components are stored
    extrapolated_sampling_vec_array = get_vector_array(extrapolated_sampling_distance_vec_vec)

    
    #get size of extrapolated autocovariance fct array
    extrapolated_autocovariance_fct_array_size = size(extrapolated_sampling_vec_array)[1:3]

    #initialize the extrapolated autocovariance fct array
    extrapolated_autocovariance_fct_array = Array{Measurements.Measurement}(undef, 
                                                        extrapolated_autocovariance_fct_array_size...)

    
    #determine the maximal sampling distances inside the single unit cell along the three axes
    max_sampling_distances = Int.( ceil.( (size_data_single_unit_cell ) ./ 2 ))

    #loop through the three coordinate axes
    for i in extrapolated_sampling_distance_vec_vec[1]
        for j in extrapolated_sampling_distance_vec_vec[2]
            for k in extrapolated_sampling_distance_vec_vec[3]

                #convert the current vector into a vector within one unit cell
                vector_within_unit_cell = ( mod.( [i,j,k] .+ max_sampling_distances, 
                                                                size_data_single_unit_cell) 
                                            .- max_sampling_distances )
                
                #if the z component of that vector is negative, it needs to be flipped, because
                #due to its mirror symmetry, the autocovariance function was only calculated for 
                #positive z components
                if vector_within_unit_cell[3] < 0
                    vector_within_unit_cell .*= (-1)
                end

                #find the entries of this vector in the sampling_distance_vec_vec of the single unit cell
                current_indices_within_unit_cell = findfirst.( .==( vector_within_unit_cell ), 
                                        autocovariance_fct_direction_dict["sampling_distance_vec_vec"]  )

                #save the corresponding autocovariance function values to the array of the extrapolated 
                #autocovariance function
                extrapolated_autocovariance_fct_array[i + extrapolated_sampling_distance_vec_vec[1][end] + 1,
                                                        j + extrapolated_sampling_distance_vec_vec[2][end] + 1,
                                                        k + 1] = (
                    autocovariance_fct_direction_dict["autocovariance_fct_array"][current_indices_within_unit_cell...]
                    )

            end
        end
    end

    #create dict to save
    extrapolated_autocovariance_fct_dict = Dict("sampling_vec_array" => extrapolated_sampling_vec_array,
                    "sampling_distance_vec_vec" => extrapolated_sampling_distance_vec_vec,
                    "autocovariance_fct_array" => extrapolated_autocovariance_fct_array,
                    "nr_measurements_per_direction" 
                            => autocovariance_fct_direction_dict["nr_measurements_per_direction"],
                    "voxel_edge_length" => autocovariance_fct_direction_dict["voxel_edge_length"],
                    "label" => autocovariance_fct_direction_dict["label"])

    #save results if desired
    if save_result
        save_dict_to_h5(copy(extrapolated_autocovariance_fct_dict);
                        save_path=save_path*"_autocovariance_fct_direction.h5")

    end

    return extrapolated_autocovariance_fct_dict

end


"""
extract autocovariance function vector along given direction
"""
function get_autocovariance_fct_along_direction_vec(direction_vec::Vector, 
                                            autocovariance_fct_array::Array)

    #if z coordinate is negative, mirror the vector, because autocovariance function
    #was only measured for positive z coordinates
    if direction_vec[3] < 0
        direction_vec = .- direction_vec
    end

    #get size of autocovariance_fct_array
    size_autocovariance_fct_array = size(autocovariance_fct_array)

    #set current position to indices of autocovariance_fct_array where autocovariance function
    #at sampling vector [0,0,0] was measured
    current_position = [ Int.( (size_autocovariance_fct_array[1:2] .+ 1) ./2 ) ..., 1]

    #initialize autocovariance fct vec
    autocovariance_fct_along_direction_vec = []

    #starting at center, walk through autocovariance_fct_array to get the desired entries
    while prod( current_position .<= size_autocovariance_fct_array  )*prod( current_position .>= [1,1,1]  ) == 1

        #add autocovariance function value to vector
        push!(autocovariance_fct_along_direction_vec, autocovariance_fct_array[Int.( round.( current_position ) ) ...] )

        #go one step along given direction
        current_position = current_position .+ direction_vec 
    end

    return autocovariance_fct_along_direction_vec

end


"""
get vector of sampling distances along the given direction
"""
function get_sampling_distance_along_direction_vec(direction_vec, 
                autocovariance_fct_along_direction_vec::Vector)
    
    #determine geometrical length of direction vector
    sampling_distance = sqrt(sum( direction_vec .^2 ))

    #get vector of sampling distances
    sampling_distance_vec = (collect(0:length(autocovariance_fct_along_direction_vec)-1) 
                                .* sampling_distance)

    return sampling_distance_vec

end


"""
get spectral density 
"""
function get_spectral_density_along_direction(wavenumber::Float64, 
                            sampling_distance_vec::Vector, 
                            autocovariance_fct_along_direction_vec::Vector)

    #determine sampling distance
    sampling_distance = sampling_distance_vec[2] - sampling_distance_vec[1]

    #calculate fourier transform for given wavenumber
    spectral_density = (1/sampling_distance 
                            * ( autocovariance_fct_along_direction_vec[1] 
                                + 2 * sum( autocovariance_fct_along_direction_vec[2:end] 
                                                .* cos.( wavenumber .* sampling_distance_vec[2:end] )  ) 
                                )
                        )

    return spectral_density

end



"""
Calculate spectral density from autocovariance fct along a given direction.
Somehow this function does not work yet, maybe because the calculation of exponentials for
Measurement type data simply takes too long.
"""
function get_spectral_density_along_direction_by_wavenumber_vec(structure_dict::Dict,
                direction_vec::Vector;
                nr_measurements_per_direction::Int64 = 1000,
                sampling_distance_vec_vec = get_sampling_distance_vec_vec(structure_dict["size_data"]),
                sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
                autocovariance_fct_dict = get_autocovariance_fct_by_sampling_vec_array(structure_dict;
                                            sampling_vec_array = sampling_vec_array,
                                            nr_measurements_per_direction = nr_measurements_per_direction),
                save_result = false,
                save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\sample_name_direction",
                voxel_edge_length = nothing,
                label = nothing)
    
    #extract autocovariance function vector along given direction
    autocovariance_fct_along_direction_vec = get_autocovariance_fct_along_direction_vec(direction_vec, 
                                                                autocovariance_fct_dict["autocovariance_fct_array"])

    #get sampling distances along given direction
    sampling_distance_along_direction_vec = get_sampling_distance_along_direction_vec(direction_vec, 
                                                                autocovariance_fct_along_direction_vec)

    #get vector of wavenumbers where spectral density will be calculated
    wavenumber_vec = get_wavenumber_vec(sampling_distance_along_direction_vec)

    #initialize spectral density vector
    spectral_density_vec = Vector{Measurements.Measurement}(undef, length(wavenumber_vec))

    #for each wavenumber, determine spectral density
    for i in eachindex(wavenumber_vec)
        spectral_density_vec[i] = get_spectral_density_along_direction(wavenumber_vec[i], 
                                            sampling_distance_along_direction_vec, 
                                            autocovariance_fct_along_direction_vec)

    end

    #create dict to save
    spectral_density_dict = Dict("wavenumber_vec" => wavenumber_vec,
                        "spectral_density_vec" => spectral_density_vec,
                        "nr_measurements_per_direction" => nr_measurements_per_direction,
                        "direction_vec" => direction_vec,
                        "voxel_edge_length" => autocovariance_fct_dict["voxel_edge_length"],
                        "label" => autocovariance_fct_dict["label"])

    #if desired, adjust voxel edge length and label
    spectral_density_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    #save results if desired
    if save_result
        save_dict_to_h5(copy(spectral_density_dict);
                        save_path=save_path*"_spectral_density_direction.h5")

    end

    return spectral_density_dict
    
end


"""
determine spectral density as a function of wavevector.
Unfortunately this function takes very long (couple of hours for a single wavevector)
"""
function get_spectral_density(wavevector::Vector{Float64}, 
                                sampling_distance_vec_vec::Vector, 
                                autocovariance_fct_array::Array)

    #initialize complex spectral density
    spectral_density = 0 + 0im

    #create an autocovariance function array that is mirrored in the first two dimensions
    autocovariance_fct_array_mirrored = reverse( reverse(autocovariance_fct_array, dims=1), dims=2) 

    for i in eachindex(sampling_distance_vec_vec[1])
        for j in eachindex(sampling_distance_vec_vec[2])

            #initialize complex spectral density sum along z direction
            spectral_density_z_component_sum = 0 + 0im

            for k in eachindex(sampling_distance_vec_vec[3])[2:end]

                spectral_density_z_component_sum += (exp(im*wavevector[3]*sampling_distance_vec_vec[3][k])
                                                            *autocovariance_fct_array_mirrored[i,j,k]
                                                    + exp(-im*wavevector[3]*sampling_distance_vec_vec[3][k])
                                                            *autocovariance_fct_array[i,j,k]   )
                                    
            end

            #add sum along z direction and other terms to sum along x and y directions
            spectral_density += (exp(-im*wavevector[1]*sampling_distance_vec_vec[1][i]) 
                                    * exp(-im*wavevector[2]*sampling_distance_vec_vec[2][j]) 
                                    * ( autocovariance_fct_array[i,j,1] + spectral_density_z_component_sum )
                                    )

            println("sampling distance "*string(sampling_distance_vec_vec[2][j])*" along y done")
        end

        println("sampling distance "*string(sampling_distance_vec_vec[1][i])*" along x done")
    end

    #multiply by the inverse of the sampling volume
    spectral_density *= 1/( (sampling_distance_vec_vec[1][2] - sampling_distance_vec_vec[1][1])
                            * (sampling_distance_vec_vec[2][2] - sampling_distance_vec_vec[2][1])
                            * (sampling_distance_vec_vec[3][2] - sampling_distance_vec_vec[3][1]) )


    return spectral_density
end



"""
determine spectral density as a function of wavevector.
Unfortunately in the moment this function takes forever (weeks?).
Instead, I think I need to use the FFT
"""
function get_spectral_density_by_wavevector_array(
                structure_dict::Dict;
                nr_measurements_per_direction::Int64 = 1000,
                sampling_distance_vec_vec = get_sampling_distance_vec_vec(structure_dict["size_data"]),
                sampling_vec_array = get_vector_array(sampling_distance_vec_vec),
                autocovariance_fct_dict = get_autocovariance_fct_by_sampling_vec_array(structure_dict;
                                            sampling_vec_array = sampling_vec_array,
                                            nr_measurements_per_direction = nr_measurements_per_direction),
                save_result = false,
                save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\sample_name",
                voxel_edge_length = nothing,
                label = nothing)

    #correct the sampling distance vec along the third dimension, where due to the mirror
    #symmetry of the autocovariance fct only positive z values where considered
    corrected_sampling_distance_vec_vec = [
                            autocovariance_fct_dict["sampling_distance_vec_vec"][1], 
                            autocovariance_fct_dict["sampling_distance_vec_vec"][2],
                            vcat(.- reverse(autocovariance_fct_dict["sampling_distance_vec_vec"][3][2:end]),
                                autocovariance_fct_dict["sampling_distance_vec_vec"][3]) ]
    
    #get vector of sampled wavenumbers along the three coordinate directions
    wavenumber_vec_vec = []

    for i in 1:3

        push!(wavenumber_vec_vec, get_wavenumber_vec(corrected_sampling_distance_vec_vec[i]) )

    end

    #create array of wavevectors
    wavevector_array = get_vector_array(wavenumber_vec_vec)

    #determine size of spectral density array
    spectral_density_array_size = size(wavevector_array)[1:3]

    #initialize spectral density array 
    spectral_density_array = Array{Complex{Measurements.Measurement}}(undef, spectral_density_array_size...)

    #for each wavevector, determine the spectral density
    for i in eachindex(wavenumber_vec_vec[1])
        for j in eachindex(wavenumber_vec_vec[2])
            for k in eachindex(wavenumber_vec_vec[3])

                spectral_density_array[i,j,k] = get_spectral_density(
                                                    wavevector_array[i,j,k,:], 
                                                    autocovariance_fct_dict["sampling_distance_vec_vec"], 
                                                    autocovariance_fct_dict["autocovariance_fct_array"])

            end
        end

        println("wavenumber "*string(wavenumber_vec_vec[1][i]*" along x done") )
    end

    #create dict to save
    spectral_density_dict = Dict("wavevector_array" => wavevector_array,
                                "wavenumber_vec_vec" => wavenumber_vec_vec,
                                "spectral_density_array" => spectral_density_array,
                                "nr_measurements_per_direction" 
                                        => autocovariance_fct_dict["nr_measurements_per_direction"],
                                "voxel_edge_length" => autocovariance_fct_dict["voxel_edge_length"],
                                "label" => autocovariance_fct_dict["label"])

    #if desired, adjust voxel edge length and label
    spectral_density_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    #save results if desired
    if save_result
        save_dict_to_h5(copy(spectral_density_dict);
                        save_path=save_path*"_spectral_density_array.h5")

    end

    return spectral_density_dict

end



"""
Get vector of vectors of sampled wavenumbers
"""
function get_wavenumber_vec_vec(complete_autocovariance_fct_array_values::Array)

    #get vectors of wavenumbers along all three dimensions
    #since a real FFT is performed, the first dimension contains only positve wavenumbers
    #whereas second and third dimension contain positive and negative wavenumbers.
    #In order to bring them into a natural order, the fftshift needs to be done
    wavenumber_vec_vec = (2*pi) .* [FFTW.rfftfreq( 
                                    size(complete_autocovariance_fct_array_values)[1] ),
                                    FFTW.fftshift( FFTW.fftfreq( 
                                        size(complete_autocovariance_fct_array_values)[2] ) ),
                                    FFTW.fftshift( FFTW.fftfreq( 
                                        size(complete_autocovariance_fct_array_values)[3] ) ) ]

    #convert to Float64
    wavenumber_vec_vec_float = []

    for wavenumber_vec in wavenumber_vec_vec
        push!(wavenumber_vec_vec_float, Float64.( wavenumber_vec ))

    end

    return  wavenumber_vec_vec_float
    
end


"""
Calculate spectral density from autocovariance fct by means of Fast Fourier
Transform
"""
function get_spectral_density_array_by_fft(complete_autocovariance_fct_direction_dict::Dict;
            save_result = false,
            save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\sample_name",
            voxel_edge_length = nothing,
            label = nothing)

    #get values of autocovariance function
    complete_autocovariance_fct_array_values = Measurements.value.( 
                complete_autocovariance_fct_direction_dict["autocovariance_fct_array"] )

    #determine fourier transform of autocovariance function values
    spectral_density_array_fft_output = FFTW.rfft(complete_autocovariance_fct_array_values)

    #bring spectral densities into an order from negative to positive wavenumbers 
    spectral_density_array = FFTW.fftshift(spectral_density_array_fft_output, [2, 3])

    #get tuple of vectors of sampled wavenumbers
    wavenumber_vec_vec = get_wavenumber_vec_vec(complete_autocovariance_fct_array_values)

    #get array of sampled sampled wavevectors 
    wavevector_array = get_vector_array(wavenumber_vec_vec)

    #create dict to save
    spectral_density_dict = Dict("wavevector_array" => wavevector_array,
        "wavenumber_vec_vec" => wavenumber_vec_vec,
        "spectral_density_array" => spectral_density_array,
        "voxel_edge_length" => complete_autocovariance_fct_direction_dict["voxel_edge_length"],
        "label" => complete_autocovariance_fct_direction_dict["label"])

    #if desired, adjust voxel edge length and label
    spectral_density_dict = modify_keys_in_dict(spectral_density_dict, voxel_edge_length, label)

    #save results if desired
    if save_result
        save_dict_to_h5(copy(spectral_density_dict);
                        save_path=save_path*"_spectral_density_array.h5")

    end

    return spectral_density_dict
    
end