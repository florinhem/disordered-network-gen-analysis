"""
Functions to calculate the autocovariance function and the spectral
density for voxelized networks with periodic boundary conditions
"""


"""
Get binary data for a spatial network where the bonds are represented by
a single line of voxels without a 'finite' bond radius
"""
function get_binary_data_from_spatial_network_bonds_only(graph_dict::Dict;
    voxel_edge_length::Float64 = 0.1)

    # generate array of zeros where data will be stored in
    data_binary = zeros(Bool, Int(round(graph_dict["supercell_edge_length"] / voxel_edge_length)), 
                            Int(round(graph_dict["supercell_edge_length"] / voxel_edge_length)), 
                            Int(round(graph_dict["supercell_edge_length"] / voxel_edge_length)))

    # loop through bonds in spatial network and set those voxels to 1 that lie
    # closer to the bond than the bond radius
    for bond in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        # get the bond vector
        direction_vec = graph_dict["spatial_network"][bond...]["vector"]

        # get the bond length
        bond_length = sqrt(graph_dict["spatial_network"][bond...]["distance_squared"])

        # get the number of voxels that are needed to represent the bond
        nr_voxels = Int( round( bond_length / voxel_edge_length ) )

        # get the voxel vector
        voxel_vector = direction_vec ./ nr_voxels

        # loop through voxels and set them to 1
        for i in 1:nr_voxels

            # get the position of the voxel by accounting for periodic boundary conditions
            voxel_position = (graph_dict["spatial_network"][bond[1]]["position"] 
                            .+ i * voxel_vector 
                            .+ graph_dict["supercell_edge_length"]) .% graph_dict["supercell_edge_length"]

            # get the index of the voxel
            voxel_index = (Int.( round.( voxel_position ./ voxel_edge_length .+ (1.000001/2) ) )
                            .- 1) .% size(data_binary) .+ 1

            # set the voxel to 1
            data_binary[voxel_index[1], voxel_index[2], voxel_index[3]] = 1

        end

    end

    return data_binary

end


"""
For a network whose bonds are represented by a single line of voxels,
give volume to the bonds by setting all voxels that lie within the bond radius to 1
"""
function add_volume_to_bonds(data_binary_bonds_only::Array{Bool,3};
    bond_radius::Float64 = 0.35 ,
    voxel_edge_length::Float64 = 0.1)

    # give volume to the bonds by setting all voxels that lie within the bond radius to 1
    data_binary = zeros(Bool, size(data_binary_bonds_only)...)

    for i in 1:size(data_binary_bonds_only)[1]
        for j in 1:size(data_binary_bonds_only)[2]
            for k in 1:size(data_binary_bonds_only)[3]

                if data_binary_bonds_only[i,j,k] == 1

                    # check window around current voxel
                    for l in (-Int( round( bond_radius / voxel_edge_length ) )
                                :Int( round( bond_radius / voxel_edge_length ) ))
                        for m in (-Int( round( bond_radius / voxel_edge_length ) )
                                    :Int( round( bond_radius / voxel_edge_length ) ))
                            for n in (-Int( round( bond_radius / voxel_edge_length ) )
                                        :Int( round( bond_radius / voxel_edge_length ) ))

                                # check if voxel is within bond radius
                                if  sqrt(l^2 + m^2 + n^2) * voxel_edge_length <= bond_radius

                                    # get index of current voxel by accounting for periodic boundary conditions
                                    cartesian_index_pbc = ( [i+l,j+m,k+n] .+ size(data_binary) .- 1 
                                                            ) .% size(data_binary) .+ 1

                                    # set voxel to 1
                                    data_binary[cartesian_index_pbc...] = 1

                                end

                            end
                        end
                    end

                end

            end
        end
    end

    return data_binary

end


"""
load data and get all its essential properties
"""
function get_binary_data_essentials(data_binary::Array{Bool})

    # get total volume fraction
    volume_fract_tot = Statistics.mean(data_binary)

    # get size and mean edge length of data array. Since data array is expected to be roughly 
    # cubic, mean should yield a sensible window length scale
    size_data = size(data_binary)
    mean_edge_length_data = Int(round( Statistics.mean(size_data) )) 

    # get the number of dimensions of data array
    nr_dimensions_data = ndims(data_binary)

    # warn if not 3d, since some functions only work for 3d data
    if nr_dimensions_data !== 3
        @warn "Data is not 3D. Most functions thus won't work!"
    end

    return [volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data]

end


"""
Get binary data for a spatial network with a given bond radius
"""
function get_binary_data_from_spatial_network(graph_dict::Dict;
    bond_radius::Float64 = 0.35,
    voxel_edge_length::Float64 = 0.1,
    save_path::String = raw"..\structures\random_networks\binary_structures\sample_name",
    filename::String = "some_structure",
    save_result::Bool=false)

    # get binary data for for only the bonds without a finite bond radius
    data_binary_bonds_only = get_binary_data_from_spatial_network_bonds_only(graph_dict; 
        voxel_edge_length = voxel_edge_length)

    # give volume to the bonds by setting all voxels that lie within the bond radius to 1
    data_binary = add_volume_to_bonds(data_binary_bonds_only; 
        bond_radius = bond_radius, 
        voxel_edge_length = voxel_edge_length)

    # get essential information about the structure data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = get_binary_data_essentials(
        data_binary)

    # save everything in structure dictionary
    structure_dict = Dict("data_binary" => data_binary, 
                            "volume_fract_tot" => volume_fract_tot, 
                            "size_data" => size_data, 
                            "mean_edge_length_data" => graph_dict["supercell_edge_length"], 
                            "nr_dimensions_data" => graph_dict["nr_dimensions"],
                            "voxel_edge_length" => voxel_edge_length ,
                            "label" => filename,
                            "coordination_nr" => graph_dict["coordination_nr"],
                            "nr_vertices" => graph_dict["nr_vertices"],
                            "bond_radius" => bond_radius )

    # if desired, save corrected data
    if save_result
        GU.save_dict_to_h5(structure_dict, save_path*"_structure.h5")

    end

    return structure_dict

end


"""
For each graph dict in a list of run and filenames, create binary data
with a given bond radius and voxel edge length
"""
function get_binary_data_from_spatial_network_single_thread(run_and_filename_chunk,
    graph_dicts_path::String,
    save_path::String;
    print_progress::Bool = true,
    print_lock = Threads.ReentrantLock(),
    bond_radius::Float64 = 0.35,
    voxel_edge_length::Float64 = 0.1)

    # loop through files
    for run_and_filename in run_and_filename_chunk

        # print current thread id and run/filename if desired
        if print_progress
            lock(print_lock) do
                Format.printfmtln("Thread {1} is creating structure dict of file {2}",
                    Threads.threadid(), run_and_filename)
            end
        end

        # load graph dict
        graph_dict = NG.load_graph_from_h5_and_gml(graph_dicts_path*run_and_filename)

        # get binary data for structure dict
        structure_dict = get_binary_data_from_spatial_network(graph_dict;
        bond_radius = bond_radius,
        voxel_edge_length = voxel_edge_length,
        save_path = save_path*run_and_filename,
        filename = run_and_filename[7:end],
        save_result=true)
            
    end

    return
end


"""
For each graph dict in a list of run and filenames, create binary data
with a given bond radius and voxel edge length
"""
function get_binary_data_from_spatial_network_multithreading(graph_dicts_path::String,
    save_path::String;
    print_progress::Bool = true,
    print_lock = Threads.ReentrantLock(),
    nr_runs::Int64 = 5,
    bond_radius::Float64 = 0.35,
    voxel_edge_length::Float64 = 0.1)

    # vector for run and filename
    run_and_filename_vec = []

    # loop through runs
    for i in 1:nr_runs

        graph_dicts_path_current_run = graph_dicts_path*"run_"*string(i)*"/"

        # get all filenames of graph dicts
        filenames = readdir(graph_dicts_path_current_run)
        filenames_graph_dicts = filter(filename -> endswith(filename, ".gml"), filenames)
        filenames_graph_dicts = [filename_graph_dict[1:end-4] 
                            for filename_graph_dict in filenames_graph_dicts]

        # append to list of all graph dict paths
        append!(run_and_filename_vec, "run_"*string(i)*"/" .* filenames_graph_dicts)

    end

    # split filenames into chunks for multi-threading
    run_and_filename_chunks = Iterators.partition(run_and_filename_vec, 
                                    length(run_and_filename_vec) ÷ Threads.nthreads())

    # run all filename chunk in parallel in different threads
    map(run_and_filename_chunks) do run_and_filename_chunk

        Threads.@spawn get_binary_data_from_spatial_network_single_thread(run_and_filename_chunk,
        graph_dicts_path,
        save_path;
        print_progress = print_progress,
        print_lock = print_lock,
        bond_radius = bond_radius,
        voxel_edge_length = voxel_edge_length)
    end

end


"""
Get vector of vectors containing the vector components at which 
the autocovariance function will be calculated
"""
function get_sampling_indices_vec_vec(size_data::Tuple)

    # determine the maximal sampling distances along the three axes
    max_sampling_distances = Int.( floor.( (size_data ) ./ 2 ))

    # get sampling distance vec vec
    # Along one axis (z direction is chosen here) only positive directions are considered,
    # because negative ones would yield redundant information
    sampling_indices_vec_vec = [
                collect(-max_sampling_distances[1]:max_sampling_distances[1]),
                collect(-max_sampling_distances[2]:max_sampling_distances[2]),
                collect(0:max_sampling_distances[3])]

    return sampling_indices_vec_vec

end


"""
Get array vectors out of three different vectors containing the
x, y and z components
"""
function get_vector_array(component_vec_vec::Vector)

    # initialize array where vectors will be stored
    vector_array = Array{typeof(component_vec_vec[1][1])}(undef, 
                                    length.(component_vec_vec)..., 3 )

    # save vectors to array
    for i in eachindex(component_vec_vec[1])
        for j in eachindex(component_vec_vec[2])
            for k in eachindex(component_vec_vec[3])

                # save current vector
                vector_array[i,j,k, :] = [component_vec_vec[1][i],
                                        component_vec_vec[2][j],
                                        component_vec_vec[3][k]]

            end
        end
    end

    return vector_array
end


"""
Get vector of vectors of sampled wavenumbers from fast fourier transform of
complete autocovariance function array
"""
function get_wavenumber_vec_vec(autocovariance_fct_direction_dict::Dict)

    # get vectors of wavenumbers along all three dimensions
    # since a real FFT is performed, the first dimension contains only positve wavenumbers
    # whereas second and third dimension contain positive and negative wavenumbers.
    # In order to bring them into a natural order, the fftshift needs to be done
    wavenumber_vec_vec = (2*pi/autocovariance_fct_direction_dict["voxel_edge_length"]) .* [
                                    FFTW.fftshift( FFTW.fftfreq( 
                size(autocovariance_fct_direction_dict["autocovariance_fct_array"])[1] ) ),
                                    FFTW.fftshift( FFTW.fftfreq( 
                size(autocovariance_fct_direction_dict["autocovariance_fct_array"])[2] ) ),
                                    FFTW.fftshift( FFTW.fftfreq( 
                size(autocovariance_fct_direction_dict["autocovariance_fct_array"])[3] ) ) ]

    # convert to Float64
    wavenumber_vec_vec_float = []

    for wavenumber_vec in wavenumber_vec_vec
        push!(wavenumber_vec_vec_float, Float64.( wavenumber_vec ))

    end

    return  wavenumber_vec_vec_float
    
end


"""
Get the autocovariance function for 3d media with periodic
boundary conditions
"""
function get_autocovariance_fct(sampling_vec::Vector{Int64},
    structure_dict::Dict)

    # use periodic boundary conditions to shift all entries of the data array
    # along the three dimensions as given by the sampling vector
    data_binary_shifted = circshift(structure_dict["data_binary"], sampling_vec)

    # calculate the contribution to the two point prob. fct. from these coodinates
    data_binary_product = structure_dict["data_binary"] .* data_binary_shifted

    # calculate 2 point prob. function
    two_point_prob_fct = Statistics.mean( data_binary_product )

    # determine autocovariance function
    autocovariance_fct = two_point_prob_fct - structure_dict["volume_fract_tot"]^2

    return autocovariance_fct
end


"""
Determine autocovariance function as a function of sampling direction for a 3d medium
with periodic boundary conditions
"""
function get_autocovariance_fct_by_sampling_indices_array(structure_dict::Dict;
                save_result = false,
                save_path = raw"..\analysis_data\sample_name",
                print_progress = false,
                thread_nr::Int64 = 0,
                print_lock = Threads.ReentrantLock())

    # get array of sampling vectors
    sampling_indices_vec_vec = get_sampling_indices_vec_vec(structure_dict["size_data"])
    sampling_indices_array = get_vector_array(sampling_indices_vec_vec)

    # get size of autocovariance fct array
    autocovariance_fct_array_size = size(sampling_indices_array)[1:3]

    # create vector where for each sampling distance the autocovariance function and its
    # uncertainty will be stored
    autocovariance_fct_array = Array{Float64}(undef, autocovariance_fct_array_size...)

    # for each sampling distance get autocovariance function and its uncertainty
    for i in 1:autocovariance_fct_array_size[1]
        for j in 1:autocovariance_fct_array_size[2]
            for k in 1:autocovariance_fct_array_size[3]

                autocovariance_fct_array[i,j,k] = get_autocovariance_fct(sampling_indices_array[i,j,k,:],
                                        structure_dict)
            end
        end

        # calculate and print progress
        if print_progress
            progress_percentage = i/autocovariance_fct_array_size[1]*100

            lock(print_lock) do
                Format.printfmtln("Current calculation progress thread nr {1:d}: {2:.1f} %", 
                    thread_nr, progress_percentage)
            end
        end

    end

    # extend the sampling distance vec along the third dimension, where due to the mirror
    # symmetry of the autocovariance fct only positive z values where considered
    complete_sampling_indices_vec_vec = [sampling_indices_vec_vec[1], 
                                            sampling_indices_vec_vec[2],
                                            vcat( .- sampling_indices_vec_vec[3][end:-1:2],
                                                sampling_indices_vec_vec[3]) ]

    # create array of sampling vectors out of sampling distances
    complete_sampling_indices_array = get_vector_array(complete_sampling_indices_vec_vec)

    # point mirror autocovariance fct to sampling distances with negative z component
    complete_autocovariance_fct_array = cat(dims=3, 
                                        autocovariance_fct_array[end:-1:1,end:-1:1,end:-1:2], 
                                        autocovariance_fct_array)

    
    # if the voxelized structure has an even number of voxels along any direction
    # one layer of the autocovariance function is redundant and needs to be removed
    if iseven(structure_dict["size_data"][1])
        complete_sampling_indices_vec_vec[1] = complete_sampling_indices_vec_vec[1][2:end]
        complete_sampling_indices_array = complete_sampling_indices_array[2:end,:,:,:]
        complete_autocovariance_fct_array = complete_autocovariance_fct_array[2:end,:,:]
    end
    if iseven(structure_dict["size_data"][2])
        complete_sampling_indices_vec_vec[2] = complete_sampling_indices_vec_vec[2][2:end]
        complete_sampling_indices_array = complete_sampling_indices_array[:,2:end,:,:]
        complete_autocovariance_fct_array = complete_autocovariance_fct_array[:,2:end,:]
    end
    if iseven(structure_dict["size_data"][3])
        complete_sampling_indices_vec_vec[3] = complete_sampling_indices_vec_vec[3][2:end]
        complete_sampling_indices_array = complete_sampling_indices_array[:,:,2:end,:]
        complete_autocovariance_fct_array = complete_autocovariance_fct_array[:,:,2:end]
    end
^
    # create dict
    autocovariance_fct_direction_dict = Dict("sampling_index_array" => complete_sampling_indices_array,
                        "sampling_index_vec_vec" => complete_sampling_indices_vec_vec,
                        "sampling_distance_array" => complete_sampling_indices_array .* structure_dict["voxel_edge_length"] ,
                        "sampling_distance_vec_vec" => complete_sampling_indices_vec_vec .* structure_dict["voxel_edge_length"] ,
                        "autocovariance_fct_array" => complete_autocovariance_fct_array,
                        "voxel_edge_length" => structure_dict["voxel_edge_length"] ,
                        "label" => structure_dict["label"])

    # save results if desired
    if save_result

        GU.save_dict_to_h5(autocovariance_fct_direction_dict, save_path*"_autocovariance_fct_direction.h5")

    end

    return autocovariance_fct_direction_dict
    
end


"""
For each structure dict in a list of run and filenames, calculate the 
autocovariance function as a function of direction
"""
function get_autocovariance_fct_direction_from_filenames_single_thread(run_and_filename_chunk,
    structure_dicts_path::String,
    save_path::String;
    print_progress::Bool = false,
    print_lock = Threads.ReentrantLock())

    # loop through files
    for run_and_filename in run_and_filename_chunk

        # check that file is structure dict
        if endswith(run_and_filename, "_structure.h5")

            # load evolution dict
            structure_dict = GU.load_h5_dict(structure_dicts_path*run_and_filename)

            # check if autocovariance dict already exists
            if !isfile(save_path*run_and_filename[1:end-13]*"_autocovariance_fct_direction.h5")

                # print current thread id and run/filename if desired
                if print_progress
                    lock(print_lock) do
                        Format.printfmtln("Thread {1} is handling file {2}",
                            Threads.threadid(), run_and_filename)
                    end
                end

                autocovariance_fct_direction_dict = get_autocovariance_fct_by_sampling_indices_array(
                    structure_dict;
                save_result = true,
                save_path = save_path*run_and_filename[1:end-13],
                print_progress = print_progress,
                thread_nr = Threads.threadid() ,
                print_lock = print_lock)
            end
        end
    end

    return
end


"""
For all structure dicts in a folder, calculate the autocovariance function as a function
of direction using multithreading
"""
function get_autocovariance_fct_direction_from_filenames_multithreading(structure_dicts_path;
    print_progress::Bool = false,
    save_path::String = "..\analysis_data\random_networks\\",
    nr_runs::Int64 = 5,
    print_lock = Threads.ReentrantLock())

    # vector for run and filename
    run_and_filename_vec = []

    # loop through runs
    for i in 1:nr_runs

        structure_dicts_path_current_run = structure_dicts_path*"run_"*string(i)*"/"

        filenames = readdir(structure_dicts_path_current_run)
        filenames_structure_dicts = filter(filename -> endswith(filename, "_structure.h5"), filenames)

        # append to list of all structure dict paths
        append!(run_and_filename_vec, "run_"*string(i)*"/" .* filenames_structure_dicts)

    end

    # split filenames into chunks for multi-threading
    run_and_filename_chunks = Iterators.partition(run_and_filename_vec, length(run_and_filename_vec) ÷ Threads.nthreads())

    # run all filename chunk in parallel in different threads
    map(run_and_filename_chunks) do run_and_filename_chunk

        Threads.@spawn get_autocovariance_fct_direction_from_filenames_single_thread(run_and_filename_chunk,
        structure_dicts_path,
        save_path;
        print_progress = print_progress,
        print_lock = print_lock)
    end

end


"""
Calculate spectral density from autocovariance fct for a 3d medium with 
periodic boundary conditions by means of Fast Fourier Transform
"""
function get_spectral_density_by_wavevector_array_fft(structure_dict::Dict;
    save_autocovariance_fct_direction_dict::Bool = false,
    save_result::Bool = false,
    save_path::String = raw"..\analysis_data\sample_name",
    autocovariance_fct_direction_dict::Dict = get_autocovariance_fct_by_sampling_indices_array(structure_dict;
            save_result = save_autocovariance_fct_direction_dict,
            save_path = save_path)
    )

    # determine fourier transform of autocovariance function values
    spectral_density_array_fft_output = FFTW.rfft(
            autocovariance_fct_direction_dict["autocovariance_fct_array"])

    # bring spectral densities into an order from negative to positive wavenumbers 
    spectral_density_array = FFTW.fftshift(spectral_density_array_fft_output, [2, 3])

    # point mirror spectral density to wavenumbers with negative x component
    spectral_density_array = cat(dims=1, 
            conj.(spectral_density_array[end:-1:2,end:-1:1,end:-1:1]), spectral_density_array)

    # get tuple of vectors of sampled wavenumbers
    wavenumber_vec_vec = get_wavenumber_vec_vec(autocovariance_fct_direction_dict)

    # get array of sampled sampled wavevectors 
    wavevector_array = get_vector_array(wavenumber_vec_vec)

    # create dict to save
    spectral_density_dict = Dict("wavevector_array" => wavevector_array,
        "wavenumber_vec_vec" => wavenumber_vec_vec,
        "spectral_density_array" => spectral_density_array,
        "voxel_edge_length" => autocovariance_fct_direction_dict["voxel_edge_length"],
        "label" => autocovariance_fct_direction_dict["label"])

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(spectral_density_dict), save_path*"_spectral_density_array.h5")

    end

    return spectral_density_dict
    
end


"""
Calculate angle averaged spectral density from 3d array of spectral density
"""
function get_spectral_density_angle_averaged(spectral_density_dict::Dict;
    gaussian_filter::Bool = true,
    gaussian_filter_sigma_x::Float64 = 2*pi/20, 
    gaussian_filter_filtered_data_x_step_length::Float64 = 2*pi/20,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name")

    # create a dictionary for angle averaged spectral density
    # with the length squared of the index vector as the key,
    # because it is proportional to the square of the wavenumber
    spectral_density_angle_averaged_dict = Dict{Int64, Vector{Float64}}()

    # determine origin vector of cartesian coordinates
    index_vector_origin = Int.(floor.( (
                size(spectral_density_dict["spectral_density_array"]) .+ 1 ) ./ 2) )

    # loop through Cartesian Indices
    for i in CartesianIndices(spectral_density_dict["spectral_density_array"])

        # determine actual index vector
        index_vector = i.I .- index_vector_origin

        # calculate length squared of index vector
        index_vector_length_squared = sum(index_vector.^2)

        # check if key exists in dictionary
        if index_vector_length_squared in keys(spectral_density_angle_averaged_dict)

            # add spectral density to dictionary
            push!(spectral_density_angle_averaged_dict[index_vector_length_squared], 
                abs(spectral_density_dict["spectral_density_array"][i]))

        else

            # add key and spectral density to dictionary
            spectral_density_angle_averaged_dict[index_vector_length_squared] = 
                [abs(spectral_density_dict["spectral_density_array"][i])]

        end

    end
    
    # initialize wavenumber vector and spectral density vector
    unfiltered_wavenumber_vec = Vector{Float64}()
    unfiltered_spectral_density_vec = Vector{Measurements.Measurement{Float64}}()

    # get lattice constant of reciprocal lattice
    reciprocal_lattice_constant = LinearAlgebra.norm(
        spectral_density_dict["wavevector_array"][(index_vector_origin .+ [1,0,0])...,:])

    # get vector of wavenumbers and angle averaged spectral density including
    # their uncertainty
    for key in keys(spectral_density_angle_averaged_dict)

        # get wavenumber from wavevector
        push!(unfiltered_wavenumber_vec, reciprocal_lattice_constant*sqrt(key))

        # get angle averaged spectral density
        push!(unfiltered_spectral_density_vec, 
            Measurements.measurement(Statistics.mean(spectral_density_angle_averaged_dict[key]),
            Statistics.std(spectral_density_angle_averaged_dict[key])))

    end

    # sort wavenumber vector and spectral density vector
    unfiltered_spectral_density_vec = unfiltered_spectral_density_vec[sortperm(unfiltered_wavenumber_vec)]
    sort!(unfiltered_wavenumber_vec)

    # create dict to save
    spectral_density_angle_averaged_dict = Dict{String, Any}(
                                "unfiltered_wavenumber_vec" => unfiltered_wavenumber_vec, 
                                "unfiltered_spectral_density_vec" => unfiltered_spectral_density_vec,
                                "label" => spectral_density_dict["label"])

    # apply gaussian filter if desired
    if gaussian_filter
        filtered_data_x, filtered_data_y = GU.gaussian_filter_1d(unfiltered_wavenumber_vec[2:end], 
            unfiltered_spectral_density_vec[2:end]; 
            sigma_x=gaussian_filter_sigma_x, 
            filtered_data_x_step_length=gaussian_filter_filtered_data_x_step_length)

        spectral_density_angle_averaged_dict["wavenumber_vec"] = filtered_data_x
        spectral_density_angle_averaged_dict["spectral_density_vec"] = filtered_data_y
        spectral_density_angle_averaged_dict["gaussian_filter_sigma_x"] = gaussian_filter_sigma_x
    end

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(spectral_density_angle_averaged_dict),
            save_path*"_spectral_density_angle_averaged.h5")

    end
                                            
    return spectral_density_angle_averaged_dict
end


"""
Define a anisotropy metric as the normalized variance of the spectral density at the
first peak of the spectral density of the diamond lattice
"""
function get_anisotropy_metric_from_spectral_density(
    spectral_density_angle_averaged_dict::Dict,
    diamond_spectral_density_peak_wavenumber::Float64 = 4.680517,
    diamond_spectral_density_peak_std::Float64 = 521.88398)

    # find the wavenumber that is the closest to the diamond peak wavenumber
    diamond_spectral_density_peak_wavenumber_index = argmin(abs.(spectral_density_angle_averaged_dict["wavenumber_vec"] .- diamond_spectral_density_peak_wavenumber))

    # get the spectral density standard deviation around the diamond peak
    peak_spectral_density_std = Measurements.uncertainty(
        spectral_density_angle_averaged_dict["spectral_density_vec"][diamond_spectral_density_peak_wavenumber_index])

    # normalize this standard deviation by the spectral density standard deviation of the diamond peak
    anisotropy_metric_from_spectral_density = peak_spectral_density_std / diamond_spectral_density_peak_std

    return anisotropy_metric_from_spectral_density
end


"""
Calculate the scaled intersection volume of 3d spheres of given radius R and
distance to the origin r
"""
function scaled_intersection_volume_3d(r_distance_to_origin::Float64
    , sphere_radius::Float64)

    # define heaviside function
    function heaviside(t)
        0.5 * (sign(t) + 1)
     end

    # calculate scaled intersection volume
    scaled_intersection_volume = ((1- 3/4 * r_distance_to_origin/sphere_radius 
            + 1/16 * (r_distance_to_origin/sphere_radius)^3) 
        * heaviside(2*sphere_radius - r_distance_to_origin))

    return scaled_intersection_volume
end


"""
From the autocovariance function, calculate the volume fraction variance
from equation 73 in 10.1016/j.physrep.2018.03.001
"""
function get_volume_fract_variance(autocovariance_fct_direction_dict::Dict;
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name")

    # get array of distances to the origin
    r_distance_to_origin_array = sqrt.(dropdims(sum(
        autocovariance_fct_direction_dict["sampling_distance_array"].^2, 
    dims=4), dims=4))

    # get supercell edge length
    voxelized_data_edge_length = (
        LinearAlgebra.norm(autocovariance_fct_direction_dict["sampling_distance_array"][1,1,1,:]
        .- autocovariance_fct_direction_dict["sampling_distance_array"][1,1,end,:]))

    # get vector of sphere radii
    sphere_radius_vec = collect(autocovariance_fct_direction_dict["voxel_edge_length"]/2:
    autocovariance_fct_direction_dict["voxel_edge_length"]/2
    :voxelized_data_edge_length/4
        )

    # initialize vector of volume fraction variances
    volume_fract_variance_vec = Vector{Float64}(undef, length(sphere_radius_vec))

    # initialize vector of volume fraction variances times window volume
    volume_fract_variance_times_window_volume_vec = Vector{Float64}(undef, length(sphere_radius_vec))

    # for each sphere radius, calculate the volume fraction variance
    for i in eachindex(sphere_radius_vec)

        # get volume of sphere
        sphere_volume = 4/3 * pi * sphere_radius_vec[i]^3

        # calculate integral of autocovariance function multiplied
        # by the scaled intersection volume
        volume_fract_variance_times_window_volume = sum(
            autocovariance_fct_direction_dict["autocovariance_fct_array"] .* 
            scaled_intersection_volume_3d.(r_distance_to_origin_array, sphere_radius_vec[i]))

        volume_fract_variance_times_window_volume_vec[i] = volume_fract_variance_times_window_volume

        # save volume fraction variance
        volume_fract_variance_vec[i] = 1/sphere_volume * volume_fract_variance_times_window_volume
    end

    # create dict to save
    volume_fract_variance_dict = Dict{String, Any}(
        "sphere_radius_vec" => sphere_radius_vec, 
        "volume_fract_variance_vec" => volume_fract_variance_vec,
        "volume_fract_variance_times_window_volume_vec" => volume_fract_variance_times_window_volume_vec,
        "label" => autocovariance_fct_direction_dict["label"])

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(volume_fract_variance_dict),
            save_path*"_volume_fraction_variance.h5")

    end

    return volume_fract_variance_dict
end


"""
Fit a function to the angle averaged spectral density
that consists of either two gaussians or a gaussian and an exponential decay
"""
function fit_spectral_density_angle_averaged(spectral_density_angle_averaged_dict::Dict)

    # define two possible fit functions that will be used depending on the behavior
    # of the spectral density at small wavenumbers
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
    
    # get wavenumber vector and spectral density vector
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

    return fit_param #  # fit_result 
end