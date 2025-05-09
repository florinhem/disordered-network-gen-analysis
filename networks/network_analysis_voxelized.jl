"""
Functions to calculate the autocovariance function and the spectral density for
voxelized networks with periodic boundary conditions
"""


"""
Get binary data for a spatial network where the bonds are represented by a
single line of voxels without a 'finite' bond radius
"""
function get_binary_data_from_spatial_network_bonds_only(
    spatial_network::MetaGraphsNext.MetaGraph;
    voxel_edge_length::Float64 = 0.1)

    # generate array of zeros where data will be stored in
    data_binary = zeros(Bool, 
        Int(round(spatial_network[]["supercell_edge_length"] 
            / voxel_edge_length)), 
        Int(round(spatial_network[]["supercell_edge_length"]
            / voxel_edge_length)), 
        Int(round(spatial_network[]["supercell_edge_length"] 
            / voxel_edge_length)))

    # loop through bonds in spatial network and set those voxels to 1 that lie
    # closer to the bond than the bond radius
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        # get the bond vector
        direction_vec = spatial_network[bond...]["vector"]

        # get the bond length
        bond_length = sqrt(spatial_network[bond...]["distance_squared"])

        # get the number of voxels that are needed to represent the bond
        nr_voxels = Int( round( bond_length / voxel_edge_length ) )

        # get the voxel vector
        voxel_vector = direction_vec ./ nr_voxels

        # loop through voxels and set them to 1
        for i in 1:nr_voxels

            # get the position of the voxel by accounting for periodic boundary
            # conditions
            voxel_position = ((spatial_network[bond[1]]["position"] 
                    .+ i * voxel_vector 
                    .+ spatial_network[]["supercell_edge_length"]) 
                .% spatial_network[]["supercell_edge_length"])

            # get the index of the voxel
            voxel_index = (Int.( round.( voxel_position ./ voxel_edge_length 
                    .+ (1.000001/2) ) )
                .- 1) .% size(data_binary) .+ 1

            # set the voxel to 1
            data_binary[voxel_index[1], voxel_index[2], voxel_index[3]] = 1

        end

    end

    return data_binary
end


"""
For a network whose bonds are represented by a single line of voxels, give
volume to the bonds by setting all voxels that lie within the bond radius to 1
"""
function add_volume_to_bonds(
    data_binary_bonds_only::Array{Bool,3};
    bond_radius::Float64 = 0.35 ,
    voxel_edge_length::Float64 = 0.1)

    # give volume to the bonds by setting all voxels that lie within the bond
    # radius to 1
    data_binary = zeros(Bool, size(data_binary_bonds_only)...)

    for i in 1:size(data_binary_bonds_only)[1]
        for j in 1:size(data_binary_bonds_only)[2]
            for k in 1:size(data_binary_bonds_only)[3]

                if data_binary_bonds_only[i,j,k] == 1

                    # check window around current voxel
                    for l in (-Int( round( bond_radius / voxel_edge_length ) )
                        :Int( round( bond_radius / voxel_edge_length ) ))
                        for m in (-Int( round( bond_radius 
                                / voxel_edge_length ) )
                            :Int( round( bond_radius / voxel_edge_length ) ))
                            for n in (-Int( round( bond_radius 
                                    / voxel_edge_length ) )
                                :Int( round( bond_radius 
                                    / voxel_edge_length ) ))

                                # check if voxel is within bond radius
                                if  (sqrt(l^2 + m^2 + n^2) 
                                    * voxel_edge_length <= bond_radius)

                                    # get index of current voxel by accounting 
                                    # for periodic boundary conditions
                                    cartesian_index_pbc = ( [i+l,j+m,k+n] 
                                            .+ size(data_binary) .- 1 
                                        ) .% size(data_binary) .+ 1

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

    # get size and mean edge length of data array. Since data array is expected
    # to be roughly cubic, mean should yield a sensible window length scale
    size_data = size(data_binary)
    mean_edge_length_data = Int(round( Statistics.mean(size_data) )) 

    nr_dimensions_data = ndims(data_binary)

    if nr_dimensions_data !== 3
        @warn "Data is not 3D. Most functions thus won't work!"
    end

    return [volume_fract_tot, size_data, mean_edge_length_data, 
        nr_dimensions_data]
end


"""
Get binary data for a spatial network with a given bond radius
"""
function get_binary_data_from_spatial_network(
    spatial_network::MetaGraphsNext.MetaGraph;
    bond_radius::Float64 = 0.35,
    voxel_edge_length::Float64 = 0.1,
    save_path::String 
        = raw"..\structures\random_networks\binary_structures\sample_name",
    filename::String = "some_structure",
    save_result::Bool=false)

    # get binary data for for only the bonds without a finite bond radius
    data_binary = get_binary_data_from_spatial_network_bonds_only(
        spatial_network; 
        voxel_edge_length = voxel_edge_length)

    # give volume to the bonds by setting all voxels that lie within the bond
    # radius to 1
    if bond_radius > voxel_edge_length/2
        data_binary = add_volume_to_bonds(data_binary; 
            bond_radius = bond_radius, 
            voxel_edge_length = voxel_edge_length)
    end

    # get essential information about the structure data
    volume_fract_tot, size_data, mean_edge_length_data, nr_dimensions_data = (
        get_binary_data_essentials(data_binary))

    # save everything in structure dictionary
    structure_dict = Dict("data_binary" => data_binary, 
        "volume_fract_tot" => volume_fract_tot, 
        "size_data" => size_data, 
        "mean_edge_length_data" => spatial_network[]["supercell_edge_length"], 
        "nr_dimensions_data" => spatial_network[]["nr_dimensions"],
        "voxel_edge_length" => voxel_edge_length ,
        "label" => filename,
        "nr_vertices" => spatial_network[]["nr_vertices"],
        "bond_radius" => bond_radius )

    # if desired, save corrected data
    if save_result
        GU.save_dict_to_h5(structure_dict, save_path*"_structure.h5")
    end

    return structure_dict
end


"""
For each network dict in a list of run and filenames, create binary data with a
given bond radius and voxel edge length
"""
function get_binary_data_from_spatial_network_single_thread(
    run_and_filename_chunk,
    spatial_networks_path::String,
    save_path::String;
    print_progress::Bool = true,
    print_lock = Threads.ReentrantLock(),
    bond_radius::Float64 = 0.35,
    voxel_edge_length::Float64 = 0.1)

    for run_and_filename in run_and_filename_chunk

        # print current thread id and run/filename if desired
        if print_progress
            lock(print_lock) do
                Format.printfmtln("Thread {1} is creating structure dict of file {2}",
                    Threads.threadid(), run_and_filename)
            end
        end

        spatial_network = NG.load_spatial_network_from_gml(
            spatial_networks_path*run_and_filename*".gml")

        structure_dict = get_binary_data_from_spatial_network(
            spatial_network;
            bond_radius = bond_radius,
            voxel_edge_length = voxel_edge_length,
            save_path = save_path*run_and_filename,
            filename = run_and_filename[7:end],
            save_result=true)
            
    end

    return
end


"""
For each network dict in a list of run and filenames, create binary data with a
given bond radius and voxel edge length
"""
function get_binary_data_from_spatial_network_multithreading(
    spatial_networks_path::String,
    save_path::String;
    print_progress::Bool = true,
    print_lock = Threads.ReentrantLock(),
    nr_runs::Int64 = 5,
    bond_radius::Float64 = 0.35,
    voxel_edge_length::Float64 = 0.1)

    # vector for run and filename
    run_and_filename_vec = []

    for i in 1:nr_runs

        spatial_networks_path_current_run = (
            spatial_networks_path*"run_"*string(i)*"/")

        # get all filenames of networks
        filenames = readdir(spatial_networks_path_current_run)
        filenames_spatial_networks = filter(
            filename -> endswith(filename, ".gml"), filenames)
        filenames_spatial_networks = [filename_spatial_network[1:end-4] 
            for filename_spatial_network in filenames_spatial_networks]

        # append to list of all network paths
        append!(run_and_filename_vec, 
            "run_"*string(i)*"/" .* filenames_spatial_networks)

    end

    # split filenames into chunks for multi-threading
    run_and_filename_chunks = Iterators.partition(run_and_filename_vec, 
        length(run_and_filename_vec) ÷ Threads.nthreads())

    # run all filename chunk in parallel in different threads
    map(run_and_filename_chunks) do run_and_filename_chunk

        Threads.@spawn get_binary_data_from_spatial_network_single_thread(
            run_and_filename_chunk,
            spatial_networks_path,
            save_path;
            print_progress = print_progress,
            print_lock = print_lock,
            bond_radius = bond_radius,
            voxel_edge_length = voxel_edge_length)
    end

    return
end


"""
Get vector of vectors containing the vector components at which  the
autocovariance function will be calculated
"""
function get_sampling_indices_vec_vec(size_data::Tuple)

    # determine the maximal sampling distances along the three axes
    max_sampling_distances = Int.( floor.( (size_data ) ./ 2 ))

    # get sampling distance vec vec
    # Along one axis (z direction is chosen here) only positive directions are
    # considered, because negative ones would yield redundant information
    sampling_indices_vec_vec = [
                collect(-max_sampling_distances[1]:max_sampling_distances[1]),
                collect(-max_sampling_distances[2]:max_sampling_distances[2]),
                collect(0:max_sampling_distances[3])]

    return sampling_indices_vec_vec
end


"""
Get array vectors out of three different vectors containing the x, y and z
components
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
    # since a real FFT is performed, the first dimension contains only positve
    # wavenumbers whereas second and third dimension contain positive and
    # negative wavenumbers. In order to bring them into a natural order, the
    # fftshift needs to be done
    wavenumber_vec_vec = (
        2*pi/autocovariance_fct_direction_dict["voxel_edge_length"]) .* [
            FFTW.fftshift( FFTW.fftfreq( size(
                autocovariance_fct_direction_dict[
                    "autocovariance_fct_array"])[1] ) ),
            FFTW.fftshift( FFTW.fftfreq( size(
                autocovariance_fct_direction_dict[
                    "autocovariance_fct_array"])[2] ) ),
            FFTW.fftshift( FFTW.fftfreq( size(
                autocovariance_fct_direction_dict[
                    "autocovariance_fct_array"])[3] ) ) ]

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
Determine autocovariance function as a function of sampling direction for a 3d
medium with periodic boundary conditions
"""
function get_autocovariance_fct_by_sampling_indices_array(
    structure_dict::Dict;
    save_result = false,
    save_path = raw"..\analysis_data\sample_name",
    print_progress = false,
    thread_nr::Int64 = 0,
    print_lock = Threads.ReentrantLock())

    # get array of sampling vectors
    sampling_indices_vec_vec = get_sampling_indices_vec_vec(
        structure_dict["size_data"])
    sampling_indices_array = get_vector_array(sampling_indices_vec_vec)

    autocovariance_fct_array_size = size(sampling_indices_array)[1:3]

    # create vector where for each sampling distance the autocovariance
    # function and its uncertainty will be stored
    autocovariance_fct_array = Array{Float64}(undef, 
        autocovariance_fct_array_size...)

    # for each sampling distance get autocovariance function and its
    # uncertainty
    for i in 1:autocovariance_fct_array_size[1]
        for j in 1:autocovariance_fct_array_size[2]
            for k in 1:autocovariance_fct_array_size[3]

                autocovariance_fct_array[i,j,k] = get_autocovariance_fct(
                    sampling_indices_array[i,j,k,:], structure_dict)
            end
        end

        # calculate and print progress
        if print_progress
            progress_percentage = i/autocovariance_fct_array_size[1]*100

            lock(print_lock) do
                Format.printfmtln("Autocovariance function calculation progress
                    thread nr {1:d}: {2:.1f} %", 
                    thread_nr, progress_percentage)
            end
        end

    end

    # extend the sampling distance vec along the third dimension, where due to
    # the mirror symmetry of the autocovariance fct only positive z values
    # where considered
    complete_sampling_indices_vec_vec = [sampling_indices_vec_vec[1], 
        sampling_indices_vec_vec[2],
        vcat( .- sampling_indices_vec_vec[3][end:-1:2],
            sampling_indices_vec_vec[3]) ]

    # create array of sampling vectors out of sampling distances
    complete_sampling_indices_array = get_vector_array(
        complete_sampling_indices_vec_vec)

    # point mirror autocovariance fct to sampling distances with negative z
    # component
    complete_autocovariance_fct_array = cat(dims=3, 
        autocovariance_fct_array[end:-1:1,end:-1:1,end:-1:2], 
        autocovariance_fct_array)

    
    # if the voxelized structure has an even number of voxels along any
    # direction one layer of the autocovariance function is redundant and needs
    # to be removed
    if iseven(structure_dict["size_data"][1])
        complete_sampling_indices_vec_vec[1] = (
            complete_sampling_indices_vec_vec[1][2:end])
        complete_sampling_indices_array = (
            complete_sampling_indices_array[2:end,:,:,:])
        complete_autocovariance_fct_array = (
            complete_autocovariance_fct_array[2:end,:,:])
    end
    if iseven(structure_dict["size_data"][2])
        complete_sampling_indices_vec_vec[2] = (
            complete_sampling_indices_vec_vec[2][2:end])
        complete_sampling_indices_array = (
            complete_sampling_indices_array[:,2:end,:,:])
        complete_autocovariance_fct_array = (
            complete_autocovariance_fct_array[:,2:end,:])
    end
    if iseven(structure_dict["size_data"][3])
        complete_sampling_indices_vec_vec[3] = (
            complete_sampling_indices_vec_vec[3][2:end])
        complete_sampling_indices_array = (
            complete_sampling_indices_array[:,:,2:end,:])
        complete_autocovariance_fct_array = (
            complete_autocovariance_fct_array[:,:,2:end])
    end
^
    autocovariance_fct_direction_dict = Dict(
        "sampling_index_array" => complete_sampling_indices_array,
        "sampling_index_vec_vec" => complete_sampling_indices_vec_vec,
        "sampling_distance_array" => (complete_sampling_indices_array 
            .* structure_dict["voxel_edge_length"]),
        "sampling_distance_vec_vec" => (complete_sampling_indices_vec_vec 
            .* structure_dict["voxel_edge_length"]),
        "autocovariance_fct_array" => complete_autocovariance_fct_array,
        "voxel_edge_length" => structure_dict["voxel_edge_length"] ,
        "label" => structure_dict["label"])

    if save_result
        GU.save_dict_to_h5(autocovariance_fct_direction_dict, 
            save_path*"_autocovariance_fct_direction.h5")
    end

    return autocovariance_fct_direction_dict
end


"""
For each structure dict in a list of run and filenames, calculate the 
autocovariance function as a function of direction
"""
function get_autocovariance_fct_direction_from_filenames_single_thread(
    run_and_filename_chunk,
    structure_dicts_path::String,
    save_path::String;
    print_progress::Bool = false,
    print_lock = Threads.ReentrantLock())

    # loop through files
    for run_and_filename in run_and_filename_chunk

        # check that file is structure dict
        if endswith(run_and_filename, "_structure.h5")

            # load structure dict
            structure_dict = GU.load_h5_dict(
                structure_dicts_path*run_and_filename)

            # check if autocovariance dict already exists
            if !isfile(save_path*run_and_filename[1:end-13]
                *"_autocovariance_fct_direction.h5")

                # print current thread id and run/filename if desired
                if print_progress
                    lock(print_lock) do
                        Format.printfmtln("Thread {1} is handling file {2}",
                            Threads.threadid(), run_and_filename)
                    end
                end

                autocovariance_fct_direction_dict = (
                    get_autocovariance_fct_by_sampling_indices_array(
                        structure_dict;
                        save_result = true,
                        save_path = save_path*run_and_filename[1:end-13],
                        print_progress = print_progress,
                        thread_nr = Threads.threadid() ,
                        print_lock = print_lock))
            end
        end
    end

    return
end


"""
For all structure dicts in a folder, calculate the autocovariance function as a
function of direction using multithreading
"""
function get_autocovariance_fct_direction_from_filenames_multithreading(
    structure_dicts_path;
    print_progress::Bool = false,
    save_path::String = "..\analysis_data\random_networks\\",
    nr_runs::Int64 = 5,
    print_lock = Threads.ReentrantLock())

    # vector for run and filename
    run_and_filename_vec = []

    # loop through runs
    for i in 1:nr_runs

        structure_dicts_path_current_run = (
            structure_dicts_path*"run_"*string(i)*"/")

        filenames = readdir(structure_dicts_path_current_run)
        filenames_structure_dicts = filter(filename 
            -> endswith(filename, "_structure.h5"), filenames)

        # append to list of all structure dict paths
        append!(run_and_filename_vec, 
            "run_"*string(i)*"/" .* filenames_structure_dicts)

    end

    # split filenames into chunks for multi-threading
    run_and_filename_chunks = Iterators.partition(run_and_filename_vec, 
        length(run_and_filename_vec) ÷ Threads.nthreads())

    # run all filename chunk in parallel in different threads
    map(run_and_filename_chunks) do run_and_filename_chunk

        Threads.@spawn (
            get_autocovariance_fct_direction_from_filenames_single_thread(
                run_and_filename_chunk,
                structure_dicts_path,
                save_path;
                print_progress = print_progress,
                print_lock = print_lock))
    end

end


"""
Calculate spectral density from autocovariance fct for a 3d medium with 
periodic boundary conditions by means of Fast Fourier Transform
"""
function get_spectral_density_by_wavevector_array_fft(
    structure_dict::Dict;
    save_autocovariance_fct_direction_dict::Bool = false,
    save_result::Bool = false,
    save_path::String = raw"..\analysis_data\sample_name",
    autocovariance_fct_direction_dict::Dict = (
        get_autocovariance_fct_by_sampling_indices_array(
            structure_dict;
            save_result = save_autocovariance_fct_direction_dict,
            save_path = save_path)))

    # determine fourier transform of autocovariance function values
    spectral_density_array_fft_output = FFTW.rfft(
        autocovariance_fct_direction_dict["autocovariance_fct_array"])

    # bring spectral densities into an order from negative to positive
    # wavenumbers 
    spectral_density_array = FFTW.fftshift(
        spectral_density_array_fft_output, [2, 3])

    # point mirror spectral density to wavenumbers with negative x component
    spectral_density_array = cat(dims=1, 
        conj.(spectral_density_array[end:-1:2,end:-1:1,end:-1:1]), 
        spectral_density_array)

    # get tuple of vectors of sampled wavenumbers
    wavenumber_vec_vec = get_wavenumber_vec_vec(
        autocovariance_fct_direction_dict)

    # get array of sampled sampled wavevectors 
    wavevector_array = get_vector_array(wavenumber_vec_vec)

    # create dict to save
    spectral_density_dict = Dict(
        "wavevector_array" => wavevector_array,
        "wavenumber_vec_vec" => wavenumber_vec_vec,
        "spectral_density_array" => spectral_density_array,
        "voxel_edge_length" 
            => autocovariance_fct_direction_dict["voxel_edge_length"],
        "label" => autocovariance_fct_direction_dict["label"])

    # save results if desired
    if save_result
        GU.save_dict_to_h5(copy(spectral_density_dict), 
            save_path*"_spectral_density_array.h5")

    end

    return spectral_density_dict
end


"""
Calculate the scaled intersection volume of 3d spheres of given radius R and
distance to the origin r
"""
function scaled_intersection_volume_3d(
    r_distance_to_origin::Float64, 
    sphere_radius::Float64)

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
From the autocovariance function, calculate the volume fraction variance from 
equation 73 in 10.1016/j.physrep.2018.03.001
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
        LinearAlgebra.norm(autocovariance_fct_direction_dict[
            "sampling_distance_array"][1,1,1,:]
        .- autocovariance_fct_direction_dict[
            "sampling_distance_array"][1,1,end,:]))

    # get vector of sphere radii
    sphere_radius_vec = collect(
        autocovariance_fct_direction_dict["voxel_edge_length"]/2
        :autocovariance_fct_direction_dict["voxel_edge_length"]/2
        :voxelized_data_edge_length/4)

    # initialize vector of volume fraction variances
    volume_fract_variance_vec = Vector{Float64}(undef, 
        length(sphere_radius_vec))

    # initialize vector of volume fraction variances times window volume
    volume_fract_variance_times_window_volume_vec = Vector{Float64}(undef,
        length(sphere_radius_vec))

    # for each sphere radius, calculate the volume fraction variance
    for i in eachindex(sphere_radius_vec)

        # get volume of sphere
        sphere_volume = 4/3 * pi * sphere_radius_vec[i]^3

        # calculate integral of autocovariance function multiplied
        # by the scaled intersection volume
        volume_fract_variance_times_window_volume = sum(
            autocovariance_fct_direction_dict["autocovariance_fct_array"] .* 
            scaled_intersection_volume_3d.(r_distance_to_origin_array, 
            sphere_radius_vec[i]))

        volume_fract_variance_times_window_volume_vec[i] = (
            volume_fract_variance_times_window_volume)

        # save volume fraction variance
        volume_fract_variance_vec[i] = (1/sphere_volume 
            * volume_fract_variance_times_window_volume)
    end

    volume_fract_variance_dict = Dict{String, Any}(
        "sphere_radius_vec" => sphere_radius_vec, 
        "volume_fract_variance_vec" => volume_fract_variance_vec,
        "volume_fract_variance_times_window_volume_vec" 
            => volume_fract_variance_times_window_volume_vec,
        "label" => autocovariance_fct_direction_dict["label"])

    if save_result
        GU.save_dict_to_h5(copy(volume_fract_variance_dict),
            save_path*"_volume_fraction_variance.h5")
    end

    return volume_fract_variance_dict
end


"""
Get a digital sphere for a given radius
"""
function get_digital_sphere(radius)

    # Initialize an empty array to hold the voxel coordinates
    voxels = []

    # Iterate through a cube that bounds the sphere
    for x in -Int(round(radius)):Int(round(radius))
        for y in -Int(round(radius)):Int(round(radius))
            for z in -Int(round(radius)):Int(round(radius))
                # Check if the current point lies within the sphere
                if x^2 + y^2 + z^2 <= radius^2
                    push!(voxels, [x, y, z])
                end
            end
        end
    end

    return voxels
end


"""
For a given digital sphere, get a mask array of the size of the voxelized data
with only the voxels of the digital sphere set to 1. The digital sphere is
centered around the index [0,0,0] of the mask array and the array has periodic
boundary conditions
"""
function get_digital_sphere_mask(radius::Float64, size_data::Tuple)

    # get digital sphere
    digital_sphere = get_digital_sphere(radius)

    # initialize mask array
    mask_array = zeros(Bool, size_data...)

    # loop through all voxels of the digital sphere
    for voxel in digital_sphere

        # set voxel to 1
        mask_array[ ((.- [1,1,1] .+ minimum(size_data) .+ voxel
                        ).% size_data .+ [1,1,1])...
                            ] = 1

    end

    return mask_array
end


"""
For a given data size store all digital sphere masks in a dictionary
"""
function get_digital_sphere_mask_dict(size_data::Tuple;
    save_result::Bool = false,
    save_path = raw"..\analysis_data\random_networks\digital_sphere_masks\\")

    # create list of digital spheres with increasing radius
    sphere_pixel_radius_vec = collect(0.5001:0.5:minimum(
        [minimum(size_data)/2, 50]))
    digital_sphere_mask_arr = Array{Bool}(undef, size_data..., 
        length(sphere_pixel_radius_vec))

    # loop through all radii
    for i in eachindex(sphere_pixel_radius_vec)

        digital_sphere_mask_arr[:,:,:,i] = get_digital_sphere_mask(
            sphere_pixel_radius_vec[i], size_data)

    end

    digital_sphere_mask_dict = Dict{String, Any}(
        "sphere_pixel_radius_vec" => sphere_pixel_radius_vec,
        "digital_sphere_mask_arr" => digital_sphere_mask_arr,
        "size_data" => size_data)


    if save_result
        GU.save_dict_to_h5(copy(digital_sphere_mask_dict),
            save_path*"digital_sphere_mask_size_" *string(size_data)* ".h5")

    end

    return digital_sphere_mask_dict
end


"""
Calculate the pore size distribution following the method described in
10.1103/PhysRevE.100.053314, modified to work with periodic boundary conditions
and by using random sampling of voxels to speed up the calculation
"""
function get_pore_size_distribution_voxelized(
    structure_dict::Dict;
    nr_sampled_voxels::Int = 20000,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name\\",
    label = nothing,
    digital_sphere_mask_path 
        = raw"..\analysis_data\random_networks\digital_sphere_masks\\",
    print_progress::Bool = false,
    thread_nr::Int64 = 0,
    print_lock = Threads.ReentrantLock())

    # if digital sphere mask file for given data size exists, load it
    if isfile(digital_sphere_mask_path*"digital_sphere_mask_size_" 
        *string(structure_dict["size_data"])* ".h5")

        digital_sphere_mask_dict = GU.load_h5_dict(digital_sphere_mask_path
            *"digital_sphere_mask_size_" *string(structure_dict["size_data"])
            * ".h5")

        # print progress
        if print_progress
            lock(print_lock) do
                Format.printfmtln(
                    "Thread nr {1:d}: digital sphere mask dict loaded", 
                        thread_nr)
            end
            
        end

    # otherwise create and save digital sphere mask dict
    else

        # print progress
        if print_progress
            lock(print_lock) do
                Format.printfmtln(
                    "Thread nr {1:d}: generating digital sphere mask dict", 
                        thread_nr)
            end
            
        end

        digital_sphere_mask_dict = get_digital_sphere_mask_dict(
            structure_dict["size_data"];
            save_result = true,
            save_path = digital_sphere_mask_path)

    end

    # create array with pore radii
    pore_pixel_radius_array = zeros(size(structure_dict["data_binary"])...)

    # create array of ones of the size of the voxelized data
    ones_array = ones(Bool, size(structure_dict["data_binary"])...)

    # initialize counter if progress is printed
    if print_progress
        i = 1
    end

    # sample given number of voxels to speed up calculation
    sampled_coords = StatsBase.sample(CartesianIndices(
        structure_dict["data_binary"]), nr_sampled_voxels, replace=false)

    # loop through all sampled voxels using cartesian indices
    for coord in sampled_coords

        # check if voxel is in pore
        if !structure_dict["data_binary"][coord]

            # loop through all digital sphere masks
            for j in eachindex(
                digital_sphere_mask_dict["sphere_pixel_radius_vec"])

                digital_sphere_mask = digital_sphere_mask_dict[
                    "digital_sphere_mask_arr"][:,:,:,j]

                # use periodic boundary conditions to shift all entries of the
                # mask array along the three dimensions as given by the
                # sampling vector
                digital_sphere_mask_shifted = circshift(
                    digital_sphere_mask, coord.I)

                # check if all voxels of the digital sphere are in the pore
                if any(structure_dict["data_binary"][
                    digital_sphere_mask_shifted])
                    break

                else
                    # set pore pixel radius to maximum of current and previous
                    # radius
                    pore_pixel_radius_array[digital_sphere_mask_shifted] = (
                        max.(pore_pixel_radius_array[
                            digital_sphere_mask_shifted],
                        digital_sphere_mask_dict["sphere_pixel_radius_vec"][j] 
                            .*  ones_array[digital_sphere_mask_shifted]))

                end
            end
        end

        # print progress
        if print_progress
            i += 1

            # print every 100th voxel
            if i % 100 == 0
                progress_percentage = i/nr_sampled_voxels*100

                lock(print_lock) do
                    Format.printfmtln("Pore size distribution calculation
                        progress thread nr {1:d}: {2:.1f} %", 
                        thread_nr, progress_percentage)
                end
            end
        end
    end

    # shape the pore pixel radius array into a vector
    pore_pixel_radius_vec = vec(pore_pixel_radius_array)

    # filter out voxels that are not in a pore
    pore_pixel_radius_filtered_vec = pore_pixel_radius_vec[
        pore_pixel_radius_vec .> 0.0]

    # create histogram of pore pixel radii
    pixel_radius_histogram = StatsBase.fit(
        StatsBase.Histogram, pore_pixel_radius_filtered_vec, 
        digital_sphere_mask_dict["sphere_pixel_radius_vec"][1]-0.25:0.5
            :digital_sphere_mask_dict["sphere_pixel_radius_vec"][end]+0.251,
        closed=:left)

    # normalize histogram
    pixel_radius_histogram = LinearAlgebra.normalize(
        pixel_radius_histogram, mode=:probability)

    pore_size_distribution = pixel_radius_histogram.weights

    # convert pixel radii to physical radii
    pore_size_vec = (digital_sphere_mask_dict["sphere_pixel_radius_vec"] 
        .* structure_dict["voxel_edge_length"])

    pore_size_distribution_dict = Dict{String, Any}(
        "pore_size_vec" => pore_size_vec,
        "pore_size_distribution" => pore_size_distribution,
        "nr_sampled_voxels" => nr_sampled_voxels)

    # add label to dictionary if label is not nothing
    if label !== nothing
        pore_size_distribution_dict["label"] = label
    end

    if save_result
        GU.save_dict_to_h5(copy(pore_size_distribution_dict),
            save_path*"_pore_size_distribution.h5")
    end

    return pore_size_distribution_dict
end


"""
Measure the standard deviation of dihedral angles. This function might need to
be revisited when considering other coordination numbers than 4 because the
dihedral angle might have a peak around 0 which is not considered here
"""
function get_dihedral_angle_std(spatial_network::MetaGraphsNext.MetaGraph)

    # initialize vector of diehedral angles
    dihedral_angle_vec = Vector{Float64}()

    # loop through all bonds
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        # TODO: store bond[1] and bond[2] here, instead of calculating it always

        # get vector along bond
        bond_vec = spatial_network[bond...]["vector"]

        # loop through all neighbors of one vertex 
        for first_neighbor in setdiff( MetaGraphsNext.neighbor_labels(
            spatial_network, bond[1]), bond[2])

            # get vector from first neighbor to first bond vertex
            first_neighbor_to_bond_vertex_vec = ( sign(bond[1] 
                - first_neighbor)*spatial_network[first_neighbor, 
                    bond[1]]["vector"])

            # loop through all neighbors of other vertex
            for second_neighbor in setdiff( MetaGraphsNext.neighbor_labels(
                spatial_network, bond[2]), bond[1])

                # get vector from second bond vertex to second neighbor
                bond_vertex_to_second_neighbor_vec = ( sign(second_neighbor 
                    - bond[2])*spatial_network[bond[2], 
                        second_neighbor]["vector"])

                # calculate dihedral angle according to the equation given in
                # https://en.wikipedia.org/wiki/Dihedral_angle# 
                dihedral_angle = atan( LinearAlgebra.norm(bond_vec)
                        *LinearAlgebra.dot(first_neighbor_to_bond_vertex_vec,
                            LinearAlgebra.cross(bond_vec,
                                bond_vertex_to_second_neighbor_vec )),
                    LinearAlgebra.dot(
                        LinearAlgebra.cross(first_neighbor_to_bond_vertex_vec,
                            bond_vec ),
                        LinearAlgebra.cross(bond_vec,
                              bond_vertex_to_second_neighbor_vec )))

                # save dihedral angle
                push!(dihedral_angle_vec, dihedral_angle)
            end
        end
    end

    # save one peak of dihedral angle distribution
    lower_limit = 0
    #upper_limit =  2 * pi / (spatial_network[]["coordination_nr"] - 1)    

    max_coordination_nr=0
    for vertex in MetaGraphsNext.labels(spatial_network)

        # get iterator of bond combinations
        bond_combinations_iter = Combinatorics.combinations(collect(
            MetaGraphsNext.neighbor_labels(spatial_network, vertex)), 2)

        vertex_coordination_nr=length(bond_combinations_iter)

        if vertex_coordination_nr>max_coordination_nr

            max_coordination_nr=vertex_coordination_nr
            
        end
    end
    upper_limit =  2 * pi / (max_coordination_nr - 1)
    

    dihedral_angle_one_peak_vec = dihedral_angle_vec[
        (dihedral_angle_vec .> lower_limit) .& (dihedral_angle_vec .< upper_limit)]

    # determine standard deviation
    dihedral_angle_std = Statistics.std(dihedral_angle_one_peak_vec)

    return [dihedral_angle_std, dihedral_angle_vec]
end


"""
For a single voxelized structure, calculate several functions that characterize
the structure
"""
function get_all_dicts_from_voxelized_structure(
    filename::String,
    spatial_network_path::String,
    analysis_data_path::String;
    bond_radius = 0.05,
    voxel_edge_length = 0.1,
    print_progress::Bool = false,
    print_lock = Threads.ReentrantLock())

    # get network
    spatial_network = NG.load_spatial_network_from_gml(
        spatial_network_path*filename*".gml")

    # create structure dict
    structure_dict = get_binary_data_from_spatial_network(
        spatial_network;
        bond_radius = bond_radius,
        voxel_edge_length = voxel_edge_length,
        save_path = structure_dict_path*filename,
        filename = filename,
        save_result=true)

    # get autocovariance function as a function of direction
    autocovariance_fct_direction_dict = (
        get_autocovariance_fct_by_sampling_indices_array(
            structure_dict;
            save_result = true,
            save_path = analysis_data_path*filename,
            print_progress = print_progress,
            thread_nr = Threads.threadid() ,
            print_lock = print_lock))
    
    # get spectral density by wavevector array
    spectral_density_dict = get_spectral_density_by_wavevector_array_fft(
        structure_dict;
        save_autocovariance_fct_direction_dict = false,
        save_result = true,
        save_path = analysis_data_path*filename,
        autocovariance_fct_direction_dict = autocovariance_fct_direction_dict)

    # get angle averaged spectral density
    spectral_density_angle_averaged_dict = get_spectral_density_angle_averaged(
        spectral_density_dict;
        gaussian_filter = true,
        gaussian_filter_sigma_x = 2*pi/25, 
        gaussian_filter_filtered_data_x_step_length = 2*pi/25,
        save_result = true,
        save_path = analysis_data_path*filename)

    # get volume fraction variance
    volume_fract_variance_dict = get_volume_fract_variance(
        autocovariance_fct_direction_dict;
        save_result = true,
        save_path = analysis_data_path*filename)

    # get pore size distribution
    pore_size_distribution_dict = get_pore_size_distribution_voxelized(
        structure_dict;
        nr_sampled_voxels = 20000,
        save_result = true,
        save_path = analysis_data_path*filename,
        label = filename,
        digital_sphere_mask_path = analysis_data_path)

    return
end

