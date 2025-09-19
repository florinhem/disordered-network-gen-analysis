"""
these functions can be used to characterize networks by means of local order
metrics
"""

"""
Measure the standard deviation of bond lengths
"""
function get_bond_length_std(spatial_network::MetaGraphsNext.MetaGraph)
    
    nr_bonds=0
    for bond in MetaGraphsNext.edge_labels(spatial_network)
        nr_bonds+=1
    end

    bond_length_vec = Vector{Float64}(undef, nr_bonds)
    bond_count = 1

    for bond in MetaGraphsNext.edge_labels(spatial_network)

        bond_length_vec[bond_count] = sqrt(
            spatial_network[bond...]["distance_squared"]
        )

        bond_count += 1
    end

    bond_length_std = Statistics.std(bond_length_vec)
    
    return [bond_length_std, bond_length_vec]
end


"""
Measure the standard deviation of bond angles
"""
function get_bond_angle_std(spatial_network::MetaGraphsNext.MetaGraph)
    
    nr_angles=0
    for vertex in MetaGraphsNext.labels(spatial_network)

        # get iterator of bond combinations
        bond_combinations_iter = Combinatorics.combinations(collect(
            MetaGraphsNext.neighbor_labels(spatial_network, vertex)), 2)

        nr_angles+=length(bond_combinations_iter)
        
    end
    
    # initialize vector of bond angles
    bond_angle_vec = Vector{Float64}(undef, nr_angles)
    angle_count = 1

    # loop through all vertices
    for vertex in MetaGraphsNext.labels(spatial_network)

        # get iterator of bond combinations
        bond_combinations_iter = Combinatorics.combinations(collect(
            MetaGraphsNext.neighbor_labels(spatial_network, vertex)), 2)

        # loop through bond combinations
        for bond_combination in bond_combinations_iter

            # get scalar product of vectors representing the current bond
            # combination

            vector_1=spatial_network[vertex, bond_combination[1]]["vector"]
            vector_2=spatial_network[vertex, bond_combination[2]]["vector"]

            scalar_product = (sign(bond_combination[1] - vertex)
                * sign(bond_combination[2] - vertex)
                * LinearAlgebra.dot(vector_1, vector_2))

            cosinus_of_bond_angle=scalar_product/ sqrt(
                LinearAlgebra.dot(vector_1,vector_1)*
                LinearAlgebra.dot(vector_2,vector_2))
            
            if cosinus_of_bond_angle ≈ 1
                cosinus_of_bond_angle=1
            elseif cosinus_of_bond_angle ≈ -1
                cosinus_of_bond_angle=-1
            end

            # calculate bond angle
            bond_angle = acos(cosinus_of_bond_angle)

            bond_angle_vec[angle_count] = bond_angle

            angle_count += 1
        end
    end

    # determine standard deviation
    bond_angle_std = Statistics.std(bond_angle_vec)

    return [bond_angle_std, bond_angle_vec]
end

"""
Measure the Shannon entropy of dihedral angles that are binned in bins of 10
degrees
"""
function get_dihedral_angle_entropy(spatial_network::MetaGraphsNext.MetaGraph)

    # initialize vector of diehedral angles
    dihedral_angle_vec = Vector{Float64}()

    # loop through all bonds
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        # store the two vertices
        vertex1=bond[1]
        vertex2=bond[2]

        # get vector along bond
        bond_vec = spatial_network[bond...]["vector"]

        # loop through all neighbors of one vertex 
        for first_neighbor in setdiff( MetaGraphsNext.neighbor_labels(
            spatial_network, vertex1), vertex2)

            # get vector from first neighbor to first bond vertex
            first_neighbor_to_bond_vertex_vec = ( sign(vertex1 
                - first_neighbor)*spatial_network[first_neighbor, 
                vertex1]["vector"])

            # loop through all neighbors of other vertex
            for second_neighbor in setdiff( MetaGraphsNext.neighbor_labels(
                spatial_network, vertex2), vertex1)

                # get vector from second bond vertex to second neighbor
                bond_vertex_to_second_neighbor_vec = ( sign(second_neighbor 
                    - vertex2)*spatial_network[vertex2, 
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

    # consider only dihedral angles in the range [0, pi]
    dihedral_angle_vec = filter( x -> x >= 0 && x <= pi, dihedral_angle_vec)

    # determine histogram of dihedral angles
    histogram = StatsBase.fit(StatsBase.Histogram, dihedral_angle_vec, 
        nbins=18)
    #dihedral_angle_bin_edges = histogram.edges
    dihedral_angle_bin_weights = histogram.weights

    # normalize histogram weights to get probability distribution
    dihedral_angle_bin_weights /= sum(dihedral_angle_bin_weights)
    
    # calculate normalized Shannon entropy
    dihedral_angle_entropy = -1/log(length(dihedral_angle_bin_weights)) *sum(
        dihedral_angle_bin_weights .* log.(dihedral_angle_bin_weights .+ 1e-12))
    
    return dihedral_angle_entropy
end


"""
Get the 12 vertices of an icosahedron. The vertices are normalized to unit
length.
"""
function icosahedron_vertices()
    # golden ratio
    φ = (1 + sqrt(5)) / 2

    vertex_vectors = [
        [-1,  φ,  0], [ 1,  φ,  0], [-1, -φ,  0], [ 1, -φ,  0],
        [ 0, -1,  φ], [ 0,  1,  φ], [ 0, -1, -φ], [ 0,  1, -φ],
        [ φ,  0, -1], [ φ,  0,  1], [-φ,  0, -1], [-φ,  0,  1]]

    # rotate all vertex vectors about the x-axis, such that the the vectors 
    # [ 0,  1,  φ] and [ 0, -1, -φ] lie on the z-axis
    rotation_angle = atan(1, φ)
    rotation_matrix = [1 0 0; 0 cos(rotation_angle) -sin(rotation_angle); 
        0 sin(rotation_angle) cos(rotation_angle)]
    for i in 1:length(vertex_vectors)
        vertex_vectors[i] = rotation_matrix * vertex_vectors[i]
    end

    normalized_vertex_vectors = [LinearAlgebra.normalize(v) 
        for v in vertex_vectors]

    return normalized_vertex_vectors
end


"""
Get the faces of the icosahedron. Each face is represented by a tuple of vertex
indices. The vertices are indexed from 1 to 12 in the order they are defined in
the icosahedron_vertices function.
"""
function icosahedron_faces()
    
    return [(1,12,6), (1,6,2), (1,2,8), (1,8,11), (1,11,12), (2,6,10), 
        (6,12,5), (12,11,3), (11,8,7), (8,2,9), (4,10,5), (4,5,3), (4,3,7), 
        (4,7,9), (4,9,10), (5,10,6), (3,5,12), (7,3,11), (9,7,8), (10,9,2)]
end


"""
Subdivide one triangle face into 4 smaller triangles
"""
function subdivide_face(v1, v2, v3)
    # Get vertices in between the corners of the triangle
    m12 = LinearAlgebra.normalize((v1 + v2) / 2)
    m23 = LinearAlgebra.normalize((v2 + v3) / 2)
    m31 = LinearAlgebra.normalize((v3 + v1) / 2)

    return [(v1, m12, m31), (v2, m23, m12), (v3, m31, m23), (m12, m23, m31)]
end

"""
Get the centers of all subdivided icosahedron faces
"""
function icosahedron_subdivided_face_centers()
    
    vertices = icosahedron_vertices()
    faces = icosahedron_faces()

    subdivided_face_centers = Vector{Vector{Float64}}()
    for (i1, i2, i3) in faces
        v1, v2, v3 = vertices[i1], vertices[i2], vertices[i3]
        for (a, b, c) in subdivide_face(v1, v2, v3)
            face_center = LinearAlgebra.normalize((a + b + c) / 3)
            push!(subdivided_face_centers, face_center)
        end
    end
    return subdivided_face_centers
end

# Assign each vector to the closest bin (by cosine similarity)
function bin_vectors(vectors, bins)
    counts = zeros(Int, length(bins))
    for v in vectors
        vn = normalize(v)
        dot_products = [dot(vn, b) for b in bins]
        idx = argmax(dot_products)
        counts[idx] += 1
    end
    return counts
end


"""
Get an entropy-based anisotropy metric from the orientations of the bonds of a
network. To this end, the bond vectors are transformed to spherical coordinates 
and the angles are binned in bins that correspond to the faces of an 
icosahedron. The entropy is then calculated from the histogram of the angles.
"""
function get_bond_orientation_entropy(
    spatial_network::MetaGraphsNext.MetaGraph)

    # get the centers of the subdivided icosahedron faces to use as bins for
    # the bonds of the network
    subdivided_face_centers = icosahedron_subdivided_face_centers()

    # consider only bins with a positive z components
    subdivided_face_centers = filter(v -> v[3] > 0, subdivided_face_centers)
    
    # check if there are 40 bins. If not, print an error message
    if length(subdivided_face_centers) != 40
        @error "The number of bins is not 40. 
            Please check the implementation of the function."
    end

    # get the bond vectors of the network
    bond_vectors = Vector{Vector{Float64}}()
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        bond_vector = spatial_network[bond...]["vector"]

        # if the z component is negative, flip the vector
        if bond_vector[3] < 0
            bond_vector = -bond_vector
        end

        # save the normalized bond vector
        push!(bond_vectors, LinearAlgebra.normalize(bond_vector)) 
    end

    # bin the bond vectors to the subdivided faces of the icosahedron
    bin_counts = zeros(Int, length(subdivided_face_centers))
    for v in bond_vectors
        dot_products = [LinearAlgebra.dot(v, b) 
            for b in subdivided_face_centers]
        closest_bin_index = argmax(dot_products)
        bin_counts[closest_bin_index] += 1
    end

    # normalize the bin counts to get a probability distribution
    bin_counts = bin_counts / sum(bin_counts)

    # calculate the normalized Shannon entropy of the bin counts
    bond_orientation_entropy = (- 1/log(length(bin_counts))
        * sum(bin_counts .* log.(bin_counts .+ 1e-12)))

    return bond_orientation_entropy
end


"""
Get the mean and standard deviation of the coordination number of a network.
"""
function get_coordination_nr_statistics(
    spatial_network::MetaGraphsNext.MetaGraph)

    # check if the network has the key "coordination_nr_vec"
    if haskey(spatial_network[], "coordination_nr_vec")
        coordination_nr_vec = spatial_network[]["coordination_nr_vec"]

    # if not, create a vector of the coordination number for each vertex
    else
        coordination_nr_vec = Vector{Int64}(
            undef, spatial_network[]["nr_vertices"])
        for vertex in MetaGraphsNext.labels(spatial_network)
            coordination_nr_vec[vertex] = length(
                MetaGraphsNext.neighbor_labels(spatial_network, vertex))
        end
    end
    
    
    # calculate mean and standard deviation of the coordination number
    coordination_nr_mean = Statistics.mean(coordination_nr_vec)
    coordination_nr_std = Statistics.std(coordination_nr_vec)

    return [coordination_nr_mean, coordination_nr_std]
end


"""
Get dictionary of q_lm averaged over bonds to neighbors for a single vertex
(also called Steinhardt local bond order parameters)  as in equation 2.1 of 
10.1103/PhysRevB.28.784 which is not rotationally invariant
"""
function get_q_lm_averaged_bonds_to_neighbors_dict_single_vertex(
    spatial_network::MetaGraphsNext.MetaGraph, 
    central_vertex::Int64,
    l_max::Int64)

    # initialize vector for spherical harmonics of all neighbors
    y_spherical_harmonic_arr_vec = Vector{SphericalHarmonics.SHVector{
            ComplexF64, 
            Vector{ComplexF64}, 
            Tuple{SphericalHarmonics.ML{SphericalHarmonics.ZeroTo{false}, 
                SphericalHarmonics.FullRange{true}}}
        }}(undef, spatial_network[central_vertex]["coordination_nr"])

    # neighbor counter
    neighbor_count = 1

    # loop through bonds to neighbors
    for neighbor in MetaGraphsNext.neighbor_labels(
                        spatial_network, central_vertex)

        # get vector from vertex to neighbor
        vertex_to_neighbor_vec = (sign(neighbor - central_vertex) 
                * spatial_network[central_vertex, neighbor]["vector"] )

        # get vector's spherical coordinates
        r_length, theta, phi = convert_cartesian_to_spherical(
            vertex_to_neighbor_vec)

        # get array of spherical harmonics
        y_spherical_harmonic_arr_vec[neighbor_count] = (
            SphericalHarmonics.computeYlm(theta, phi; lmax = l_max))

        neighbor_count += 1
    end

    # average over bonds to neighbors

    # initialize dict of Steinhardt order parameters
    q_lm_averaged_bonds_to_neighbors_dict = Dict{Tuple{Int64, Int64}, 
        Complex{Float64}}()

    # loop through values of l
    for l in 0:l_max

        # loop through values of m
        for m in -l:l

            # initialize current average steinhardt order parameter
            q_lm = 0.0

            # average over neighbors
            for neighbor_count in 1:spatial_network[central_vertex]["coordination_nr"]

                q_lm += (1/spatial_network[central_vertex]["coordination_nr"]
                    * 
                    y_spherical_harmonic_arr_vec[neighbor_count][(l,m)]
                )

            end

            # save current average steinhardt order parameter
            q_lm_averaged_bonds_to_neighbors_dict[(l,m)] = q_lm
        end

    end

    return q_lm_averaged_bonds_to_neighbors_dict
end


"""
Get vector of mean values of q_l (rotationally invariant Steinhardt local bond 
order parameters) for a single vertex and for all parameters l up to l_max 
where l is the index of the spherical harmonic Y_{lm}. The equations are taken 
from references 10.1103/PhysRevB.28.784 and 10.1063/1.2977970
"""
function get_q_l_averaged_single_vertex_dict(
    spatial_network::MetaGraphsNext.MetaGraph,
    central_vertex::Int64,
    l_max::Int64)

    # get average steinhardt order parameter dict for current vertex this 
    # quantity is nor yet rotationally invariant
    q_lm_averaged_bonds_to_neighbors_dict = (
        get_q_lm_averaged_bonds_to_neighbors_dict_single_vertex(
            spatial_network, central_vertex, l_max))

    # initialize dict of Steinhardt order parameters 
    q_l_averaged_single_vertex_dict = Dict{Int64, Float64}()

    # loop through values of l
    for l in 0:l_max

        # initialize sum over m
        q_lm_squared_sum = 0.0

        # loop through values of m
        for m in -l:l

            # add absolute value squared of q_lm to sum
            q_lm_squared_sum += abs2(
                q_lm_averaged_bonds_to_neighbors_dict[(l,m)])
        end

        # calculate steinhardt order parameter for current l
        q_l_averaged_single_vertex_dict[l] = sqrt( 
            4*pi / (2*l + 1) * q_lm_squared_sum )

    end

    return q_l_averaged_single_vertex_dict
end


"""
Get vector of mean values of q_l (rotationally invariant Steinhardt local bond
order parameters) for the entire network and for all parameters l up to l_max 
where l is the index of the spherical harmonic Y_{lm}.
"""
function get_q_l_total_network_mean_dict(
    spatial_network::MetaGraphsNext.MetaGraph,
    l_max::Int64)

    # initialize dictionary of q_l averaged over entire network with all values
    # set to 0
    q_l_total_network_mean_arr = Array{Float64}(undef, 
        spatial_network[]["nr_vertices"], l_max+1)


    # loop through vertices
    for vertex in MetaGraphsNext.labels(spatial_network)

        # get vector of steinhardt order parameters for current vertex
        q_l_averaged_single_vertex_dict = (
            get_q_l_averaged_single_vertex_dict(spatial_network, vertex, 
            l_max))

        # for each l, add current vertex' contribution to sum of all vertices
        for l in 0:l_max
            q_l_total_network_mean_arr[vertex, l+1] = (
                q_l_averaged_single_vertex_dict[l])
    
        end
    end

    # calculate mean and standard deviation of q_l for the entire network
    q_l_total_network_mean_dict = Dict{Int64, 
        Measurements.Measurement{Float64}}()

    # calculate mean and standard deviation for each l individually
    for l in 0:l_max
        q_l_total_network_mean_dict[l] = Measurements.measurement(
            Statistics.mean(q_l_total_network_mean_arr[:, l+1]),
            Statistics.std(q_l_total_network_mean_arr[:, l+1])
        )

    end

    return q_l_total_network_mean_dict
end


"""
Get the minimal distance between a point and a bond of a spatial network by
taking into account periodic boundary conditions.
"""
function get_minimal_distance_to_bond(
    spatial_network::MetaGraphsNext.MetaGraph, 
    point::Vector{Float64}, bond::Tuple{Int64, Int64})

    vertex_1 = spatial_network[bond[1]]["position"]
    vertex_2 = spatial_network[bond[2]]["position"]

    bond_vector_normalized = LinearAlgebra.normalize(
        spatial_network[bond...]["vector"])

    bond_length = spatial_network[bond...]["distance_squared"]

    # get vector from vertex 1 to point taking into account periodic boundary
    # conditions
    vertex_1_to_point_vector = NG.get_distance_vector_pbc(
        vertex_1, point, spatial_network[]["supercell_edge_length"])

    # get the projection of the vector from vertex 1 to point onto the bond
    # vector
    t = LinearAlgebra.dot(vertex_1_to_point_vector, bond_vector_normalized)

    # limit the projection to the length of the bond
    t = clamp(t, 0, bond_length)

    # get the point on the bond closest to the point taking into account 
    # periodic boundary conditions
    closest_point_on_bond = mod.(vertex_1 .+ (t .* bond_vector_normalized), 
        spatial_network[]["supercell_edge_length"])

    # get vector from point to closest point on bond taking into account
    # periodic boundary conditions
    point_to_closest_point_on_bond_vector = NG.get_distance_vector_pbc(
        point, closest_point_on_bond, spatial_network[]["supercell_edge_length"])

    # get the distance between point and closest point on bond
    distance = LinearAlgebra.norm(point_to_closest_point_on_bond_vector)

    return distance
end


"""
Calculate the pore size distribution of the mathematical network with
infinitely thin bonds following the method described in
10.1103/PhysRevE.100.053314, to work with periodic boundary conditions
"""
function get_pore_size_distribution(
    spatial_network::MetaGraphsNext.MetaGraph;
    sampling_grid_size = 0.2,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing,
    digital_sphere_mask_path 
        = raw"..\analysis_data\random_networks\digital_sphere_masks\\",
    print_progress::Bool = false,
    thread_nr::Int64 = 0,
    print_lock = Threads.ReentrantLock())

    # determine the actual sampling grid size such that the grid is uniform
    # across the periodic boundary conditions
    sampling_grid_size = (spatial_network[]["supercell_edge_length"]
        / round(spatial_network[]["supercell_edge_length"]
        / sampling_grid_size))

    # get the vector of grid points along one direction
    grid_points = collect(
        0:sampling_grid_size:spatial_network[]["supercell_edge_length"])
    nr_grid_points_per_direction = length(grid_points)
    
    # create array of radii of the spheres centered at the grid points
    sphere_radii = Array{Float64}(undef, nr_grid_points_per_direction, 
        nr_grid_points_per_direction, nr_grid_points_per_direction)

    if print_progress
        coord_count = 0
        nr_coords = nr_grid_points_per_direction^3
        lock(print_lock) do
            Format.printfmtln("Pore size distribution calculation started 
                thread nr {1:d}. Nr grid points: {2:d}", thread_nr, nr_coords)
        end
    end


    # loop through sampling grid points
    for coord in CartesianIndices(sphere_radii)

        grid_point = [grid_points[coord[1]], grid_points[coord[2]], 
            grid_points[coord[3]]]

        # start from the maximal possible pore radius
        sphere_radius = spatial_network[]["supercell_edge_length"]/2

        # loop through all bonds of the network
        for bond in MetaGraphsNext.edge_labels(spatial_network)

            # get minimal distance between grid point and bond
            min_distance = get_minimal_distance_to_bond(
                spatial_network, grid_point, bond)

            # update the pore radius if the distance is smaller than the
            # current pore radius
            if min_distance < sphere_radius
                sphere_radius = min_distance
            end
        end

        # save the pore radius for the current grid point
        sphere_radii[coord] = sphere_radius

        # print progress
        if print_progress
            coord_count += 1

            # print every 100th voxel
            if coord_count % 1000 == 0
                progress_percentage = coord_count/nr_coords*100

                lock(print_lock) do
                    Format.printfmtln("Pore size distribution calculation
                        progress thread nr {1:d}: {2:.1f} %", 
                        thread_nr, progress_percentage)
                end
            end
        end
    end

    # get digital sphere mask to determine the largest pore radius for each
    # grid point
    grid_array_size = (nr_grid_points_per_direction, 
        nr_grid_points_per_direction, nr_grid_points_per_direction)

    # if digital sphere mask file for given data size exists, load it
    if isfile(digital_sphere_mask_path*"digital_sphere_mask_size_" 
        *string(grid_array_size)* ".h5")

        digital_sphere_mask_dict = GU.load_h5_dict(digital_sphere_mask_path
            *"digital_sphere_mask_size_" *string(grid_array_size)
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
            grid_array_size;
            save_result = true,
            save_path = digital_sphere_mask_path)
    end

    # get the vector of radii of the digital sphere mask in units of pixels or 
    # grid indices
    digital_sphere_mask_pixel_radius_vec = (
        digital_sphere_mask_dict["sphere_pixel_radius_vec"])
    # convert this vector to units of the equilibrium bond length
    pore_size_vec = (digital_sphere_mask_pixel_radius_vec 
        * sampling_grid_size)

    # create array of pore radii, which, for each grid point, contains the
    # radius of the largest pore, that this point can be part of
    pore_radii = zeros(nr_grid_points_per_direction, 
        nr_grid_points_per_direction, nr_grid_points_per_direction)

    # create array of ones of the size of the grid
    ones_array = ones(Bool, size(pore_radii)...)

    for coord in CartesianIndices(pore_radii)

        # only consider grid points that are not too close to any bond
        if sphere_radii[coord] > 1/4 * sampling_grid_size

            # get the sphere mask whose radius is closest to the current grid 
            # point's sphere radius
            digital_sphere_mask = (digital_sphere_mask_dict[
                "digital_sphere_mask_arr"][:,:,:,argmin(
                abs.(pore_size_vec .- sphere_radii[coord]))])

            # use periodic boundary conditions to shift all entries of the
            # mask array to center the mask at the current grid point
            digital_sphere_mask_shifted = circshift(
                digital_sphere_mask, coord.I)

            # set pore radius of surrounding grid points to the maximum of 
            # their current pore radius and the radius of the digital sphere 
            # mask
            pore_radii[digital_sphere_mask_shifted] = (
                max.(pore_radii[digital_sphere_mask_shifted],
                    sphere_radii[coord] 
                    .*  ones_array[digital_sphere_mask_shifted]))
        end
    end

    # shape the pore pixel radius array into a vector
    pore_radius_vec = vec(pore_radii)

    # filter out voxels that are not in a pore
    pore_radius_filtered_vec = pore_radius_vec[
        pore_radius_vec .> 0.0]

    # create histogram of pore radii
    radius_histogram = StatsBase.fit(
        StatsBase.Histogram, pore_radius_filtered_vec, 
        pore_size_vec[1]-0.25*sampling_grid_size:sampling_grid_size/2
            :pore_size_vec[end]+0.251*sampling_grid_size,
        closed=:left)

    # normalize histogram
    radius_histogram = LinearAlgebra.normalize(
        radius_histogram, mode=:probability)

    pore_size_distribution = radius_histogram.weights

    pore_size_distribution_dict = Dict{String, Any}(
        "pore_size_vec" => pore_size_vec,
        "pore_size_distribution" => pore_size_distribution,
        "sampling_grid_size" => sampling_grid_size)

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
Get the second moment of the pore size distribution which, according to 
10.1103/PhysRevE.104.014127 is an estimate of the critical pore radius squared
(although they use another definition of the pore size distribution)
"""
function get_pore_size_distribution_second_moment(
    pore_size_distribution_dict::Dict)

    pore_size_step_length = (pore_size_distribution_dict["pore_size_vec"][2] 
        - pore_size_distribution_dict["pore_size_vec"][1])

    # get second moment of pore size distribution
    pore_size_distribution_second_moment = pore_size_step_length * sum(
        pore_size_distribution_dict["pore_size_vec"].^2 
        .* pore_size_distribution_dict["pore_size_distribution"])

    return pore_size_distribution_second_moment
end


"""
Get an estimate of the critical pore radius from the second moment of the
pore size distribution according to 10.1103/PhysRevE.104.014127
"""
function get_critical_pore_radius(pore_size_distribution_dict::Dict)

    # get second moment of pore size distribution
    pore_size_distribution_second_moment = (
        get_pore_size_distribution_second_moment(
            pore_size_distribution_dict))

    # estimate critical pore radius as square root of second moment
    critical_pore_radius = sqrt(pore_size_distribution_second_moment)

    return critical_pore_radius
end


"""
Check if a cycle contains a given bond. The bond is represented as a tuple of
two vertex indices. The cycle is represented as a vector of vertex indices.
"""
function cycle_contains_bond(cycle::Vector{Int64}, bond::Tuple{Int64, Int64})
    
    for i in 1:length(cycle)-1
        if (cycle[i], cycle[i+1]) == bond || (cycle[i+1], cycle[i]) == bond
            return true
        end
    end
    return false
end

"""
Check if each bond is contained in at least one cycle. The bond is represented
as a tuple of two vertex indices. The cycle is represented as a vector of
vertex indices.
"""
function each_bond_in_at_least_one_cycle(
    spatial_network::MetaGraphsNext.MetaGraph,
    simple_cycles::Vector{Vector{Int64}})

    for bond in MetaGraphsNext.edge_labels(spatial_network)
        bond_found = false
        for cycle in simple_cycles
            if cycle_contains_bond(cycle, bond)
                bond_found = true
                break
            end
        end

        if !bond_found
            return false
        end
    end

    return true
end


"""
Get a list of all very strong rings in a spatial network as defined in 
10.1016/0022-3093(91)90145-V and explained in 10.1016/S0927-0256(01)00256-7.
"""
function get_very_strong_rings_vec(spatial_network::MetaGraphsNext.MetaGraph;
    max_ring_size_to_check::Int64 = 12)

    # get the number of simple cycles up to the max_ring_size_to_check
    simple_cycles = MetaGraphsNext.simplecycles_limited_length(spatial_network,
        max_ring_size_to_check, 10^6)

    # filter out cycles that contain only two vertices
    filter!(cycle -> length(cycle) > 2, simple_cycles)

    # avoid double counting of cycles that are the same but have inverse vertex
    # orderings
    for cycle in simple_cycles
        cycle_inverse = vcat(cycle[1], reverse(cycle[2:end]))
        if cycle_inverse in simple_cycles
            filter!(!=(cycle_inverse), simple_cycles)
        end
    end

    # to each cycle attach the first vertex at the end of the cycle such that
    # contained bonds are dicrecly represented
    for cycle in simple_cycles
        push!(cycle, cycle[1])
    end

    # check if each bond is contained in at least one cycle
    if !each_bond_in_at_least_one_cycle(spatial_network, simple_cycles)
        @warn "Not all bonds are contained in at least one cycle. Consider
            increasing the max_ring_size_to_check parameter."
    end

    # loop through all bonds and for each bond, find all very strong rings, 
    # that are all cycles with minimal length containing the bond
    very_strong_rings_vec = Vector{Vector{Int64}}()
    for bond in MetaGraphsNext.edge_labels(spatial_network)

        # get the cycles that contain the bond
        cycles_containing_bond = filter(cycle -> cycle_contains_bond(
            cycle, bond), simple_cycles)

        if length(cycles_containing_bond) > 0
            min_cycle_length = length(cycles_containing_bond[1])
            for cycle in cycles_containing_bond
                if length(cycle) < min_cycle_length
                    min_cycle_length = length(cycle)
                end
            end

            # get only rings with minimal length containing the current bond
            very_strong_rings_vec_current_bond = filter(
                cycle -> length(cycle) <= min_cycle_length, 
                cycles_containing_bond)

            # add the cycles to the list of very strong rings if they are not
            # already in the list
            for ring in very_strong_rings_vec_current_bond
                if !(ring in very_strong_rings_vec)
                    push!(very_strong_rings_vec, ring)
                end
            end
        end
    end

    # sort the rings by their length
    sort!(very_strong_rings_vec, by = length)

    return very_strong_rings_vec
end


"""
Get ring size distribution of the very strong rings as defined in 
10.1016/0022-3093(91)90145-V and explained in 10.1016/S0927-0256(01)00256-7.
"""
function get_ring_size_distribution(spatial_network::MetaGraphsNext.MetaGraph;
    max_ring_size_to_check::Int64 = 10,
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing)

    # get the very strong rings of the network
    very_strong_rings_vec = get_very_strong_rings_vec(spatial_network, 
        max_ring_size_to_check = max_ring_size_to_check)    

    # create histogram of ring lengths
    ring_length_vec = Vector{Int64}()
    for ring in very_strong_rings_vec
        push!(ring_length_vec, length(ring)-1)
    end

    ring_size_vec = collect(3:maximum(ring_length_vec))

    histogram = StatsBase.fit(StatsBase.Histogram, ring_length_vec, 
        2.5:1:ring_size_vec[end]+0.5, closed=:left)

    # normalize histogram
    histogram = LinearAlgebra.normalize(histogram, mode=:probability)
    ring_size_distribution = histogram.weights

    # convert the vector of very strong rings into a matrix to be ably to save 
    # it in the h5 file
    very_strong_rings_vec_mat = zeros(Int64, length(very_strong_rings_vec),
        maximum(length.(very_strong_rings_vec)))
    for i in 1:length(very_strong_rings_vec)
        very_strong_rings_vec_mat[i, 1:length(very_strong_rings_vec[i])] = 
            very_strong_rings_vec[i]
    end

    ring_size_distribution_dict = Dict{String, Any}(
        "ring_size_vec" => ring_size_vec,
        "ring_size_distribution" => ring_size_distribution,
        "very_strong_rings_vec_mat" => very_strong_rings_vec_mat)

    # add label to dictionary if label is not nothing
    if label !== nothing
        ring_size_distribution_dict["label"] = label
    end

    if save_result
        GU.save_dict_to_h5(copy(ring_size_distribution_dict),
            save_path*"_ring_size_distribution.h5")
    end
    
    return ring_size_distribution_dict
end


"""
Get the mean and standard deviation of the very strong ring sizes of a network
according to the definition in 10.1016/0022-3093(91)90145-V and 
10.1016/S0927-0256(01)00256-7 .
"""
function get_ring_statistics(ring_size_distribution_dict::Dict{String, Any})

    # get the ring size vector and the ring size distribution
    ring_size_vec = ring_size_distribution_dict["ring_size_vec"]
    ring_size_distribution = (
        ring_size_distribution_dict["ring_size_distribution"])

    # calculate mean and standard deviation of the ring size
    ring_size_mean = sum(ring_size_vec .* ring_size_distribution)
    ring_size_std = sum(
        (ring_size_vec .- ring_size_mean).^2 .* ring_size_distribution)

    return [ring_size_mean, ring_size_std]
end


"""
Get the normal vector of a polygon defined by its vertices. Use Newell's method
as described in 10.1016/B978-0-08-050755-2.50052-X and on page 493 of the book
"Real-Time Collision Detection" by Christer Ericson to estimate the normal 
vector along which the projection of the polygon has maximal size.
"""
function polygon_normal(points)
    n = zeros(3)
    for i in 1:length(points)
        p1 = points[i]
        p2 = points[mod1(i+1, length(points))]
        n[1] += (p1[2] - p2[2]) * (p1[3] + p2[3])
        n[2] += (p1[3] - p2[3]) * (p1[1] + p2[1])
        n[3] += (p1[1] - p2[1]) * (p1[2] + p2[2])
    end
    return LinearAlgebra.normalize(n)
end


"""
Project points onto a plane defined by a normal vector n. The projection is
done by creating an orthonormal basis on the plane and then projecting the
points onto the basis vectors.
"""
function project_to_plane(points, n)
    # Create orthonormal basis u, v on plane perpendicular to n
    u = LinearAlgebra.normalize(
        LinearAlgebra.cross(n, abs(n[1]) < 0.9 ? [1.0, 0.0, 0.0] 
        : [0.0, 1.0, 0.0]))
    v = LinearAlgebra.normalize(LinearAlgebra.cross(n, u))
    
    # Project each point onto 2D coordinates (u, v)
    return [(LinearAlgebra.dot(p, u), LinearAlgebra.dot(p, v)) for p in points]
end


"""
Get the distance from a point to a line segment defined by two endpoints a and 
b. The distance is calculated by projecting the point onto the line defined by 
a and b
"""
function point_to_segment_distance(point::Tuple{Float64, Float64}, 
        vertex_1::Tuple{Float64, Float64}, vertex_2::Tuple{Float64, Float64})
    segment_vec = vertex_2 .- vertex_1
    vertex_1_to_point_vec = point .- vertex_1
    t = clamp(LinearAlgebra.dot(vertex_1_to_point_vec, segment_vec) 
        / LinearAlgebra.dot(segment_vec, segment_vec), 0.0, 1.0)
    projection = vertex_1 .+ t .* segment_vec
    return LinearAlgebra.norm(point .- projection)
end


"""
Get the minimum distance from a point to the edges of a polygon defined by its
vertices. The polygon is represented as a vector of tuples, where each tuple
contains the x and y coordinates of a vertex. The point is represented as a

"""
function min_distance_to_bonds(point::Tuple{Float64, Float64},
        polygon_vertices::Vector{Tuple{Float64, Float64}})

    min_dist = Inf
    nr_bonds = length(polygon_vertices)-1  
    for i in 1:nr_bonds
        vertex_1 = polygon_vertices[i]
        vertex_2 = polygon_vertices[i+1] 
        point_to_bond_distance = point_to_segment_distance(
            point, vertex_1, vertex_2)
        min_dist = min(min_dist, point_to_bond_distance)
    end
    return min_dist
end


"""
Get ring radius distribution from the very strong rings of a network. To this
end, use Newell's method as described in 10.1016/B978-0-08-050755-2.50052-X
and on page 493 of the book "Real-Time Collision Detection" by Christer Ericson
to estimate the normal vector of each ring along which the ring has maximal
size. Then, determine the maximal radius of a circle that can be placed in the
projection of the ring onto the plane defined by the normal vector.
"""
function get_ring_radius_distribution(
    spatial_network::MetaGraphsNext.MetaGraph,
    ring_size_distribution_dict::Dict{String, Any};
    save_result::Bool = false,
    save_path = raw"..\analysis_data\sample_name",
    label = nothing)

    # get the very strong rings of the network
    very_strong_rings_mat = ring_size_distribution_dict[
        "very_strong_rings_vec_mat"]
    # convert the matrix of very strong rings into a vector of vectors
    very_strong_rings_vec = Vector{Vector{Int64}}()
    for i in 1:size(very_strong_rings_mat, 1)
        push!(very_strong_rings_vec, 
            filter(v -> v != 0, very_strong_rings_mat[i, :]))
    end

    # for each ring, get the normal vector and the maximal radius of a circle
    # that can be placed in the projection of the ring onto the plane defined
    # by the normal vector
    ring_radius_vec = Vector{Float64}()

    nr_rings = length(very_strong_rings_vec)
    for ring in very_strong_rings_vec

        # get the virtural positions of the vertices of the ring with respect
        # to the first vertex of the ring considering periodic boundary
        # conditions
        vertex_positions = [spatial_network[ring[1]]["position"] ]
        for vertex in ring[2:end]
            push!(vertex_positions, 
                NG.get_virtual_position(spatial_network[ring[1]]["position"],
                spatial_network[vertex]["position"],
                spatial_network[]["supercell_edge_length"]))
        end

        # get the normal vector of the ring
        normal_vector = polygon_normal(vertex_positions[1:end-1])

        # project the ring onto the plane defined by the normal vector
        projected_ring = project_to_plane(vertex_positions, normal_vector)

        # add some noise to the projected ring points to avoid issues with
        # the Polylabel package when rings are very symmetric
        d = Distributions.Normal(0, 1e-6)
        noise_map = Dict{Tuple{Float64,Float64},Tuple{Float64,Float64}}()
        function noisy_point(p, d, noise_map)
            if !haskey(noise_map, p)
                noise_map[p] = (Random.rand(d), Random.rand(d))  
            end
            return (p[1] + noise_map[p][1], p[2] + noise_map[p][2])
        end
        projected_ring_noisy = [noisy_point(p, d, noise_map) 
            for p in projected_ring]

        # create a polygon from the projected ring
        polygon = GeometryBasics.Polygon(
            GeometryBasics.Point{2, Float64}[projected_ring_noisy...])

        # determine the pole of inaccessibility of the polygon
        # (the point in the polygon that is farthest away from the edges of the
        # polygon)
        pole_of_inaccessibility = Polylabel.polylabel(polygon, atol=0.01)

        # get the distance from the pole of inaccessibility to the edges of
        # the polygon
        min_distance = min_distance_to_bonds(pole_of_inaccessibility, 
            projected_ring)
        # save the distance as the ring radius of the ring
        push!(ring_radius_vec, min_distance)

    end

    # create histogram of ring radii
    ring_radius_bin_centers = collect(0.1:0.1:round(
        maximum(ring_radius_vec), digits = 1))

    histogram = StatsBase.fit(StatsBase.Histogram, ring_radius_vec,
        ring_radius_bin_centers[1]-0.05:0.1
        :ring_radius_bin_centers[end]+0.0501, closed=:left)
    
    # normalize histogram
    histogram = LinearAlgebra.normalize(histogram, mode=:probability)

    ring_radius_distribution = histogram.weights

    ring_radius_distribution_dict = Dict{String, Any}(
        "ring_radius_vec" => ring_radius_bin_centers,
        "ring_radius_distribution" => ring_radius_distribution)

    # add label to dictionary if label is not nothing
    if label !== nothing
        ring_radius_distribution_dict["label"] = label
    end

    if save_result
        GU.save_dict_to_h5(copy(ring_radius_distribution_dict),
            save_path*"_ring_radius_distribution.h5")
    end
    
    return ring_radius_distribution_dict
end


"""
Get the mean and standard deviation of the radii of the very strong rings of a
network according to the definition in 10.1016/0022-3093(91)90145-V and
10.1016/S0927-0256(01)00256-7 . To get the radius of a ring in 3d space, the
ring is projected onto a plane using Newell's method as described in
10.1016/B978-0-08-050755-2.50052-X and on page 493 of the book "Real-Time
Collision Detection" by Christer Ericson. The radius is then calculated as
the distance from the pole of inaccessibility of the polygon to the edges
of the polygon.
"""
function get_ring_radius_statistics(
    ring_radius_distribution_dict::Dict{String, Any})

    # get the ring size vector and the ring size distribution
    ring_radius_vec = ring_radius_distribution_dict["ring_radius_vec"]
    ring_radius_distribution = (
        ring_radius_distribution_dict["ring_radius_distribution"])

    # calculate mean and standard deviation of the ring size
    ring_radius_mean = sum(ring_radius_vec .* ring_radius_distribution)
    ring_radius_std = sum(
        (ring_radius_vec .- ring_radius_mean).^2 .* ring_radius_distribution)

    return [ring_radius_mean, ring_radius_std]
end