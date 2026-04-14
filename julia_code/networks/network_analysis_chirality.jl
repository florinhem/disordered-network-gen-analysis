"""
these functions can be used to analyze the chirality of spatial networks
"""


"""
rotate a spatial network using a given rotation matrix
"""
function rotate_spatial_network(
    spatial_network::MetaGraphsNext.MetaGraph,
    rotation_matrix::Array{Float64,2})

    new_spatial_network = deepcopy(spatial_network)

    # rotate all vertices of the spatial network
    for vertex in MetaGraphsNext.labels(spatial_network)
        new_spatial_network[vertex]["position"] = (rotation_matrix 
            * (spatial_network[vertex]["position"] 
                .- spatial_network[]["supercell_edge_length"] ./ 2) 
            .+ spatial_network[]["supercell_edge_length"] ./ 2 )
    end

    # rotate all bond vectors of the spatial network
    for bond in MetaGraphsNext.edge_labels(spatial_network)
        new_spatial_network[bond...]["vector"] = (
            rotation_matrix * spatial_network[bond...]["vector"])
    end

    return new_spatial_network
end


"""
get the 8 spatial networks that result from the C3 rotations about the 8 
possible (111) axes of a cubic unit cell
(table 5.1  and fig. 5.2 in 10.1007/978-3-540-32899-5)
"""
function get_spatial_networks_c3_rotations(
    spatial_network::MetaGraphsNext.MetaGraph)

    # define the 8 (111) axes of a cubic unit cell (in positive and negative 
    # directions)
    axes = [[1.0, 1.0, 1.0], 
            [1.0, 1.0, -1.0], 
            [1.0, -1.0, 1.0], 
            [-1.0, 1.0, 1.0], 
            [1.0, -1.0, -1.0], 
            [-1.0, 1.0, -1.0], 
            [-1.0, -1.0, 1.0], 
            [-1.0, -1.0, -1.0]]

    # get the rotation matrices
    rotation_matrices = [NG.rotation_matrix(axis, 120.0) for axis in axes]

    spatial_networks = [
        rotate_spatial_network(spatial_network, rotation_matrix) 
        for rotation_matrix in rotation_matrices]
    
    return spatial_networks
end


"""
get the 3 spatial networks that result from the C2 rotations about the 3 
possible (100) axes of a cubic unit cell
(table 5.1  and fig. 5.2 in 10.1007/978-3-540-32899-5)
"""
function get_spatial_networks_c2_rotations(
    spatial_network::MetaGraphsNext.MetaGraph)

    # define the 3 (100) axes of a cubic unit cell (in positive directions)
    axes = [[1.0, 0.0, 0.0], 
            [0.0, 1.0, 0.0], 
            [0.0, 0.0, 1.0]]

    # get the rotation matrices
    rotation_matrices = [NG.rotation_matrix(axis, 180.0) for axis in axes]

    spatial_networks = [
        rotate_spatial_network(spatial_network, rotation_matrix) 
        for rotation_matrix in rotation_matrices]

    return spatial_networks
end


"""
get the 6 spatial networks that result from the C2' rotations about the 6
possible (110) axes of a cubic unit cell
(table 5.1  and fig. 5.2 in 10.1007/978-3-540-32899-5)
"""
function get_spatial_networks_c2_prime_rotations(
    spatial_network::MetaGraphsNext.MetaGraph)

    # define the 6 (110) axes of a cubic unit cell (in positive directions)
    axes = [[1.0, 1.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, -1.0, 0.0],
            [-1.0, 0.0, 1.0],
            [0.0, 1.0, -1.0]]

    # get the rotation matrices
    rotation_matrices = [NG.rotation_matrix(axis, 180.0) for axis in axes]

    spatial_networks = [
        rotate_spatial_network(spatial_network, rotation_matrix) 
        for rotation_matrix in rotation_matrices]

    return spatial_networks
end


"""
get the 6 spatial networks that result from the C4 rotations about the 6
possible (100) axes of a cubic unit cell
(table 5.1  and fig. 5.2 in 10.1007/978-3-540-32899-5)
"""
function get_spatial_networks_c4_rotations(
    spatial_network::MetaGraphsNext.MetaGraph)

    # define the 6 (100) axes of a cubic unit cell (in positive and negative 
    # directions)
    axes = [[1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [-1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0]]

    # get the rotation matrices
    rotation_matrices = [NG.rotation_matrix(axis, 90.0) for axis in axes]

    spatial_networks = [
        rotate_spatial_network(spatial_network, rotation_matrix) 
        for rotation_matrix in rotation_matrices]

    return spatial_networks
end


"""
get the positions of points decorating all bonds of a spatial network
"""
function get_decorated_spatial_network_points(
    spatial_network::MetaGraphsNext.MetaGraph;
    points_per_bond::Int64 = 3,
    periodic_boundary_conditions::Bool = true,
    exclude_layer_thickness::Float64 = 0.0,)

    # get minimal and maximal vertex coordinates along all three axes in case
    # of non-periodic boundary conditions
    if !periodic_boundary_conditions
        min_vertex_coords, max_vertex_coords = get_min_max_vertex_coords(
        spatial_network)
    end

    supercell_edge_length = spatial_network[]["supercell_edge_length"]

    points = Vector{Vector{Float64}}()

    for vertex in MetaGraphsNext.labels(spatial_network)
        if periodic_boundary_conditions
            push!(points, spatial_network[vertex]["position"])
        else
            # get the positions of the vertex
            vertex_pos = spatial_network[vertex]["position"]

            # check if the vertex is within the excluded layer thickness
            if (all(vertex_pos 
                    .> (min_vertex_coords .+ exclude_layer_thickness))
                && all(vertex_pos 
                    .< (max_vertex_coords .- exclude_layer_thickness)))

                push!(points, spatial_network[vertex]["position"])
            end
        end
    end

    if periodic_boundary_conditions 
        considered_bonds = MetaGraphsNext.edge_labels(spatial_network)
    else
        # get all considered bonds
        considered_bonds = get_considered_bonds(
            spatial_network;
            periodic_boundary_conditions=periodic_boundary_conditions,
            exclude_layer_thickness=exclude_layer_thickness)
    end

    for bond in considered_bonds
        start_position = spatial_network[bond[1]]["position"]
        bond_vector = spatial_network[bond...]["vector"]

        for i in 1:points_per_bond-1
            point = (start_position .+ (i / points_per_bond) .* bond_vector)

            # shift the point inside the supercell using PBCs
            for dim in 1:3
                if point[dim] < 0.0
                    point[dim] += supercell_edge_length
                elseif point[dim] >= supercell_edge_length
                    point[dim] -= supercell_edge_length
                end
            end
            push!(points, point)
        end
    end

    return points
end


"""
for all symmetry operations of the chiral octahedral group O (table 5.1 
in 10.1007/978-3-540-32899-5), get the points decorating the spatial networks 
that result from applying these symmetry operations to a given spatial network
"""
function get_decorated_spatial_networks_chiral_o_group(
    spatial_network::MetaGraphsNext.MetaGraph;
    points_per_bond::Int64 = 3,)

    all_spatial_network_points = Vector{Vector{Vector{Float64}}}()

    # identity E
    push!(all_spatial_network_points, 
        get_decorated_spatial_network_points(
            spatial_network;
            points_per_bond = points_per_bond))

    # C3 rotations about (111) axes
    spatial_networks_c3 = get_spatial_networks_c3_rotations(spatial_network)
    for spatial_network_c3 in spatial_networks_c3
        push!(all_spatial_network_points, 
            get_decorated_spatial_network_points(
                spatial_network_c3;
                points_per_bond = points_per_bond))
    end

    # C2 rotations about (100) axes
    spatial_networks_c2 = get_spatial_networks_c2_rotations(spatial_network)
    for spatial_network_c2 in spatial_networks_c2
        push!(all_spatial_network_points, 
            get_decorated_spatial_network_points(
                spatial_network_c2;
                points_per_bond = points_per_bond))
    end

    # C2' rotations about (110) axes
    spatial_networks_c2_prime = get_spatial_networks_c2_prime_rotations(
        spatial_network)
    for spatial_network_c2_prime in spatial_networks_c2_prime
        push!(all_spatial_network_points, 
            get_decorated_spatial_network_points(
                spatial_network_c2_prime;
                points_per_bond = points_per_bond))
    end

    # C4 rotations about (100) axes
    spatial_networks_c4 = get_spatial_networks_c4_rotations(spatial_network)
    for spatial_network_c4 in spatial_networks_c4
        push!(all_spatial_network_points, 
            get_decorated_spatial_network_points(
                spatial_network_c4;
                points_per_bond = points_per_bond))
    end

    return all_spatial_network_points
end


"""
get the overlap of two spatial networks decorated with 3d Gaussian functions
"""
function get_overlap_spatial_networks(
    spatial_network_1_points::Vector{Vector{Float64}},
    spatial_network_2_points::Vector{Vector{Float64}},
    supercell_edge_length::Float64;
    sigma::Float64 = 1/6)

    overlap = 0.0
    for point_1 in spatial_network_1_points
        for point_2 in spatial_network_2_points
            distance_squared = NG.get_distance_sq_pbc(
                point_1, point_2, supercell_edge_length)
            overlap += exp(-distance_squared / (4.0 * sigma^2))
        end
    end

    overlap *=  pi^(3/2) * (sigma)^3

    return overlap
end


"""
get the self overlap of the projection of the spatial network onto an 
irreducible representation of the chiral octahedral group O with given
dimension and characters
(table 5.1  in 10.1007/978-3-540-32899-5)
"""
function get_irrep_projection_self_overlap(
    group_order::Int64,
    irrep_dimension::Int64,
    characters_by_classes::Vector{Int64},
    total_state_self_overlap::Float64,
    all_spatial_network_points::Vector{Vector{Vector{Float64}}},
    supercell_edge_length::Float64;
    sigma::Float64 = 1/6)

    characters_by_operations = vcat(
        fill(characters_by_classes[1], 1), 
        fill(characters_by_classes[2], 8), 
        fill(characters_by_classes[3], 3), 
        fill(characters_by_classes[4], 6), 
        fill(characters_by_classes[5], 6))

    # loop through all symmetry operations of the group
    irrep_projection_self_overlap = 0.0

    for (operation_1, character_1) in enumerate(characters_by_operations)

        # add the self overlap of the total state multiplied by the square of
        # the current character to the projected self overlap
        irrep_projection_self_overlap += character_1^2 * total_state_self_overlap

        # get the spatial network points for the current operation
        spatial_network_1_points = all_spatial_network_points[operation_1]

        # loop through all other operations and add the overlap of the current
        # operation with the other operation multiplied by the product of their 
        # characters to the projected self overlap
        for operation_2 in operation_1+1:length(characters_by_operations)

            character_2 = characters_by_operations[operation_2]

            spatial_network_2_points = all_spatial_network_points[operation_2]

            overlap = get_overlap_spatial_networks(
                spatial_network_1_points, spatial_network_2_points, 
                supercell_edge_length; sigma = sigma)

            irrep_projection_self_overlap += 2.0 * character_1 * character_2 * overlap
        end
    end

    irrep_projection_self_overlap *= (irrep_dimension / group_order)^2

    return irrep_projection_self_overlap
end


function get_all_irrep_projection_self_overlaps(
    spatial_network::MetaGraphsNext.MetaGraph;
    points_per_bond::Int64 = 3,
    sigma::Float64 = 1/6,
    print_progress::Bool = false,
    thread_nr::Int64 = 0,
    print_lock = Threads.ReentrantLock())

    all_spatial_network_points = get_decorated_spatial_networks_chiral_o_group(
        spatial_network,
        points_per_bond = points_per_bond)

    total_state_self_overlap = get_overlap_spatial_networks(
        all_spatial_network_points[1], all_spatial_network_points[1], 
        spatial_network[]["supercell_edge_length"]; sigma = sigma)

    group_order = 24

    # get the irrep dimensions for the irreps A1, A2, E, T1, T2
    irrep_dimensions = [1, 1, 2, 3, 3]

    # get the characters for the classes E, 8C3, 3C2, 6C2', 6C4 for all irreps
    character_table = [[1, 1, 1, 1, 1],
                       [1, 1, 1, -1, -1],
                       [2, -1, 2, 0, 0],
                       [3, 0, -1, -1, 1],
                       [3, 0, -1, 1, -1]]

    # get the projected self overlaps for all irreps
    irrep_projection_self_overlaps = Vector{Float64}()

    for (irrep_nr, irrep_dimension) in enumerate(irrep_dimensions)
        characters_by_classes = character_table[irrep_nr]

        irrep_projection_self_overlap = get_irrep_projection_self_overlap(
            group_order, irrep_dimension, characters_by_classes, 
            total_state_self_overlap, all_spatial_network_points, 
            spatial_network[]["supercell_edge_length"]; sigma = sigma)
        push!(irrep_projection_self_overlaps, irrep_projection_self_overlap)

        if print_progress
            lock(print_lock) do
                Format.printfmtln("Chirality metric for thread nr {1:d}: 
                irrep projection nr {2:d}/5 finished", thread_nr, irrep_nr)
            end
        end
    end

    # normalize the projected self overlaps by the total state self overlap
    irrep_contributions = [
        irrep_projection_self_overlap / total_state_self_overlap 
        for irrep_projection_self_overlap in irrep_projection_self_overlaps]

    return irrep_contributions
end


"""
calculate the Hausdorff chirality measure (HCM) of a spatial network using the
Hausdorff distance with periodic boundary conditions (PBCs) between the 
original spatial network and the enantiomer obtained by point inversion of the 
original spatial network
"""
function get_hausdorff_dist(
    set1::AbstractMatrix, set2::AbstractMatrix, supercell_edge_length::Real,
    periodic_boundary_conditions::Bool)
    
    if periodic_boundary_conditions
        # Define the metric for periodic boundary conditions
        metric = Distances.PeriodicEuclidean(
            [supercell_edge_length, supercell_edge_length, 
            supercell_edge_length])
    else
        # Define the standard Euclidean metric for non-periodic boundary 
        # conditions
        metric = Distances.Euclidean()
    end
    
    # calculate pairwise distances under Minimum Image Convention
    distances_mat = Distances.pairwise(metric, set1, set2)
    
    # Directed Hausdorff distance
    h1 = maximum(minimum(distances_mat, dims=2)) 
    h2 = maximum(minimum(distances_mat, dims=1))
    
    return max(h1, h2)
end


"""
calculate the Hausdorff chirality measure (HCM) of a spatial network using the
Hausdorff distance with periodic boundary conditions (PBCs) between the
original spatial network and the enantiomer obtained by point inversion of the
original spatial network, where the HCM is defined as the minimum Hausdorff
distance between the original spatial network and the enantiomer obtained by
point inversion of the original spatial network after applying any possible
rotation to the enantiomer, normalized by the diameter of the structure
https://pubs.acs.org/doi/abs/10.1021/ja00041a016
"""
function get_hausdorff_chirality(
    spatial_network::MetaGraphsNext.MetaGraph;
    points_per_bond::Int64 = 3,
    exclude_layer_thickness::Float64 = 0.0,
    periodic_boundary_conditions::Bool = true,)

    supercell_edge_length = spatial_network[]["supercell_edge_length"]

    spatial_network_points = get_decorated_spatial_network_points(
        spatial_network;
        points_per_bond = points_per_bond,
        periodic_boundary_conditions=periodic_boundary_conditions,
        exclude_layer_thickness=exclude_layer_thickness)

    # shift the points to the center of the supercell for easier handling of 
    # PBCs
    for point in spatial_network_points
        for dim in 1:3
            point[dim] -= supercell_edge_length / 2
        end
    end

    # convert the list of points to a 3xN matrix
    points = hcat(spatial_network_points...)
    
    inversion_matrix = LinearAlgebra.Diagonal([-1.0, -1.0, -1.0])
    inverted_points = inversion_matrix * points

    # find the rotation and translation that minimize the Hausdorff distance 
    # between the original points and the inverted points under PBCs
    function objective(params)
        angles = params[1:3]
        translation = params[4:6]

        rotation_matrix = Rotations.RotXYZ(angles[1], angles[2], angles[3])
        transformed_points = (
            Array(rotation_matrix) * inverted_points) .+ translation

        hausdorff_dist = get_hausdorff_dist(
            points, transformed_points, supercell_edge_length, 
            periodic_boundary_conditions)
        return hausdorff_dist
    end
    
    if periodic_boundary_conditions

        # calculate the diameter of the structure as the maximum distance between
        # any two points in the original spatial network. This equals half the
        # diagonal of the cubic supercell 
        diameter = supercell_edge_length * sqrt(3) / 2

        # Start at zero rotation and zero translation
        initial_params = zeros(6)

        # Optimization using Nelder-Mead
        res = Optim.optimize(objective, initial_params, Optim.NelderMead(), 
                             Optim.Options(iterations=2000, g_tol=1e-6))
        min_hausdorff = Optim.minimum(res)
        
    else
        # get the minimal and maximal vertex coords along the three axes
        min_vertex_coords, max_vertex_coords = get_min_max_vertex_coords(
            spatial_network)

        # get the diameter from these extremal vertex coordinates
        diameter = sqrt(sum((max_vertex_coords .- min_vertex_coords).^2)) / 2

        # Define six different initial guesses that rotate the structure to
        # orientations +-z, +-y, +-x and do not apply any translation
        guesses = [
             [ 0,       0,      0,       0, 0, 0 ],
             [ π,       0,      0,       0, 0, 0 ],
             [ 0,  -π/2,      0,       0, 0, 0 ],
             [ 0,   π/2,      0,       0, 0, 0 ],
             [ π/2,     0,      0,       0, 0, 0 ],
             [ -π/2,    0,      0,       0, 0, 0 ]
            ]

        # map the optimization over the list of guesses
        results = map(guesses) do initial_params
            Optim.optimize(objective, initial_params, Optim.NelderMead(), 
                           Optim.Options(iterations=2000, g_tol=1e-6))
        end

        # find the best optimization result object
        best_result_object = argmin(res -> Optim.minimum(res), results)
        min_hausdorff = Optim.minimum(best_result_object)
    end

    # to get the HCM, we normalize the minimum Hausdorff distance by the 
    # diameter of the structure
    hcm = min_hausdorff / diameter
    
    return hcm
end


"""
Finds the optimal inversion center that minimizes the continuous chirality
measure CCM as defined in 10.1021/ja00106a053. We choose the normalization of
10.1103/PhysRevB.110.174112 that should yield CCM values in approximately the
range [0, 1]
"""
function get_continuous_chirality_measure(
    spatial_network::MetaGraphsNext.MetaGraph;
    points_per_bond::Int64 = 3,
    exclude_layer_thickness::Float64 = 0.0,
    periodic_boundary_conditions::Bool = true,)

    supercell_edge_length = spatial_network[]["supercell_edge_length"]

    spatial_network_points = get_decorated_spatial_network_points(
        spatial_network;
        points_per_bond = points_per_bond,
        periodic_boundary_conditions=periodic_boundary_conditions,
        exclude_layer_thickness=exclude_layer_thickness)

    # shift the points to the center of the supercell for easier handling of 
    # PBCs
    for point in spatial_network_points
        for dim in 1:3
            point[dim] -= supercell_edge_length / 2
        end
    end

    nr_points = length(spatial_network_points)
    
    # Convert Vector of Vectors to a 3xN Matrix
    points_mat = hcat(spatial_network_points...)
    
    # calculate the inverted coordinates at the origin
    inverted_base = -1.0 .* points_mat
    
    # We allocate these once to prevent the Garbage Collector from 
    # freezing the optimization loop.
    transformed_points = zeros(Float64, 3, nr_points)
    dist_sq_mat = zeros(Float64, nr_points, nr_points)
    
    # Optimization objective
    function objective(translation_params::Vector{Float64})
        
        # Update transformed points in-place
        @inbounds for j in 1:nr_points
            for d in 1:3
                # Apply translation and wrap into [0, L]
                val = inverted_base[d, j] + translation_params[d]
                transformed_points[d, j] = mod(val, supercell_edge_length)
            end
        end
        
        # get the squared distance matrix considering periodic boundary
        # conditions
        @inbounds for j in 1:nr_points
            for i in 1:nr_points
                dsq = 0.0
                for d in 1:3
                    # Absolute distance along dimension d
                    delta = abs(points_mat[d, i] - transformed_points[d, j])
                    
                    if periodic_boundary_conditions
                        # Apply Minimum Image Convention (PBC)
                        delta = min(delta, supercell_edge_length - delta)
                    end
                    
                    dsq += delta^2
                end
                dist_sq_mat[i, j] = dsq
            end
        end
        
        # Solve the exact 1-to-1 assignment problem
        assignment, cost = Hungarian.hungarian(dist_sq_mat)
        
        # Calculate the sum of squared distances to the achiral midpoint
        total_dist_to_achiral_sq = cost / 4.0
        
        # Normalize by the number of points
        return total_dist_to_achiral_sq / nr_points
    end

    # Set up the optimization
    initial_guess = zeros(Float64, 3)
    
    # We restrict the iterations slightly since Nelder-Mead with a global 
    # assignment problem can still be heavy for very large nr_points.
    opt_options = Optim.Options(iterations=200, show_trace=false)
    
    res = Optim.optimize(
        objective, 
        initial_guess, 
        Optim.NelderMead(), 
        opt_options
    )
    
    # Extract continuous chirality measure
    min_ccm = Optim.minimum(res)
    
    return min_ccm
end