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
    points_per_bond::Int64 = 3,)

    supercell_edge_length = spatial_network[]["supercell_edge_length"]

    points = Vector{Vector{Float64}}()

    for vertex in MetaGraphsNext.labels(spatial_network)
        push!(points, spatial_network[vertex]["position"])
    end

    for bond in MetaGraphsNext.edge_labels(spatial_network)
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
            distance_vector = NG.get_distance_vector_pbc(
                point_1, point_2, supercell_edge_length)
            distance_squared = LinearAlgebra.norm(distance_vector)^2
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
function hausdorff_dist_pbc(
    set1::AbstractMatrix, set2::AbstractMatrix, supercell_edge_length::Real)
    
    metric = Distances.PeriodicEuclidean(
        [supercell_edge_length, supercell_edge_length, supercell_edge_length])
    
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
    points_per_bond::Int64 = 3,)

    supercell_edge_length = spatial_network[]["supercell_edge_length"]

    spatial_network_points = get_decorated_spatial_network_points(
        spatial_network;
        points_per_bond = points_per_bond,)

    # convert the list of points to a 3xN matrix
    points = hcat(spatial_network_points...)
    
    inversion_matrix = LinearAlgebra.Diagonal([-1.0, -1.0, -1.0])
    inverted_points = inversion_matrix * points
    
    # calculate the diameter of the structure as the maximum distance between
    # any two points in the original spatial network. This equals half the
    # diagonal of the cubic supercell 
    diameter = supercell_edge_length * sqrt(3) / 2

    # find the rotation and translation that minimize the Hausdorff distance 
    # between the original points and the inverted points under PBCs
    function objective(params)
        angles = params[1:3]
        translation = params[4:6]
        
        rotation_matrix = Rotations.RotXYZ(angles[1], angles[2], angles[3])
        
        transformed_points = (
            Array(rotation_matrix) * inverted_points) .+ translation
        
        hausdorff_dist = hausdorff_dist_pbc(
            points, transformed_points, supercell_edge_length)
        return hausdorff_dist
    end
    
    # Start at zero rotation and zero translation
    initial_params = zeros(6)
    
    # Optimization using Nelder-Mead
    res = Optim.optimize(objective, initial_params, Optim.NelderMead(), 
                         Optim.Options(iterations=2000, g_tol=1e-6))
    
    min_hausdorff = Optim.minimum(res)
    
    # to get the HCM, we normalize the minimum Hausdorff distance by the 
    # diameter of the structure
    hcm = min_hausdorff / diameter
    
    return hcm
end

