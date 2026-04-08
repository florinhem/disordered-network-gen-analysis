"""
These functions are used to calculate energies in spatial networks
"""


"""
Get bond stretching energy for a given bond
"""
function local_bond_stretching_energy_keating(
    spatial_network::MetaGraphsNext.MetaGraph,
    bond::Tuple{Int64, Int64})

    bond_stretching_energy = (
        3/16 * (spatial_network[bond...]["distance_squared"] - 1 )^2 ) 
    
    return bond_stretching_energy
end


"""
Get bond bending energy for a given vertex
"""
function local_bond_bending_energy_keating(
    spatial_network::MetaGraphsNext.MetaGraph,
    vertex_label::Int64)

    neighbor_label_vec::Vector{Int64} = collect(MetaGraphsNext.neighbor_labels(
        spatial_network, 
        vertex_label))

    vertex_coordination_nr=length(neighbor_label_vec)

    theta_ground_state=spatial_network[]["theta_ground_state"]

    bond_bending_sum::Float64 = 0.0

    @inbounds @simd for j in 1:(vertex_coordination_nr-1)
        sign1::Int64=sign(neighbor_label_vec[j] - vertex_label)
        vector_j::Vector{Float64}=
            spatial_network[vertex_label, neighbor_label_vec[j]]["vector"]
        
        @inbounds @simd for k in j+1:vertex_coordination_nr
            sign2::Int64=sign1*sign(neighbor_label_vec[k] - vertex_label)
            vector_k::Vector{Float64}=
                spatial_network[vertex_label, neighbor_label_vec[k]]["vector"]
            bond_bending_sum += (sign2*LinearAlgebra.dot(vector_j,vector_k) 
                - cosd(theta_ground_state) )^2 
            
        end
    end

    bond_bending_energy::Float64 = (
        3/8 * spatial_network[]["bond_bending_const"] * bond_bending_sum)

    return bond_bending_energy
end


"""
Get torsional energy for a given bond
"""
function local_torsional_energy_keating(
        spatial_network::MetaGraphsNext.MetaGraph,
        bond::Tuple{Int64, Int64})

    # since for a bond (i,j), always i<j, the vector by definition points from 
    # i to j and the sign is positive
    bond_vector = spatial_network[bond...]["vector"]

    # get the normalized bond vector
    bond_vector_normalized = bond_vector / sqrt(
        spatial_network[bond...]["distance_squared"])

    neighbor_label_vec_1::Vector{Int64} = [x for x in 
        MetaGraphsNext.neighbor_labels(spatial_network, bond[1]) 
        if x != bond[2]]

    neighbor_label_vec_2::Vector{Int64} = [x for x in 
        MetaGraphsNext.neighbor_labels(spatial_network, bond[2]) 
        if x != bond[1]]

    vertex_1_coordination_nr = spatial_network[bond[1]]["coordination_nr"]
    vertex_2_coordination_nr = spatial_network[bond[2]]["coordination_nr"]

    # get the number of minima in the torsional potential which is the least 
    # common multiple of the coordination numbers of the two vertices minus one 
    # since the bond itself is not included
    l = lcm(vertex_1_coordination_nr-1, vertex_2_coordination_nr-1)

    # delta phi is the favoured dihedral angle
    delta_phi = spatial_network[]["delta_phi"]

    torsional_sum::Float64 = 0.0

    @inbounds @simd for j in 1:(vertex_1_coordination_nr-1)
        sign1::Int64=sign(bond[1] - neighbor_label_vec_1[j])
        vector_j::Vector{Float64} = (sign1 .*
            spatial_network[bond[1], neighbor_label_vec_1[j]]["vector"])

        # first normal vector
        n1 = LinearAlgebra.cross(vector_j, bond_vector)

        @inbounds @simd for k in 1:(vertex_2_coordination_nr-1)
            sign2::Int64=sign(neighbor_label_vec_2[k] - bond[2])
            vector_k::Vector{Float64} = (sign2 .*
                spatial_network[bond[2], neighbor_label_vec_2[k]]["vector"])
            
            # second normal vector
            n2 = LinearAlgebra.cross(bond_vector, vector_k)

            # compute dihedral angle phi robustly
            x = LinearAlgebra.dot(n1, n2)
            y = LinearAlgebra.dot(bond_vector_normalized, 
                LinearAlgebra.cross(n1, n2))
            phi = atan(y, x)

            # torsional energy
            torsional_sum += 1 - cos(l * (phi - delta_phi*pi/180))
        end
    end
    
    torsional_energy::Float64 = (
        1/2 * spatial_network[]["torsional_const"] * torsional_sum)

    return torsional_energy
end


"""
Calculate the bending energy of a spatial network
"""
function get_bending_energy_keating(spatial_network::MetaGraphsNext.MetaGraph)

    bending_energy = 0.0

    for vertex in MetaGraphsNext.labels(spatial_network)
        bending_energy += local_bond_bending_energy_keating(
            spatial_network, vertex)
    end

    return bending_energy
end


"""
Calculate the stretching energy of a spatial network
"""
function get_stretching_energy_keating(spatial_network::MetaGraphsNext.MetaGraph)

    stretching_energy = 0.0

    for bond in MetaGraphsNext.edge_labels(spatial_network)
        stretching_energy += local_bond_stretching_energy_keating(
            spatial_network, bond)
    end

    return stretching_energy
end


"""
Calculate the torsional energy of a spatial network
"""
function get_torsional_energy_keating(spatial_network::MetaGraphsNext.MetaGraph)

    torsional_energy = 0.0

    for bond in MetaGraphsNext.edge_labels(spatial_network)
        torsional_energy += local_torsional_energy_keating(
            spatial_network, bond)
    end

    return torsional_energy
end


"""
Calculate the total energy of a spatial network
"""
function get_total_energy_keating(spatial_network::MetaGraphsNext.MetaGraph)

    total_energy = 0.0

    total_energy+=get_bending_energy_keating(
        spatial_network::MetaGraphsNext.MetaGraph)

    total_energy+=get_stretching_energy_keating(
        spatial_network::MetaGraphsNext.MetaGraph)

    if (haskey(spatial_network[], "torsional_const") 
        && spatial_network[]["torsional_const"] > 0.0)
        total_energy+=get_torsional_energy_keating(
            spatial_network::MetaGraphsNext.MetaGraph)
    end

    return total_energy
end


"""
Get energy for a cluster of vertices whose vertices and bonds are stored in the
    respective dictionary
"""
function get_cluster_energy(spatial_network, cluster_dict)

    # initialize cluster energy
    cluster_energy = 0.0

    # get vector of all cluster vertices
    all_cluster_vertices_vec = vcat(
        cluster_dict["cluster_vertices_to_move_vec"], 
        cluster_dict["cluster_vertices_outer_shell_vec"]) 

    for vertex in all_cluster_vertices_vec
        cluster_energy += local_bond_bending_energy_keating(
            spatial_network, vertex)
    end

    # loop through all bonds inside the cluster
    for bond in cluster_dict["cluster_bonds_inside_vec"]
        cluster_energy += local_bond_stretching_energy_keating(
            spatial_network, bond)
    end

    # loop through all bonds on the edge of the cluster which only contribute
    # half
    for bond in cluster_dict["cluster_bonds_edge_vec"]
        cluster_energy += 1/2 * local_bond_stretching_energy_keating(
            spatial_network, bond)
    end

    # get torsional energy if desired
    if (haskey(spatial_network[], "torsional_const") 
        && spatial_network[]["torsional_const"] > 0.0)

        for bond in cluster_dict["cluster_bonds_inside_vec"]
            cluster_energy += local_torsional_energy_keating(
                spatial_network, bond)
        end
    end

    return cluster_energy
end


"""
Calculate the negative Keating force (-f) on a given vertex
which is the gradient of its local Keating energy from
the distances and vectors that have already been calculated
"""
function gradient_keating_efficient(
    spatial_network::MetaGraphsNext.MetaGraph,
    central_vertex::Int64)

    # initialize gradient
    gradient = zeros(spatial_network[]["nr_dimensions"])

    # get vector of neighbors
    neighbor_vec = collect(
        MetaGraphsNext.neighbor_labels(spatial_network, central_vertex))

    vertex_coordination_nr=length(neighbor_vec)

    theta_ground_state=spatial_network[]["theta_ground_state"]

    # loop through central vertex's neighbors
    for j in 1:vertex_coordination_nr
        
        # get vector pointing from central vertex to neighbor j
        distance_vector_j = (float(sign(neighbor_vec[j] - central_vertex))
            * spatial_network[central_vertex, neighbor_vec[j]]["vector"])

        # get bond stretching term
        bond_stretching_term = ( - 3/4 * ( 
            spatial_network[central_vertex, neighbor_vec[j]][
                "distance_squared"] - 1 ) ) .* distance_vector_j

        # get bond bending term
        bond_bending_sum = zeros(spatial_network[]["nr_dimensions"])

        for k in j+1:vertex_coordination_nr
            
            # get vector pointing from central vertex to neighbor k
            distance_vector_k = (float(sign(neighbor_vec[k] - central_vertex))
                * spatial_network[central_vertex, neighbor_vec[k]]["vector"])

            bond_bending_sum .-= ( (
                3/4 * spatial_network[]["bond_bending_const"]  
                * ( LinearAlgebra.dot( distance_vector_j, distance_vector_k ) 
                - cosd(theta_ground_state) ) )  
            .* ( distance_vector_j .+ distance_vector_k )) 
            
        end

        # get bond bending terms due to next to nearest neighbors
        neighbor_bond_bending_sum = zeros(spatial_network[]["nr_dimensions"])

        # get neighbors of current neighbor excluding the central vertex
        neighbors_neighbor_vec = setdiff( collect(
            MetaGraphsNext.neighbor_labels(spatial_network, neighbor_vec[j])
                ), central_vertex )

        neighbor_vertex_coordination_nr=length(neighbors_neighbor_vec)

        for l in 1:neighbor_vertex_coordination_nr-1
            
            # get vector pointing from next-to-nearest-neighbor l to neighbor j
            distance_vector_l = (
                    float(sign(neighbor_vec[j] - neighbors_neighbor_vec[l]))
                * spatial_network[neighbors_neighbor_vec[l], 
                    neighbor_vec[j]]["vector"])

            neighbor_bond_bending_sum .-= (
                3/4 * spatial_network[]["bond_bending_const"]  
                * ( LinearAlgebra.dot( distance_vector_j, distance_vector_l )  
                    - cosd(theta_ground_state) )    
                .* distance_vector_l )
            
        end

        # sum bond stretching and bending terms
        gradient .+= (bond_stretching_term 
            .+ bond_bending_sum
            .+ neighbor_bond_bending_sum)

    end

    return gradient
end



"""
Calculate the Hessian matrix for a given vertex from
the distances and vectors that have already been calculated
"""
function hessian_keating_efficient(
    spatial_network::MetaGraphsNext.MetaGraph,
    central_vertex::Int64)

    nr_dimensions=spatial_network[]["nr_dimensions"]

    # initialize hessian
    hessian = zeros(nr_dimensions,nr_dimensions)
    
    # get vector of neighbors
    neighbor_vec = collect(
        MetaGraphsNext.neighbor_labels(spatial_network, central_vertex))

    vertex_coordination_nr=length(neighbor_vec)

    theta_ground_state=spatial_network[]["theta_ground_state"]

    bond_bending_const=spatial_network[]["bond_bending_const"]

    for a in 1:nr_dimensions

        # since the hessian is symmetric, it only needs to be calculated for
        # the lower left half
        for b in a:nr_dimensions

            for j in 1:vertex_coordination_nr

                # get vector pointing from central vertex to neighbor j
                distance_vector_j = (float(sign(neighbor_vec[j] 
                        - central_vertex))
                    * spatial_network[central_vertex, neighbor_vec[j]][
                        "vector"])
        
                # get bond stretching term
                bond_stretching_term = (
                        3/2 * distance_vector_j[b] * distance_vector_j[a]
                    + ==(a,b) * 3/4 * (
                    spatial_network[central_vertex, neighbor_vec[j]][
                        "distance_squared"] - 1 ) )
        
                # get bond bending term 
                bond_bending_sum = 0.0
        
                for k in j+1:vertex_coordination_nr

                    # get vector pointing from central vertex to neighbor k
                    distance_vector_k = (
                            float(sign(neighbor_vec[k] - central_vertex))
                        * spatial_network[central_vertex, neighbor_vec[k]][
                            "vector"])
        
                    bond_bending_sum += ( 
                            3/4 * bond_bending_const  
                        *( (distance_vector_j .+ distance_vector_k)[b]
                            *(distance_vector_j .+ distance_vector_k)[a]
                            + ==(a,b) * 2 *( LinearAlgebra.dot(
                                    distance_vector_j, distance_vector_k ) 
                                    - cosd(theta_ground_state) )))
                    
                end

                # get bond bending terms due to next to nearest neighbors
                neighbor_bond_bending_sum = 0.0

                # get neighbors of current neighbor excluding the central
                # vertex
                neighbors_neighbor_vec = setdiff( collect(
                    MetaGraphsNext.neighbor_labels(spatial_network, 
                        neighbor_vec[j]) ), central_vertex )
        
                for l in 1:vertex_coordination_nr-1

                    # get vector pointing from next-to-nearest-neighbor l to
                    # neighbor j
                    distance_vector_l = (float(sign(neighbor_vec[j] 
                            - neighbors_neighbor_vec[l]))
                        * spatial_network[neighbors_neighbor_vec[l],
                        neighbor_vec[j]]["vector"])
        
                    neighbor_bond_bending_sum += (  
                        3/4 * bond_bending_const  
                        * distance_vector_l[b] * distance_vector_l[a] )
                    
                end
        
                # sum bond stretching and bending terms
                hessian[a,b] += (bond_stretching_term 
                    + bond_bending_sum 
                    + neighbor_bond_bending_sum)
        
            end
            
        end

        # use symmetry of hessian and copy entries from lower left to upper
        # right half
        for b in 1:a-1
            hessian[a,b] = hessian[b,a]
        end
    end

    return LinearAlgebra.Symmetric(hessian)
end



"""
Get relaxation weight of a single vertex, which contributes to the Metropolis
acceptance probability when thermal fluctuations are included
"""
function get_vertex_relaxation_weight(
    excited_spatial_network::MetaGraphsNext.MetaGraph, 
    relaxed_spatial_network::MetaGraphsNext.MetaGraph, 
    vertex::Int64, 
    temperature)

    # get displacement vector 
    displacement_vec = (excited_spatial_network[vertex]["position"]
        .- relaxed_spatial_network[vertex]["position"])

    # get hessian matrix of current vertex at relaxed position
    hessian = hessian_keating_efficient(relaxed_spatial_network, vertex)

    # get vector of sigmas
    sigma_vec = sqrt.(temperature ./ LinearAlgebra.diag(hessian) ) 

    # get vector of gauss functions
    gauss_fct_vec = ( (1/sqrt(2*pi)) ./ sigma_vec
        .* exp.( .- displacement_vec.^2 ./ (2 .* (sigma_vec .^2 ) ) ) )

    # get relaxation weight which is the product of the gauss fct for all
    # dimensions
    vertex_relaxation_weight = prod(gauss_fct_vec)

    return vertex_relaxation_weight
end


"""
Get total weight corresponding to thermal fluctuations when relaxing or 
exciting a cluster. This weight contributes to the Metropolis acceptance
probability when thermal fluctuations are included
"""
function get_cluster_fluctuation_weight(
    excited_spatial_network::MetaGraphsNext.MetaGraph, 
    relaxed_spatial_network::MetaGraphsNext.MetaGraph, 
    cluster_dict::Dict,
    temperature)

    # initialize cluster relaxation weight
    cluster_relaxation_weight = 1

    # loop through all vertices that are allowed to move
    for vertex in cluster_dict["cluster_vertices_to_move_vec"]

        # get relaxation weight for current vertex
        vertex_relaxation_weight = get_vertex_relaxation_weight(
            excited_spatial_network, 
            relaxed_spatial_network, 
            vertex, 
            temperature)

        cluster_relaxation_weight *= vertex_relaxation_weight
    end

    return cluster_relaxation_weight
end


"""
Determine the Keating energy from the vertex coordinates with the option of
including a non-Keating torsional energy
"""
function cluster_energy_optim(
    vertices_to_move_coord_vec::AbstractVector{T},
    vertices_outer_shell_coord_vec::Vector{Float64},
    vertices_bonds_edge_coord_vec::Vector{Float64},
    vertices_bonds_edge_vec::Vector{Int64},
    bonds_connected_to_vertex_dict::Dict{Int64, Vector{Tuple{Int64, Int64}}},
    bond_vertex_neighbor_dict::Dict{Tuple{Int64, Int64}, 
        Tuple{Vector{Int64}, Vector{Int64}}},
    bond_l_dict::Dict{Tuple{Int64, Int64}, Int64},
    cluster_dict::Dict,
    spatial_network::MetaGraphsNext.MetaGraph) where T 

    all_cluster_vertices_vec = vcat(
        cluster_dict["cluster_vertices_to_move_vec"], 
        cluster_dict["cluster_vertices_outer_shell_vec"]) 
    all_cluster_bonds_vec = vcat(
        cluster_dict["cluster_bonds_inside_vec"], 
        cluster_dict["cluster_bonds_edge_vec"])
    
    # create a dictionary for the vertex coordinates 
    vertex_coord_dict = Dict{Int64, StaticArrays.SVector{3, T}}()
    for (index, vertex) in enumerate(
            cluster_dict["cluster_vertices_to_move_vec"])
        idx = (index - 1) * 3
        vertex_coord_dict[vertex] = StaticArrays.SVector{3, T}(
            vertices_to_move_coord_vec[idx + 1],
            vertices_to_move_coord_vec[idx + 2],
            vertices_to_move_coord_vec[idx + 3]
        )
    end
    for (index, vertex) in enumerate(
            cluster_dict["cluster_vertices_outer_shell_vec"])
        idx = (index - 1) * 3
        vertex_coord_dict[vertex] = StaticArrays.SVector{3, T}(
            vertices_outer_shell_coord_vec[idx + 1],
            vertices_outer_shell_coord_vec[idx + 2],
            vertices_outer_shell_coord_vec[idx + 3]
        )
    end
    for (index, vertex) in enumerate(vertices_bonds_edge_vec)
        idx = (index - 1) * 3
        vertex_coord_dict[vertex] = StaticArrays.SVector{3, T}(
            vertices_bonds_edge_coord_vec[idx + 1],
            vertices_bonds_edge_coord_vec[idx + 2],
            vertices_bonds_edge_coord_vec[idx + 3]
        )
    end

    # create a dictionary for the bond vectors and their squared distances to 
    # avoid redundant calculations in the bending and torsional energy 
    # calculations
    bond_vector_dict = Dict{Tuple{Int64, Int64}, StaticArrays.SVector{3, T}}()
    bond_distance_squared_dict = Dict{Tuple{Int64, Int64}, T}()
    
    for bond in all_cluster_bonds_vec
        bond_vector = StaticArrays.SVector{3, T}(
            get_distance_vector_pbc(
                vertex_coord_dict[bond[1]],
                vertex_coord_dict[bond[2]],
                spatial_network[]["supercell_edge_length"]
            )
        )
        bond_vector_dict[bond] = bond_vector
        bond_distance_squared_dict[bond] = LinearAlgebra.norm(bond_vector)^2
    end

    # calculate bond bending energy
    bond_bending_sum = zero(T)

    for vertex in all_cluster_vertices_vec
        vertex_bonds_vec = bonds_connected_to_vertex_dict[vertex] 
        
        for (index1, bond1) in enumerate(vertex_bonds_vec)
            sign1 = bond1[2] == vertex ? -1 : 1
            for index2 in index1+1:length(vertex_bonds_vec)
                bond2 = vertex_bonds_vec[index2]
                sign2 = (bond2[2] == vertex ? -1 : 1) * sign1
                
                bond_bending_sum += (sign2 * LinearAlgebra.dot(
                        bond_vector_dict[bond1], bond_vector_dict[bond2]) 
                    - cosd(spatial_network[]["theta_ground_state"]) )^2 
            end
        end
    end

    bond_bending_energy = (
        3/8 * spatial_network[]["bond_bending_const"] * bond_bending_sum)

    # determine bond stretching energy
    bond_stretching_energy = 0.0
    for bond in cluster_dict["cluster_bonds_inside_vec"]
        bond_stretching_energy += (
            3/16 * (bond_distance_squared_dict[bond] - 1)^2 )
    end

    # bonds on the edge of the cluster only contribute half to the energy
    for bond in cluster_dict["cluster_bonds_edge_vec"]
        bond_stretching_energy += (
            3/32 * (bond_distance_squared_dict[bond] - 1)^2 )
    end

    cluster_energy = bond_bending_energy + bond_stretching_energy
    
    # torsional energy if required
    if (haskey(spatial_network[], "torsional_const") 
            && spatial_network[]["torsional_const"] > 0.0)
        torsional_sum = zero(T)
        delta_phi_rad = spatial_network[]["delta_phi"] * pi / 180
        
        for bond in cluster_dict["cluster_bonds_inside_vec"]
            bond_vector = bond_vector_dict[bond]
            bond_vector_normalized = (bond_vector 
                / sqrt(bond_distance_squared_dict[bond]))

            vertex_1_coordination_nr = spatial_network[bond[1]][
                "coordination_nr"]
            vertex_2_coordination_nr = spatial_network[bond[2]][
                "coordination_nr"]

            neighbor_label_vec_1 = bond_vertex_neighbor_dict[bond][1]
            neighbor_label_vec_2 = bond_vertex_neighbor_dict[bond][2]
            l = bond_l_dict[bond]

            @inbounds @simd for j in 1:(vertex_1_coordination_nr-1)
                neighbor_label_j = neighbor_label_vec_1[j]
                sign1::Int64=sign(bond[1] - neighbor_label_j)
                bond_j_key = (min(bond[1], neighbor_label_j), 
                    max(bond[1], neighbor_label_j))
                vector_j = (
                    sign1 .* bond_vector_dict[bond_j_key])

                # first normal vector
                n1 = LinearAlgebra.cross(vector_j, bond_vector)

                @inbounds @simd for k in 1:(vertex_2_coordination_nr-1)
                    neighbor_label_k = neighbor_label_vec_2[k]
                    sign2::Int64=sign(neighbor_label_k - bond[2])
                    bond_k_key = (min(bond[2], neighbor_label_k), 
                         max(bond[2], neighbor_label_k))
                    vector_k = (
                        sign2 .* bond_vector_dict[bond_k_key])

                    # second normal vector
                    n2 = LinearAlgebra.cross(bond_vector, vector_k)

                    # compute dihedral angle phi robustly
                    x = LinearAlgebra.dot(n1, n2)
                    y = LinearAlgebra.dot(bond_vector_normalized, 
                        LinearAlgebra.cross(n1, n2))
                    phi = atan(y, x)

                    # torsional energy
                    torsional_sum += 1 - cos(l * (phi - delta_phi_rad))
                end
            end
        end

        torsional_energy = (
            1/2 * spatial_network[]["torsional_const"] * torsional_sum)

        cluster_energy += torsional_energy
    end
    
    return cluster_energy
end
