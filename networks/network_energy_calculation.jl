"""
These functions are used to calculate energies in spatial networks
"""


"""
Get bond stretching energy for a given bond
"""
function local_bond_stretching_energy_keating(graph_dict::Dict, bond::Tuple{Int64, Int64})

    #get bond stretching energy
    bond_stretching_energy = (3/16 * ( 
            graph_dict["spatial_network"][bond...]["distance_squared"] - 1 
                                            )^2 ) 
    
    return bond_stretching_energy

end


"""
Get bond bending energy for a given vertex
"""
function local_bond_bending_energy_keating(graph_dict::Dict, vertex_label::Int64)

    #get vector of neighbor labels
    neighbor_label_vec = collect(MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], 
                                vertex_label))

    #initialize bond bending sum
    bond_bending_sum = 0

    #loop through all bond combinations
    for j in 1:graph_dict["coordination_nr"]

        for k in j+1:graph_dict["coordination_nr"]

            bond_bending_sum += (  LinearAlgebra.dot( sign(neighbor_label_vec[j] - vertex_label) .* 
                        graph_dict["spatial_network"][vertex_label, neighbor_label_vec[j]]["vector"], 
                        sign(neighbor_label_vec[k] - vertex_label) .* 
                        graph_dict["spatial_network"][vertex_label, neighbor_label_vec[k]]["vector"]
                                 ) + 1/3 )^2 
            
        end

    end

    #calculate bond bending energy
    bond_bending_energy =  3/8 * graph_dict["bond_bending_const"] * bond_bending_sum

    return bond_bending_energy

end


"""
Calculate the total energy of a spatial network
"""
function get_total_energy_keating(graph_dict::Dict)

    total_energy = 0

    #loop through all vertices and sum bond bending energies
    for vertex in MetaGraphsNext.labels(graph_dict["spatial_network"])

        total_energy += local_bond_bending_energy_keating(graph_dict, vertex)
    end

    #loop through all bonds and sum bond stretching energies
    for bond in MetaGraphsNext.edge_labels(graph_dict["spatial_network"])

        total_energy += local_bond_stretching_energy_keating(graph_dict, bond)
    end

    return total_energy
end


"""
Get energy for a cluster of vertices whose vertices and bonds
are stored in the respective dictionary
"""
function get_cluster_energy(graph_dict, cluster_dict)

    #initialize cluster energy
    cluster_energy = 0

    #get vector of all cluster vertices
    all_cluster_vertices_vec = vcat(cluster_dict["cluster_vertices_to_move_vec"], 
                            cluster_dict["cluster_vertices_outer_shell_vec"]) 

    #loop through all vertices and sum bond bending energies
    for vertex in all_cluster_vertices_vec

        cluster_energy += local_bond_bending_energy_keating(graph_dict, vertex)
    end

    #loop through all bonds inside the cluster
    for bond in cluster_dict["cluster_bonds_inside_vec"]

        cluster_energy += local_bond_stretching_energy_keating(
                                                            graph_dict, bond)
    end

    #loop through all bonds on the edge of the cluster which only
    #contributy half
    #loop through all bonds inside the cluster
    for bond in cluster_dict["cluster_bonds_edge_vec"]

        cluster_energy += 1/2 * local_bond_stretching_energy_keating(
                                                            graph_dict, bond)
    end

    return cluster_energy

end


"""
Calculate the local Keating energy for a given vertex
"""
function local_energy_keating(vertex_label::Int64, graph_dict::Dict)

    #get bond bending energy
    local_energy = local_bond_bending_energy_keating(graph_dict, 
                                                                vertex_label)

    #sum bond stretching energy contributions by considering that each bond
    #is shared by two vertices
    for neighbor in MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"],
                                                        vertex_label)

        local_energy += 1/2 * local_bond_stretching_energy_keating(graph_dict, 
                                                        (vertex_label, neighbor))

    end
    
    return local_energy

end


"""
Calculate the contibution of a single vertex position x to the total
Keating energy from the position of its neighbors and next to nearest neighbors
"""
function energy_from_position_keating(x::Vector, 
                        graph_dict::Dict,
                        neighbor_positions_mat::Matrix{Float64},
                        next_neighbor_positions_arr::Array{Float64})

    local_energy = 0

    for j in 1:graph_dict["coordination_nr"]

        #get vector pointing from central vertex to neighbor
        distance_vector_j = neighbor_positions_mat[:,j] .- x

        #get bond stretching term
        bond_stretching_term = (3/16 * ( LinearAlgebra.norm(distance_vector_j)^2 - 1 )^2 ) 

        #get bond bending term
        bond_bending_sum = 0

        for k in j+1:graph_dict["coordination_nr"]

            bond_bending_sum += ( 3/8 * graph_dict["bond_bending_const"]  
                * ( LinearAlgebra.dot( distance_vector_j, 
                                    (neighbor_positions_mat[:,k] .- x) ) + 1/3 )^2 )
            
        end

        #get bond bending terms due to next to nearest neighbors
        neighbor_bond_bending_sum = 0

        for l in 1:graph_dict["coordination_nr"]-1

            neighbor_bond_bending_sum += ( 3/8 * graph_dict["bond_bending_const"]  
                * ( LinearAlgebra.dot( distance_vector_j, 
                    (neighbor_positions_mat[:,j] .- next_neighbor_positions_arr[j,:,l]) ) 
                + 1/3 )^2 )
            
        end

        #sum bond stretching and bending terms
        local_energy += bond_stretching_term + bond_bending_sum + neighbor_bond_bending_sum

    end
    
    return local_energy

end


"""
Calculate the negative Keating force (-f) on a given vertex
which is the gradient of its local Keating energy from
the distances and vectors that have already been calculated
"""
function gradient_keating_efficient(graph_dict::Dict, central_vertex::Int64)
    
    #initialize gradient
    gradient = zeros(graph_dict["nr_dimensions"])

    #get vector of neighbors
    neighbor_vec = collect(
        MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], central_vertex))

    #loop through central vertex's neighbors
    for j in 1:graph_dict["coordination_nr"]

        #get vector pointing from central vertex to neighbor j
        distance_vector_j = (sign(neighbor_vec[j] - central_vertex)
            * graph_dict["spatial_network"][central_vertex, neighbor_vec[j]]["vector"])

        #get bond stretching term
        bond_stretching_term = ( - 3/4 * ( 
            graph_dict["spatial_network"][central_vertex, neighbor_vec[j]]["distance_squared"] - 1 
            ) ) .* distance_vector_j

        #get bond bending term
        bond_bending_sum = zeros(graph_dict["nr_dimensions"])

        for k in j+1:graph_dict["coordination_nr"]

            #get vector pointing from central vertex to neighbor k
            distance_vector_k = (sign(neighbor_vec[k] - central_vertex)
                * graph_dict["spatial_network"][central_vertex, neighbor_vec[k]]["vector"])

            bond_bending_sum .-= ( ( 3/4 * graph_dict["bond_bending_const"]  
                * ( LinearAlgebra.dot( distance_vector_j, distance_vector_k ) + 1/3 ) )
            .* ( distance_vector_j .+ distance_vector_k )) 
            
        end

        #get bond bending terms due to next to nearest neighbors
        neighbor_bond_bending_sum = zeros(graph_dict["nr_dimensions"])

        #get neighbors of current neighbor excluding the central vertex
        neighbors_neighbor_vec = setdiff( collect(
            MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], neighbor_vec[j])
                ), central_vertex )

        for l in 1:graph_dict["coordination_nr"]-1

            #get vector pointing from next-to-nearest-neighbor l to neighbor j
            distance_vector_l = (sign(neighbor_vec[j] - neighbors_neighbor_vec[l])
            * graph_dict["spatial_network"][neighbors_neighbor_vec[l], neighbor_vec[j]]["vector"])

            neighbor_bond_bending_sum .-= ( 3/4 * graph_dict["bond_bending_const"]  
                * ( LinearAlgebra.dot( distance_vector_j, distance_vector_l )  + 1/3 ) 
                .* distance_vector_l )
            
        end

        #sum bond stretching and bending terms
        gradient .+= bond_stretching_term .+ bond_bending_sum .+ neighbor_bond_bending_sum

    end

    return gradient
end



"""
Calculate the Hessian matrix for a given vertex from
    the distances and vectors that have already been calculated
"""
function hessian_keating_efficient(graph_dict::Dict, central_vertex::Int64)

    #initialize hessian
    hessian = zeros(graph_dict["nr_dimensions"],graph_dict["nr_dimensions"])
    
    #get vector of neighbors
    neighbor_vec = collect(
        MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], central_vertex))

    for a in 1:graph_dict["nr_dimensions"]

        #since the hessian is symmetric, it only needs to be calculated for the lower left half
        for b in a:graph_dict["nr_dimensions"]

            for j in 1:graph_dict["coordination_nr"]

                #get vector pointing from central vertex to neighbor j
                distance_vector_j = (sign(neighbor_vec[j] - central_vertex)
                * graph_dict["spatial_network"][central_vertex, neighbor_vec[j]]["vector"])
        
                #get bond stretching term
                bond_stretching_term = (3/2 * distance_vector_j[b] * distance_vector_j[a]
                    + ==(a,b) * 3/4 * (
                    graph_dict["spatial_network"][central_vertex, neighbor_vec[j]]["distance_squared"] 
                    - 1 ) )
        
                #get bond bending term 
                bond_bending_sum = 0
        
                for k in j+1:graph_dict["coordination_nr"]

                    #get vector pointing from central vertex to neighbor k
                    distance_vector_k = (sign(neighbor_vec[k] - central_vertex)
                    * graph_dict["spatial_network"][central_vertex, neighbor_vec[k]]["vector"])
        
                    bond_bending_sum += (  3/4 * graph_dict["bond_bending_const"]  
                        *( (distance_vector_j .+ distance_vector_k)[b]
                            *(distance_vector_j .+ distance_vector_k)[a]
                            + ==(a,b) * 2 *( LinearAlgebra.dot( distance_vector_j, distance_vector_k ) 
                                        + 1/3 ) )
                        )
                    
                end

                #get bond bending terms due to next to nearest neighbors
                neighbor_bond_bending_sum = 0

                #get neighbors of current neighbor excluding the central vertex
                neighbors_neighbor_vec = setdiff( collect(
                    MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"], neighbor_vec[j])
                        ), central_vertex )
        
                for l in 1:graph_dict["coordination_nr"]-1

                    #get vector pointing from next-to-nearest-neighbor l to neighbor j
                    distance_vector_l = (sign(neighbor_vec[j] - neighbors_neighbor_vec[l])
                    * graph_dict["spatial_network"][neighbors_neighbor_vec[l], neighbor_vec[j]]["vector"])
        
                    neighbor_bond_bending_sum += (  3/4 * graph_dict["bond_bending_const"]  
                        * distance_vector_l[b] * distance_vector_l[a] 
                        )
                    
                end
        
                #sum bond stretching and bending terms
                hessian[a,b] += bond_stretching_term + bond_bending_sum + neighbor_bond_bending_sum
        
            end
            
        end

        #use symmetry of hessian and copy entries from lower left to upper right half
        for b in 1:a-1
            hessian[a,b] = hessian[b,a]
        end
    end

    return LinearAlgebra.Symmetric(hessian)
end



"""
Calculate the negative Keating force (-f) on a given vertex
which is the gradient of its local Keating energy.
x is the position vector of the vertex
"""
function gradient_keating!(gradient::Vector, x::Vector{Float64}, 
                        graph_dict::Dict,
                        neighbor_positions_mat::Matrix{Float64},
                        next_neighbor_positions_arr::Array{Float64})

    #reset gradient (for some reason the optimization does not
    #work if this is done vectorized)
    gradient[:] = zeros(graph_dict["nr_dimensions"])

    for j in 1:graph_dict["coordination_nr"]

        #get vector pointing from central vertex to neighbor
        distance_vector_j = neighbor_positions_mat[:,j] .- x

        #get bond stretching term
        bond_stretching_term = ( - 3/4 * ( LinearAlgebra.norm(distance_vector_j)^2 - 1 ) 
                                ) .* distance_vector_j

        #get bond bending term
        bond_bending_sum = zeros(graph_dict["nr_dimensions"])

        for k in j+1:graph_dict["coordination_nr"]

            bond_bending_sum .+= ( ( 3/4 * graph_dict["bond_bending_const"]  
                * ( LinearAlgebra.dot( distance_vector_j, 
                                    neighbor_positions_mat[:,k] .- x ) + 1/3 ) )
            .* ( 2 .* x .- ( neighbor_positions_mat[:,j] .+ neighbor_positions_mat[:,k] ) )
            )
            
        end

        #get bond bending terms due to next to nearest neighbors
        neighbor_bond_bending_sum = zeros(graph_dict["nr_dimensions"])

        for l in 1:graph_dict["coordination_nr"]-1

            neighbor_bond_bending_sum .+= ( 3/4 * graph_dict["bond_bending_const"]  
                * ( LinearAlgebra.dot( distance_vector_j, 
                    (neighbor_positions_mat[:,j] .- next_neighbor_positions_arr[j,:,l]) ) 
                + 1/3 )
            .* (next_neighbor_positions_arr[j,:,l] .- neighbor_positions_mat[:,j]) )
            
        end

        #sum bond stretching and bending terms
        gradient .+= bond_stretching_term .+ bond_bending_sum .+ neighbor_bond_bending_sum

    end

end
    

"""
Calculate the Hessian matrix for a given vertex and Keating energy
which is the a matrix of second derivatives of its energy
"""
function hessian_keating!(hessian::Matrix{Float64}, x::Vector{Float64}, 
                        graph_dict::Dict, 
                        neighbor_positions_mat::Matrix{Float64},
                        next_neighbor_positions_arr::Array{Float64})

    #get hessian
    hessian[:,:] = zeros(graph_dict["nr_dimensions"],graph_dict["nr_dimensions"])

    for a in 1:graph_dict["nr_dimensions"]

        #since the hessian is symmetric, it only needs to be calculated for the lower left half
        for b in a:graph_dict["nr_dimensions"]

            for j in 1:graph_dict["coordination_nr"]
        
                #get bond stretching term
                bond_stretching_term = (3/2 * (x[b]-neighbor_positions_mat[b,j] ) 
                                            * (x[a]-neighbor_positions_mat[a,j] )
                        + ==(a,b) * 3/4 * ( LinearAlgebra.norm( 
                                neighbor_positions_mat[:,j] .- x )^2 - 1 ) )
        
                #get bond bending term 
                bond_bending_sum = 0
        
                for k in j+1:graph_dict["coordination_nr"]
        
                    bond_bending_sum += (  3/4 * graph_dict["bond_bending_const"]  
                        *( (2*x[b] - ( neighbor_positions_mat[b,j] + neighbor_positions_mat[b,k] ))
                            *(2*x[a] - ( neighbor_positions_mat[a,j] + neighbor_positions_mat[a,k] ))
                            + ==(a,b) * 2 *( LinearAlgebra.dot( neighbor_positions_mat[:,j] .- x, 
                                                            neighbor_positions_mat[:,k] .- x ) 
                                        + 1/3 ) )
                        )
                    
                end

                #get bond bending terms due to next to nearest neighbors
                neighbor_bond_bending_sum = 0
        
                for l in 1:graph_dict["coordination_nr"]-1
        
                    neighbor_bond_bending_sum += (  3/4 * graph_dict["bond_bending_const"]  
                        *( next_neighbor_positions_arr[j,b,l] - neighbor_positions_mat[b,j] )
                        *( next_neighbor_positions_arr[j,a,l] - neighbor_positions_mat[a,j] )
                        )
                    
                end
        
                #sum bond stretching and bending terms
                hessian[a,b] += bond_stretching_term + bond_bending_sum + neighbor_bond_bending_sum
        
            end
            
        end

        #use symmetry of hessian and copy entries from lower left to upper right half
        for b in 1:a-1
            hessian[a,b] = hessian[b,a]
        end
    end
end


"""
Get relaxation weight of a single vertex, which contributes to the
Metropolis acceptance probability when thermal fluctuations are
included
"""
function get_vertex_relaxation_weight(excited_graph_dict::Dict, 
    relaxed_graph_dict::Dict, vertex::Int64, temperature::Real)

    #get displacement vector 
    displacement_vec = (excited_graph_dict["spatial_network"][vertex]["position"]
                    .- relaxed_graph_dict["spatial_network"][vertex]["position"])

    #get hessian matrix of current vertex at relaxed position
    hessian = hessian_keating_efficient(relaxed_graph_dict, vertex)

    #get vector of sigmas
    sigma_vec = sqrt.(temperature ./ LinearAlgebra.diag(hessian) ) 

    #get vector of gauss functions
    gauss_fct_vec = ( (1/sqrt(2*pi)) ./ sigma_vec
                    .* exp.( .- displacement_vec.^2 ./ (2 .* (sigma_vec .^2 ) ) ) )

    #get relaxation weight which is the product of the gauss fct for
    #all dimensions
    vertex_relaxation_weight = prod(gauss_fct_vec)

    return vertex_relaxation_weight
end


"""
Get total weight corresponding to thermal fluctuations when relaxing or 
exciting a cluster. This weight contributes to the Metropolis acceptance
probability when thermal fluctuations are included
"""
function get_cluster_fluctuation_weight(excited_graph_dict::Dict, 
                                    relaxed_graph_dict::Dict, 
                                    cluster_dict::Dict,
                                    temperature::Real)

    #initialize cluster relaxation weight
    cluster_relaxation_weight = 1

    #loop through all vertices that are allowed to move
    for vertex in cluster_dict["cluster_vertices_to_move_vec"]

        #get relaxation weight for current vertex
        vertex_relaxation_weight = get_vertex_relaxation_weight(excited_graph_dict, 
                                                        relaxed_graph_dict, vertex, temperature)

        cluster_relaxation_weight *= vertex_relaxation_weight
    end

    return cluster_relaxation_weight
end

