"""
These functions are used in different processes when working with
the network graphs
"""


"""
Calculate vector pointing from position a to position b 
under periodic boundary conditions
"""
function get_distance_vector_pbc(position_a::Vector,
                                position_b::Vector,
                                supercell_edge_length::Real )

    #get vector pointing from position a to b without considering boundary conditions
    distance_vector_without_pbc = (position_b .- position_a)

    #modify the vector according to boundary conditions
    distance_vector = ( ( abs.(distance_vector_without_pbc) 
                    .< (supercell_edge_length/2) ) 
                .* distance_vector_without_pbc
            .+ ( abs.(distance_vector_without_pbc .+ supercell_edge_length) 
                    .< (supercell_edge_length/2) ) 
                .* (distance_vector_without_pbc .+ supercell_edge_length)
            .+ ( abs.(distance_vector_without_pbc .- supercell_edge_length) 
                    .< (supercell_edge_length/2) ) 
                .* (distance_vector_without_pbc .- supercell_edge_length)
            )    

    return distance_vector
    
end


"""
Calculate virtual position of an atom relative to a central atom
by placing it outside of the supercell if periodic boundary conditions
have to be taken into account
"""
function get_virtual_position(central_atom_position::Vector{Float64},
                                other_atom_position::Vector{Float64},
                                supercell_edge_length::Real )

    #get vector pointing from central atom to neighbor without considering 
    #boundary conditions
    distance_vector_without_pbc = (other_atom_position .- central_atom_position)

    #get the other atoms virtual position according to boundary conditions
    virtual_other_atom_position = ( ( abs.(distance_vector_without_pbc) 
                    .< (supercell_edge_length/2) ) 
                .* other_atom_position
            .+ ( abs.(distance_vector_without_pbc .+ supercell_edge_length) 
                    .< (supercell_edge_length/2) ) 
                .* (other_atom_position .+ supercell_edge_length)
            .+ ( abs.(distance_vector_without_pbc .- supercell_edge_length) 
                    .< (supercell_edge_length/2) ) 
                .* (other_atom_position .- supercell_edge_length)
            )    

    return virtual_other_atom_position
    
end


"""
Get matrix with the coordinates of an atoms neighbors
by taking periodic boundary conditions into account
"""
function get_neighbor_positions_mat(graph_dict::Dict, central_atom::Int64;
                                    exclude_atoms::Vector = [])

    #get central atom's position
    central_atom_position = graph_dict["spatial_network"][central_atom]["position"]

    #create matrix to store neighbors coordinates
    neighbor_positions_mat = Matrix{Float64}(undef, 
                            graph_dict["nr_dimensions"],
                            graph_dict["coordination_nr"]-length(exclude_atoms))
    
    #save coordinates to matrix and array
    current_neighbor = 1

    for neighbor in MetaGraphsNext.neighbor_labels(
        graph_dict["spatial_network"], central_atom)

        if !in(neighbor, exclude_atoms)

            #get neighbor's virtual coordinates which might be outside of the 
            #supercell if periodic boundary conditions play a role
            neighbor_positions_mat[:,current_neighbor] = get_virtual_position(
                            central_atom_position,
                            graph_dict["spatial_network"][neighbor]["position"],
                            graph_dict["supercell_edge_length"] )

            current_neighbor += 1
        end

    end

    return neighbor_positions_mat

end


"""
Get array with the coordinates of an atoms next to nearest neighbors
by taking periodic boundary conditions into account
"""
function get_next_neighbor_positions_arr(graph_dict::Dict, central_atom::Int64)

    #get central atom's position
    central_atom_position = graph_dict["spatial_network"][central_atom]["position"]

    #get central atoms neighbors 
    neighbor_vec = collect(MetaGraphsNext.neighbor_labels(
                                graph_dict["spatial_network"], central_atom))

    #create array to store next to nearest neighbors coordinates
    #The first array index labels the number of the direct neighbor
    next_neighbor_positions_arr = Array{Float64}(undef, 
                                                graph_dict["coordination_nr"],
                                                graph_dict["nr_dimensions"],
                                                graph_dict["coordination_nr"]-1)
    
    #loop through central atoms neighbors
    for i in 1:graph_dict["coordination_nr"]

        current_next_neighbor = 1

        #loop through the current neighbor's neighbors
        for next_neighbor in MetaGraphsNext.neighbor_labels(
                                        graph_dict["spatial_network"], neighbor_vec[i])

            if next_neighbor !== central_atom

                #get next neighbor's virtual coordinates which might be outside of the 
                #supercell if periodic boundary conditions play a role
                next_neighbor_positions_arr[i,:,current_next_neighbor] = get_virtual_position(
                            central_atom_position,
                            graph_dict["spatial_network"][next_neighbor]["position"],
                            graph_dict["supercell_edge_length"] )

                current_next_neighbor += 1
            end

        end

    end

    return next_neighbor_positions_arr

end



"""
get all bonds inside and on the edge of cluster
"""
function get_cluster_bonds_vec(graph_dict::Dict,
                                cluster_atoms_to_move_vec::Vector{Int64},
                                cluster_atoms_outer_shell_vec::Vector{Int64})

    #initialize vectors for bonds
    cluster_bonds_inside_vec = []
    cluster_bonds_edge_vec = []

    #get vector of all cluster atoms
    all_cluster_atoms_vec = vcat(cluster_atoms_to_move_vec, 
                                cluster_atoms_outer_shell_vec)

    #loop through all cluster atoms
    for cluster_atom in all_cluster_atoms_vec

        #loop through all neighbors of current atom
        for neighbor in MetaGraphsNext.neighbor_labels(graph_dict["spatial_network"],
                                                            cluster_atom)

            #add current bond to the respective vector if the it is not stored yet
            if neighbor in all_cluster_atoms_vec

                if cluster_atom < neighbor
                    push!(cluster_bonds_inside_vec, (cluster_atom, neighbor) )
                end

            else
                push!(cluster_bonds_edge_vec, Tuple(sort([cluster_atom, neighbor])) )
            
            end
    
        end

    end

    return [cluster_bonds_inside_vec, cluster_bonds_edge_vec]
end



"""
Get a dictionary of vectors containing all atoms neighboring
the central atoms up to the given shell
"""
function get_cluster_in_shells_dict(graph_dict::Dict, 
                                    central_atoms::Tuple; 
                                    shell_nr::Int64 = 5)

    #initialize dictionary for all neighbors sorted by shells
    cluster_in_shells_dict = Dict(0 => copy(collect(central_atoms)) )

    #initialize vector for cluster atoms which will be allowed to move
    cluster_atoms_to_move_vec = collect(central_atoms)

    #initialize vector for atoms in the outer shell which will not be
    #allowed to move
    cluster_atoms_outer_shell_vec = Vector{Int64}()

    #loop through neighbor shells
    for current_shell in 1:shell_nr

        #initialize vector of atoms in current shell
        current_shell_atoms_vec = Vector{Int64}()

        #loop through atoms of lower shell
        for lower_shell_atom in cluster_in_shells_dict[current_shell-1]

            #loop through neighbors of current atom
            for neighbor in MetaGraphsNext.neighbor_labels(
                            graph_dict["spatial_network"], lower_shell_atom)

                #if not in outer shell, save as atom to move
                if current_shell < shell_nr

                    #save current atom if it was not considered before
                    if !(neighbor in cluster_atoms_to_move_vec)

                        push!(current_shell_atoms_vec, neighbor)
                        push!(cluster_atoms_to_move_vec, neighbor)
                    end

                #if in outer shell, the atom is not allowed to move
                else
                    if (!(neighbor in cluster_atoms_to_move_vec) 
                            && !(neighbor in cluster_atoms_outer_shell_vec))

                        push!(current_shell_atoms_vec, neighbor)
                        push!(cluster_atoms_outer_shell_vec, neighbor)
                    end
                end

            end

        end

        #save to cluster in shells dict
        cluster_in_shells_dict[current_shell] = current_shell_atoms_vec

    end

    #get all bonds inside and on the edge of cluster
    cluster_bonds_inside_vec, cluster_bonds_edge_vec = get_cluster_bonds_vec(
                                                graph_dict,
                                                cluster_atoms_to_move_vec,
                                                cluster_atoms_outer_shell_vec
                                                )

    #collect cluster information into dictionary
    cluster_dict = Dict("cluster_in_shells_dict" => cluster_in_shells_dict, 
            "cluster_atoms_to_move_vec" => cluster_atoms_to_move_vec, 
            "cluster_atoms_outer_shell_vec" => cluster_atoms_outer_shell_vec,
            "cluster_bonds_inside_vec" => cluster_bonds_inside_vec, 
            "cluster_bonds_edge_vec" => cluster_bonds_edge_vec
            )

    #add cluster energy to dictionary
    cluster_dict["cluster_energy"] = get_cluster_energy(graph_dict, cluster_dict)
    cluster_dict["cluster_energy_up_to_date"] = true

    return cluster_dict

end


"""
Pick a random bond that has not been declined since the last accepted move
"""
function get_random_bond(graph_dict::Dict; declined_bonds = [], seed = Nothing)

    #set seed if desired
    if seed !== Nothing
        Random.seed!(seed)
    end

    #determine nr of bonds
    nr_bonds = graph_dict["nr_atoms"] * graph_dict["coordination_nr"] / 2 

    #check if all bonds have been attempted already
    if length(declined_bonds) == nr_bonds
            @warn "All bonds have been attempted without success"
        random_bond = []

    #if the list of declined bonds is already very long
    #pick one of the remaining ones
    elseif length(declined_bonds) > nr_bonds/2
        all_bonds_vec = collect(
                MetaGraphsNext.edge_labels(graph_dict["spatial_network"]))

        random_bond = rand(all_bonds_vec)

    #otherwise get random bond without listing all bonds
    else

        #pick a random atom
        atom_1 = rand(1:graph_dict["nr_atoms"])

        #pick a random neighbor
        atom_2 = collect(MetaGraphsNext.neighbor_labels(
                            graph_dict["spatial_network"], atom_1)
                            )[rand(1:graph_dict["coordination_nr"])]

        #create bond
        random_bond = Tuple(sort([atom_1, atom_2]))

        #find new bond if current one was already declined
        if random_bond in declined_bonds
            random_bond = get_random_bond(graph_dict; declined_bonds = declined_bonds)
        end
    end

    return random_bond
end


"""
Check whether all vertices in network have the correct coordination number
"""
function get_incorrectly_coordinated_atoms(graph_dict::Dict)

    incorrectly_coordinated_atoms = []

    #for each atom, check whether it has the correct coordination nr
    for atom in MetaGraphsNext.labels(graph_dict["spatial_network"])
        if (length( collect( MetaGraphsNext.neighbor_labels(
                                graph_dict["spatial_network"], atom) ) ) 
            !== graph_dict["coordination_nr"])

            push!(incorrectly_coordinated_atoms, atom)
        end
    end

    return incorrectly_coordinated_atoms
    
end


"""
This function relaxes a given cluster in two ways, first using the
Newton method which is slower but more accurate, and then more efficiently
but less accurate using the method from 10.1142/S0217984987000065. The two
relaxation methods are compared by plotting the evolution of atomic positions
and cluster energies
"""
function compare_relaxation_methods(original_graph_dict,
    central_cluster_atoms,
    filename;
    nr_max_relaxation_cycles = 25,
    shell_nr::Int64 = 4,
    relaxation_overshoot_factor_r::Real = 1.5,
    relaxation_optimization_parameter_l::Real = 1,
    save_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\plots\random_networks\\" )

    #initialize arrays for atomic positions and cluster energy as a
    #function of relaxation cycle
    atomic_position_arr = Array{Float64}(undef, nr_max_relaxation_cycles+1, 3, 2)
    cluster_energy_arr = Array{Float64}(undef, nr_max_relaxation_cycles+1, 2)

    #loop through relaxation methods
    relax_efficiently_vec = [false, true]

    for i in 1:2

        #reset graph dict to original one
        graph_dict = deepcopy(original_graph_dict)

        #get the cluster dict
        cluster_dict = get_cluster_in_shells_dict(
                                            graph_dict, 
                                            central_cluster_atoms; 
                                            shell_nr = shell_nr)

        #store initial atomic position and cluster energy
        atomic_position_arr[1,:,i] = graph_dict["spatial_network"][central_cluster_atoms[1]]["position"]
        cluster_energy_arr[1,i] = cluster_dict["cluster_energy"]

        #perform relaxation cycles
        for relaxation_cycle in 2:nr_max_relaxation_cycles+1
            graph_dict, cluster_dict = relax_cluster_one_cycle_keating!(graph_dict, 
            cluster_dict; 
            relax_efficiently = relax_efficiently_vec[i],
            relaxation_overshoot_factor_r = relaxation_overshoot_factor_r,
            relaxation_optimization_parameter_l = relaxation_optimization_parameter_l,
            update_cluster_energy = true )

            #keep track of atomic position and cluster energy
            atomic_position_arr[relaxation_cycle,:,i] = graph_dict["spatial_network"][central_cluster_atoms[1]]["position"]
            cluster_energy_arr[relaxation_cycle,i] = cluster_dict["cluster_energy"]

        end
    end

    #plot evolution of atomic position and cluster energy
    Plots.plot(xlabel= Latex.L"x",
    ylabel=Latex.L"y" ,
    legend = true, dpi=250)
    Plots.plot!(atomic_position_arr[:,1,1], atomic_position_arr[:,2,1], label = "Newton", markershape = :auto)
    Plots.plot!(atomic_position_arr[:,1,2], atomic_position_arr[:,2,2], label = "efficient", markershape = :auto)

    Plots.savefig(save_path*filename*"_x_y_pos.png")

    Plots.plot(xlabel= Latex.L"x",
    ylabel=Latex.L"z" ,
    legend = true, dpi=250)
    Plots.plot!(atomic_position_arr[:,1,1], atomic_position_arr[:,3,1], label = "Newton", markershape = :auto)
    Plots.plot!(atomic_position_arr[:,1,2], atomic_position_arr[:,3,2], label = "efficient", markershape = :auto)

    Plots.savefig(save_path*filename*"_x_z_pos.png")

    Plots.plot(xlabel= "relaxation cycle",
    ylabel="cluster energy" ,
    legend = true, dpi=250)
    Plots.plot!(collect(0:nr_max_relaxation_cycles), cluster_energy_arr[:,1], label = "Newton")
    Plots.plot!(collect(0:nr_max_relaxation_cycles), cluster_energy_arr[:,2], label = "efficient")

    Plots.savefig(save_path*filename*"_cluster_energy.png")

    return [atomic_position_arr, cluster_energy_arr]
end


"""
Convert cartesian to spherical coordinates
"""
function convert_cartesian_to_spherical(cartesian_vec::Vector)

    #check if vector is 3d
    if length(cartesian_vec) !== 3
        @error "conversion to spherical coordinates only implemented in 3d"
        return []
    end

    #calculate r, theta and phi 
    r_length = LinearAlgebra.norm(cartesian_vec)
    theta = acos(cartesian_vec[3]/r_length)
    phi = acos(cartesian_vec[1]/(r_length*sin(theta)) )

    return [r_length, theta, phi]
end
    