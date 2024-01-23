"""
These functions generate geometric graphs
"""


"""
generate matrix of vertex positions in a simple cubic lattice, 
where each column is a position vector
"""
function get_simple_cubic_vertex_position_mat(nr_vertices_per_dimension::Int64, 
                                            nr_dimensions::Int64)

    #generate empty matrix for vertexic positions
    vertex_position_mat = Matrix{Float64}(undef, nr_dimensions, nr_vertices_per_dimension^nr_dimensions)

    #loop through Cartesian Indices
    for i in CartesianIndices(vertex_position_mat)

        #this sophisticated equation gives the right entries such that all vertices have 
        #a distance from the supercell boundary of half an equilibrium bond length
        vertex_position_mat[i] =  ( ( ceil( i[2] / nr_vertices_per_dimension^(i[1]-1) ) + 1 )
                                    %nr_vertices_per_dimension + 0.5  ) 

    end

    return vertex_position_mat
end


"""
generate a simple cubic network using the graphs package
"""
function get_simple_cubic_network(nr_vertices; 
                                    nr_dimensions::Int64 = 3 )

    #calculate actual nr vertices from given nr vertices such that it can build 
    #a simple cubic lattice inside a cubic box without defects
    nr_vertices_per_dimension = Int(round( nr_vertices^(1/nr_dimensions) ))
    nr_vertices = nr_vertices_per_dimension^nr_dimensions

    #calculate edge length of supercell
    supercell_edge_length = nr_vertices_per_dimension

    #get matrix of vertex positions, where each column is a position vector
    vertex_position_mat = get_simple_cubic_vertex_position_mat(nr_vertices_per_dimension, 
                                                            nr_dimensions)

    #generate a graph by connecting all vertices of specified vertexic positions that are closer
    #to each other than the distance cutoff
    #p=2 is the Euclidean distance
    original_graph, edge_length_vec = Graphs.euclidean_graph(vertex_position_mat, 
                                            L= supercell_edge_length,
                                            p=2, 
                                            cutoff=1.1,
                                            bc=:periodic)

    #create dictionary out of original graph and its properties
    original_graph_dict = Dict("original_graph" => original_graph,
                    "edge_length_vec" => edge_length_vec,
                    "coordination_nr" => 2*nr_dimensions,
                    "nr_vertices" => nr_vertices,
                    "nr_dimensions" => nr_dimensions,
                    "supercell_edge_length" => supercell_edge_length,
                    "vertex_position_mat" => vertex_position_mat
                    )
    
    return original_graph_dict

end


"""
generate matrix of vertex positions in the cubic diamond structure,
where each column is a position vector.
Inside a unit cells, the vertexic positions in units of the nearest
neighbor distance are
1/sqrt(3) .* [(0,0,0), (0,2,2), (2,0,2), (2,2,0),
    (3,3,3), (3,1,1), (1,3,1), (1,1,3)]
"""
function get_diamond_vertex_position_mat(nr_unit_cells_per_dimension::Int64, 
                                    nr_vertices,
                                    edge_length_unit_cell)

    #generate empty matrix for vertexic positions
    vertex_position_mat = Matrix{Float64}(undef, 3, nr_vertices)

    #set the coordinates inside a unit cell in units of the equilibrium bond length
    coordinates_inside_unit_cell_vec = ( ( 1/sqrt(3) ) 
                                                .* [[0,0,0], [0,2,2], [2,0,2], [2,2,0],
                                                    [3,3,3], [3,1,1], [1,3,1], [1,1,3]] )

    #shift all coordinates, such that none lie on the edge of the supercell
    coordinate_shift_vector =  ( 1/(2*sqrt(3)) ) .* [1,1,1]

    #set counter of current vertex
    current_vertex_nr = 1

    #loop through all three dimensions
    for i in 0:nr_unit_cells_per_dimension-1
        for j in 0:nr_unit_cells_per_dimension-1
            for k in 0:nr_unit_cells_per_dimension-1

                for nr_vertex_inside_unit_cell in 1:8

                    #calculate position of the current vertex
                    vertex_position_mat[:, current_vertex_nr] = ( [i,j,k] .* edge_length_unit_cell
                                    .+ coordinates_inside_unit_cell_vec[nr_vertex_inside_unit_cell]
                                    .+ coordinate_shift_vector   )

                    #increase vertex counter
                    current_vertex_nr += 1

                end

            end
        end
    end

    return vertex_position_mat
end


"""
generate a diamond network using the graphs package.
This algorithm is based on the information that the unit cell contains 8 vertices
"""
function get_diamond_network(nr_vertices )

    #determine the edge length of a unit cell
    edge_length_unit_cell = 4/sqrt(3)

    #calculate the actual nr vertices, given that we require a 
    #cubic supercell and using the fact that the unit cell contains 8 vertices 
    nr_unit_cells_per_dimension = max(1, Int(round( (nr_vertices/8)^(1/3) )) )
    nr_vertices = 8 * nr_unit_cells_per_dimension^3

    #calculate edge length of supercell
    supercell_edge_length = nr_unit_cells_per_dimension*edge_length_unit_cell


    #get matrix of vertex positions, where each column is a position vector
    vertex_position_mat = get_diamond_vertex_position_mat(nr_unit_cells_per_dimension, 
                                                            nr_vertices,
                                                            edge_length_unit_cell)

    #generate a graph by connecting all vertices of specified vertexic positions that are closer
    #to each other than the distance cutoff
    #p=2 is the Euclidean distance
    original_graph, edge_length_vec = Graphs.euclidean_graph(vertex_position_mat, 
                                            L= supercell_edge_length,
                                            p=2, 
                                            cutoff=1.1,
                                            bc=:periodic)

    #create dictionary out of original graph and its properties
    original_graph_dict = Dict("original_graph" => original_graph,
                    "edge_length_vec" => edge_length_vec,
                    "coordination_nr" => 4,
                    "nr_vertices" => nr_vertices,
                    "nr_dimensions" => 3,
                    "supercell_edge_length" => supercell_edge_length,
                    "vertex_position_mat" => vertex_position_mat
                    )
    
    return original_graph_dict

end


"""
add information about vertexic positions and edge vectors to the original graph
"""
function convert_original_graph_to_spatial_network( original_graph_dict::Dict )

    #create an empty network graph where vertexic positions and edge vectors will be stored
    spatial_network = MetaGraphsNext.MetaGraph(Graphs.Graph(); 
                                        label_type = Int64,
                                        vertex_data_type = Dict{String, Any},
                                        edge_data_type = Dict{String, Any} )

    #label each vertex by its code integer and assign it its position vector
    for vertex in Graphs.vertices(original_graph_dict["original_graph"])

        spatial_network[vertex] = Dict( "position" => original_graph_dict["vertex_position_mat"][:,vertex] )

    end

    #label each edge by the vertices it connects and assign to it the vector where it points

    #get the nr and vector of original edges
    nr_edges = Graphs.ne(original_graph_dict["original_graph"])
    original_edges_vec = collect(Graphs.edges(original_graph_dict["original_graph"]))
    
    #loop through orinal edges to get edge descriptions
    for edge_nr in 1:nr_edges

        #get source and target of edge
        source = Graphs.src(original_edges_vec[edge_nr])
        target = Graphs.dst(original_edges_vec[edge_nr])

        #calculate vector from source to target considering periodic boundary conditions
        edge_vector = get_distance_vector_pbc(original_graph_dict["vertex_position_mat"][:,source],
                                            original_graph_dict["vertex_position_mat"][:,target],
                                            original_graph_dict["supercell_edge_length"] )

        #save edge vector and its length
        spatial_network[source, target] =  Dict("vector" => edge_vector, 
        "distance_squared" => (original_graph_dict["edge_length_vec"][original_edges_vec[edge_nr]])^2 )

    end
    
    #create dictionary out of graph and its properties
    graph_dict = Dict("spatial_network" => spatial_network,
                    "coordination_nr" => original_graph_dict["coordination_nr"],
                    "nr_vertices" => original_graph_dict["nr_vertices"],
                    "nr_dimensions" => original_graph_dict["nr_dimensions"],
                    "supercell_edge_length" => original_graph_dict["supercell_edge_length"]
                    )

    return graph_dict


end


"""
create a network graph representing the given network structure
"""
function get_periodic_network(evolution_dict)

    #depending on the network structure, create an original graph that does not
    #contain vertexic positions and bond information

    #simple cubic which is defined for any dimensionality
    if cmp(evolution_dict["network_type"], "simple cubic") == 0

        original_graph_dict = get_simple_cubic_network(evolution_dict["nr_vertices"];
                                    nr_dimensions = evolution_dict["nr_dimensions"]) 

    #diamond which is only defined for 3d
    elseif cmp(evolution_dict["network_type"], "diamond") == 0

        if evolution_dict["nr_dimensions"] == 3

            original_graph_dict = get_diamond_network(evolution_dict["nr_vertices"] ) 
        else
            @error "The diamond network is only defined in 3d."
        end

    else
        @error "Only simple cubic and diamond networks are implemented so far."

    end

    #convert original graph into a network graph that contains positional information
    graph_dict = convert_original_graph_to_spatial_network( original_graph_dict )

    #set bond bending constant
    graph_dict["bond_bending_const"] = evolution_dict["bond_bending_const"]

    
    #thermally excite network if desired
    if evolution_dict["thermal_fluctuations"]
        graph_dict["total_energy_up_to_date"] = false

        graph_dict = excite_entire_network!(graph_dict,
            evolution_dict;
            relax_first = false,
            update_total_energy = true)

    #otherwise just get total energy
    else

        #get total energy
        graph_dict["total_energy"] = get_total_energy_keating(graph_dict)
        graph_dict["total_energy_up_to_date"] = true

    end

    return graph_dict

end


"""
Create network of random vertex positions and connections
"""
function get_poisson_random_network(evolution_dict::Dict)

    #diamond which is only defined for 3d
    if cmp(evolution_dict["network_type"], "diamond") == 0

        if evolution_dict["nr_dimensions"] == 3

            original_graph_dict = get_diamond_network(evolution_dict["nr_vertices"] ) 
        else
            @error "The diamond network is only defined in 3d."
        end

    else
        @error "Only simple cubic and diamond networks are implemented so far."

    end

    #create an empty network graph where vertexic positions and edge vectors will be stored
    spatial_network = MetaGraphsNext.MetaGraph(Graphs.Graph(); 
                                        label_type = Int64,
                                        vertex_data_type = Dict{String, Any},
                                        edge_data_type = Dict{String, Any} )

    #label each vertex by its code integer and assign it its position vector
    for vertex in Graphs.vertices(original_graph_dict["original_graph"])

        spatial_network[vertex] = Dict( "position" => rand(Float64, (3)) .* original_graph_dict["supercell_edge_length"] )

    end

    #label each edge by the vertices it connects and assign to it the vector where it points

    #get the nr and vector of original edges
    nr_edges = Graphs.ne(original_graph_dict["original_graph"])
    original_edges_vec = collect(Graphs.edges(original_graph_dict["original_graph"]))
    
    #loop through orinal edges to get edge descriptions
    for edge_nr in 1:nr_edges

        #get source and target of edge
        source = Graphs.src(original_edges_vec[edge_nr])
        target = Graphs.dst(original_edges_vec[edge_nr])

        #calculate vector from source to target considering periodic boundary conditions
        edge_vector = get_distance_vector_pbc(spatial_network[source]["position"],
                                            spatial_network[target]["position"],
                                            original_graph_dict["supercell_edge_length"] )

        #save edge vector and its length
        spatial_network[source, target] =  Dict("vector" => edge_vector, 
        "distance_squared" => LinearAlgebra.norm(edge_vector)^2 )

    end
    
    #create dictionary out of graph and its properties
    graph_dict = Dict("spatial_network" => spatial_network,
                    "coordination_nr" => original_graph_dict["coordination_nr"],
                    "nr_vertices" => original_graph_dict["nr_vertices"],
                    "nr_dimensions" => original_graph_dict["nr_dimensions"],
                    "supercell_edge_length" => original_graph_dict["supercell_edge_length"]
                    )


    #set bond bending constant
    graph_dict["bond_bending_const"] = evolution_dict["bond_bending_const"]

    
    #thermally excite network if desired
    if evolution_dict["thermal_fluctuations"]
        graph_dict["total_energy_up_to_date"] = false

        graph_dict = excite_entire_network!(graph_dict,
            evolution_dict;
            relax_first = false,
            update_total_energy = true)

    #otherwise just get total energy
    else

        #get total energy
        graph_dict["total_energy"] = get_total_energy_keating(graph_dict)
        graph_dict["total_energy_up_to_date"] = true

    end

    return graph_dict
    
end