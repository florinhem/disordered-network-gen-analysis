"""
These functions generate spatial networks
"""


"""
generate matrix of vertex positions in a primitive cubic lattice, 
where each column is a position vector
"""
function get_primitive_cubic_vertex_position_mat(
    nr_vertices_per_dimension::Int64, 
    nr_dimensions::Int64)

    vertex_position_mat = Matrix{Float64}(
        undef, 
        nr_dimensions, 
        nr_vertices_per_dimension^nr_dimensions)

    for i in CartesianIndices(vertex_position_mat)

        # this sophisticated equation gives the right entries such that all 
        # vertices have a distance from the supercell boundary of half an
        # equilibrium bond length
        vertex_position_mat[i] =  (
            ( ceil( i[2] / nr_vertices_per_dimension^(i[1]-1) ) + 1 )
            %nr_vertices_per_dimension + 0.5  ) 

    end

    return vertex_position_mat
end


"""
generate a primitive cubic network using the graphs package
"""
function get_primitive_cubic_network(
    nr_vertices::Int64; 
    nr_dimensions::Int64 = 3 )

    # calculate actual nr vertices from given nr vertices such that it can
    # build a primitive cubic lattice inside a cubic box without defects
    nr_vertices_per_dimension = Int(round( nr_vertices^(1/nr_dimensions) ))
    nr_vertices = nr_vertices_per_dimension^nr_dimensions

    # calculate edge length of supercell
    supercell_edge_length = nr_vertices_per_dimension

    # get matrix of vertex positions, where each column is a position vector
    vertex_position_mat = get_primitive_cubic_vertex_position_mat(
        nr_vertices_per_dimension, 
        nr_dimensions)

    # generate a graph by connecting all vertices of specified vertex positions
    # that are closer to each other than the distance cutoff
    # p=2 is the Euclidean distance
    original_graph, edge_length_vec = Graphs.euclidean_graph(
        vertex_position_mat, 
        L= supercell_edge_length,
        p=2, 
        cutoff=1.1,
        bc=:periodic)

    original_spatial_network = Dict("original_graph" => original_graph,
                    "edge_length_vec" => edge_length_vec,
                    "coordination_nr_vec" => fill(2*nr_dimensions,nr_vertices),
                    "nr_vertices" => nr_vertices,
                    "nr_dimensions" => nr_dimensions,
                    "supercell_edge_length" => supercell_edge_length,
                    "vertex_position_mat" => vertex_position_mat
                    )
    
    return original_spatial_network
end


"""
generate matrix of vertex positions in the bcc structure,
where each column is a position vector.
Inside a unit cell, the vertex positions in units of the nearest
neighbor distance are (1/(2*sqrt(3))) .* [[1, 1, 1], [3, 3, 3]]
"""
function get_bcc_vertex_position_mat(
    nr_unit_cells_per_dimension::Int64, 
    nr_vertices::Int64,
    edge_length_unit_cell)

    vertex_position_mat = Matrix{Float64}(undef, 3, nr_vertices)

    # set the coordinates inside a unit cell in units of the equilibrium bond
    # length
    coordinates_inside_unit_cell_vec = (
        (1/(2*sqrt(3))) .* [[1, 1, 1], [3, 3, 3]] )

    current_vertex_nr = 1

    for i in 0:nr_unit_cells_per_dimension-1
        for j in 0:nr_unit_cells_per_dimension-1
            for k in 0:nr_unit_cells_per_dimension-1

                for nr_vertex_inside_unit_cell in 1:2

                    vertex_position_mat[:, current_vertex_nr] = ( 
                        [i,j,k] .* edge_length_unit_cell
                        .+ coordinates_inside_unit_cell_vec[
                            nr_vertex_inside_unit_cell]  )

                    current_vertex_nr += 1

                end
            end
        end
    end

    return vertex_position_mat
end


"""
generate a bcc network using the graphs package. This algorithm is based on the
information that the unit cell contains 2
"""
function get_bcc_network(nr_vertices)

    edge_length_unit_cell = 2/sqrt(3)

    # calculate the actual nr vertices, given that we require a 
    # cubic supercell and using the fact that the unit cell contains 2 vertices 
    nr_unit_cells_per_dimension = max(1, Int(round( (nr_vertices/2)^(1/3) )) )
    nr_vertices = 2 * nr_unit_cells_per_dimension^3

    supercell_edge_length = nr_unit_cells_per_dimension*edge_length_unit_cell

    # get matrix of vertex positions, where each column is a position vector
    vertex_position_mat = get_bcc_vertex_position_mat(
        nr_unit_cells_per_dimension, 
        nr_vertices,
        edge_length_unit_cell)

    # generate a graph by connecting all vertices of specified vertex positions
    # that are closer to each other than the distance cutoff
    # p=2 is the Euclidean distance
    original_graph, edge_length_vec = Graphs.euclidean_graph(
        vertex_position_mat, 
        L= supercell_edge_length,
        p=2, 
        cutoff=1.05,
        bc=:periodic)

    original_spatial_network = Dict("original_graph" => original_graph,
                    "edge_length_vec" => edge_length_vec,
                    "coordination_nr_vec" => fill(8,nr_vertices),
                    "nr_vertices" => nr_vertices,
                    "nr_dimensions" => 3,
                    "supercell_edge_length" => supercell_edge_length,
                    "vertex_position_mat" => vertex_position_mat
                    )
    
    return original_spatial_network
end


"""
generate matrix of vertex positions in the fcc structure,
where each column is a position vector.
Inside a unit cell, the vertex positions in units of the nearest
neighbor distance are
(sqrt(2)/4) .* [[1, 1, 1], [1, 3, 3], [3, 1, 3], [3, 3, 1]]
"""
function get_fcc_vertex_position_mat(
    nr_unit_cells_per_dimension::Int64, 
    nr_vertices::Int64,
    edge_length_unit_cell)

    # generate empty matrix for vertex positions
    vertex_position_mat = Matrix{Float64}(undef, 3, nr_vertices)

    # set the coordinates inside a unit cell in units of the equilibrium bond
    # length
    coordinates_inside_unit_cell_vec = (
        (sqrt(2)/4) .* [[1, 1, 1], [1, 3, 3], [3, 1, 3], [3, 3, 1]] )

    current_vertex_nr = 1

    # loop through all three dimensions
    for i in 0:nr_unit_cells_per_dimension-1
        for j in 0:nr_unit_cells_per_dimension-1
            for k in 0:nr_unit_cells_per_dimension-1

                for nr_vertex_inside_unit_cell in 1:4

                    # calculate position of the current vertex
                    vertex_position_mat[:, current_vertex_nr] = ( 
                        [i,j,k] .* edge_length_unit_cell
                        .+ coordinates_inside_unit_cell_vec[
                            nr_vertex_inside_unit_cell]  )

                    # increase vertex counter
                    current_vertex_nr += 1

                end

            end
        end
    end

    return vertex_position_mat
end


"""
generate a fcc network using the graphs package. This algorithm is based on the
information that the unit cell contains 4 vertices
"""
function get_fcc_network(nr_vertices )

    # determine the edge length of a unit cell
    edge_length_unit_cell = sqrt(2)

    # calculate the actual nr vertices, given that we require a 
    # cubic supercell and using the fact that the unit cell contains 4 vertices 
    nr_unit_cells_per_dimension = max(1, Int(round( (nr_vertices/4)^(1/3) )) )
    nr_vertices = 4 * nr_unit_cells_per_dimension^3

    # calculate edge length of supercell
    supercell_edge_length = nr_unit_cells_per_dimension*edge_length_unit_cell


    # get matrix of vertex positions, where each column is a position vector
    vertex_position_mat = get_fcc_vertex_position_mat(
        nr_unit_cells_per_dimension, 
        nr_vertices,
        edge_length_unit_cell)

    # generate a graph by connecting all vertices of specified vertex positions
    # that are closer to each other than the distance cutoff
    # p=2 is the Euclidean distance
    original_graph, edge_length_vec = Graphs.euclidean_graph(
        vertex_position_mat, 
        L= supercell_edge_length,
        p=2, 
        cutoff=1.05,
        bc=:periodic)

    # create dictionary out of original graph and its properties
    original_spatial_network = Dict("original_graph" => original_graph,
                    "edge_length_vec" => edge_length_vec,
                    "coordination_nr_vec" => fill(12,nr_vertices),
                    "nr_vertices" => nr_vertices,
                    "nr_dimensions" => 3,
                    "supercell_edge_length" => supercell_edge_length,
                    "vertex_position_mat" => vertex_position_mat
                    )
    
    return original_spatial_network
end


"""
generate matrix of vertex positions in the cubic diamond structure,
where each column is a position vector.
Inside a unit cell, the vertex positions in units of the nearest
neighbor distance are
1/sqrt(3) .* [[0,0,0], [0,2,2], [2,0,2], [2,2,0],
    [3,3,3], [3,1,1], [1,3,1], [1,1,3]]
"""
function get_diamond_vertex_position_mat(
    nr_unit_cells_per_dimension::Int64, 
    nr_vertices,
    edge_length_unit_cell)

    # generate empty matrix for vertex positions
    vertex_position_mat = Matrix{Float64}(undef, 3, nr_vertices)

    # set the coordinates inside a unit cell in units of the equilibrium bond
    # length
    coordinates_inside_unit_cell_vec = (
        ( 1/sqrt(3) ) .* [[0,0,0], [0,2,2], [2,0,2], [2,2,0],
                          [3,3,3], [3,1,1], [1,3,1], [1,1,3]] )

    # shift all coordinates, such that none lie on the edge of the supercell
    coordinate_shift_vector =  ( 1/(2*sqrt(3)) ) .* [1,1,1]

    current_vertex_nr = 1

    for i in 0:nr_unit_cells_per_dimension-1
        for j in 0:nr_unit_cells_per_dimension-1
            for k in 0:nr_unit_cells_per_dimension-1

                for nr_vertex_inside_unit_cell in 1:8

                    vertex_position_mat[:, current_vertex_nr] = ( 
                        [i,j,k] .* edge_length_unit_cell
                        .+ coordinates_inside_unit_cell_vec[
                            nr_vertex_inside_unit_cell]
                        .+ coordinate_shift_vector   )

                    current_vertex_nr += 1

                end

            end
        end
    end

    return vertex_position_mat
end


"""
generate a diamond network using the graphs package. This algorithm is based on
the information that the unit cell contains 8 vertices
"""
function get_diamond_network(nr_vertices)
    
    edge_length_unit_cell = 4/sqrt(3)

    # calculate the actual nr vertices, given that we require a 
    # cubic supercell and using the fact that the unit cell contains 8 vertices 
    nr_unit_cells_per_dimension = max(1, Int(round( (nr_vertices/8)^(1/3) )) )
    nr_vertices = 8 * nr_unit_cells_per_dimension^3

    supercell_edge_length = nr_unit_cells_per_dimension*edge_length_unit_cell


    # get matrix of vertex positions, where each column is a position vector
    vertex_position_mat = get_diamond_vertex_position_mat(
        nr_unit_cells_per_dimension, 
        nr_vertices,
        edge_length_unit_cell)

    # generate a graph by connecting all vertices of specified vertex positions
    # that are closer to each other than the distance cutoff
    # p=2 is the Euclidean distance
    original_graph, edge_length_vec = Graphs.euclidean_graph(
        vertex_position_mat, 
        L= supercell_edge_length,
        p=2, 
        cutoff=1.1,
        bc=:periodic)

    original_spatial_network = Dict("original_graph" => original_graph,
                    "edge_length_vec" => edge_length_vec,
                    "coordination_nr_vec" => fill(4,nr_vertices),
                    "nr_vertices" => nr_vertices,
                    "nr_dimensions" => 3,
                    "supercell_edge_length" => supercell_edge_length,
                    "vertex_position_mat" => vertex_position_mat
                    )
    
    return original_spatial_network
end



"""
generate matrix of vertex positions in the cubic gyroid structure,
where each column is a position vector.
Inside a unit cell, the vertex positions in units of the nearest
neighbor distance are
(1/8) .* [[7, 5, 1], [5, 3, 1], [3, 3, 3], [1, 5, 3], 
          [1, 7, 5], [3, 1, 5], [5, 1, 7], [7, 7, 7]]
"""
function get_gyroid_vertex_position_mat(
    nr_unit_cells_per_dimension::Int64, 
    nr_vertices,
    edge_length_unit_cell)

    # generate empty matrix for vertex positions
    vertex_position_mat = Matrix{Float64}(undef, 3, nr_vertices)

    # set the coordinates inside a unit cell in units of the equilibrium bond
    # length
    coordinates_inside_unit_cell_vec = (
        ((sqrt(8)/8)) .*    [[7, 5, 1], [5, 3, 1], [3, 3, 3], [1, 5, 3], 
                     [1, 7, 5], [3, 1, 5], [5, 1, 7], [7, 7, 7]] )

    current_vertex_nr = 1

    for i in 0:nr_unit_cells_per_dimension-1
        for j in 0:nr_unit_cells_per_dimension-1
            for k in 0:nr_unit_cells_per_dimension-1

                for nr_vertex_inside_unit_cell in 1:8   #TODO Change 8 to eachindex

                    vertex_position_mat[:, current_vertex_nr] = ( 
                        [i,j,k] .* edge_length_unit_cell
                        .+ coordinates_inside_unit_cell_vec[
                            nr_vertex_inside_unit_cell]
                        )

                    current_vertex_nr += 1

                end

            end
        end
    end

    return vertex_position_mat
end


"""
generate a gyroid network using the graphs package. This algorithm is based on
the information that the unit cell contains 8 vertices
"""
function get_gyroid_network(nr_vertices)
    
    edge_length_unit_cell = sqrt(8) #TODO: Ask if this is correct

    # calculate the actual nr vertices, given that we require a 
    # cubic supercell and using the fact that the unit cell contains 8 vertices 
    nr_unit_cells_per_dimension = max(1, Int(round( (nr_vertices/8)^(1/3) )) )
    nr_vertices = 8 * nr_unit_cells_per_dimension^3

    supercell_edge_length = nr_unit_cells_per_dimension*edge_length_unit_cell


    # get matrix of vertex positions, where each column is a position vector
    vertex_position_mat = get_gyroid_vertex_position_mat(
        nr_unit_cells_per_dimension, 
        nr_vertices,
        edge_length_unit_cell)

    # generate a graph by connecting all vertices of specified vertex positions
    # that are closer to each other than the distance cutoff
    # p=2 is the Euclidean distance
    original_graph, edge_length_vec = Graphs.euclidean_graph(
        vertex_position_mat, 
        L= supercell_edge_length,
        p=2, 
        cutoff=1.1,
        bc=:periodic)

    original_spatial_network = Dict("original_graph" => original_graph,
                    "edge_length_vec" => edge_length_vec,
                    "coordination_nr_vec" => fill(3,nr_vertices),
                    "nr_vertices" => nr_vertices,
                    "nr_dimensions" => 3,
                    "supercell_edge_length" => supercell_edge_length,
                    "vertex_position_mat" => vertex_position_mat
                    )
    
    return original_spatial_network
end




"""
generate matrix of vertex positions in the cubic ctn structure,
where each column is a position vector.
Inside a unit cell, the vertex positions are given
"""
function get_ctn_vertex_position_mat(
    nr_unit_cells_per_dimension::Int64, 
    nr_vertices,
    edge_length_unit_cell)

    # generate empty matrix for vertex positions
    vertex_position_mat = Matrix{Float64}(undef, 3, nr_vertices)

    # set the coordinates inside a unit cell in units of the equilibrium bond
    # length
    coordinates_inside_unit_cell_vec = 
    1/(sqrt((3/8-0.2083)^2+(0-0.2083)^2+(1/4-0.2083)^2)) .*
    ([
        [0,0.25,0.375],
        [0,0.75,0.125],
        [0.0417000000000001,0.4583000000000002,0.5416999999999998],
        [0.0417000000000001,0.5416999999999996,0.9583000000000002],
        [0.1249999999999998,0,0.7499999999999998],
        [0.2083000000000002,0.2083000000000002,0.2083000000000002],
        [0.2083000000000002,0.7916999999999998,0.2916999999999998],
        [0.25,0.375,0],
        [0.25,0.625,0.5],
        [0.2916999999999998,0.2083000000000002,0.7916999999999998],
        [0.2916999999999998,0.7916999999999998,0.7083000000000002],
        [0.375,0,0.25],
        [0.4583000000000002,0.4583000000000002,0.4583000000000002],
        [0.4582999999999999,0.5417000000000001,0.0416999999999998],
        [0.5000000000000002,0.25,0.625],
        [0.4999999999999998,0.75,0.8750000000000002],
        [0.5417000000000001,0.0416999999999996,0.4583000000000002],
        [0.5417000000000001,0.9583000000000002,0.0416999999999998],
        [0.625,0.5000000000000002,0.2499999999999998],
        [0.7083000000000002,0.2916999999999998,0.7916999999999998],
        [0.7083000000000002,0.7083000000000002,0.7083000000000002],
        [0.75,0.125,0],
        [0.75,0.8750000000000002,0.4999999999999998],
        [0.7916999999999998,0.2916999999999998,0.2083000000000002],
        [0.7916999999999998,0.7083000000000002,0.2916999999999998],
        [0.8750000000000002,0.4999999999999998,0.75],
        [0.9583000000000002,0.0417000000000001,0.5416999999999998],
        [0.9583000000000002,0.9583000000000002,0.9583000000000002],
    ] )

    current_vertex_nr = 1

    for i in 0:nr_unit_cells_per_dimension-1
        for j in 0:nr_unit_cells_per_dimension-1
            for k in 0:nr_unit_cells_per_dimension-1

                for nr_vertex_inside_unit_cell in eachindex(coordinates_inside_unit_cell_vec)

                    vertex_position_mat[:, current_vertex_nr] = ( 
                        [i,j,k] .* edge_length_unit_cell
                        .+ coordinates_inside_unit_cell_vec[
                            nr_vertex_inside_unit_cell]
                        )

                    current_vertex_nr += 1

                end

            end
        end
    end

    return vertex_position_mat
end


"""
generate a ctn network using the graphs package. This algorithm is based on
the information that the unit cell contains 28 vertices
"""
function get_ctn_network(nr_vertices)
    
    edge_length_unit_cell = 1/(sqrt((3/8-0.2083)^2+(0-0.2083)^2+(1/4-0.2083)^2)) #TODO Check this

    # calculate the actual nr vertices, given that we require a 
    # cubic supercell and using the fact that the unit cell contains 28 vertices 
    nr_unit_cells_per_dimension = max(1, Int(round( (nr_vertices/28)^(1/3) )) )
    nr_vertices = 28 * nr_unit_cells_per_dimension^3

    supercell_edge_length = nr_unit_cells_per_dimension*edge_length_unit_cell


    # get matrix of vertex positions, where each column is a position vector
    vertex_position_mat = get_ctn_vertex_position_mat(
        nr_unit_cells_per_dimension, 
        nr_vertices,
        edge_length_unit_cell)

    # generate a graph by connecting all vertices of specified vertex positions
    # that are closer to each other than the distance cutoff
    # p=2 is the Euclidean distance
    original_graph, edge_length_vec = Graphs.euclidean_graph(
        vertex_position_mat, 
        L= supercell_edge_length,
        p=2, 
        cutoff=1.1,
        bc=:periodic)

        println(vertex_position_mat)
        println(size(vertex_position_mat))

    coordination_nr_vec::Vector{Int64}=fill(-1,size(vertex_position_mat,2))
   
    for vertex in Graphs.vertices(original_graph)
        nr_vertex=0
        for edge in Graphs.neighbors(original_graph,vertex)
            #println("edge, $edge")
            nr_vertex+=1
        end
        coordination_nr_vec[vertex]=nr_vertex
    end

    println("coordination_nr_vec, $coordination_nr_vec")

    original_spatial_network = Dict("original_graph" => original_graph,
                    "edge_length_vec" => edge_length_vec,
                    "coordination_nr_vec" => coordination_nr_vec,
                    "nr_vertices" => nr_vertices,
                    "nr_dimensions" => 3,
                    "supercell_edge_length" => supercell_edge_length,
                    "vertex_position_mat" => vertex_position_mat
                    )
    
    return original_spatial_network
end






"""
add information about vertex positions and edge vectors to the original graph
"""
function convert_original_graph_to_spatial_network(
    original_spatial_network::Dict)

    # create an empty network graph where vertex positions and edge vectors
    # will be stored
    spatial_network = MetaGraphsNext.MetaGraph(
        Graphs.Graph(); 
        label_type = Int64,
        vertex_data_type = Dict{String, Any},
        edge_data_type = Dict{String, Union{Float64, Vector{Float64}}},
        graph_data = Dict{String, Any}(
            "nr_vertices" => original_spatial_network["nr_vertices"],
            "nr_dimensions" => original_spatial_network["nr_dimensions"],
            "supercell_edge_length" 
                => original_spatial_network["supercell_edge_length"])
        )

    # label each vertex by its code integer and assign it its position vector
    for vertex in Graphs.vertices(original_spatial_network["original_graph"])
        
        spatial_network[vertex] = Dict{String, Any}( 
            "position" => original_spatial_network["vertex_position_mat"][
                :,vertex],
            "coordination_nr" => original_spatial_network["coordination_nr_vec"][vertex]
        )

    end

    # label each edge by the vertices it connects and assign to it the vector
    # where it points

    # get the nr and vector of original edges
    nr_edges = Graphs.ne(original_spatial_network["original_graph"])
    original_edges_vec = collect(
        Graphs.edges(original_spatial_network["original_graph"]))
    
    for edge_nr in 1:nr_edges

        source = Graphs.src(original_edges_vec[edge_nr])
        target = Graphs.dst(original_edges_vec[edge_nr])

        # calculate vector from source to target considering periodic boundary
        # conditions
        edge_vector = get_distance_vector_pbc(
            original_spatial_network["vertex_position_mat"][:,source],
            original_spatial_network["vertex_position_mat"][:,target],
            original_spatial_network["supercell_edge_length"] )

        spatial_network[source, target] =  Dict("vector" => edge_vector, 
        "distance_squared" => (
            original_spatial_network["edge_length_vec"][
                original_edges_vec[edge_nr]] )^2 )

    end

    return spatial_network
end


"""
create a network graph representing the given network structure
"""
function get_periodic_network(evolution_dict)

    # depending on the network structure, create an original graph that does not
    # contain vertex positions and bond information

    # primitive cubic which is defined for any dimensionality
    if cmp(evolution_dict["network_type"], "primitive cubic") == 0

        original_spatial_network = get_primitive_cubic_network(
            evolution_dict["nr_vertices"];
            nr_dimensions = evolution_dict["nr_dimensions"]) 

    # bcc which is only defined for 3d
    elseif cmp(evolution_dict["network_type"], "bcc") == 0

        if evolution_dict["nr_dimensions"] == 3

            original_spatial_network = get_bcc_network(
                evolution_dict["nr_vertices"] ) 
        else
            @error "The bcc network is only defined in 3d."
        end

    # fcc which is only defined for 3d
    elseif cmp(evolution_dict["network_type"], "fcc") == 0

        if evolution_dict["nr_dimensions"] == 3

            original_spatial_network = get_fcc_network(
                evolution_dict["nr_vertices"] ) 
        else
            @error "The fcc network is only defined in 3d."
        end

    # diamond which is only defined for 3d
    elseif cmp(evolution_dict["network_type"], "diamond") == 0

        if evolution_dict["nr_dimensions"] == 3

            original_spatial_network = get_diamond_network(
                evolution_dict["nr_vertices"] ) 
        else
            @error "The diamond network is only defined in 3d."
        end

    elseif cmp(evolution_dict["network_type"], "gyroid") == 0

        if evolution_dict["nr_dimensions"] == 3

            original_spatial_network = get_gyroid_network(
                evolution_dict["nr_vertices"] ) 
        else
            @error "The gyroid network is only defined in 3d."
        end

    elseif cmp(evolution_dict["network_type"], "ctn") == 0

        if evolution_dict["nr_dimensions"] == 3

            original_spatial_network = get_ctn_network(
                evolution_dict["nr_vertices"] ) 
        else
            @error "The ctn network is only defined in 3d."
        end

    else
        @error "Only primitive cubic and diamond networks are implemented so
            far."

    end

    # convert original graph into a network graph that contains positional
    # information
    spatial_network = convert_original_graph_to_spatial_network(
        original_spatial_network)

    spatial_network[]["bond_bending_const"] = evolution_dict[
        "bond_bending_const"]
    
    spatial_network[]["theta_ground_state"] = evolution_dict[
        "theta_ground_state"]

    # thermally excite network if desired
    if evolution_dict["thermal_fluctuations"]
        spatial_network[]["total_energy_up_to_date"] = false

        spatial_network = excite_entire_network!(spatial_network,
            evolution_dict;
            relax_first = false,
            update_total_energy = true)

    # otherwise just get total energy
    else
        spatial_network[]["total_energy"] = get_total_energy_keating(
            spatial_network)
        spatial_network[]["total_energy_up_to_date"] = true

    end

    return spatial_network
end


"""
Create network of random vertex positions and connections
"""
function get_poisson_random_network(evolution_dict::Dict)

    # diamond which is only defined for 3d
    if cmp(evolution_dict["network_type"], "diamond") == 0

        if evolution_dict["nr_dimensions"] == 3

            original_spatial_network = get_diamond_network(
                evolution_dict["nr_vertices"] ) 
        else
            @error "The diamond network is only defined in 3d."
        end

    else
        @error "Only diamond networks are implemented so far."

    end

    # create an empty network graph where vertex positions and edge vectors
    # will be stored
    spatial_network = MetaGraphsNext.MetaGraph(
        Graphs.Graph(); 
        label_type = Int64,
        vertex_data_type = Dict{String, Any},
        edge_data_type = Dict{String, Any},
        graph_data = Dict{String, Any}(
            "nr_vertices" => original_spatial_network["nr_vertices"],
            "nr_dimensions" => original_spatial_network["nr_dimensions"],
            "supercell_edge_length" 
                => original_spatial_network["supercell_edge_length"])
        )

    # label each vertex by its code integer and assign it its position vector
    for vertex in Graphs.vertices(original_spatial_network["original_graph"])

        spatial_network[vertex] = Dict( 
            "position" => rand(Float64, (3)) 
                .* original_spatial_network["supercell_edge_length"],
            "coordination_nr" => original_spatial_network["coordination_nr_vec"][vertex]
        )

    end

    # label each edge by the vertices it connects and assign to it the vector
    # where it points

    nr_edges = Graphs.ne(original_spatial_network["original_graph"])
    original_edges_vec = collect(
        Graphs.edges(original_spatial_network["original_graph"]))
    
    # loop through orinal edges to get edge descriptions
    for edge_nr in 1:nr_edges

        # get source and target of edge
        source = Graphs.src(original_edges_vec[edge_nr])
        target = Graphs.dst(original_edges_vec[edge_nr])

        edge_vector = get_distance_vector_pbc(
            spatial_network[source]["position"],
            spatial_network[target]["position"],
            original_spatial_network["supercell_edge_length"] )

        spatial_network[source, target] =  Dict("vector" => edge_vector, 
        "distance_squared" => LinearAlgebra.norm(edge_vector)^2 )

    end

    spatial_network[]["bond_bending_const"] = evolution_dict[
        "bond_bending_const"]

    spatial_network[]["theta_ground_state"] = evolution_dict[
        "theta_ground_state"]
    
    # thermally excite network if desired
    if evolution_dict["thermal_fluctuations"]
        spatial_network[]["total_energy_up_to_date"] = false

        spatial_network = excite_entire_network!(spatial_network,
            evolution_dict;
            relax_first = false,
            update_total_energy = true)

    # otherwise just get total energy
    else
        spatial_network[]["total_energy"] = get_total_energy_keating(
            spatial_network)
        spatial_network[]["total_energy_up_to_date"] = true

    end

    return spatial_network
end