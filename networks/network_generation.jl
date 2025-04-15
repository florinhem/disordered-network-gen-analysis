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
function get_diamond_network_old(nr_vertices)
    
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

    println("vertex_position_mat, $vertex_position_mat")

    # generate a graph by connecting all vertices of specified vertex positions
    # that are closer to each other than the distance cutoff
    # p=2 is the Euclidean distance
    original_graph, edge_length_vec = Graphs.euclidean_graph(
        vertex_position_mat, 
        L= supercell_edge_length,
        p=2, 
        cutoff=1.1,
        bc=:periodic)

    println("original_graph, $original_graph")
    println("edge_length_vec, $edge_length_vec")

    return
    

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
returns the rotation matrix for a certain axis and certain angle
"""
function rotation_matrix(axis::Vector{Float64}, angle::Float64)
    axis = axis / LinearAlgebra.norm(axis)
    angle=angle/360*2*pi
    cos_angle = LinearAlgebra.cos(angle)
    sin_angle = LinearAlgebra.sin(angle)
    one_minus_cos = 1 - cos_angle
    x, y, z = axis

    return [
        cos_angle + x^2 * one_minus_cos  x * y * one_minus_cos - z * sin_angle  x * z * one_minus_cos + y * sin_angle;
        y * x * one_minus_cos + z * sin_angle  cos_angle + y^2 * one_minus_cos  y * z * one_minus_cos - x * sin_angle;
        z * x * one_minus_cos - y * sin_angle  z * y * one_minus_cos + x * sin_angle  cos_angle + z^2 * one_minus_cos
    ]
end

"""
makes a copy of edges and rotates n_fold times around a certain axis with
a certain shift. After that we translate the edges.
"""
function copy_and_rotate_and_translate(
    edges,
    axis::Vector{Float64}, 
    shift::Vector{Float64}, 
    translate::Vector{Float64}, 
    n_fold::Int64)

    new_edges=deepcopy(edges)
    current_edge=length(edges)+1
    for n in 1:(n_fold-1)
        angle=360/n_fold*n
        rot_matrix=rotation_matrix(axis,angle)

        for (edge_nr,(vertex_position_start,vertex_position_end,length)) in edges
            new_vertex_position_start=rot_matrix * (vertex_position_start .- shift) .+ shift .+ translate
            new_vertex_position_end=rot_matrix * (vertex_position_end .- shift) .+ shift .+ translate
            new_edges[current_edge]=(new_vertex_position_start,new_vertex_position_end,length)
            current_edge+=1
        end
    end

    #println(new_edges)
    return new_edges
end


"""
makes a copy of edges and applies a glide plane (mirror at a plane and then
translate along the mirror)
"""
function copy_and_glide(
    edges,
    translate::Vector{Float64},
    normal_vector_to_mirror_plane::Vector{Float64}
)
    # Normalize the normal vector
    normal_vector_to_mirror_plane = normal_vector_to_mirror_plane / LinearAlgebra.norm(normal_vector_to_mirror_plane)
    
    new_edges = deepcopy(edges)
    current_edge = length(edges) + 1

    for (edge_nr, (vertex_position_start, vertex_position_end, length)) in edges
        # Mirror the start and end positions
        mirrored_start = vertex_position_start - 2 * LinearAlgebra.dot(vertex_position_start, normal_vector_to_mirror_plane) * normal_vector_to_mirror_plane
        mirrored_end = vertex_position_end - 2 * LinearAlgebra.dot(vertex_position_end, normal_vector_to_mirror_plane) * normal_vector_to_mirror_plane
        
        # Translate the mirrored positions
        new_vertex_position_start = mirrored_start .+ translate
        new_vertex_position_end = mirrored_end .+ translate
        
        # Add the new edge to the dictionary
        new_edges[current_edge] = (new_vertex_position_start, new_vertex_position_end, length)
        current_edge += 1
    end
    
    return new_edges
end



"""
makes a copy of edges and translates all edges once in translate direction
"""
function copy_and_translate(
    edges, 
    translate::Vector{Float64})

    new_edges=deepcopy(edges)
    current_edge=length(edges)+1
    for (edge_nr,(vertex_position_start,vertex_position_end,length)) in edges
        new_vertex_position_start=vertex_position_start .+ translate
        new_vertex_position_end=vertex_position_end .+ translate
        new_edges[current_edge]=(new_vertex_position_start,new_vertex_position_end,length)
        current_edge+=1
    end
    
    #println(new_edges)
    return new_edges
end


"""
makes a copy of edges and translates n times in the translation direction
"""
function copy_and_translate_n_times(
    edges, 
    translate,
    n)

    new_edges=deepcopy(edges)
    current_edge=length(edges)+1
    

    for i in 1:(n-1)
        for (edge_nr,(vertex_position_start,vertex_position_end,length)) in edges 
            new_vertex_position_start=vertex_position_start .+ translate*i
            new_vertex_position_end=vertex_position_end .+ translate*i
            new_edges[current_edge]=(new_vertex_position_start,new_vertex_position_end,length)
           
            current_edge+=1
        end
    end
    
    return new_edges
end

"""
makes a copy of edges and translates with a vector in all 3 dimensions a 
certain number of times.
"""
function array_3D(
    edges, 
    translate_x,
    nbr_x,
    translate_y,
    nbr_y, 
    translate_z,
    nbr_z)
    
    edges = copy_and_translate_n_times(edges,translate_x,nbr_x)
    edges = copy_and_translate_n_times(edges,translate_y,nbr_y)
    edges = copy_and_translate_n_times(edges,translate_z,nbr_z)
    
    return edges
end

"""
deletes copy of edges. 
Compares if the edge starts and the edge ends is the same
for 2 different edges 
OR 
Compares if the (edge start and edge end) and (edge end 
and edge start) are the same. 

If its the case, then we have one edge that is 
doubled and can be deleted.
"""
function delete_copys(edges)

    epsilon=0.001
    new_edges=deepcopy(edges)

    for (edge_nr_1,(vertex_position_start_1,vertex_position_end_1,length_1)) in edges
        for (edge_nr_2,(vertex_position_start_2,vertex_position_end_2,length_2)) in edges
            if edge_nr_1<edge_nr_2
                if (length_1-length_2)<epsilon
                    if ((LinearAlgebra.norm(vertex_position_start_1-vertex_position_start_2)<epsilon &&
                       LinearAlgebra.norm(vertex_position_end_1-vertex_position_end_2)<epsilon) 
                       ||
                       (LinearAlgebra.norm(vertex_position_start_1-vertex_position_end_2)<epsilon &&
                       LinearAlgebra.norm(vertex_position_end_1-vertex_position_start_2)<epsilon))

                       delete!(new_edges, edge_nr_2)
                    end
                end
            end
        end
    end
    return new_edges
end 

"""
returns the distance for a periodic boundary condition
"""
function distance_periodic_bc(position_1, position_2, L)
    Δ = abs.(position_1-position_2)
    Δ = min.(L .- Δ, Δ)
    distance=LinearAlgebra.norm(Δ)

    return distance
end

"""
returns all vertex positions as a dictionary 
"""
function get_vertex_positions_dict(edges, epsilon, L)
    
    vertex_positions_dict=Dict()
    next_vertex_nr=1
    
    for (edge_nr,(vertex_position_start,vertex_position_end,length)) in edges 
        add_start_to_vertex_positions_dict=true
        add_end_to_vertex_positions_dict=true
        for (vertex_nr, vertex_position) in vertex_positions_dict
            
            if (distance_periodic_bc(vertex_position_start,vertex_position,L)<epsilon)
                add_start_to_vertex_positions_dict=false
            end

            if (distance_periodic_bc(vertex_position_end,vertex_position,L)<epsilon)
                add_end_to_vertex_positions_dict=false
            end
        end

        if add_start_to_vertex_positions_dict==true
            vertex_positions_dict[next_vertex_nr]=vertex_position_start
            next_vertex_nr+=1
        end

        if add_end_to_vertex_positions_dict==true
            vertex_positions_dict[next_vertex_nr]=vertex_position_end
            next_vertex_nr+=1
        end
    end

    return vertex_positions_dict
end

"""
returns a matrix for the dictionary
"""
function get_mat_from_dict(vertex_positions_dict)

    nr_vertices = length(vertex_positions_dict)
    nr_dimensions = length(vertex_positions_dict[1])
    vertex_position_mat = Matrix{Float64}(undef, nr_dimensions, nr_vertices)

    for (vertex_nr, vertex_position) in vertex_positions_dict
        vertex_position_mat[:, vertex_nr] = vertex_position
    end

    return vertex_position_mat
end

"""
returns the edges (starting position, ending position, length, 
starting vertex, ending vertex)
"""
function get_edges_with_vertex(edges, vertex_positions_dict, epsilon, L)

    edges_with_vertex=Dict()

    for (edge_nr,(vertex_position_start,vertex_position_end,length)) in edges
        vertex_start=nothing
        vertex_end=nothing
        for (vertex_nr, vertex_position) in vertex_positions_dict
            if vertex_start===nothing && distance_periodic_bc(vertex_position_start, vertex_position, L)<epsilon
                vertex_start=vertex_nr
            end
            if vertex_end===nothing && distance_periodic_bc(vertex_position_end, vertex_position, L)<epsilon
                vertex_end=vertex_nr
            end
        end
        if vertex_start<vertex_end
            edges_with_vertex[edge_nr]=(vertex_position_start,vertex_position_end,length,vertex_start, vertex_end)
        else
            edges_with_vertex[edge_nr]=(vertex_position_end,vertex_position_start,length,vertex_end, vertex_start)
        end
    end

    return edges_with_vertex
end

"""
returns the original graph that has all vertices and all edges
"""
function create_graph(edges_with_vertex, nr_vertices)
    
    original_graph=Graphs.SimpleGraph(nr_vertices)

    for (edge_nr,(vertex_position_start,vertex_position_end,length,vertex_start,vertex_end)) in edges_with_vertex
        Graphs.add_edge!(original_graph, vertex_start, vertex_end)
    end
    
    return original_graph
end

function get_edge_length_vec(original_graph, edges_with_vertex)

    edge_length_vec = Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Float64}()

    for edge in Graphs.edges(original_graph)
        v1 = Graphs.src(edge)
        v2 = Graphs.dst(edge)

        for (edge_nr,(vertex_position_start,vertex_position_end,length,vertex_start,vertex_end)) in edges_with_vertex
            if v1==vertex_start && v2==vertex_end
                edge_length_vec[edge]=length
            end
        end 
    end

    return edge_length_vec
end

"""
returns a dictionary with edges (starting position, ending position, length)
that are now in the block spanning from (0,0,0) to (size, size, size)
"""
function fold_to_block(edges, size)

    old_edges=deepcopy(edges)
    new_edges=Dict()
    current_edge=1

    for (edge_nr,(vertex_position_start,vertex_position_end,length)) in old_edges
        new_vertex_position_start=mod.(vertex_position_start,size)
        new_vertex_position_end=mod.(vertex_position_end,size)
        new_edges[current_edge]=(new_vertex_position_start,new_vertex_position_end,length)
        current_edge+=1
    end

    return new_edges
end

"""
returns a dictionary with edges (starting position, ending position, length)
that fullfill the condition to be in the cube spanning from 
(x_min,y_min,z_min) to (x_max,y_max,z_max)
"""
function filter_to_unitcell(
    edges,
    x_min,
    x_max,
    y_min,
    y_max,
    z_min,
    z_max)

    new_edges=deepcopy(edges)

    for (edge_nr,(vertex_position_start,vertex_position_end,length)) in edges
           
        x,y,z=vertex_position_start
        x2,y2,z2=vertex_position_end
        if (x_min>x || x_max<x || y_min>y || y_max<y || z_min>z || z_max<z) &&
           (x_min>x2 || x_max<x2 || y_min>y2 || y_max<y2 || z_min>z2 || z_max<z2)
            
           delete!(new_edges, edge_nr)
        end
    end

    return new_edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the diamond structure with the symmetry operations
"""
function get_edges_dia(edge_length_unit_cell)
    # define the edges with the help of rcsr.net
    # look at the spacegroup name and find the spacegroup number
    edges = Dict(1 => ([0.0,0.0,0.0] .* edge_length_unit_cell, [1/4,1/4,1/4] .* edge_length_unit_cell, sqrt(3)/4 .* edge_length_unit_cell))
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,[0.0,0.0,1.0],[0.0,1/4,0.0].* edge_length_unit_cell,[0.0,0.0,1/2] .* edge_length_unit_cell,2)
    edges = copy_and_rotate_and_translate(edges,[0.0,1.0,0.0],[1/4,0.0,0.0].* edge_length_unit_cell,[0.0,1/2,0.0] .* edge_length_unit_cell,2)
    edges = copy_and_rotate_and_translate(edges,[1.0,1.0,1.0],[0.0,0.0,0.0].* edge_length_unit_cell,[0.0,0.0,0.0] .* edge_length_unit_cell,3)
    edges = copy_and_rotate_and_translate(edges,[1.0,1.0,0.0],[0.0,-1/4,3/8].* edge_length_unit_cell,[1/2,1/2,0.0] .* edge_length_unit_cell,2)

    edges = copy_and_translate(edges, [0,1/2,1/2].* edge_length_unit_cell)
    edges = copy_and_translate(edges, [1/2,0,1/2].* edge_length_unit_cell)
    edges = copy_and_translate(edges, [1/2,1/2,0].* edge_length_unit_cell)

    # clean unnecessairy edges that are repeated
    edges = delete_copys(edges)

    return edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the srs (gyroid) structure with the symmetry operations
"""
function get_edges_srs(edge_length_unit_cell)
    # define the edges with the help of rcsr.net
    # look at the spacegroup name and find the spacegroup number
    edges = Dict(
        1 => ([1/8, 1/8, 1/8] .* edge_length_unit_cell, [-1/8, 3/8, 1/8] .* edge_length_unit_cell, sqrt(2)/4 .* edge_length_unit_cell),
        )
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,[0.0,0.0,1.0],[1/4,0.0,0.0].* edge_length_unit_cell,[0.0,0.0,1/2] .* edge_length_unit_cell,2)
    edges = copy_and_rotate_and_translate(edges,[0.0,1.0,0.0],[0.0,0.0,1/4].* edge_length_unit_cell,[0.0,1/2,0.0] .* edge_length_unit_cell,2)
    edges = copy_and_rotate_and_translate(edges,[1.0,1.0,1.0],[0.0,0.0,0.0].* edge_length_unit_cell,[0.0,0.0,0.0] .* edge_length_unit_cell,3)
    edges = copy_and_rotate_and_translate(edges,[1.0,1.0,0.0],[0.0,-1/4,1/8].* edge_length_unit_cell,[1/2,1/2,0.0] .* edge_length_unit_cell,2)

    edges = copy_and_translate(edges, [1/2,1/2,1/2].* edge_length_unit_cell)

    # clean unnecessairy edges that are repeated
    edges = delete_copys(edges)

    return edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the srd structure with the symmetry operations
"""
function get_edges_srd(edge_length_unit_cell)
    # define the edges with the help of rcsr.net
    # look at the spacegroup name and find the spacegroup number
    edges = Dict(
        1 => ([0.0, 1/4, 1/2] .* edge_length_unit_cell, [0.0, 3/4, 1/2] .* edge_length_unit_cell, 1/2 .* edge_length_unit_cell),
        2 => ([1/4, 1/4, 1/4] .* edge_length_unit_cell, [0.0, 1/4, 1/2] .* edge_length_unit_cell, sqrt(2)/4 .* edge_length_unit_cell)
        )
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,[0.0,0.0,1.0],[0.0,0.0,0.0].* edge_length_unit_cell,[0.0,0.0,0.0] .* edge_length_unit_cell,2)
    edges = copy_and_rotate_and_translate(edges,[0.0,1.0,0.0],[0.0,0.0,0.0].* edge_length_unit_cell,[0.0,0.0,0.0] .* edge_length_unit_cell,2)
    edges = copy_and_rotate_and_translate(edges,[1.0,1.0,1.0],[0.0,0.0,0.0].* edge_length_unit_cell,[0.0,0.0,0.0] .* edge_length_unit_cell,3)
    edges = copy_and_rotate_and_translate(edges,[1.0,1.0,0.0],[0.0,0.0,1/4].* edge_length_unit_cell,[1/2,1/2,0.0] .* edge_length_unit_cell,2)

    # clean unnecessairy edges that are repeated
    edges = delete_copys(edges)

    return edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the ctn structure with the symmetry operations
"""
function get_edges_ctn(edge_length_unit_cell)
    # define the edges with the help of rcsr.net
    # look at the spacegroup name and find the spacegroup number

    # it is possible to change the x value for the first vertex V1. V2 is fixed.
    x=0.2082
    V1=[x, x, x]
    V2=[3/8, 0.0, 1/4]

    # the difference between V1 and V2 gives us the length
    D1=V1 .- V2
    L1=LinearAlgebra.norm(D1)
    println(1/L1)
    
    edges = Dict(
        1 => (V1 .* edge_length_unit_cell, V2 .* edge_length_unit_cell, L1 .* edge_length_unit_cell)
    )
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,[0.0,0.0,1.0],[1/4,0.0,0.0].* edge_length_unit_cell,[0.0,0.0,1/2] .* edge_length_unit_cell,2)
    edges = copy_and_rotate_and_translate(edges,[0.0,1.0,0.0],[0.0,0.0,1/4].* edge_length_unit_cell,[0.0,1/2,0.0] .* edge_length_unit_cell,2)
    edges = copy_and_rotate_and_translate(edges,[1.0,1.0,1.0],[0.0,0.0,0.0].* edge_length_unit_cell,[0.0,0.0,0.0] .* edge_length_unit_cell,3)
    edges = copy_and_glide(edges,[1/4,1/4,1/4].* edge_length_unit_cell,[1.0,-1.0,0.0])

    edges = copy_and_translate(edges, [1/2,1/2,1/2].* edge_length_unit_cell)

    # clean unnecessairy edges that are repeated
    edges = delete_copys(edges)

    println("get_edges_ctn finished")
    return edges
end


"""
returns a vector for the coordination numbers of all vertices of the 
original graph
"""
function get_coordination_nr_vec(original_graph)
    coordination_nr_vec::Vector{Int64}=fill(-1,length(Graphs.vertices(original_graph)))

    for vertex in Graphs.vertices(original_graph)
        nr_vertex=0
        for edge in Graphs.neighbors(original_graph,vertex)
            nr_vertex+=1
        end

        coordination_nr_vec[vertex]=nr_vertex
    end

    return coordination_nr_vec
end


"""
returns the original space graph. It checks which network we want to generate.
Then it calculates for the unitcell and supercell the edges 
(number, starting vertex, ending vertex and length) of the network.
"""
function get_network(nr_vertices, network_name)
    if cmp(network_name , "dia") == 0       #diamond
        println("get_network, dia")
        nr_dimensions = 3
        edge_length_unit_cell = 4/sqrt(3)
        nr_vertices_per_unit_cell = 8
    elseif cmp(network_name , "srs") == 0   #gyroid
        println("get_network, srs")
        nr_dimensions = 3
        edge_length_unit_cell = 4/sqrt(2)
        nr_vertices_per_unit_cell = 8
    elseif cmp(network_name , "srd") == 0   #srd
        println("get_network, srd")
        nr_dimensions = 3
        edge_length_unit_cell = 2.3075 #analytically calculated when E_length minimal
        nr_vertices_per_unit_cell = 10
    elseif cmp(network_name , "ctn") == 0   #ctn
        println("get_network, ctn")
        nr_dimensions = 3
        edge_length_unit_cell = 3.7033 #analytically calculated when E_length minimal
        nr_vertices_per_unit_cell = 28
    else
        @error "Only dia, srd, srs, ctn are implemented, $network_name not."
    end


    # calculate the actual nr vertices, given that we require a 
    # cubic supercell and using the fact that the unit cell contains 10 vertices 
    nr_unit_cells_per_dimension = max(1, Int(round( (nr_vertices/nr_vertices_per_unit_cell)^(1/3) )) )
    nr_vertices = nr_vertices_per_unit_cell * nr_unit_cells_per_dimension^3
    supercell_edge_length = nr_unit_cells_per_dimension*edge_length_unit_cell

    if cmp(network_name , "dia") == 0
        edges = get_edges_dia(edge_length_unit_cell)
    elseif cmp(network_name , "srd") == 0
        edges = get_edges_srd(edge_length_unit_cell)
    elseif cmp(network_name , "srs") == 0
        edges = get_edges_srs(edge_length_unit_cell)
    elseif cmp(network_name , "ctn") == 0
        edges = get_edges_ctn(edge_length_unit_cell)
    end
  
    # copy all edges L times in x,y,z direction to get all edges 
    L=5
    epsilon=0.001
    x=[1.0,0.0,0.0].* edge_length_unit_cell
    y=[0.0,1.0,0.0].* edge_length_unit_cell
    z=[0.0,0.0,1.0].* edge_length_unit_cell

    edges = array_3D(edges,x,L,y,L,z,L)
    edges = delete_copys(edges)

    # look at the cube spanning from (1,1,1) and (2,2,2), only take these edges
    delta=0.01
    filter_min=(2.0-delta).* edge_length_unit_cell
    filter_max=(3.0+delta).* edge_length_unit_cell
    edges = filter_to_unitcell(edges, filter_min, filter_max, filter_min, filter_max, filter_min, filter_max)

    # copy the unitcell in all 3 dimensions to get the actual network
    edges = array_3D(edges,x,nr_unit_cells_per_dimension,y,nr_unit_cells_per_dimension,z,nr_unit_cells_per_dimension)
    edges = delete_copys(edges)

    # fold all edges back into the supercell
    edges = fold_to_block(edges, supercell_edge_length)
    edges = delete_copys(edges)
  
    # get the positions, edges, nr_vertices, original_graph and edges lengths
    vertex_positions_dict=get_vertex_positions_dict(edges, epsilon, supercell_edge_length)
    vertex_position_mat=get_mat_from_dict(vertex_positions_dict)
    edges_with_vertex=get_edges_with_vertex(edges, vertex_positions_dict, epsilon, supercell_edge_length)
    nr_vertices=length(vertex_positions_dict)
    original_graph=create_graph(edges_with_vertex, nr_vertices)
    edge_length_vec=get_edge_length_vec(original_graph, edges_with_vertex)

    # this if else can be removed, if you are sure that the network is correct
    println("edge_length_vec:, $(length(edge_length_vec))")
    if cmp(network_name , "dia") == 0
        println("nr of edges in unitcell should be 16=?=$(length(edge_length_vec)/(nr_unit_cells_per_dimension^3))")
    elseif cmp(network_name , "srd") == 0
        println("nr of edges in unitcell should be 18=?=$(length(edge_length_vec)/(nr_unit_cells_per_dimension^3))")     
    elseif cmp(network_name , "srs") == 0
        println("nr of edges in unitcell should be 12=?=$(length(edge_length_vec)/(nr_unit_cells_per_dimension^3))") 
    elseif cmp(network_name , "ctn") == 0
        println("nr of edges in unitcell should be 48=?=$(length(edge_length_vec)/(nr_unit_cells_per_dimension^3))")      
    end

    # calculate the coordination number
    coordination_nr_vec=get_coordination_nr_vec(original_graph)
    println("coordination_nr_vec, $coordination_nr_vec")

    count_3 = count(x -> x == 3, coordination_nr_vec)
    count_4 = count(x -> x == 4, coordination_nr_vec)
    # this if else can be removed, if you are sure that the network is correct
    if cmp(network_name , "dia") == 0
        println("CN_all=$(length(coordination_nr_vec))=?=CN4=$(count_4)")
    elseif cmp(network_name , "srs") == 0
        println("CN_all=$(length(coordination_nr_vec))=?=CN3=$(count_3)")
    elseif cmp(network_name , "srd") == 0
        println("CN3/CN4=4/6=0.666=?=$(count_3/count_4)")
    elseif cmp(network_name , "ctn") == 0
        println("CN3/CN4=16/12=1.333=?=$(count_3/count_4)")
    end

    # create original spatial network
    original_spatial_network = Dict(
        "original_graph" => original_graph,
        "edge_length_vec" => edge_length_vec,
        "coordination_nr_vec" => coordination_nr_vec,
        "nr_vertices" => nr_vertices,
        "nr_dimensions" => nr_dimensions,
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

    original_spatial_network = get_network( #*#
        evolution_dict["nr_vertices"],
        evolution_dict["network_type"]
    ) 

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