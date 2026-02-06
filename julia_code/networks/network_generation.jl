"""
These functions generate spatial networks
"""

"""
returns the rotation matrix for a certain axis and certain angle
"""
function rotation_matrix(
    axis::Vector{Float64}, 
    angle::Float64)

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
            new_vertex_position_start=(rot_matrix * 
                (vertex_position_start .- shift) .+ shift .+ translate)
            new_vertex_position_end=(rot_matrix * 
                (vertex_position_end .- shift) .+ shift .+ translate)
            new_edges[current_edge]=(
                new_vertex_position_start,new_vertex_position_end,length)
            current_edge+=1
        end
    end

    return new_edges
end


"""
rotates all edges n_fold times around a certain axis
"""
function rotate(
    edges,
    axis::Vector{Float64},
    n_fold::Int64)

    new_edges=Dict{Int, Tuple{Vector{Float64}, Vector{Float64}, Float64}}()
    current_edge=1
    for n in 1:(n_fold-1)
        angle=360/n_fold*n
        rot_matrix=rotation_matrix(axis,angle)

        for (edge_nr,(vertex_position_start,vertex_position_end,length)) in edges
            new_vertex_position_start=(rot_matrix * vertex_position_start)
            new_vertex_position_end=(rot_matrix * vertex_position_end)
            new_edges[current_edge]=(
                new_vertex_position_start,new_vertex_position_end,length)
            current_edge+=1
        end
    end

    return new_edges
end


"""
makes a copy of edges and applies a glide plane 
(mirror at a plane and then translate along the mirror)
"""
function copy_and_glide(
    edges,
    translate::Vector{Float64},
    normal_vector_to_mirror_plane::Vector{Float64})

    # Normalize the normal vector
    normal_vector_to_mirror_plane = (normal_vector_to_mirror_plane / 
        LinearAlgebra.norm(normal_vector_to_mirror_plane))
    
    new_edges = deepcopy(edges)
    current_edge = length(edges) + 1

    for (edge_nr, (vertex_position_start, vertex_position_end, length)) in edges
        # Mirror the start and end positions
        mirrored_start = (vertex_position_start - 
            2 * LinearAlgebra.dot(vertex_position_start, 
                normal_vector_to_mirror_plane) *
            normal_vector_to_mirror_plane)
        mirrored_end = (vertex_position_end - 
            2 * LinearAlgebra.dot(vertex_position_end, 
                normal_vector_to_mirror_plane) * 
            normal_vector_to_mirror_plane)
        
        # Translate the mirrored positions
        new_vertex_position_start = mirrored_start .+ translate
        new_vertex_position_end = mirrored_end .+ translate
        
        # Add the new edge to the dictionary
        new_edges[current_edge] = (
            new_vertex_position_start, new_vertex_position_end, length)
        current_edge += 1
    end
    
    return new_edges
end


"""
makes a copy of edges and applies a point invertion 
(mirror all points to the opposite side of the invertion point)
"""
function copy_and_invert(
    edges,
    invertion_point::Vector{Float64})

    new_edges = deepcopy(edges)
    current_edge = length(edges) + 1

    for (edge_nr, (vertex_position_start, vertex_position_end, length)) in edges
        # Perform point inversion for start and end positions
        inverted_start = 2 .* invertion_point .- vertex_position_start
        inverted_end = 2 .* invertion_point .- vertex_position_end

        # Add the new inverted edge to the dictionary
        new_edges[current_edge] = (inverted_start, inverted_end, length)
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
        new_edges[current_edge]=(
            new_vertex_position_start,new_vertex_position_end,length)
        current_edge+=1
    end
    
    return new_edges
end


"""
makes a copy of edges and translates n times in the translation direction
"""
function copy_and_translate_n_times(
    edges, 
    translate::Vector{Float64},
    n)

    new_edges=deepcopy(edges)
    current_edge=length(edges)+1

    for i in 1:(n-1)
        for (edge_nr,(vertex_position_start,vertex_position_end,length)) in edges 
            new_vertex_position_start=vertex_position_start .+ translate*i
            new_vertex_position_end=vertex_position_end .+ translate*i
            new_edges[current_edge]=(
                new_vertex_position_start,new_vertex_position_end,length)
           
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

    # change the labels of the edges dictionary to ascending numbers from 1 to
    # length(edges)
    edges = Dict(1:length(edges) .=> collect(values(edges)))

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
function delete_copys(edges, epsilon)

    new_edges=deepcopy(edges)

    for (edge_nr_1,(vertex_position_start_1,vertex_position_end_1,length_1)) in edges
        for (edge_nr_2,(vertex_position_start_2,vertex_position_end_2,length_2)) in edges
            if edge_nr_1<edge_nr_2
                if (length_1-length_2)<epsilon
                    if ((LinearAlgebra.norm(
                        vertex_position_start_1-vertex_position_start_2)<epsilon 
                        &&
                       LinearAlgebra.norm(
                        vertex_position_end_1-vertex_position_end_2)<epsilon) 
                       ||
                       (LinearAlgebra.norm(
                        vertex_position_start_1-vertex_position_end_2)<epsilon 
                        &&
                       LinearAlgebra.norm(
                        vertex_position_end_1-vertex_position_start_2)<epsilon))

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
            
            if (distance_periodic_bc(
                vertex_position_start,vertex_position,L)<epsilon)

                add_start_to_vertex_positions_dict=false
            end

            if (distance_periodic_bc(
                vertex_position_end,vertex_position,L)<epsilon)

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
            if vertex_start===nothing && distance_periodic_bc(
                vertex_position_start, vertex_position, L)<epsilon

                vertex_start=vertex_nr
            end
            if vertex_end===nothing && distance_periodic_bc(
                vertex_position_end, vertex_position, L)<epsilon
                
                vertex_end=vertex_nr
            end
        end
        if vertex_start<vertex_end
            edges_with_vertex[edge_nr]=(
                vertex_position_start,vertex_position_end,length,
                vertex_start,vertex_end)
        else
            edges_with_vertex[edge_nr]=(
                vertex_position_end,vertex_position_start,length,
                vertex_end,vertex_start)
        end
    end

    return edges_with_vertex
end

"""
returns the original graph that has all vertices and all edges
"""
function create_graph(edges_with_vertex, nr_vertices)
    
    original_graph=Graphs.SimpleGraph(nr_vertices)

    for (edge_nr,(vertex_position_start,vertex_position_end,
        length,vertex_start,vertex_end)) in edges_with_vertex
        
        Graphs.add_edge!(original_graph, vertex_start, vertex_end)
    end
    
    return original_graph
end

"""
returns a vector containing the length of each edge
"""
function get_edge_length_vec(original_graph, edges_with_vertex)

    edge_length_vec = Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Float64}()

    for edge in Graphs.edges(original_graph)
        v1 = Graphs.src(edge)
        v2 = Graphs.dst(edge)

        for (edge_nr,(vertex_position_start,vertex_position_end,
            length,vertex_start,vertex_end)) in edges_with_vertex

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
        new_edges[current_edge]=(
            new_vertex_position_start,new_vertex_position_end,length)
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
that shifted all the edges (the center of the edge) into the unitcell 
"""
function shift_edge_middle_into_unitcell(edges, unitcell_min, unitcell_max)

    unitcell_length=unitcell_max-unitcell_min
    new_edges=Dict()
    current_edge=1

    for (edge_nr,(vertex_position_start,vertex_position_end,length)) in edges

        # Calculate how to translate the middle of the edge into the unitcell
        #=
        edge_middle_position=vertex_position_start .+ ( 
            (vertex_position_end .+ vertex_position_start) ./ 2 )
            =#
        edge_middle_position=(vertex_position_start .+ vertex_position_end) ./ 2
        new_edge_middle_position=(
            mod.(edge_middle_position, unitcell_length) .+ unitcell_min)
        translation=new_edge_middle_position .- edge_middle_position

        # Shift the starting and ending position of the edge
        new_vertex_position_start=vertex_position_start .+ translation
        new_vertex_position_end=vertex_position_end .+ translation

        new_edges[current_edge]=(
            new_vertex_position_start,new_vertex_position_end,length)
        current_edge+=1
    end

    return new_edges
end

"""
returns a dictionary with edges (starting position, ending position, length)
that has now scaled up all the points and length with the factor 
(usually the unitcell length) 
"""
function scale(edges, factor)
    new_edges=Dict()
    current_edge=1

    for (edge_nr,(vertex_position_start,vertex_position_end,length)) in edges
        new_edges[current_edge]=(vertex_position_start .* 
            factor,vertex_position_end .* factor,length * factor)
        current_edge+=1
    end

    return new_edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the diamond structure with the symmetry operations
"""
function get_edges_dia()
    # define the edges with the help of rcsr.net
    # look at the space group name: Fd-3m
    # and find the space group number: 227

    # define the position of vertex V1 and edge E1. 
    # V2 is twice the distance from V1 to E1. 
    V1=[1/8,1/8,1/8]
    E1=[0.0,0.0,0.0]

    # the difference between V1 and V2 gives us the length
    D1=2 .* (V1 .- E1)
    V2=V1 .- D1
    L1=LinearAlgebra.norm(D1)
    
    edges = Dict(1 => (V1, V2, L1))
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,
        [0.0,0.0,1.0],[3/8,1/8,0.0],[0.0,0.0,1/2],2)
    edges = copy_and_rotate_and_translate(edges,
        [0.0,1.0,0.0],[1/8,0.0,3/8],[0.0,1/2,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],3)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,0.0],[0.0,-1/4,1/4],[1/2,1/2,0.0],2)
    edges = copy_and_invert(edges,[0.0,0.0,0.0])

    edges = copy_and_translate(edges, [0.0,1/2,1/2])
    edges = copy_and_translate(edges, [1/2,0.0,1/2])
    edges = copy_and_translate(edges, [1/2,1/2,0.0])

    return edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the srs (gyroid) structure with the symmetry operations
"""
function get_edges_srs()
    # define the edges with the help of rcsr.net
    # look at the space group name: I4(1)32
    # and find the space group number: 214

    # define the position of vertex V1 and edge E1. 
    # V2 is twice the distance from V1 to E1. 
    V1=[1/8,1/8,1/8]
    E1=[0.0,1/4,1/8]

    # the difference between V1 and V2 gives us the length
    D1=2 .* (V1 .- E1)
    V2=V1 .- D1
    L1=LinearAlgebra.norm(D1)
    
    edges = Dict(1 => (V1, V2, L1))

    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,
        [0.0,0.0,1.0],[1/4,0.0,0.0],[0.0,0.0,1/2],2)
    edges = copy_and_rotate_and_translate(edges,
        [0.0,1.0,0.0],[0.0,0.0,1/4],[0.0,1/2,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],3)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,0.0],[0.0,-1/4,1/8],[1/2,1/2,0.0],2)

    edges = copy_and_translate(edges, [1/2,1/2,1/2])

    return edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the srd structure with the symmetry operations
"""
function get_edges_srd()
    # define the edges with the help of rcsr.net
    # look at the space group name: P4(2)32
    # and find the space group number: 213

    # define the position of vertex V1 and edge E1. 
    # V3 is twice the distance from V1 to E1. 
    V1=[0.0,1/4,1/2]
    E1=[0.0,1/2,1/2]

    # the difference between V1 and V2 gives us the length
    D1=2 .* (V1 .- E1)
    V3=V1 .- D1
    L1=LinearAlgebra.norm(D1)

    # define the position of vertex V2 and edge E2. 
    # V4 is twice the distance from V2 to E2. 
    x=1/8
    V2=[1/4,1/4,1/4]
    E2=[x,1/4,1/2-x]

    # the difference between V1 and V2 gives us the length
    D2=2 .* (V2 .- E2)
    V4=V2 .- D2
    L2=LinearAlgebra.norm(D2)

    edges = Dict(
        1 => (V1, V3, L1),
        2 => (V2, V4, L2)
        )
    
    println("edges, $edges")
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,
        [0.0,0.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [0.0,1.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],3)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,0.0],[0.0,0.0,1/4],[1/2,1/2,0.0],2)
    
    return edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the ctn structure with the symmetry operations
"""
function get_edges_ctn()
    # define the edges with the help of rcsr.net
    # look at the space group name: I-43d
    # and find the space group number: 220

    # it is possible to change the x value for the first vertex V1. V2 is fixed.
    x=0.2082
    V1=[x, x, x]
    V2=[3/8, 0.0, 1/4]

    # the difference between V1 and V2 gives us the length
    D1=V1 .- V2
    L1=LinearAlgebra.norm(D1)
    
    edges = Dict(1 => (V1, V2, L1))
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,
        [0.0,0.0,1.0],[1/4,0.0,0.0],[0.0,0.0,1/2],2)
    edges = copy_and_rotate_and_translate(edges,
        [0.0,1.0,0.0],[0.0,0.0,1/4],[0.0,1/2,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],3)
    edges = copy_and_glide(edges,[1/4,1/4,1/4],[1.0,-1.0,0.0])
    
    edges = copy_and_translate(edges, [1/2,1/2,1/2])

    return edges
end

"""
returns a dictionary with edges (starting position, ending position, length)
for the pto structure with the symmetry operations
"""
function get_edges_pto()
    # define the edges with the help of rcsr.net
    # look at the space group name: Pm-3n
    # and find the space group number: 223

    # V1 & V2 are fixed.
    V1=[1/4, 1/4, 1/4]
    V2=[1/4, 0.0, 1/2]

    # the difference between V1 and V2 gives us the length
    D1=V1 .- V2
    L1=LinearAlgebra.norm(D1)
    
    edges = Dict(1 => (V1, V2, L1))
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,
        [0.0,0.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [0.0,1.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],3)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,0.0],[0.0,0.0,1/4],[1/2,1/2,0.0],2)
    edges = copy_and_invert(edges, [0.0,0.0,0.0])
    
    return edges
end

function get_edges_pto_half_fill()
    # define the edges with the help of rcsr.net
    # look at the space group name: Pm-3n
    # and find the space group number: 223

    # V1 & V2 are fixed.
    V1=[1/4, 1/4, 1/4]
    V2=[1/4, 0.0, 1/2]

    # the difference between V1 and V2 gives us the length
    D1=V1 .- V2
    L1=LinearAlgebra.norm(D1)
    
    edges = Dict(1 => (V1, V2, L1))
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,
        [0.0,0.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [0.0,1.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],3)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,0.0],[0.0,0.0,1/4],[1/2,1/2,0.0],2)
    edges = copy_and_invert(edges, [0.0,0.0,0.0])

    edges[length(edges)+1] = ([1/2, 1/2, 1/2], [1/4, 1/4, 1/4], sqrt(3)/4)
    edges[length(edges)+1] = ([1/2, 1/2, 1/2], [1/4, 3/4, 3/4], sqrt(3)/4)
    edges[length(edges)+1] = ([1/2, 1/2, 1/2], [3/4, 1/4, 3/4], sqrt(3)/4)
    edges[length(edges)+1] = ([1/2, 1/2, 1/2], [3/4, 3/4, 1/4], sqrt(3)/4)
    
    return edges
end


function get_edges_pto_fill()
    # define the edges with the help of rcsr.net
    # look at the space group name: Pm-3n
    # and find the space group number: 223

    # V1 & V2 are fixed.
    V1=[1/4, 1/4, 1/4]
    V2=[1/4, 0.0, 1/2]

    # the difference between V1 and V2 gives us the length
    D1=V1 .- V2
    L1=LinearAlgebra.norm(D1)
    
    edges = Dict(1 => (V1, V2, L1))
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,
        [0.0,0.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [0.0,1.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],3)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,0.0],[0.0,0.0,1/4],[1/2,1/2,0.0],2)
    edges = copy_and_invert(edges, [0.0,0.0,0.0])

    edges[length(edges)+1] = ([1/2, 1/2, 1/2], [1/4, 1/4, 1/4], sqrt(3)/4)
    edges[length(edges)+1] = ([1/2, 1/2, 1/2], [1/4, 3/4, 3/4], sqrt(3)/4)
    edges[length(edges)+1] = ([1/2, 1/2, 1/2], [3/4, 1/4, 3/4], sqrt(3)/4)
    edges[length(edges)+1] = ([1/2, 1/2, 1/2], [3/4, 3/4, 1/4], sqrt(3)/4)

    edges[length(edges)+1] = ([0.9999,0.9999,0.9999], [3/4, 3/4, 3/4], sqrt(3)/4)
    edges[length(edges)+1] = ([0.9999,0.0,0.0], [3/4, 1/4, 1/4], sqrt(3)/4)
    edges[length(edges)+1] = ([0.0,0.9999,0.0], [1/4, 3/4, 1/4], sqrt(3)/4)
    edges[length(edges)+1] = ([0.0,0.0,0.9999], [1/4, 1/4, 3/4], sqrt(3)/4)
    
    return edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the lcs structure with the symmetry operations
"""
function get_edges_lcs()
    # define the edges with the help of rcsr.net
    # look at the space group name: Ia-3d
    # and find the space group number: 230

    # it is possible to change the x value for the first edge E1. V1 is fixed.
    x=7/16

    # we change E1=[x, 7/8, 3/4-x] to E1=[x, -1/8, 3/4-x] 
    # to have simpler (closer to origin (0,0,0)) values
    E1=[x, -1/8, 3/4-x] 
    V1=[3/8, 0.0, 1/4]

    # calculate the V2 position
    V2 = 2*E1 .- V1

    # the difference between V1 and V2 gives us the length
    D1=V1 .- V2
    L1=LinearAlgebra.norm(D1)

    edges = Dict(1 => (V1, V2, L1))
    
    # symmetry operations for space group number with the help of the book:
    # "International Tables for Crystallography"
    edges = copy_and_rotate_and_translate(edges,
        [0.0,0.0,1.0],[1/4,0.0,0.0],[0.0,0.0,1/2],2)
    edges = copy_and_rotate_and_translate(edges,
        [0.0,1.0,0.0],[0.0,0.0,1/4],[0.0,1/2,0.0],2)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,1.0],[0.0,0.0,0.0],[0.0,0.0,0.0],3)
    edges = copy_and_rotate_and_translate(edges,
        [1.0,1.0,0.0],[0.0,-1/4,1/8],[1/2,1/2,0.0],2)
    edges = copy_and_invert(edges, [0.0,0.0,0.0])
    edges = copy_and_translate(edges, [1/2,1/2,1/2])

    return edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the bcu_cn_5_6_7_8 structure which is the bcu network with bonds removed
to achieve coordination numbers of 5, 6, 7, and 8 as in the disordered networks
of 10.1002/adfm.202302720. Here, we place the edges by hand, because the
network has no apparent symmetry.
"""
function get_edges_bcu_cn_5_6_7_8()

    L1 = 1/(2*sqrt(3/2))

    ax = [1/8, 0, 0]
    ay = [0, 1/8, 0]
    az = [0, 0, 1/8]

    edges = Dict(
        1 => (ax .+ ay .+ az, 3 .* (ax .+ ay .+ az), L1),
        2 => ( 5 .* ax .+ ay .+ az, 3 .* (ax .+ ay .+ az), L1),
        3 => ( 5 .* ax .+ ay .+ az, 7 .* ax .+ 3 .* (ay .+ az), L1),
        4 => ( 9 .* ax .+ ay .+ az, 7 .* ax .+ 3 .* (ay .+ az), L1),

        5 => (ax .+ 5 .* ay .+ az, 3 .* (ax .+ ay .+ az), L1),
        6 => (5 .* ax .+ 5 .* ay .+ az, 3 .* (ax .+ ay .+ az), L1),
        7 => (5 .* ax .+ 5 .* ay .+ az, 7 .* ax .+ 3 .* (ay .+ az), L1),

        8 => (ax .+ 5 .* ay .+ az, 3 .* ax .+ 7 .* ay .+ 3 .* az, L1),
        9 => (5 .* ax .+ 5 .* ay .+ az, 3 .* ax .+ 7 .* ay .+ 3 .* az, L1),
        10 => (9 .* ax .+ 5 .* ay .+ az,  7 .* ax .+ 7 .* ay .+ 3 .* az, L1),

        11 => (ax .+ 9 .* ay .+ az, 3 .* ax .+ 7 .* ay .+ 3 .* az, L1),
        12 => (5 .* ax .+ 9 .* ay .+ az, 3 .* ax .+ 7 .* ay .+ 3 .* az, L1),
        13 => (5 .* ax .+ 9 .* ay .+ az, 7 .* ax .+ 7 .* ay .+ 3 .* az, L1),
        14 => (9 .* ax .+ 9 .* ay .+ az, 7 .* ax .+ 7 .* ay .+ 3 .* az, L1),

        15 => (ax .+ ay .+ 5 .* az, 3 .* (ax .+ ay .+ az), L1),
        16 => (5 .* ax .+ ay .+ 5 .* az, 3 .* (ax .+ ay .+ az), L1),
        17 => (9 .* ax .+ ay .+ 5 .* az, 7 .* ax .+ 3 .* (ay .+ az), L1),

        18 => (ax .+ 5 .* ay .+ 5 .* az, 3 .* (ax .+ ay .+ az), L1),
        19 => (5 .* ax .+ 5 .* ay .+ 5 .* az, 3 .* (ax .+ ay .+ az), L1),
        20 => (5 .* ax .+ 5 .* ay .+ 5 .* az, 7 .* ax .+ 3 .* (ay .+ az), L1),
        21 => (9 .* ax .+ 5 .* ay .+ 5 .* az, 7 .* ax .+ 3 .* (ay .+ az), L1),

        22 => (3 .* ax .+ 7 .* ay .+ 3 .* az, 5 .* (ax .+ ay .+ az), L1),
        23 => (7 .* ax .+ 7 .* ay .+ 3 .* az, 5 .* (ax .+ ay .+ az), L1),
        24 => (7 .* ax .+ 7 .* ay .+ 3 .* az, 9 .* ax .+ 5 .* (ay .+ az), L1),

        25 => (3 .* ax .+ 7 .* ay .+ 3 .* az, 5 .* ax .+ 9 .* ay .+ 5 .* az, L1),
        26 => (7 .* ax .+ 7 .* ay .+ 3 .* az, 5 .* ax .+ 9 .* ay .+ 5 .* az, L1),
        27 => (7 .* ax .+ 7 .* ay .+ 3 .* az, 9 .* ax .+ 9 .* ay .+ 5 .* az, L1),

        28 => (ax .+ ay .+ 5 .* az, 3 .* ax .+ 3 .* ay .+ 7 .* az, L1),
        29 => (5 .* ax .+ ay .+ 5 .* az, 3 .* ax .+ 3 .* ay .+ 7 .* az, L1),
        30 => (5 .* ax .+ ay .+ 5 .* az, 7 .* ax .+ 3 .* ay .+ 7 .* az, L1),
        31 => (9 .* ax .+ ay .+ 5 .* az, 7 .* ax .+ 3 .* ay .+ 7 .* az, L1),

        32 => (ax .+ 5 .* ay .+ 5 .* az, 3 .* ax .+ 3 .* ay .+ 7 .* az, L1),
        33 => (5 .* ax .+ 5 .* ay .+ 5 .* az, 3 .* ax .+ 3 .* ay .+ 7 .* az, L1),
        34 => (5 .* ax .+ 5 .* ay .+ 5 .* az, 7 .* ax .+ 3 .* ay .+ 7 .* az, L1),
        35 => (9 .* ax .+ 5 .* ay .+ 5 .* az, 7 .* ax .+ 3 .* ay .+ 7 .* az, L1),

        36 => (ax .+ 5 .* ay .+ 5 .* az, 3 .* ax .+ 7 .* ay .+ 7 .* az, L1),
        37 => (5 .* ax .+ 5 .* ay .+ 5 .* az, 7 .* ax .+ 7 .* ay .+ 7 .* az, L1),

        38 => (ax .+ 9 .* ay .+ 5 .* az, 3 .* ax .+ 7 .* ay .+ 7 .* az, L1),
        39 => (5 .* ax .+ 9 .* ay .+ 5 .* az, 3 .* ax .+ 7 .* ay .+ 7 .* az, L1),
        40 => (9 .* ax .+ 9 .* ay .+ 5 .* az, 7 .* ax .+ 7 .* ay .+ 7 .* az, L1),

        41 => (3 .* ax .+ 3 .* ay .+ 7 .* az, ax .+ ay .+ 9 .* az, L1),
        42 => (7 .* ax .+ 3 .* ay .+ 7 .* az, 9 .* ax .+ ay .+ 9 .* az, L1),

        43 => (3 .* ax .+ 3 .* ay .+ 7 .* az, 5 .* ax .+ 5 .* ay .+ 9 .* az, L1),
        44 => (7 .* ax .+ 3 .* ay .+ 7 .* az, 5 .* ax .+ 5 .* ay .+ 9 .* az, L1),
        45 => (7 .* ax .+ 3 .* ay .+ 7 .* az, 9 .* ax .+ 5 .* ay .+ 9 .* az, L1),

        46 => (3 .* ax .+ 7 .* ay .+ 7 .* az, 1 .* ax .+ 5 .* ay .+ 9 .* az, L1),
        47 => (3 .* ax .+ 7 .* ay .+ 7 .* az, 5 .* ax .+ 5 .* ay .+ 9 .* az, L1),
        48 => (7 .* ax .+ 7 .* ay .+ 7 .* az, 5 .* ax .+ 5 .* ay .+ 9 .* az, L1),
        49 => (7 .* ax .+ 7 .* ay .+ 7 .* az, 9 .* ax .+ 5 .* ay .+ 9 .* az, L1),

        50 => (3 .* ax .+ 7 .* ay .+ 7 .* az, 1 .* ax .+ 9 .* ay .+ 9 .* az, L1),
        51 => (3 .* ax .+ 7 .* ay .+ 7 .* az, 5 .* ax .+ 9 .* ay .+ 9 .* az, L1),
        52 => (7 .* ax .+ 7 .* ay .+ 7 .* az, 5 .* ax .+ 9 .* ay .+ 9 .* az, L1),
    )
    
    return edges
end


"""
returns a dictionary with edges (starting position, ending position, length)
for the pcu_cn_4_5_6 structure which is the pcu network with bonds removed
to achieve coordination numbers of 4, 5, and 6 as in the disordered network
of 10.1016/j.mtadv.2024.100524. Here, we place the edges by hand, because the
network has no apparent symmetry.
"""
function get_edges_pcu_cn_4_5_6()

    L1 = 1/2

    ax = [1/4, 0, 0]
    ay = [0, 1/4, 0]
    az = [0, 0, 1/4]

    edges = Dict(
        1 => (ax .+ ay .+ az, ax .+ ay .+ 3 .* az, L1),
        2 => (ax .+ ay .+ az, ax .+ ay .- az, L1),
        3 => (ax .+ ay .+ az, 3 .* ax .+ ay .+ az, L1),
        4 => (ax .+ ay .+ az, ax .+ 3 .* ay .+ az, L1),
        5 => (ax .+ ay .+ az, ax .- ay .+ az, L1),

        6 => (3 .* ax .+ ay .+ az, 3 .* ax .+ ay .+ 3 .* az, L1),
        7 => (3 .* ax .+ ay .+ az, 3 .* ax .+ ay .- az, L1),
        8 => (3 .* ax .+ ay .+ az, 3 .* ax .+ 3 .* ay .+ az, L1),
        9 => (3 .* ax .+ ay .+ az, 3 .* ax .- ay .+ az, L1),

        10 => (ax .+ 3 .* ay .+ az, 3 .* ax .+ 3 .* ay .+ az, L1),
        11 => (ax .+ 3 .* ay .+ az, .- ax .+ 3 .* ay .+ az, L1),

        12 => (3 .* ax .+ 3 .* ay .+ az, 3 .* ax .+ 3 .* ay .+ 3 .* az, L1),
        13 => (3 .* ax .+ 3 .* ay .+ az, 3 .* ax .+ 3 .* ay .- az, L1),

        14 => (ax .+ ay .+ 3 .* az, ax .+ 3 .* ay .+ 3 .* az, L1),
        15 => (ax .+ ay .+ 3 .* az, ax .- ay .+ 3 .* az, L1),
        16 => (ax .+ ay .+ 3 .* az, .- ax .+ ay .+ 3 .* az, L1),

        17 => (3 .* ax .+ ay .+ 3 .* az, 3 .* ax .+ 3 .* ay .+ 3 .* az, L1),
        18 => (3 .* ax .+ ay .+ 3 .* az, 3 .* ax .- ay .+ 3 .* az, L1),

        19 => (ax .+ 3 .* ay .+ 3 .* az, 3 .* ax .+ 3 .* ay .+ 3 .* az, L1),
        20 => (ax .+ 3 .* ay .+ 3 .* az, .- ax .+ 3 .* ay .+ 3 .* az, L1),
    )

    return edges
end


"""
returns a vector for the coordination numbers of all vertices of the 
original graph
"""
function get_coordination_nr_vec(original_graph)
    coordination_nr_vec::Vector{Int64}=fill(-1,length(
        Graphs.vertices(original_graph)))

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
get the original graph from the edges by copying and shifting them into one
unit cell. Then, we scale the network up from a unitcell lenght of 1 to 
edge_length_unit_cell. Finally, we create the original spatial network as a 
dictionary with all these values.
"""
function get_original_graph_from_edges(
    edges::Dict{Int, Tuple{Vector{Float64}, Vector{Float64}, Float64}},
    edge_length_unit_cell,
    nr_vertices::Int64,
    nr_vertices_per_unit_cell::Int64,
    nr_dimensions::Int64 = 3,
    delta::Float64 = 0.01,
    epsilon::Float64 = 0.001,
    unitcell_min=0,
    unitcell_max=1)

    # calculate the actual nr vertices, given that we require a 
    # cubic supercell and using the fact that the unit cell contains a certain
    # number of vertices 
    nr_unit_cells_per_dimension = max(
        1, Int(round( (nr_vertices/nr_vertices_per_unit_cell)^(1/3) )) )
    nr_vertices = nr_vertices_per_unit_cell * nr_unit_cells_per_dimension^3
    supercell_edge_length = nr_unit_cells_per_dimension*edge_length_unit_cell

    # we delete the edges that are repeating
    edges = delete_copys(edges, epsilon)

    # we shift all edges into the unitcell. We look at the center point of 
    # each vertex
    edges=shift_edge_middle_into_unitcell(edges, unitcell_min, unitcell_max)
    edges = delete_copys(edges, epsilon)

    # define the size of the box such that we can connect all edges into one
    # big network
    L=5
    x=[1.0,0.0,0.0]
    y=[0.0,1.0,0.0]
    z=[0.0,0.0,1.0]

    # copy all edges L times in x,y,z direction to get all edges that connect
    # to all the other edges
    edges = array_3D(edges,x,L,y,L,z,L)

    # look at the cube spanning from (2,2,2) and (3,3,3), only take these edges
    filter_min=(2.0-delta)
    filter_max=(3.0+delta)
    edges = filter_to_unitcell(edges, 
        filter_min, filter_max, 
        filter_min, filter_max, 
        filter_min, filter_max)

    # copy the unitcell in all 3 dimensions to get the actual network
    edges = array_3D(edges,
        x,nr_unit_cells_per_dimension,
        y,nr_unit_cells_per_dimension,
        z,nr_unit_cells_per_dimension)

    # fold all edges back into the supercell
    edges = fold_to_block(edges, nr_unit_cells_per_dimension)
    
    edges = delete_copys(edges, epsilon)
    
    # scale the network up from a unitcell lenght of 1 to edge_length_unit_cell
    edges=scale(edges,edge_length_unit_cell)

    # get the positions, edges, nr_vertices, original_graph and edges lengths
    vertex_positions_dict=get_vertex_positions_dict(edges, 
        epsilon, supercell_edge_length)
    vertex_position_mat=get_mat_from_dict(vertex_positions_dict)
    edges_with_vertex=get_edges_with_vertex(edges, vertex_positions_dict, 
        epsilon, supercell_edge_length)
    nr_vertices=length(vertex_positions_dict)
    original_graph=create_graph(edges_with_vertex, nr_vertices)

    # metrics to check that the network is correct
    edge_length_vec=get_edge_length_vec(original_graph, edges_with_vertex)
    coordination_nr_vec=get_coordination_nr_vec(original_graph)

    # create original spatial network
    original_graph = Dict(
        "original_graph" => original_graph,
        "edge_length_vec" => edge_length_vec,
        "coordination_nr_vec" => coordination_nr_vec,
        "nr_vertices" => nr_vertices,
        "nr_dimensions" => nr_dimensions,
        "supercell_edge_length" => supercell_edge_length,
        "vertex_position_mat" => vertex_position_mat
        )

    return original_graph
end


"""
returns the original space graph. It checks which network we want to generate.
Then it calculates for the unitcell and supercell the edges 
(number, starting vertex, ending vertex and length) of the network.
"""
function get_network(nr_vertices, network_name)
    if cmp(network_name , "dia") == 0 || cmp(network_name , "diamond") == 0
        nr_dimensions = 3
        edge_length_unit_cell = 4/sqrt(3)
        nr_vertices_per_unit_cell = 8
        nr_edges_per_unit_cell=16
        edges = get_edges_dia()
    elseif cmp(network_name , "srs") == 0
        nr_dimensions = 3
        edge_length_unit_cell = 4/sqrt(2)
        nr_vertices_per_unit_cell = 8
        nr_edges_per_unit_cell=12
        edges = get_edges_srs()
    elseif cmp(network_name , "srd") == 0
        nr_dimensions = 3
        # we calculated numerically the size of the unitcell,
        # such that the bond length energy is minimal. This value is close
        # but not exactly the same as the weighted lengths of the different
        # coordination number 3 and 4.
        edge_length_unit_cell = 2.3075 
        nr_vertices_per_unit_cell = 10
        nr_edges_per_unit_cell=18
        edges = get_edges_srd()
    elseif cmp(network_name , "ctn") == 0
        nr_dimensions = 3
        edge_length_unit_cell = 3.7033
        nr_vertices_per_unit_cell = 28
        nr_edges_per_unit_cell=48
        edges = get_edges_ctn()
    elseif cmp(network_name , "pto") == 0
        nr_dimensions = 3
        edge_length_unit_cell = 2.8284
        nr_vertices_per_unit_cell = 14
        nr_edges_per_unit_cell=24
        edges = get_edges_pto()
    elseif cmp(network_name , "pto_half_fill") == 0
        nr_dimensions = 3
        edge_length_unit_cell = 2.8284
        nr_vertices_per_unit_cell = 15
        nr_edges_per_unit_cell=28
        edges = get_edges_pto_half_fill()
    elseif cmp(network_name , "pto_fill") == 0
        nr_dimensions = 3
        edge_length_unit_cell = 2.8284
        nr_vertices_per_unit_cell = 16
        nr_edges_per_unit_cell=32
        edges = get_edges_pto_half_fill()
    elseif cmp(network_name , "lcs") == 0
        nr_dimensions = 3
        edge_length_unit_cell = 3.2660
        nr_vertices_per_unit_cell = 24
        nr_edges_per_unit_cell=48
        edges = get_edges_lcs()
    elseif cmp(network_name , "bcu_cn_5_6_7_8") == 0
        nr_dimensions = 3
        edge_length_unit_cell = 2*sqrt(3/2)
        nr_vertices_per_unit_cell = 16
        nr_edges_per_unit_cell=52
        edges = get_edges_bcu_cn_5_6_7_8()
    elseif cmp(network_name , "pcu_cn_4_5_6") == 0
        nr_dimensions = 3
        edge_length_unit_cell = 2
        nr_vertices_per_unit_cell = 8
        nr_edges_per_unit_cell=20
        edges = get_edges_pcu_cn_4_5_6()
    else
        @error ("$network_name not implemented.")
    end

    original_graph = get_original_graph_from_edges(
        edges,
        edge_length_unit_cell,
        nr_vertices,
        nr_vertices_per_unit_cell,
        nr_dimensions)
    
    return original_graph
end


"""
add information about vertex positions and edge vectors to the original graph
"""
function convert_original_graph_to_spatial_network(
    original_graph::Dict)

    # create an empty network graph where vertex positions and edge vectors
    # will be stored
    spatial_network = MetaGraphsNext.MetaGraph(
        Graphs.Graph(); 
        label_type = Int64,
        vertex_data_type = Dict{String, Any},
        edge_data_type = Dict{String, Union{Float64, Vector{Float64}}},
        graph_data = Dict{String, Any}(
            "nr_vertices" => original_graph["nr_vertices"],
            "nr_dimensions" => original_graph["nr_dimensions"],
            "supercell_edge_length" 
                => original_graph["supercell_edge_length"])
        )

    # label each vertex by its code integer and assign it its position vector
    for vertex in Graphs.vertices(original_graph["original_graph"])
        
        spatial_network[vertex] = Dict{String, Any}( 
            "position" => original_graph[
                "vertex_position_mat"][:,vertex],
            "coordination_nr" => original_graph[
                "coordination_nr_vec"][vertex]
        )

    end

    # label each edge by the vertices it connects and assign to it the vector
    # where it points

    # get the nr and vector of original edges
    nr_edges = Graphs.ne(original_graph["original_graph"])
    original_edges_vec = collect(
        Graphs.edges(original_graph["original_graph"]))
    
    for edge_nr in 1:nr_edges

        source = Graphs.src(original_edges_vec[edge_nr])
        target = Graphs.dst(original_edges_vec[edge_nr])

        # calculate vector from source to target considering periodic boundary
        # conditions
        edge_vector = get_distance_vector_pbc(
            original_graph["vertex_position_mat"][:,source],
            original_graph["vertex_position_mat"][:,target],
            original_graph["supercell_edge_length"] )

        spatial_network[source, target] =  Dict("vector" => edge_vector, 
        "distance_squared" => (
            original_graph["edge_length_vec"][
                original_edges_vec[edge_nr]] )^2 )

    end

    return spatial_network
end


"""
create a network graph representing the given network structure
"""
function get_periodic_network(evolution_dict)

    original_graph = get_network( #*#
        evolution_dict["nr_vertices"],
        evolution_dict["network_type"]
    ) 

    # convert original graph into a network graph that contains positional
    # information
    spatial_network = convert_original_graph_to_spatial_network(
        original_graph)

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
    if (cmp(evolution_dict["network_type"] , "dia") == 0 
        || cmp(evolution_dict["network_type"], "diamond") == 0)

        if evolution_dict["nr_dimensions"] == 3

            original_graph = get_diamond_network(
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
            "nr_vertices" => original_graph["nr_vertices"],
            "nr_dimensions" => original_graph["nr_dimensions"],
            "supercell_edge_length" 
                => original_graph["supercell_edge_length"])
        )

    # label each vertex by its code integer and assign it its position vector
    for vertex in Graphs.vertices(original_graph["original_graph"])

        spatial_network[vertex] = Dict( 
            "position" => rand(Float64, (3)) 
                .* original_graph["supercell_edge_length"],
            "coordination_nr" => original_graph[
                "coordination_nr_vec"][vertex]
        )

    end

    # label each edge by the vertices it connects and assign to it the vector
    # where it points

    nr_edges = Graphs.ne(original_graph["original_graph"])
    original_edges_vec = collect(
        Graphs.edges(original_graph["original_graph"]))
    
    # loop through orinal edges to get edge descriptions
    for edge_nr in 1:nr_edges

        # get source and target of edge
        source = Graphs.src(original_edges_vec[edge_nr])
        target = Graphs.dst(original_edges_vec[edge_nr])

        edge_vector = get_distance_vector_pbc(
            spatial_network[source]["position"],
            spatial_network[target]["position"],
            original_graph["supercell_edge_length"] )

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
