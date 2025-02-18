import Plots
import MetaGraphsNext
import LightGraphs

# Point between
vertex_list = [
    [7, 7, 3],
    [1, 7, 5],
    [1, 3, 5],
    [3, 3, 3],
    [5, 1, 3],
    [5, 1, 7],
    [3, 7, 7],
    [3, 5, 1],
    [7, 5, 1],
    [7, 3, 7],
]

vertex_list=vertex_list ./ 8

println("vertex_list, $vertex_list")

vertex_mod_array = []

for vertex in vertex_list
    vertex_mod = mod.(vertex, [1, 1, 1])
    push!(vertex_mod_array, vertex_mod)
end

println("vertex_mod_array, $vertex_mod_array")

# Extract x, y, z coordinates
x = [v[1] for v in vertex_mod_array]
y = [v[2] for v in vertex_mod_array]
z = [v[3] for v in vertex_mod_array]

# Create 3D plot
Plots.plotlyjs()

Plots.scatter3d(x, y, z, markersize=5, marker=:sphere, label="Vertices")
Plots.xlims!(0,1)
Plots.ylims!(0,1)
Plots.zlims!(0,1)
Plots.xlabel!("X-axis")
Plots.ylabel!("Y-axis")
Plots.zlabel!("Z-axis")
Plots.title!("3D Plot of Vertex List")






nr_vertices=10
vertex_position_mat=vertex_mod_array
original_graph = MetaGraphsNext.MetaGraph(LightGraphs.SimpleGraph(nr_vertices); label_type=Int)

# add edges based on custom logic
for i in 1:nr_vertices
    for j in i+1:nr_vertices
        # Add your custom logic to decide if an edge should be added
        # For example, you can add an edge if the distance between vertices is less than a certain value
        if norm(vertex_position_mat[:, i] - vertex_position_mat[:, j]) < 0.6
            MetaGraphsNext.add_edge!(original_graph, i, j)
        end

        #coordination_nr_vec
    end
end

edge_length_vec = [norm(vertex_position_mat[:, e.src] - vertex_position_mat[:, e.dst]) for e in edges(original_graph)]


#Plot

# create a line plot for each bond
figure = Plots.plot()
for bond in MetaGraphsNext.edge_labels(original_graph)
    pos_1 = original_graph[bond[1]]
    pos_2 = original_graph[bond[2]]
    Plots.plot!([pos_1[1], pos_2[1]], 
    [pos_1[2], pos_2[2]], 
    [pos_1[3], pos_2[3]], 
    type="scatter3d", mode="lines", color="black", showlegend=false)
end

Plots.gr()
return figure