import Plots
Plots.plotlyjs()

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

#=
[6, 6, 2],
    [0, 6, 4],
    [0, 2, 4],
    [2, 2, 2],
    [4, 0, 2],
    [4, 0, 6],
    [2, 6, 6],
    [2, 4, 0],
    [6, 4, 0],
    [6, 2, 6],
    =#

#=
[3, 3, 1],
    [0, 3, 2],
    [0, 1, 2],
    [1, 1, 1],
    [2, 0, 1],
    [2, 0, 3],
    [1, 3, 3],
    [1, 2, 0],
    [3, 2, 0],
    [3, 1, 3],
    =#

#=
[0.75, 0.75, 1.25],
    [1, 0.75, 1.5],
    [1, 1.25, 1.5],
    [1.25, 1.25, 1.25],
    [1.5, 1, 1.25],
    [1.5, 1, 0.75],
    [1.25, 0.75, 0.75],
    [1.25, 1.5, 1],
    [0.75, 1.5, 1],
    [0.75, 1.25, 0.75],
    =#

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
Plots.scatter3d(x, y, z, markersize=5, marker=:sphere, label="Vertices")
Plots.xlims!(0,1)
Plots.ylims!(0,1)
Plots.zlims!(0,1)
Plots.xlabel!("X-axis")
Plots.ylabel!("Y-axis")
Plots.zlabel!("Z-axis")
Plots.title!("3D Plot of Vertex List")

# Save the plot as an image file
#Plots.savefig("vertex_list_plot.png")