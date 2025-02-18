import Plots
plotlyjs()

# Point between
vertex_list = [
    [7, 5, 1],
    [5, 3, 1],
    [3, 3, 3],
    [1, 5, 3],
    [1, 7, 5],
    [3, 1, 5],
    [5, 1, 7],
    [7, 7, 7],
]

vertex_list=vertex_list ./ 8

#=
[7/8, 5/8, 1/8],
    [5/8, 3/8, 1/8],
    [3/8, 3/8, 3/8],
    [1/8, 5/8, 3/8],
    [1/8, 7/8, 5/8],
    [3/8, 1/8, 5/8],
    [5/8, 1/8, 7/8],
    [7/8, 7/8, 7/8],
=#

#=[0.875, 0.625, 1.125],
    [0.625, 0.375, 1.125],
    [0.375, 0.375, 1.375],
    [0.125, 0.625, 1.375],
    [0.125, 0.875, 1.625],
    [0.375, 1.125, 1.625],
    [0.625, 1.125, 1.875],
    [0.875, 0.875, 1.875],=#

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