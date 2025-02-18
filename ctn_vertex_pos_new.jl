import Plots
Plots.plotlyjs()

function cubePlot()
    min=0
    max=1
    vertices = [
       (min, min, min),
       (max,min,min),
       (max,max,min),
       (min,max,min),
       (min,min,max),
       (max,min,max),
       (max,max,max),
       (min,max,max)
    ]

    # Define the edges of the cube
    edges = [
        (1, 2), (2, 3), (3, 4), (4, 1), # Bottom face
        (5, 6), (6, 7), (7, 8), (8, 5), # Top face
        (1, 5), (2, 6), (3, 7), (4, 8)  # Vertical edges
    ]

    # Extract x, y, z coordinates
    x = [v[1] for v in vertices]
    y = [v[2] for v in vertices]
    z = [v[3] for v in vertices]

    #println(x)

    # Plot the cube
    Plots.plot3d(x, y, z, seriestype = :scatter, markersize = 5, label = "Vertices")
    for (i, j) in edges
        Plots.plot3d!([x[i], x[j]], [y[i], y[j]], [z[i], z[j]], seriestype = :line, label = "")
    end

    display(Plots.plot3d)
end

# Point between
vertex_list = [
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
]

vertex_mod_array = []

for vertex in vertex_list
    vertex_mod = mod.(vertex, [1, 1, 1])
    push!(vertex_mod_array, vertex_mod)
end

vertex_mod_array=vertex_list
println("vertex_mod_array, $vertex_mod_array")

# Extract x, y, z coordinates
x = [v[1] for v in vertex_mod_array]
y = [v[2] for v in vertex_mod_array]
z = [v[3] for v in vertex_mod_array]

cubePlot()

# Create 3D plot
Plots.scatter3d!(x, y, z, markersize=1, marker=:sphere, label="Vertices")
a=0
b=1
Plots.xlims!(a,b)
Plots.ylims!(a,b)
Plots.zlims!(a,b)
Plots.xlabel!("X-axis")
Plots.ylabel!("Y-axis")
Plots.zlabel!("Z-axis")
Plots.title!("3D Plot of Vertex List")

# Save the plot as an image file
#Plots.savefig("vertex_list_plot.png")