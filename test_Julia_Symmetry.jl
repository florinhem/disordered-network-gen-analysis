using Graphs
using GraphPlot3D
using Colors

# Example arrays of vertices and edges
vertices = 1:5
edges = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)]

# Create an empty graph with the specified number of vertices
g = SimpleGraph(length(vertices))

# Add edges to the graph
for (u, v) in edges
    add_edge!(g, u, v)
end

# Define the 3D positions of each vertex
positions = Dict(
    1 => (0.0, 0.0, 0.0),
    2 => (1.0, 0.0, 0.0), 
    3 => (1.0, 1.0, 0.0),
    4 => (0.0, 1.0, 0.0),
    5 => (0.5, 0.5, 1.0)
)

# Convert positions to a format suitable for gplot3d
x = [positions[v][1] for v in vertices]
y = [positions[v][2] for v in vertices]
z = [positions[v][3] for v in vertices]

# Plot the graph with the specified 3D positions
gplot3d(g, layout=(x, y, z), node_color=colorant"blue", node_labels=1:nv(g))



#=
using Graphs
using GraphPlot
using Colors

# Example arrays of vertices and edges
vertices = 1:5
edges = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)]

# Create an empty graph with the specified number of vertices
g = SimpleGraph(length(vertices))

# Add edges to the graph
for (u, v) in edges
    add_edge!(g, u, v)
end

# Define the positions of each vertex
positions = Dict(
    1 => (0.0, 0.0),
    2 => (1.0, 0.0),
    3 => (1.0, 1.0),
    4 => (0.0, 1.0),
    5 => (0.5, 0.5)
)

# Convert positions to a format suitable for gplot
x = [positions[v][1] for v in vertices]
y = [positions[v][2] for v in vertices]

# Plot the graph with the specified positions
gplot(g, layout=(x, y), node_color=colorant"blue", node_labels=1:nv(g))
=#



#=using Graphs
using GraphPlot

# Example arrays of vertices and edges
vertices = 1:5
edges = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 1)]

# Create an empty graph with the specified number of vertices
g = Graphs.Graph(length(vertices))

# Add edges to the graph
for (u, v) in edges
    add_edge!(g, u, v)
end

# Print the graph to verify
println(g)
gplot(g)
=#



#=
using LightGraphs
using GraphPlot
using Colors

# Create a simple graph with 2 points and a connection between them
g = SimpleGraph(2)
add_edge!(g, 1, 2)

# Function to apply rotation symmetry operation and add edges
function apply_rotation_symmetry!(g, angle)
    n = nv(g)  # number of vertices
    for i in 1:n
        for j in i+1:n
            if has_edge(g, i, j)
                # Example rotation operation: add edges to the rotated vertices
                new_i = mod(i + angle, n) + 1
                new_j = mod(j + angle, n) + 1
                add_edge!(g, new_i, new_j)
            end
        end
    end
end

# Function to apply translation symmetry operation and add edges
function apply_translation_symmetry!(g, translation)
    n = nv(g)  # number of vertices
    for i in 1:n
        for j in i+1:n
            if has_edge(g, i, j)
                # Example translation operation: add edges to the translated vertices
                new_i = mod(i + translation, n) + 1
                new_j = mod(j + translation, n) + 1
                add_edge!(g, new_i, new_j)
            end
        end
    end
end

# Apply rotation symmetry with angle 1
apply_rotation_symmetry!(g, 1)

# Apply translation symmetry with translation 1
apply_translation_symmetry!(g, 1)

# Plot the graph
gplot(g, layout=spring_layout, node_color=colorant"blue", node_labels=1:nv(g))
=#