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