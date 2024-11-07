using Plots

# Sample data
x = 1:10
y = rand(10)
z = rand(10)

# Create a scatter plot with a colorbar at the top
scatter = Plots.scatter(x, y, z, zcolor=z, size=(400, 300), label="", colorbar=true)

# Adjust the layout to position the colorbar at the top
Plots.plot!(scatter, layout=(2, 1), size=(600, 400), colorbar=:left)