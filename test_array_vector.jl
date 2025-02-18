# Initialize an empty vector to hold 3D arrays
next_neighbor_positions_arr = Vector{Array{Float64, 2}}()

# Example of populating the vector with 3D arrays of varying third dimension
for i in 1:4  # 4 vertices
    n = rand(2:7)  # Randomly choose a value between 2 and 7 for the third dimension
    arr = Array{Float64}(undef, 3, n)  # Create a 3D array with dimensions (3, n)
    
    # Populate the array with some values (e.g., random values)
    for j in 1:3
        for k in 1:n
            arr[j, k] = rand()
        end
    end
    
    # Push the array into the vector
    push!(next_neighbor_positions_arr, arr)
end

#=
# Print the vector of 3D arrays
for (i, arr) in enumerate(next_neighbor_positions_arr)
    println("Array for vertex $i:")
    println(arr)
end
=#
println(next_neighbor_positions_arr)