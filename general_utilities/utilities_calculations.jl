"""
Include calculation utilities
"""

"""
Apply gaussian filter of given sigma (standard deviation) to 1D data that is 
non equidistantly spaced. Calculation is inspired by
https://stackoverflow.com/questions/24143320/gaussian-sum-filter-for-irregular-spaced-points
"""
function gaussian_filter_1d(
    data_x::Vector{Float64}, 
    data_y::Vector; 
    sigma_x::Float64 = 2*pi/10, 
    filtered_data_x_step_length::Float64 = 2*pi/10)
    
    # create a vector of new sampling points
    filtered_data_x = collect(minimum(data_x)+sigma_x
        :filtered_data_x_step_length
        :maximum(data_x)-sigma_x) 

    # create a vector to store the filtered data
    filtered_data_y = zeros(typeof(data_y[1]), length(filtered_data_x))
    
    # loop through each data point
    for i in eachindex(filtered_data_x)

        # get vector of distances to the data points
        delta_x = data_x .- filtered_data_x[i]

        # determine data window around the current data point
        # that is used for the calculation of the filtered data point
        window_lower_index = findfirst(x -> x > -3*sigma_x, delta_x)
        window_upper_index = findlast(x -> x < 3*sigma_x, delta_x)

        # calculate unnormalized weigths for each data point
        weights = exp.(-0.5 * delta_x[window_lower_index:window_upper_index] 
            .^ 2 / sigma_x^2)

        # normalize the weights
        weights = weights ./ sum(weights)

        # calculate the filtered data point
        filtered_data_y[i] = LinearAlgebra.dot(
            data_y[window_lower_index:window_upper_index], weights)
        
    end
    
    return filtered_data_x, filtered_data_y
end