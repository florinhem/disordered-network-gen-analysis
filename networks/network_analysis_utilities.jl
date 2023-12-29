"""
these functions provide utilities for analyzing networks
"""


"""
Average y(x) over different graph realizations. The y_arr contains
values for the same graph and different x along the 1 dimension
and for the same x and different graphs along the 2 dimension
"""
function get_multiple_graph_average_vec(
    y_arr::Array{Real})

    #get the number of x values, for which y(x) was evaluated
    nr_x_values = length(y_arr[:,1])

    #initialize vector for average values of y
    y_average_vec = Vector{Measurements.Measurement{Float64}}(undef, nr_x_values)

    #loop through index corresponding to different x values and average
    #over different realizations
    for x_index in 1:nr_x_values

        y_average_vec[x_index] = Measurements.measurement(
            Statistics.mean(y_arr[x_index,:]),
            Statistics.std(y_arr[x_index,:])
        )
    end

    return y_average_vec
end


"""
Average the quantity y over different graph realizations.
The y_vec contains this quantity for different graphs
"""
function get_multiple_graph_average(y_vec::Vector{Real})

    #calculate mean and standard deviation of quantity
    y_average = Measurements.measurement(
            Statistics.mean(y_vec),
            Statistics.std(y_vec)
        )

    return y_average
end