


function main2()
    # For different T_max and beta, you would need to fit polynomials for each combination
    # Example data (replace with your actual data)
    T_max_values = [0.1, 0.2]
    beta_values = [0.285, 0.3]
    theta_values = [110, 130, 150, 170, 180]
    
    bond_angle_std_values = [
        [0.1, 0.2, 0.15, 0.25, 0.3],  # T_max=0.1, beta=0.285
        [0.12, 0.22, 0.18, 0.28, 0.32],  # T_max=0.1, beta=0.3
        [0.11, 0.21, 0.16, 0.26, 0.31],  # T_max=0.2, beta=0.285
        [0.13, 0.23, 0.19, 0.29, 0.33]   # T_max=0.2, beta=0.3
    ]
    bond_length_std_values = [
        [0.05, 0.1, 0.08, 0.12, 0.15],  # T_max=0.1, beta=0.285
        [0.06, 0.11, 0.09, 0.13, 0.16],  # T_max=0.1, beta=0.3
        [0.055, 0.105, 0.085, 0.125, 0.155],  # T_max=0.2, beta=0.285
        [0.065, 0.115, 0.095, 0.135, 0.165]   # T_max=0.2, beta=0.3
    ]

    # Combine bond_angle_std_values and bond_length_std_values into a single array
    combined_data_3d = [hcat(bond_angle_std_values[i], bond_length_std_values[i]) for i in 1:length(bond_angle_std_values)]
    println("combined_data_3d, $combined_data_3d")

    # Fit polynomials for each combination
    poly_combined_3d = [Polynomials.fit(theta_values, combined_data_3d[i][:, j], 2) for i in 1:length(combined_data_3d) for j in 1:2]
    println("poly_combined_3d, $poly_combined_3d")

    # Interpolated values for different T_max and beta
    T_max_interp = 0.15
    beta_interp = 0.29
    theta_interp = 140

    # Find the closest indices for T_max and beta
    T_max_idx = findall(x -> x >= T_max_interp, T_max_values)[1]
    beta_idx = findall(x -> x >= beta_interp, beta_values)[1]
    index = (T_max_idx - 1) * length(beta_values) + beta_idx

    bond_angle_std_interp_3d = poly_combined_3d[(index - 1) * 2 + 1](theta_interp)
    bond_length_std_interp_3d = poly_combined_3d[(index - 1) * 2 + 2](theta_interp)

    println("bond_angle_std_interp_3d, $bond_angle_std_interp_3d")
    println("bond_length_std_interp_3d, $bond_length_std_interp_3d")

    println("Interpolated bond_angle_std at T_max=$T_max_interp, β=$beta_interp, θ=$theta_interp: $bond_angle_std_interp_3d")
    println("Interpolated bond_length_std at T_max=$T_max_interp, β=$beta_interp, θ=$theta_interp: $bond_length_std_interp_3d")
end

#main1()
main2()





###




#=
import Polynomials



function main1()
    theta = [110, 130, 150, 170, 180]
    bond_angle_std = [0.1, 0.2, 0.15, 0.25, 0.3]
    bond_length_std = [0.05, 0.1, 0.08, 0.12, 0.15]

    # Fit quadratic polynomials to the data
    poly_angle = Polynomials.fit(theta, bond_angle_std, 2)
    poly_length = Polynomials.fit(theta, bond_length_std, 2)

    # Interpolated values (example)
    theta_interp = 140
    bond_angle_std_interp = poly_angle(theta_interp)
    bond_length_std_interp = poly_length(theta_interp)

    println("Interpolated bond_angle_std at θ=$theta_interp: $bond_angle_std_interp")
    println("Interpolated bond_length_std at θ=$theta_interp: $bond_length_std_interp")
end


function main2()
    # For different T_max and beta, you would need to fit polynomials for each combination
# Example data (replace with your actual data)
T_max_values = [0.1, 0.2]
beta_values = [0.285, 0.3]
theta_values = [110, 130, 150, 170, 180]
bond_angle_std_values = [
    [0.1, 0.2, 0.15, 0.25, 0.3],  # T_max=0.1, beta=0.285
    [0.12, 0.22, 0.18, 0.28, 0.32],  # T_max=0.1, beta=0.3
    [0.11, 0.21, 0.16, 0.26, 0.31],  # T_max=0.2, beta=0.285
    [0.13, 0.23, 0.19, 0.29, 0.33]   # T_max=0.2, beta=0.3
]
bond_length_std_values = [
    [0.05, 0.1, 0.08, 0.12, 0.15],  # T_max=0.1, beta=0.285
    [0.06, 0.11, 0.09, 0.13, 0.16],  # T_max=0.1, beta=0.3
    [0.055, 0.105, 0.085, 0.125, 0.155],  # T_max=0.2, beta=0.285
    [0.065, 0.115, 0.095, 0.135, 0.165]   # T_max=0.2, beta=0.3
]

# Fit polynomials for each combination
poly_angle_3d = [Polynomials.fit(theta_values, bond_angle_std_values[i], 2) for i in 1:length(T_max_values)*length(beta_values)]
poly_length_3d = [Polynomials.fit(theta_values, bond_length_std_values[i], 2) for i in 1:length(T_max_values)*length(beta_values)]

# Interpolated values for different T_max and beta
T_max_interp = 0.15
beta_interp = 0.29
theta_interp = 140

# Find the closest indices for T_max and beta
T_max_idx = findall(x -> x >= T_max_interp, T_max_values)[1]
beta_idx = findall(x -> x >= beta_interp, beta_values)[1]
index = (T_max_idx - 1) * length(beta_values) + beta_idx

bond_angle_std_interp_3d = poly_angle_3d[index](theta_interp)
bond_length_std_interp_3d = poly_length_3d[index](theta_interp)

println("Interpolated bond_angle_std at T_max=$T_max_interp, β=$beta_interp, θ=$theta_interp: $bond_angle_std_interp_3d")
println("Interpolated bond_length_std at T_max=$T_max_interp, β=$beta_interp, θ=$theta_interp: $bond_length_std_interp_3d")
end

#main1()
main2()
=#


#=
function main1()
    theta = [110, 130, 150, 170, 180]
    bond_angle_std = [0.1, 0.2, 0.15, 0.25, 0.3]
    bond_length_std = [0.05, 0.1, 0.08, 0.12, 0.15]

    # Combine bond_angle_std and bond_length_std into a single array
    combined_data = hcat(bond_angle_std, bond_length_std)

    # Fit quadratic polynomials to the combined data
    poly_combined = [Polynomials.fit(theta, combined_data[:, i], 2) for i in 1:size(combined_data, 2)]

    # Interpolated values (example)
    theta_interp = 140
    bond_angle_std_interp = poly_combined[1](theta_interp)
    bond_length_std_interp = poly_combined[2](theta_interp)

    println("Interpolated bond_angle_std at θ=$theta_interp: $bond_angle_std_interp")
    println("Interpolated bond_length_std at θ=$theta_interp: $bond_length_std_interp")
end
=#





###

#=
import Polynomials

# Define the number of values for T_max, beta, and theta
T_max_values = Float64[1.0, 2.0, 3.0, 4.0, 5.0]  # Example T_max values
beta_values = Float64[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]  # Example beta values
theta_values = Float64[10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]  # Example theta values

# Sample angle and length values corresponding to T_max, beta, and theta
# This should contain the outputs as per the combinations of input values
angle_values = rand(length(T_max_values) * length(beta_values) * length(theta_values)) * 360  # Random angles
length_values = rand(length(T_max_values) * length(beta_values) * length(theta_values)) * 10  # Random lengths

# For actual data, replace the above angle_values, and length_values with your real output values.

function interpolate(T_max_interpolation, beta_interpolation, theta_interpolation)
    # Create a 3D grid for the T_max, beta, and theta combinations
    results = zeros(length(T_max_interpolation), length(beta_interpolation), length(theta_interpolation), 2)
    
    for i in 1:length(T_max_values)
        for j in 1:length(beta_values)
            for k in 1:length(theta_values)
                idx = (i - 1) * length(beta_values) * length(theta_values) + (j - 1) * length(theta_values) + k
                results[i, j, k, 1] = angle_values[idx]  # First output is angle
                results[i, j, k, 2] = length_values[idx] # Second output is length
            end
        end
    end

    # Create a Polynomial interpolation
    poly_angle = Polynomials.Polynomial(
        [results[i, j, k, 1] for i in 1:length(T_max_values), j in 1:length(beta_values), k in 1:length(theta_values)]
    )
    
    poly_length = Polynomials.Polynomial(
        [results[i, j, k, 2] for i in 1:length(T_max_values), j in 1:length(beta_values), k in 1:length(theta_values)]
    )

    return poly_angle, poly_length
end

# Execute interpolation
poly_angle, poly_length = interpolate(T_max_values, beta_values, theta_values)

# Test the interpolation with example inputs
T_max_input = 3.0
beta_input = 0.4
theta_input = 30.0

interpolated_angle = poly_angle(T_max_input, beta_input, theta_input)
interpolated_length = poly_length(T_max_input, beta_input, theta_input)

println("Interpolated Angle: ", interpolated_angle)
println("Interpolated Length: ", interpolated_length)
=#