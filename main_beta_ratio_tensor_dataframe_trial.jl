# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots
import Colors
import Glob
import DataFrames
import LaTeXStrings
import IterTools
import Combinatorics
import Statistics
import Polynomials
import LinearAlgebra
using StatsPlots

# Function to fit the data and predict bond_angle_std and bond_length_std
function fit_and_predict_1(df::DataFrames.DataFrame, T_max_interp, beta_interp, theta_interp)
    # Separate input and output columns
    input_data = df[:, [:MaxT, :Beta, :Theta]]
    output_data = df[:, [:BondAngleStd, :BondLengthStd]]

    # Fit quadratic polynomials for bond_angle_std and bond_length_std
    poly_bond_angle = Polynomials.fit(input_data[:, 3], output_data[:, 1], 2)
    poly_bond_length = Polynomials.fit(input_data[:, 3], output_data[:, 2], 2)

    println("poly_bond_angle, $poly_bond_angle")
    println("poly_bond_length, $poly_bond_length")

    # Predict the bond_angle_std and bond_length_std for given T_max, beta, and theta
    bond_angle_std_interp = poly_bond_angle(theta_interp)
    bond_length_std_interp = poly_bond_length(theta_interp)

    return bond_angle_std_interp, bond_length_std_interp
end

function fit_and_predict(df::DataFrames.DataFrame, T_max_interp, beta_interp, theta_interp)
    # Separate input and output columns
    input_data = df[:, [:MaxT, :Beta, :Theta]]
    output_data = df[:, [:BondAngleStd, :BondLengthStd]]

    # Create the design matrix for a quadratic fit
    X = hcat(ones(size(input_data, 1)), input_data[:, 1], input_data[:, 2], input_data[:, 3], input_data[:, 1].^2, input_data[:, 2].^2, input_data[:, 3].^2, input_data[:, 1].*input_data[:, 2], input_data[:, 1].*input_data[:, 3], input_data[:, 2].*input_data[:, 3])
    #println("X:")
    #display(X)
    # Fit the model for bond_angle_std
    coeffs_bond_angle = X \ output_data[:, 1]
    #println("coeffs_bond_angle, $coeffs_bond_angle")

    # Fit the model for bond_length_std
    coeffs_bond_length = X \ output_data[:, 2]

    # Create the design matrix for the interpolation point
    X_interp = [1, T_max_interp, beta_interp, theta_interp, T_max_interp^2, beta_interp^2, theta_interp^2, T_max_interp*beta_interp, T_max_interp*theta_interp, beta_interp*theta_interp]
    #println("X_interp, $X_interp")
    # Predict the bond_angle_std and bond_length_std for given T_max, beta, and theta
    bond_angle_std_interp = LinearAlgebra.dot(X_interp, coeffs_bond_angle)
    bond_length_std_interp = LinearAlgebra.dot(X_interp, coeffs_bond_length)

    return bond_angle_std_interp, bond_length_std_interp
end


function get_data_and_predict(;
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature_array,
    T_max_interp,
    beta_interp,
    theta_interp,
    save_path,
    filename_start,
    plot_save_path,
    plot_filename_start)

    # test before we begin
    @assert length(nr_vertices_array)>=1
    @assert length(maximal_temperature_array)>=1
    @assert length(bond_bending_const_array)>=1
    @assert length(temperature_gradient_array)>=1
    @assert length(nr_monte_carlo_steps_per_temperature_array)>=1
    @assert length(theta_ground_state_array)>=1
    @assert length(nr_trials_per_temperature_array)>=1

    # store array with all the paths in the directory to check if we have
    # this .h5 and .gml to be able to plot it.
    path_array=Glob.glob(filename_start*"*",save_path)

    # for data storage 
    data::Vector{Vector{Float64}}=[]

    for k in eachindex(nr_vertices_array)

        nr_vertices=nr_vertices_array[k]

        for j in eachindex(maximal_temperature_array)

            maximal_temperature=maximal_temperature_array[j]

            for m in eachindex(bond_bending_const_array)

                bond_bending_const=bond_bending_const_array[m]

                for n in eachindex(temperature_gradient_array)

                    temperature_gradient=temperature_gradient_array[n]
                    
                    for o in eachindex(nr_monte_carlo_steps_per_temperature_array)

                        nr_monte_carlo_steps_per_temperature=nr_monte_carlo_steps_per_temperature_array[o]
                    
                        for p in eachindex(theta_ground_state_array)
                            
                            theta_ground_state=theta_ground_state_array[p]

                            for i in eachindex(nr_trials_per_temperature_array)

                                trial=nr_trials_per_temperature_array[i]
                                
                                filename = (filename_start
                                    *"_N="*"$nr_vertices"
                                    *"_T="*"$maximal_temperature"
                                    *"_Beta="*"$bond_bending_const"
                                    *"_GradT="*"$temperature_gradient"
                                    *"_StepsPerT="*"$nr_monte_carlo_steps_per_temperature"
                                    *"_Theta_GS="*"$theta_ground_state"
                                    *"_Trial="*"$trial"
                                    )

                                total_path=save_path*filename

                                if(total_path*".gml" in path_array)
                                    println(filename)

                                    spatial_network=NG.load_spatial_network_from_gml(total_path*".gml")
                                    #Energies
                                    stretching_energy=NG.get_stretching_energy_keating(spatial_network)
                                    bending_energy=NG.get_bending_energy_keating(spatial_network)
                                    total_energy=NG.get_total_energy_keating(spatial_network)
                                    @assert stretching_energy+bending_energy===total_energy
                                    ratio_energy=stretching_energy/bending_energy
                                    #Std of bond
                                    bond_length_std, bond_length_vec = NA.get_bond_length_std(spatial_network)
                                    bond_angle_std, bond_angle_vec = NA.get_bond_angle_std(spatial_network)
                                    #Accepted moves
                                    evolution_dict = GU.load_h5_dict(total_path*"_evolution.h5")
                                    accepted_moves=sum(evolution_dict["move_accepted_vec"])

                                    #=display(bending_energy)
                                    display(stretching_energy)
                                    display(total_energy)
                                    display(ratio_energy)=#

                                    data=push!(data,
                                        [
                                        nr_vertices,
                                        maximal_temperature,
                                        bond_bending_const,
                                        temperature_gradient,
                                        nr_monte_carlo_steps_per_temperature,
                                        theta_ground_state,
                                        trial,
                                        bond_length_std,
                                        bond_angle_std,
                                        accepted_moves,
                                        stretching_energy,
                                        bending_energy,
                                        total_energy,
                                        ratio_energy
                                        ]
                                    )
                        
                                else
                                    println("file not in directory")
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    
    #println(data)
    #println(typeof(data))

    # data cleaning

    matrix=mapreduce(permutedims, vcat, data)

    df=DataFrames.DataFrame(matrix,:auto)

    DataFrames.rename!(df, 
        [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13, :x14] 
        .=>  [:NrVertices, :MaxT, :Beta, :GradT, :MCsteps, :Theta, :Trial, :BondLengthStd, 
            :BondAngleStd, :AccMoves, :StretchEnergy, :BendEnergy ,
            :TotalEnergy, :RatioEnergy])

    println("df:")
    display(df)

    df_filtered=filter(row -> row.AccMoves > 250, df)

    println("df_filtered:")
    display(df_filtered)
    
    bond_angle_std_interp, bond_length_std_interp = fit_and_predict(df_filtered, T_max_interp, beta_interp, theta_interp)
    println("Interpolated bond_angle_std: $bond_angle_std_interp")
    println("Interpolated bond_length_std: $bond_length_std_interp")

    

end

get_data_and_predict(
    nr_vertices_array=[216],
    maximal_temperature_array=[0.1,0.125,0.15,0.175,0.2], 
    bond_bending_const_array=[0.2,0.25,0.3,0.35,0.4],
    temperature_gradient_array=[0.1],
    nr_monte_carlo_steps_per_temperature_array=[0.01],
    theta_ground_state_array=[100.0,120.0,160.0,180.0],                 
    nr_trials_per_temperature_array=[1,2,3,4,5],
    T_max_interp = 0.15,
    beta_interp = 0.3,
    theta_interp = 140,
    save_path = raw".\simulations\multiple_trials\\", 
    filename_start="m_t_",
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_t_1_"
)
