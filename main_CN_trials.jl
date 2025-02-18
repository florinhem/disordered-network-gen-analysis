
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
using StatsPlots

function scatter_plot_for_mulitple_gml(;
    nr_vertices_array,
    maximal_temperature_array,
    bond_bending_const_array,
    temperature_gradient_array,
    nr_monte_carlo_steps_per_temperature_array,
    theta_ground_state_array,
    nr_trials_per_temperature_array,
    theta_compare,                                   
    theta_compare2,
    ratio_small,
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

                            stretching_energy_array=[]
                            bending_energy_array=[]
                            total_energy_array=[]
                            ratio_energy_array=[]
                            bond_length_std_array=[]
                            bond_angle_std_array=[]
                            accepted_moves_array=[]

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

                                    append!(stretching_energy_array,stretching_energy)
                                    append!(bending_energy_array,bending_energy)
                                    append!(total_energy_array,total_energy)
                                    append!(ratio_energy_array,ratio_energy)
                                    append!(bond_length_std_array,bond_length_std)
                                    append!(bond_angle_std_array, bond_angle_std)
                                    append!(accepted_moves_array,accepted_moves)
                                    
                        
                                else
                                    println("file not in directory")
                                end
                            end

                            


                            if(stretching_energy_array!=[])
                                display("stretching_energy_array, $stretching_energy_array")
                                display("bending_energy_array, $bending_energy_array")
                                display("total_energy_array, $total_energy_array")
                                display("ratio_energy_array, $ratio_energy_array")
                                display("bond_length_std_array, $bond_length_std_array")
                                display("bond_angle_std_array, $bond_angle_std_array")
                                display("accepted_moves_array, $accepted_moves_array")

                                stretching_energy_mean=Statistics.mean(stretching_energy_array)
                                bending_energy_mean=Statistics.mean(bending_energy_array)
                                total_energy_mean=Statistics.mean(total_energy_array)
                                ratio_energy_mean=Statistics.mean(ratio_energy_array)
                                bond_length_std_mean=Statistics.mean(bond_length_std_array)
                                bond_angle_std_mean=Statistics.mean(bond_angle_std_array)
                                accepted_moves_mean=Statistics.mean(accepted_moves_array)

                                display("stretching_energy_mean, $stretching_energy_mean")

                                data=push!(data,
                                    [
                                        nr_vertices,
                                        maximal_temperature,
                                        bond_bending_const,
                                        temperature_gradient,
                                        nr_monte_carlo_steps_per_temperature,
                                        theta_ground_state,
                                        bond_length_std_mean,
                                        bond_angle_std_mean,
                                        accepted_moves_mean,
                                        stretching_energy_mean,
                                        bending_energy_mean,
                                        total_energy_mean,
                                        ratio_energy_mean
                                    ]
                                )
                            end

                        end
                    end
                end
            end
        end
    end

    
    println(data)
    #println(typeof(data))

    # data cleaning

    matrix=mapreduce(permutedims, vcat, data)

    df=DataFrames.DataFrame(matrix,:auto)

    DataFrames.rename!(df, 
        [:x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13] 
        .=>  [:NrVertices, :MaxT, :Beta, :GradT, :MCsteps, :Theta, :BondLenghtStd, 
            :BondAngleStd, :AccMoves, :StretchEnergy, :BendEnergy ,
            :TotalEnergy, :RatioEnergy])

    df.AccMovesLg10 = log10.(df.AccMoves .+ 1)
    #println(df)

    #Cleaning
    #df=filter(row -> row.AccMoves > 250, df)
    
    display(df)

    #-----------------------
    # Delta
    #-----------------------



    Load_M=zeros(size(df,1),size(df,1),11)
    for i in 1:size(df,1)
        for j in 1:size(df,1)
            Load_M[i,j,1]=df[i,2]
            Load_M[i,j,2]=df[j,2]
            Load_M[i,j,3]=df[i,3]
            Load_M[i,j,4]=df[j,3]
            Load_M[i,j,5]=df[i,6]
            Load_M[i,j,6]=df[j,6]
            Load_M[i,j,7]=df[i,10]
            Load_M[i,j,8]=df[j,10]
            Load_M[i,j,9]=df[i,11]/df[i,3]
            Load_M[i,j,10]=df[j,11]/df[j,3]

            deltaSigma=sqrt((1-df[i,7]/df[j,7])^2+(1-df[i,8]/df[j,8])^2)
            Load_M[i,j,11]=deltaSigma
            if  i===j || df[i,6]===df[j,6]
                Load_M[i,j,11]=Inf64
            end
            
        end
    end
    #display(Load_M)



    Round_M=zeros(size(df,1),size(df,1),11)
    for i in 1:size(Load_M,1)
        for j in 1:size(Load_M,2)
            for k in 1:size(Load_M,3)
                Round_M[i,j,k]=round(Load_M[i,j,k],digits=3)
            end
        end
    end
    #display(Round_M)

   


    d=Round_M[:,:,11] |> vec |> sort
    d=filter(!isinf,d)

    #display(d)




    number_of_smallest=ceil(Int,ratio_small*size(d,1))
    max=maximum(first(d, number_of_smallest), dims = 1)[1]
    #display("max: $max")




    #=Countsb1b2=zeros(size(b1,1),size(b2,1))
    PossibleCountsb1b2=zeros(size(b1,1),size(b2,1))
    Countst1t2=zeros(size(t1,1),size(t2,1))
    PossibleCountst1t2=zeros(size(t1,1),size(t2,1))=#
    #println(typeof(number_of_smallest))
    ListRatioBBEE=zeros(number_of_smallest)
    count=1
    for i in 1:size(Round_M,1)
        for j in 1:size(Round_M,2)
            T1=Round_M[i,j,1]
            T2=Round_M[i,j,2]
            B1=Round_M[i,j,3]
            B2=Round_M[i,j,4]
            Theta1=Round_M[i,j,5]
            Theta2=Round_M[i,j,6]
            
            StrechE1=Round_M[i,j,7]
            StrechE2=Round_M[i,j,8]
            BendE1=Round_M[i,j,9]
            BendE2=Round_M[i,j,10]
            D=Round_M[i,j,11]
            #display("$B1, $B2, $Theta1, $Theta2, $StrechE1, $StrechE2, $BendE1, $BendE2, $D ")

            if D<=max   #if D>100 && D<200 #Extreme case 

                if Theta1===theta_compare && Theta2===theta_compare2
                    println()
                    println("Compare 2 networks:")

                    println("T1: $T1, T2: $T2")
                    println("B1: $B1, B2: $B2")
                    println("Theta1: $Theta1, Theta2: $Theta2")
                    println("StrechE1: $StrechE1, StrechE2: $StrechE2")
                    println("BendE1: $BendE1, BendE2: $BendE2")
                    println("$(round(B2/B1,digits=4)) =?= $(round(BendE1/BendE2,digits=4))")
                    println("$(round(B2/B1*BendE2/BendE1,digits=4)) =?= 1")
                    ListRatioBBEE[count]=round(B2/B1*BendE2/BendE1,digits=4)
                    count+=1

                elseif Theta2===theta_compare && Theta1===theta_compare2
                    println()
                    println("Compare two networks:")

                    println("T2: $T2, T1: $T1")
                    println("B2: $B2, B1: $B1")
                    println("Theta2: $Theta2, Theta1: $Theta1")
                    println("StrechE2: $StrechE2, StrechE1: $StrechE1")
                    println("BendE2: $BendE2, BendE1: $BendE1")
                    println("$(round(B1/B2,digits=4)) =?= $(round(BendE2/BendE1,digits=4))")
                    println("$(round(B1/B2*BendE1/BendE2,digits=4)) =?= 1")
                    ListRatioBBEE[count]=round(B1/B2*BendE1/BendE2,digits=4)
                    count+=1
                else
                    println("Err")
                end
            end
           
            
        end
    end
    println("List:")
    println(ListRatioBBEE)
    println("Mean:")
    println(Statistics.mean(ListRatioBBEE))
    println("Std:")
    println(Statistics.std(ListRatioBBEE))
    #display(Countsb1b2)
    #display(PossibleCountsb1b2)
    #display(Countst1t2)
    #display(PossibleCountst1t2)


    #=

    b1=Round_M[:,:,1] |> unique |> sort
    b1=filter(!iszero,b1)
    b2=Round_M[:,:,2] |> unique |> sort
    b2=filter(!iszero,b2)
    t1=Round_M[:,:,3] |> unique |> sort
    t1=filter(!iszero,t1)
    t2=Round_M[:,:,4] |> unique |> sort
    t2=filter(!iszero,t2)
    d=Round_M[:,:,5] |> vec |> sort
    d=filter(!isinf,d)
    theta1=Round_M[:,:,6] |> unique |> sort
    theta1=filter(!iszero,theta1)
    theta2=Round_M[:,:,7] |> unique |> sort
    theta2=filter(!iszero,theta2)

    #display(b1)
    #display(b2)
    #display(t1)
    #display(t2)
    #display(d)
    #display(theta1)
    #display(theta2)

    
    
    number_of_smallest=ceil(Int,ratio_small*size(d,1))
    max=maximum(first(d, number_of_smallest), dims = 1)[1]
    #display("max: $max")
    


    Countsb1b2=zeros(size(b1,1),size(b2,1))
    PossibleCountsb1b2=zeros(size(b1,1),size(b2,1))
    Countst1t2=zeros(size(t1,1),size(t2,1))
    PossibleCountst1t2=zeros(size(t1,1),size(t2,1))
    for i in 1:size(Delta2MatrixRound,1)
        for j in 1:size(Delta2MatrixRound,2)
            B1=Delta2MatrixRound[i,j,1]
            B2=Delta2MatrixRound[i,j,2]
            T1=Delta2MatrixRound[i,j,3]
            T2=Delta2MatrixRound[i,j,4]
            D=Delta2MatrixRound[i,j,5]
            Theta1=Delta2MatrixRound[i,j,6]
            Theta2=Delta2MatrixRound[i,j,7]
            #display("$B1, $B2, $T1, $T2, $D")
            #display(B)
            #display(D)
           
            for k in 1:size(b1,1)
                for m in 1:size(b2,1)
                    if B1===b1[k] && B2===b2[m] && Theta1===theta_compare2 && Theta2===theta_compare
                        if D<=max
                            #display("$B1, $B2")
                            Countsb1b2[k,m]+=1
                        end
                        if D!=Inf64
                            #display("$T,$B,$D")
                            PossibleCountsb1b2[k,m]+=1
                        end
                    end
                end
            end

            for k in 1:size(t1,1)
                for m in 1:size(t2,1)
                    if T1===t1[k] && T2===t2[m] && Theta1===theta_compare2 && Theta2===theta_compare && B2-B1===0.1
                        if D<=max
                            display("$B1, $B2")
                            Countst1t2[k,m]+=1
                        end
                        if D!=Inf64
                            #display("$T,$B,$D")
                            PossibleCountst1t2[k,m]+=1
                        end
                    end
                end
            end
        end
    end
    #display(Countsb1b2)
    #display(PossibleCountsb1b2)
    #display(Countst1t2)
    #display(PossibleCountst1t2)

    NormalizedCountsb1b2=zeros(size(b1,1),size(b2,1))
    for i in 1:size(b1,1)
        for j in 1:size(b2,1)
            NormalizedCountsb1b2[i,j]=Countsb1b2[i,j]/PossibleCountsb1b2[i,j]
        end
    end

    #display(NormalizedCountsb1b2)

    NormalizedCountst1t2=zeros(size(t1,1),size(t2,1))
    for i in 1:size(t1,1)
        for j in 1:size(t2,1)
            NormalizedCountst1t2[i,j]=Countst1t2[i,j]/PossibleCountst1t2[i,j]
        end
    end

    #display(NormalizedCountst1t2)


    Delta2Sort=Delta2MatrixRound[Delta2MatrixRound[:,:,5].<=max,:]
    
    =#
end


scatter_plot_for_mulitple_gml(
    nr_vertices_array=[64],             #[216],
    maximal_temperature_array=[0.1,0.125,0.15,0.175,0.2],
    bond_bending_const_array=[0.15,0.2,0.25,0.3,0.35],
    temperature_gradient_array=[0.1],     
    nr_monte_carlo_steps_per_temperature_array=[0.01],    
    theta_ground_state_array=[100.0,110.0,],
    nr_trials_per_temperature=5,
    theta_compare=110.0,                                   
    theta_compare2=180.0,
    ratio_small=0.025,
    save_path = raw".\simulations\multiple_parameters\\",
    filename_start = "m_a3_CN",    
    plot_save_path = raw".\simulations\analysis_plot\\",
    plot_filename_start = "m_r_a1_"
)
