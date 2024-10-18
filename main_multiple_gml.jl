
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import MetaGraphsNext
import Graphs
import Plots

function save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array,
    maximal_temperature_array,
    nr_trials_per_temperature,
    bond_bending_const)
    
    
    for nr_vertices in nr_vertices_array

        for maximal_temperature in maximal_temperature_array

            for i in 1:nr_trials_per_temperature
                
                println("$nr_vertices"*", "*"$maximal_temperature"*", "*"$i")

                evolution_dict = NA.get_evolution_dict(;nr_vertices = nr_vertices, network_type="diamond", bond_bending_const=bond_bending_const, min_ring_size=3)
                spatial_network = NG.get_periodic_network(evolution_dict)
            
                temperature_vec, nr_monte_carlo_steps_per_temperature_vec = 
                    NA.get_temperature_sequence_heating_cooling_gradient(
                        maximal_temperature;
                        temperature_gradient = 10.0, 
                        nr_monte_carlo_steps_per_temperature = 0.01,
                        quench = false)

                evolution_dict["temperature_vec"] = temperature_vec
                evolution_dict["nr_monte_carlo_steps_per_temperature_vec"] = nr_monte_carlo_steps_per_temperature_vec

                spatial_network, total_energy_vec, move_accepted_vec = NG.evolve_network_temperature_sequence!(
                    spatial_network,
                    evolution_dict;
                    print_progress = true,
                    print_every_nr_attempted_bond_switches = 10)

                save_path = raw".\my_networks\multiple_gml\\"
                filename = ("multiple_gml_6"
                    *"_N="*"$nr_vertices"
                    *"_T="*"$maximal_temperature"
                    *"_trial="*"$i"
                    *"_beta="*"$bond_bending_const")

            

                NG.save_spatial_network_to_gml(
                    spatial_network,
                    filename;
                    evolution_dict = evolution_dict,
                    save_path = save_path)
                

            end
        end
    end
end


save_multiple_N_T_trials_beta_gml(;
    nr_vertices_array=[216,512],
    maximal_temperature_array=[0.0,0.5,1.0,1.5],
    nr_trials_per_temperature=1,
    bond_bending_const=0.285
)