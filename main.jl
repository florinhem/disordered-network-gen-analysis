
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU
import Plots
import LaTeXStrings as Latex # to display latex symbols in plot labels
import NaNStatistics

# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105


dict_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\random_networks\without_ring_size_limitation\\"


start_temperature = 2
end_temperature = 1
temperature_decrease_per_monte_carlo_step = 0.5
nr_monte_carlo_steps_per_temperature = 0.01
quench = true 

nr_monte_carlo_steps_during_temperature_decrease = ((start_temperature - end_temperature)
    /temperature_decrease_per_monte_carlo_step)

temperature_vec = (start_temperature 
    .- temperature_decrease_per_monte_carlo_step 
        .* collect(0
            :nr_monte_carlo_steps_per_temperature
            :nr_monte_carlo_steps_during_temperature_decrease))

# create temperature sequence
nr_monte_carlo_steps_per_temperature_vec = vcat([2], ones(length(temperature_vec)-1) .* nr_monte_carlo_steps_per_temperature )

if quench
    if end_temperature == 0
        nr_monte_carlo_steps_per_temperature_vec[end] = 50
    else
        push!(temperature_vec, 0)
        push!(nr_monte_carlo_steps_per_temperature_vec, 50)
    end
end
