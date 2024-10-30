using HDF5

path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\data_old\test\216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_evolution.h5"

c = h5open(path, "r") do file
    read(file)
end

#%%

c = h5open(path, "r") do file
    names(c)
end

#%%

using HDF5
f = h5open(path)
lvl1keys = keys(f)
for key in lvl1keys
    println(key)
    println(f[key])
end

#data = [f[joinpath(key, "rest/of/path")] for key in lvl1keys]

import .GeneralUtilities as GU

dict=GU.load_h5_dict(path)

Keys=keys(dict)

for k in keys(dict)
    println(k)
end

#%%

for i in dict 
    println(i)
end

#%%

include("structure_analysis_modules.jl")
import .GeneralUtilities as GU

#path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\data_old\test\216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_evolution.h5"
path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\multiple_p_6_N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=110.0_GradT=0.1_StepsPerT=0.01_evolution.h5"

function pretty_print(d::Dict, pre=1)
    for (k,v) in d
        if k in ["temperature_vec","total_energy_vec","move_accepted_vec","nr_monte_carlo_steps_per_temperature_vec"]
            println("nothing")
        else
            if typeof(v) <: Dict
                s = "$(repr(k)) => "
                println(join(fill(" ", pre)) * s)
                pretty_print(v, pre+1+length(s))
            else
                println(join(fill(" ", pre)) * "$(repr(k)) => $(repr(v))")
            end
        end
    end
end


dict=GU.load_h5_dict(path)

pretty_print(dict)

dict["theta_ground_state"]