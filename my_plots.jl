
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import LaTeXStrings as Latex


path = raw"C:\Users\HemmannF\switchdrive\presentations\material\\"

fontsize=18

Plots.gr()
Plots.default(grid=false, 
legend = true, 
dpi=250,
xtickfontsize=fontsize,
ytickfontsize=fontsize,
xguidefontsize=fontsize,
yguidefontsize=fontsize,
legendfontsize=fontsize,
bottom_margin = 3Plots.mm,
linewidth=3, 
thickness_scaling = 1)

# get structure factor dicts

dict_path_1 = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\random_networks\1000_vertices_T_1_quenched_structure_factor_isotrope.h5"
dict_path_4 = raw"C:\Users\HemmannF\switchdrive\structure_analysis\analysis_data\random_networks\1000_vertices_T_4_quenched_structure_factor_isotrope.h5"

structure_factor_dict_1 = GU.load_h5_dict(dict_path_1)
structure_factor_dict_4 = GU.load_h5_dict(dict_path_4)

# get effective hyperuniformity parameter and fit parameters for T=1 and T=4
hyperuniformity_parameter_1, polynomial_fit_1 = NA.get_hyperuniformity_metric(structure_factor_dict_1)
hyperuniformity_parameter_4, polynomial_fit_4 = NA.get_hyperuniformity_metric(structure_factor_dict_4)

# plot structure factor
x_vec = collect(0:10/200:10)
fit_1_vec = polynomial_fit_1.(x_vec)
fit_4_vec = polynomial_fit_4.(x_vec)

Plots.plot(structure_factor_dict_4["wavenumber_vec"], structure_factor_dict_4["structure_factor_vec"], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=4" )
Plots.plot!(structure_factor_dict_1["wavenumber_vec"], structure_factor_dict_1["structure_factor_vec"], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=1" )
Plots.plot!(x_vec, fit_4_vec, linecolor="sienna", ls=:dash, label = "fit "*Latex.L"kT=4")
Plots.plot!(x_vec, fit_1_vec, linecolor="cyan2", ls=:dash, label = "fit "*Latex.L"kT=1")
Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,22.5), ylims=(0,20), size = (500, 600), bottom_margin = 0Plots.mm)

Plots.savefig(path*"structure_factor_T_1_4_fits.png")


Plots.plot(structure_factor_dict_4["wavenumber_vec"], structure_factor_dict_4["structure_factor_vec"], linecolor=Plots.palette(:tab10)[2], label = Latex.L"kT=4" )
Plots.plot!(structure_factor_dict_1["wavenumber_vec"], structure_factor_dict_1["structure_factor_vec"], linecolor=Plots.palette(:tab10)[1], label = Latex.L"kT=1" )
Plots.plot!(xlabel="wavenumber", ylabel = "structure factor", xlims=(0,22.5), ylims=(0,20), size = (500, 600), bottom_margin = 0Plots.mm)

Plots.savefig(path*"structure_factor_T_1_4_stretched.png")