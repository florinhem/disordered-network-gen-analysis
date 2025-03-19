
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import LaTeXStrings as Latex

fontsize=16

Plots.gr()
#Plots.plotlyjs()
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
thickness_scaling = 1,
framestyle = :box)

# functions to have pi ticks
function pitick(start, stop, denom; mode=:text)
    a = Int(cld(start, 2*π/denom))
    b = Int(fld(stop, 2*π/denom))
    tick = range(a*2*π/denom, b*2*π/denom; step=2*π/denom)
    ticklabel = piticklabel.( 2 .* (a:b) .// denom, Val(mode))
    tick, ticklabel
end

function piticklabel(x::Rational, ::Val{:text})
    iszero(x) && return "0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return S * N * "π"
    S * N * "π\\" * repr(d)
end

function piticklabel(x::Rational, ::Val{:latex})
    iszero(x) && return Latex.L"0"
    S = x < 0 ? "-" : ""
    n, d = abs(numerator(x)), denominator(x)
    N = n == 1 ? "" : repr(n)
    d == 1 && return Latex.L"%$S%$N\pi"
    Latex.L"%$S\frac{%$N\pi}{%$d}"
end


"""
    get_tickslogscale(lims; skiplog=false)
Return a tuple (ticks, ticklabels) for the axis limit `lims`
where multiples of 10 are major ticks with label and minor ticks have no label
skiplog argument should be set to true if `lims` is already in log scale.
"""
function get_tickslogscale(lims::Tuple{T, T}; skiplog::Bool=false) where {T<:AbstractFloat}
    mags = if skiplog
        # if the limits are already in log scale
        floor.(lims)
    else
        floor.(log10.(lims))
    end
    rlims = if skiplog; 10 .^(lims) else lims end

    total_tickvalues = []
    total_ticknames = []

    rgs = range(mags..., step=1)
    for (i, m) in enumerate(rgs)
        if m >= 0
            tickvalues = range(Int(10^m), Int(10^(m+1)); step=Int(10^m))
            ticknames  = vcat([string(round(Int, 10^(m)))],
                              ["" for i in 2:9],
                              [string(round(Int, 10^(m+1)))])
        else
            tickvalues = range(10^m, 10^(m+1); step=10^m)
            ticknames  = vcat([string(10^(m))], ["" for i in 2:9], [string(10^(m+1))])
        end

        if i==1
            # lower bound
            indexlb = findlast(x->x<rlims[1], tickvalues)
            if isnothing(indexlb); indexlb=1 end
        else
            indexlb = 1
        end
        if i==length(rgs)
            # higher bound
            indexhb = findfirst(x->x>rlims[2], tickvalues)
            if isnothing(indexhb); indexhb=10 end
        else
            # do not take the last index if not the last magnitude
            indexhb = 9
        end

        total_tickvalues = vcat(total_tickvalues, tickvalues[indexlb:indexhb])
        total_ticknames = vcat(total_ticknames, ticknames[indexlb:indexhb])
    end
    return (total_tickvalues, total_ticknames)
end

"""
    fancylogscale!(p; forcex=false, forcey=false)
Transform the ticks to log scale for the axis with scale=:log10.
forcex and forcey can be set to true to force the transformation
if the variable is already expressed in log10 units.
"""
function fancylogscale!(p::Plots.Subplot; forcex::Bool=false, forcey::Bool=false)
    kwargs = Dict()
    for (ax, force, lims) in zip((:x, :y), (forcex, forcey), (Plots.xlims, Plots.ylims))
        axis = Symbol("$(ax)axis")
        ticks = Symbol("$(ax)ticks")

        if force || p.attr[axis][:scale] == :log10
            # Get limits of the plot and convert to Float
            ls = float.(lims(p))
            ts = if force
                (vals, labs) = get_tickslogscale(ls; skiplog=true)
                (log10.(vals), labs)
            else
                get_tickslogscale(ls)
            end
            kwargs[ticks] = ts
        end
    end

    if length(kwargs) > 0
        Plots.plot!(p; kwargs...)
    end
    p
end
fancylogscale!(p::Plots.Plot; kwargs...) = (fancylogscale!(p.subplots[1]; kwargs...); return p)
fancylogscale!(; kwargs...) = fancylogscale!(Plots.plot!(); kwargs...)


plots_save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\random_networks\216_vertices_heat_cool_bond_bending_comparison\\"

bond_bending_vec = [0.06, 0.135, 0.21, 0.285, 0.36, 0.435, 0.51]

nr_vertices_vec = [216]

markershapes = [:utriangle, :pentagon, :diamond,  :circle, :star5, :rect, :dtriangle]

markersizes = [4]

color_palette = Plots.palette(:tab10)

order_metrics_names = ["bond_length_std_vec", "bond_angle_std_vec", "dihedral_angle_std_vec", "cluster_metric_vec", 
    "pore_size_distribution_second_moment_vec",
    "anisotropy_metric_from_structure_factor_vec", "anisotropy_metric_from_spectral_density_vec"]


order_metrics_labels = ["Bond length std", "Bond angle std", "Dihedral angle std", "Cluster metric",
    "Pore size dist. 2nd m.",
    "Anisotropy from s. f.", "Anisotropy from s. d."]


order_metrics_dict_arr = Array{Dict{Any, Any}, 2}(undef, length(bond_bending_vec), length(nr_vertices_vec))

for j in eachindex(bond_bending_vec) 

    # loop through folders and append all order metrics to the order_metrics_dict
    for k in eachindex(nr_vertices_vec)

        analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\random_networks\\"*string(nr_vertices_vec[k])*"_vertices_bond_bending_"*string(bond_bending_vec[j])*"\\"

        order_metrics_dict = Dict()

        for i in 1:5

            current_analysis_data_path = analysis_data_path*"run_"*string(i)*"\\"

            current_order_metrics_dict = GU.load_h5_dict(current_analysis_data_path*"all_order_metrics.h5")

            for (key, value) in current_order_metrics_dict
                if haskey(order_metrics_dict, key)
                    order_metrics_dict[key] = vcat(order_metrics_dict[key], value)
                else
                    order_metrics_dict[key] = value
                end
            end

        end

        # sort all vectors in order of the total keating energy
        for order_metric_name in order_metrics_names
            order_metrics_dict[order_metric_name] = order_metrics_dict[order_metric_name][sortperm(order_metrics_dict["total_keating_energy_vec"])]
        end
        order_metrics_dict["filenames_vec"] = order_metrics_dict["filenames_vec"][sortperm(order_metrics_dict["total_keating_energy_vec"])]
        sort!(order_metrics_dict["total_keating_energy_vec"])

        order_metrics_dict_arr[j, k] = order_metrics_dict

    end

end


Plots.scatter()
#Plots.plot()

for j in reverse(eachindex(bond_bending_vec) )
    for k in eachindex(nr_vertices_vec)
        # get the temperature from the filtered filenames
        pattern = r"T_([0-9\.]+)"
        extracted_numbers = [match(pattern, s).captures[1] for s in order_metrics_dict_arr[j, k]["filenames_vec"]]
        temperatures = parse.(Float64, extracted_numbers)

        min_temp = 0.08
        max_temp = 0.26
        normalized_temperatures = (temperatures .- min_temp) ./ (max_temp - min_temp)
        colormap = Plots.cgrad(:roma, rev = true, scale = :exp)
        mapped_colors = [colormap[normalized_temperature] for normalized_temperature in  normalized_temperatures]

        Plots.scatter!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, order_metrics_dict_arr[j, k]["bond_length_std_vec"], markershape = markershapes[j], color = mapped_colors, markersize=markersizes[k], alpha = 0.5, label = false)
        #Plots.plot!(order_metrics_dict_arr[j, k]["bond_angle_std_vec"] *180/pi, order_metrics_dict_arr[j, k]["bond_length_std_vec"], seriestype = :scatter, markershape = markershapes[j],  color=mapped_colors, label = Latex.L"\beta = "*string(bond_bending_vec[j])*", "*string(nr_vertices_vec[k])*" vertices", legend=:topleft, markersize=markersizes[k], alpha = 0.5)
        #Plots.plot!(x, y, zcolor = z, seriestype = :scatter, markersize=5, label = "points")zcolor = temperatures,

        
        Plots.scatter!([], [], color = "black", markershape = markershapes[j], label = Latex.L"\beta = "*Latex.latexstring(bond_bending_vec[j]), markersize=markersizes[k], alpha = 0.5, legend = true)

    end

end

Plots.scatter!([15.3, 13.6, 13.5], [0.045, 0.047, 0.066], color = ["#7BA8D6","#4549A9","#8C66EC"], markershape = [:diamond, :circle, :square], label = false, markersize=10, alpha = 1)

min_temp = 0.08
max_temp = 0.26
Plots.scatter!(size=(470,400), legend = true, xlabel = Latex.L" \sigma_\mathrm{angle} / °", ylabel = Latex.L" \sigma_\mathrm{length} / d", rightmargin=5Plots.mm, colorbar=true, clim=(min_temp, max_temp), color=Plots.cgrad(:roma, rev = true, scale = :exp), xlims=(-4, 23), ylims=(0.01, 0.097))
#Plots.plot!(size=(550,400), legend = false, xlabel = "Bond angle st. d. / °", ylabel = "Bond length st. d. / "*Latex.L"d", rightmargin=5Plots.mm, clim=(min_temp, max_temp), color_palette=Plots.cgrad(:roma, rev = true, scale = :exp))

Plots.savefig(plots_save_path*"bond_bending_stretching_highlighted.png")