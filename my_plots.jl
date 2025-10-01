
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

import Plots
import LaTeXStrings as Latex
import Measurements
import GLMakie

import CSV
import DataFrames
import Statistics

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
bottom_margin = 1Plots.mm,
linewidth=3, 
thickness_scaling = 1,
framestyle = :box,
fontfamily="DejaVu Sans")

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

function load_table_as_dict(filename::String)
    df = CSV.File(filename; delim='\t') |> DataFrames.DataFrame
    return Dict(col => collect(df[!, col]) for col in names(df))
end

save_filename = "dia_nr_vertices_216_beta_0.2500_t_max_0.1787_t_gradient_0.0477"
network_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\local_relaxation\targeted\shell_nr_4\run_1\\"

analysis_data_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\structure_factor_comparison\\"
digital_sphere_mask_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\analysis_data\crystals\digital_sphere_masks\\"

plot_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\plots\local_relaxation\structure_factor_comparison\\"

maximal_wavevector_int = 5

spatial_network = NG.load_spatial_network_from_gml(
    network_path*save_filename*".gml")

structure_factor_dict = NA.get_structure_factor_by_wavevector_array(
    spatial_network;
    consider_bonds = false,
    maximal_wavevector_int = maximal_wavevector_int,
    periodic_boundary_conditions = true,
    wavevector_array_positive_z = NA.get_wavevector_array_positive_z(spatial_network; 
        maximal_wavevector_int=maximal_wavevector_int,
            periodic_boundary_conditions=true),
    save_result = false,
    save_path = analysis_data_path*save_filename,
    label = nothing,
    print_progress = true,
    thread_nr = 0,
    print_lock = Threads.ReentrantLock())

structure_factor_angle_averaged_dict = NA.get_structure_factor_angle_averaged(
        structure_factor_dict;
        consider_bonds = false,
        gaussian_filter = true,
        gaussian_filter_sigma_x = 2*pi/25, 
        gaussian_filter_filtered_data_x_step_length = 2*pi/25,
        save_result = false,
        save_path = analysis_data_path*save_filename,
        label = nothing)

# convert to a spatial_network_without periodic boundaries
spatial_network_no_pbc = NA.convert_periodic_to_non_periodic(spatial_network)

structure_factor_no_pbc_dict = NA.get_structure_factor_by_wavevector_array(
    spatial_network_no_pbc;
    consider_bonds = false,
    maximal_wavevector_int = maximal_wavevector_int,
    periodic_boundary_conditions = false,
    wavevector_array_positive_z = NA.get_wavevector_array_positive_z(spatial_network_no_pbc; 
        maximal_wavevector_int=maximal_wavevector_int,
            periodic_boundary_conditions=false),
    save_result = false,
    save_path = analysis_data_path*save_filename,
    label = nothing,
    print_progress = true,
    thread_nr = 0,
    print_lock = Threads.ReentrantLock())

structure_factor_no_pbc_angle_averaged_dict = NA.get_structure_factor_angle_averaged(
        structure_factor_no_pbc_dict;
        consider_bonds = false,
        gaussian_filter = true,
        gaussian_filter_sigma_x = 2*pi/25, 
        gaussian_filter_filtered_data_x_step_length = 2*pi/25,
        save_result = false,
        save_path = analysis_data_path*save_filename,
        label = nothing)

Plots.plot(structure_factor_angle_averaged_dict["wavenumber_vec"], Measurements.value.(structure_factor_angle_averaged_dict["structure_factor_vec"]), ribbon =  Measurements.uncertainty.(structure_factor_angle_averaged_dict["structure_factor_vec"]), label="with PBC")
Plots.plot!(structure_factor_no_pbc_angle_averaged_dict["wavenumber_vec"], Measurements.value.(structure_factor_no_pbc_angle_averaged_dict["structure_factor_vec"]), ribbon =  Measurements.uncertainty.(structure_factor_no_pbc_angle_averaged_dict["structure_factor_vec"]), label="without PBC")
Plots.xlabel!("Wavenumber / (1/d)")
Plots.ylabel!("Structure factor")
Plots.savefig(plot_path*save_filename*"_structure_factor_comparison.png")
