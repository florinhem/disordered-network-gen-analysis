
"""
This module contains all functions with general utitilities like
IO functions
"""
module GeneralUtilities

# import the necessary packages
import Measurements    # for handling data with uncertainty and error propagation
import FileIO   # for file loading and saving
import LinearAlgebra    # to perform linear algebra calculations like the dot product

# set path to files
if Sys.iswindows()
    load_path = raw"..\code_photonic_structures\general_utilities\\"
else
    load_path = raw"../code_photonic_structures/general_utilities/"
end

# include IO utilities
include(load_path*"utilities_io.jl")

# include calculation utilities
include(load_path*"utilities_calculations.jl")

end


"""
This module combines all functions to statistically analyze binary structure data
"""
module BinaryDataAnalysis

# import the necessary packages
import ..GeneralUtilities as GU  # for general utitilities like IO functions
import FileIO   # for file loading and saving
import Images # for image conversion, in our case from color to grayscale
import Statistics   # for statistical operations like mean()
import Measurements    # for handling data with uncertainty and error propagation
import FFTW     # to calculate the Fast Fourier Transform
import Plots    # for plotting
import GLMakie  # Makie backend for plotting
import LaTeXStrings as Latex # to display latex symbols in plot labels
import Format   # for python-like string formatting 
import MetaGraphsNext   # to convert spatial network to binary data

# set path to files
if Sys.iswindows()
    load_path = raw"..\code_photonic_structures\binary_structures\\"
else
    load_path = raw"../code_photonic_structures/binary_structures/"
end

# include functions to load binary 3d data, convert and correct it and to extract
# some basic measures
include(load_path*"binary_structure_io_basics.jl")

# include functions to generate binary 3d data from nodal surfaces
include(load_path*"binary_structure_generation.jl")

# include functions to calculate local volume
# fraction variance by means of different measuring windows 
include(load_path*"binary_structure_volume_fraction_variance.jl")

# include functions to calculate the autocovariance function and the spectral
# density for 3d data assuming statistical isotropy and only working with 
# absolute values of distance and momentum
include(load_path*"binary_structure_spectral_density_isotrope.jl")

# include functions to calculate the autocovariance function and the spectral
# density for 3d data as a function of distance vector and momentum vector
include(load_path*"binary_structure_spectral_density_direction.jl")

# include functions that are used to plot different statistical measures of 
# 3d binary structures
include(load_path*"binary_structure_plotting.jl")

end


"""
This module combines all functions to generate network structures
"""
module NetworkGeneration

# import the necessary packages
import ..GeneralUtilities as GU  # for general utitilities like IO functions
import Graphs   # for simple graphs
import MetaGraphsNext   # to deal with graphs with labelled vertices and edges
import Combinatorics    # mainly used to get all possible combinations of bonds
import LinearAlgebra    # to perform linear algebra calculations like the dot product
import Optim    # for optimization such as relaxation of individual vertices
import GLMakie  # Makie backend for plotting
import GraphMakie   # additions to GLMakie for graph plotting
import FileIO   # for file loading and saving
import Random   # to access more features for random number generation like seeds
import Statistics   # for statistical operations like mean()
import GeometryBasics   # to create meshes out of networks
import Plots    # for plotting
import LaTeXStrings as Latex # to display latex symbols in plot labels
import Format   # for python-like string formatting 

# set path to files
if Sys.iswindows()
    load_path = raw"..\code_photonic_structures\networks\\"
else
    load_path = raw"../code_photonic_structures/networks/"
end

# these functions are utilities for network generation and modification 
include(load_path*"network_utilities.jl")

# these functions generate networks
include(load_path*"network_generation.jl")

# these functions calculate the energy of networks
include(load_path*"network_energy_calculation.jl")

# these functions modify networks
include(load_path*"network_modification.jl")

# functions to plot networks
include(load_path*"network_plotting.jl")

# functions for loading and saving graph data
include(load_path*"network_io.jl")

end



"""
This module combines all functions to analyze network structures
"""
module NetworkAnalysis

# import the necessary packages
import ..GeneralUtilities as GU  # for general utitilities like IO functions
import ..NetworkGeneration as NG    # to generate network structures
import MetaGraphsNext   # to deal with graphs with labelled vertices and edges
import Optim    # for optimization such as relaxation of individual vertices
import LinearAlgebra    # to perform linear algebra calculations like the dot product
import Statistics   # for statistical operations like mean()
import Measurements    # for handling data with uncertainty and error propagation
import Polynomials  # to curve fit polynomials
import Peaks    # to locate peaks of a function
import Combinatorics    # mainly used to get all possible combinations of bonds
import SphericalHarmonics   # to calculate spherical harmonics
import Plots    # for plotting
import LaTeXStrings as Latex # to display latex symbols in plot labels
import Format   # for python-like string formatting 

# set path to files
if Sys.iswindows()
    load_path = raw"..\code_photonic_structures\networks\\"
else
    load_path = raw"../code_photonic_structures/networks/"
end

# these functions can be used to characterize networks
# by means of local order parameters
include(load_path*"network_analysis_local.jl")

# these functions can be used to characterize networks
# by means of order parameters measuring correlations
include(load_path*"network_analysis_correlations.jl")

# these functions provide utilities for analyzing networks
include(load_path*"network_analysis_utilities.jl")

# these functions are used to evolve and eventually analyze networks
include(load_path*"network_analysis_and_evolution.jl")

# functions to plot network analysis results
include(load_path*"network_analysis_plotting.jl")

end


