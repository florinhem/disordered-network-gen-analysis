

"""
This module combines all functions to statistically analyze binary structure data
"""
module BinaryDataAnalysis

#import the necessary packages
import FileIO   #for file loading and saving
import Images #for image conversion, in our case from color to grayscale
import Statistics   #for statistical operations like mean()
import Measurements    #for handling data with uncertainty and error propagationy
import FFTW     #to calculate the Fast Fourier Transform
import Plots    #for plotting
import LaTeXStrings as Latex #to display latex symbols in plot labels
import Formatting as Fmt    #for python-like string formatting

#set path to files
load_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\code_photonic_structures\binary_structes\\"

#include functions to load binary 3d data, convert and correct it and to extract
#some basic measures
include(load_path*"binary_structure_io_basics.jl")

#include functions to generate binary 3d data from nodal surfaces
include(load_path*"binary_structure_generation.jl")

#include functions to calculate local volume
#fraction variance by means of different measuring windows 
include(load_path*"binary_structure_volume_fraction_variance.jl")

#include functions to calculate the autocovariance function and the spectral
#density for 3d data assuming statistical isotropy and only working with 
#absolute values of distance and momentum
include(load_path*"binary_structure_spectral_density_isotrope.jl")

#include functions to calculate the autocovariance function and the spectral
#density for 3d data as a function of distance vector and momentum vector
include(load_path*"binary_structure_spectral_density_direction.jl")

#include functions that are used to plot different statistical measures of 
#3d binary structures
include(load_path*"binary_structure_plotting.jl")

end


"""
This module combines all functions to generate network structures
"""
module NetworkGeneration

#import the necessary packages
import Graphs   #for simple graphs
import MetaGraphsNext   #to deal with graphs with labelled vertices and edges
import Combinatorics    #mainly used to get all possible combinations of bonds
import LinearAlgebra    #to perform linear algebra calculations like the dot product
import Optim    #for optimization such as relaxation of individual vertices
import GLMakie  #Makie backend for plotting
import GraphMakie   #additions to GLMakie for graph plotting
import FileIO   #for file loading and saving
import Random   #to access more features for random number generation like seeds
import Statistics   #for statistical operations like mean()
import CoordinateTransformations    #to convert between coordinate systems
import SphericalHarmonics   #to calculate spherical harmonics
import GeometryBasics   #to create meshes out of networks
import Plots    #for plotting
import LaTeXStrings as Latex #to display latex symbols in plot labels
import Formatting as Fmt    #for python-like string formatting

#set path to files
load_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\code_photonic_structures\networks\\"

#these functions are utilities for network generation and modification 
include(load_path*"network_utilities.jl")

#these functions generate networks
include(load_path*"network_generation.jl")

#these functions calculate the energy of networks
include(load_path*"network_energy_calculation.jl")

#these functions modify networks
include(load_path*"network_modification.jl")

#these functions can be used to characterize networks
#by means of local order parameters
include(load_path*"network_characterization_local.jl")

#these functions can be used to characterize networks
#by means of order parameters measuring correlations
include(load_path*"network_characterization_correlations.jl")

#functions for graph plotting
include(load_path*"network_plotting.jl")

#functions for loading and saving graph data
include(load_path*"network_io.jl")

end



"""
This module combines all functions to analyze network structures
"""
module NetworkAnalysis

#import the necessary packages
import MetaGraphsNext   #to deal with graphs with labelled vertices and edges
import LinearAlgebra    #to perform linear algebra calculations like the dot product

#set path to files
load_path = raw"C:\Users\HemmannF\switchdrive\structure_analysis\code_photonic_structures\networks\\"

#these functions can be used to characterize networks
#by means of local order parameters
include(load_path*"network_characterization_local.jl")

#these functions can be used to characterize networks
#by means of order parameters measuring correlations
include(load_path*"network_characterization_correlations.jl")

end


