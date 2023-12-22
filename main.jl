
#include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

#import my module that contains all functions for the analysis of binary structure data
import .BinaryDataAnalysis as BDA
import Hankel

data_path_raw = raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\3Dvolumes\Blue_SI.tif"


structure_dict = BDA.get_structure_dict_from_colorscale(data_path_raw; 
    voxel_size=(10,12,10), 
    label = "P. c. mirabilis blue",
    save_result=true, 
    save_path=raw"C:\Users\HemmannF\switchdrive\structure_analysis\structures\pachy\pachy_blue")
