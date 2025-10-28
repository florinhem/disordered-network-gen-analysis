
# include file where structure analysis modules are stored
include("structure_analysis_modules_no_plotting.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

network_path = raw"..\structures\biological\networks\pachy\\"
analysis_data_path = raw"..\analysis_data\biological\networks\pachy\\"



filename = "pachy_red"
    
NA.get_order_metrics(filename,
    network_path,
    analysis_data_path;
    hyperuniformity_min_wavenumber_to_consider=pi/2,
    save_result = true)