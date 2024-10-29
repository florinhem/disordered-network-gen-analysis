include("structure_analysis_modules.jl")
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU

#path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\data_old\test\216_vertices_T_0.1_heat_cool_0.1_per_mc_quenched_evolution.h5"
path=raw"C:\Users\GlauserV\OneDrive - Université de Fribourg\Anlagen\AMI\Projekt\GitFlorin\code_photonic_structures\simulations\multiple_parameters\\"

filename=raw"multiple_p_quench_false_theta_array_N=216_T=0.1_Trial=1_Beta=0.285_Theta_GS=180.0_GradT=0.1_StepsPerT=0.01"

#spatial_network = NG.load_spatial_network_from_gml(path*filename)

evolution_dict = GU.load_h5_dict(path*filename*"_evolution.h5")

println(sum(evolution_dict["move_accepted_vec"]))