
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU



# possible choices of nr_vertices for diamond: 64, 216, 512, 1000, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 1000 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

crystals = ["ctn", "dia", "lcs", "srs"]

for crystal in crystals
    spatial_networks_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\neural_network_networks\\"*crystal*"\\run_2\\"
    save_path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\neural_network_networks\for_simulation\\"*crystal*"\\run_2\\"
    # get all filenames in the spatial network path
    filenames = readdir(spatial_networks_path)
    # filter the filenames to only contain the ones that end with ".gml"
    filenames = filter(filename -> endswith(filename, ".gml"), filenames)

    # filter out those files that have already been processed
    filenames = filter(filename -> !isfile(joinpath(save_path, replace(filename, ".gml" => "_for_sim.gml"))), filenames)

    # cut away the ".gml" extension from the filenames
    filenames = map(filename -> replace(filename, r"\.gml$" => ""), filenames)
    for filename in filenames
        spatial_network = NG.load_spatial_network_from_gml(joinpath(spatial_networks_path, filename)*".gml")
        spatial_network = NG.get_spatial_network_for_simulation!(
            spatial_network;
            vector_out_of_supercell_length = 1,
            duplicate_bonds_close_to_supercell_edge = true,
            save_result = true,
            filename = filename,
            save_path = save_path)
    end
end