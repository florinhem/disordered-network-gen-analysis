
# include file where structure analysis modules are stored
include("structure_analysis_modules.jl")

# import my module that contains all functions for the generation and analysis of networks
import .NetworkGeneration as NG
import .NetworkAnalysis as NA
import .GeneralUtilities as GU


# possible choices of nr_vertices for diamond: 64, 216, 512, 216, that is (2*n)^3 with natural nr natural

# the supercell edge lengths are 
# 216 vertices: supercell_edge_length = 11.547005383792516
# 512 vertices: supercell_edge_length = 9.237604307034013
# 216 vertices: supercell_edge_length = 6.9282032302755105
# 64 vertices: supercell_edge_length = 4.619802153517007
# which is the cube root of the number of vertices times 2/sqrt(3)

# julia --threads 23

path = raw"C:\Users\HemmannF\OneDrive - Université de Fribourg\structure_analysis\structures\random_networks\\"

function myfunc()
    count = 0
    filepath_list = []
    # Iterate through all directories and subdirectories
    for (root, dirs, files) in walkdir(path)
        for file in files
            if endswith(file, ".gml")

                joined_path = joinpath(root, file)
                spatial_network = NG.load_spatial_network_from_gml(joined_path)
                NG.save_spatial_network_to_gml(spatial_network,
                file[1:end-4];
                    save_path=root*"\\")

                #println(joinpath(root, file))  # Print the full path to the .gml file
                push!(filepath_list, joinpath(root, file))

                #print every 50th file
                count += 1
                if count % 50 == 0
                    println(count, " ", joined_path)
                end
            end
        end
    end

    return
end

myfunc()