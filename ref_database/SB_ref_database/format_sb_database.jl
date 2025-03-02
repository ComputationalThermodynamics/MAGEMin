# Julia script to transform the STIX11.jld2 file into the format of the SB_ref_database.h file
# This reads the data exported as 

using JLD2
using DataFrames

include("functions_ss.jl")

sb_ver = "sb21"

if sb_ver == "sb11"
    data = read_data("stx11_data.json")
elseif sb_ver == "sb21"
    data = read_data("stx21_data.json")
end

out = format_em(data)
print(out)
