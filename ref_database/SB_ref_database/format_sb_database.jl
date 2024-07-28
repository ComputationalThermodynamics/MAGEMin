# Julia script to transform the STIX11.jld2 file into the format of the SB_ref_database.h file

using JLD2
using DataFrames

@load "STIX11.jld2" data


sb_ver = "sb11"

tab = "    "
out = ""
out *= "/* SiO2:0 CaO:1 Al2O3:2 FeO:3 MgO:4 Na2O:5 */\n"
out *= "EM_db_sb arr_em_db_sb_$(sb_ver)[$(size(data)[1])] = {\n"

# loop through all phases
for i=1:size(data)[1]
out *= tab*"{\n"
out *= tab*tab*"\"$(data[i,:id])\", \n"

# retrieve the composition
composition = join(collect(values(data[i, :oxides])), " ")
out *= tab*tab*"{$composition},\n"

# retrieve site occupancies
fml = data[i, :fml]
n_sites = count(c -> c == '[', fml)

# retrieve multiplicities
matches = eachmatch(r"\[([^\[\]]+)\]", data[i, :fml])
contents = [m.captures[1] for m in matches]
mul = []
for i in contents
    matches = eachmatch(r"_(\d+)", i)
    numbers = [parse(Int, m.captures[1]) for m in matches]
    sum_numbers = sum(numbers)
    push!(mul, sum_numbers)
end
mul = join(collect(mul), " ")


out *= tab*tab*"{$n_sites $mul},\n"


out *= tab*"},\n"
end

out *= "}\n"
print(out)






# EM_db_sb arr_em_db_sb_sb11[47] = {
#     {
#         "fo",
#         {1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.0},
#         {-2172.59, 0.0951, 4.366},
#         {0.2333, 1.494e-06, -603.8, -1.8697},
#         {2.85e-05, 1285.0, 3.84, -0.003, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
#         {810.0,0.00182,-140.0}
#     },
# }