# Julia script to transform the STIX11.jld2 file into the format of the SB_ref_database.h file
# This reads the data exported as 

using JLD2
using DataFrames

# @load "STIX11.jld2" data

@load "STIX11.jld2" data2

data = data2
sb_ver = "sb11"
elems  = ["Si", "Ca", "Al", "Fe", "Mg", "Na"]

tab = "    "
out = ""
out *= "/* SiO2:0 CaO:1 Al2O3:2 FeO:3 MgO:4 Na2O:5 */\n"
out *= "EM_db_sb arr_em_db_sb_$(sb_ver)[$(size(data)[1])] = {\n"

# loop through all phases
for i=1:size(data)[1]
    out *= tab*"{\n"
    out *= tab*tab*"\"$(data[i,:abbrev])\", \"$(data[i,:id])\", \"$(data[i,:fml])\",\n"

    # retrieve the composition
    composition = join(collect(values(data[i, :oxides])), ", ")
    out *= tab*tab*"{$composition},\n"

    # retrieve site occupancies
    fml = data[i, :fml]
    n_sites = count(c -> c == '[', fml)

    # retrieve multiplicities and site compositions
    matches     = eachmatch(r"\[([^\[\]]+)\]", data[i, :fml])
    contents    = [m.captures[1] for m in matches]
    mul         = []
    site_cmp    = []
    for i in contents
        matches = eachmatch(r"_(\d+)", i)
        n_atoms = [parse(Int, m.captures[1]) for m in matches]

        matches = eachmatch(r"[A-Za-z]+", i)
        elements = [m.match for m in matches]

        if length(n_atoms) == 1
            push!(site_cmp, length(n_atoms))
            push!(site_cmp, findfirst(elems .== elements[1]), Float64(n_atoms[1]))
        else
            push!(site_cmp, length(n_atoms))
            for j=1:length(n_atoms)
                push!(site_cmp,  findfirst(elems .== elements[j]), Float64(n_atoms[j]))
            end
        end

        sum_atoms = sum(n_atoms)
        push!(mul, Float64(sum_atoms))
    end

    # line3   = zeros(Float64,16)
    # l_tmp   = cat(n_sites, mul, site_cmp, dims=1)

    # line3[1:length(l_tmp)]  .= l_tmp
    # line3   = join(collect(line3), ", ")

    # out    *= tab*tab*"{$line3},\n"

    # retrieve standard state properties
    line4   = join(collect(Vector(data[i, 5:14])), ", ")
    out *= tab*tab*"{$line4},\n"

    out *= tab*tab*"{$(data[i,:aSm]),$(data[i,:pd]),$(data[i,:td])},\n"

    out *= tab*"},\n"
end

out *= "};\n"

print(out)
