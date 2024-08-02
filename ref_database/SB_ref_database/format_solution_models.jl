# Julia script to transform the STIX11.jld2 file into the format of the SB_ref_database.h file
# This reads the data exported as 


struct ModelJSON
    name            :: String
    abbrev          :: String
    endmembers      :: Dict{String, Vector{String}}
    margules        :: Dict{String, Float64}
    van_laar        :: Vector{Float64}
end

using JSON3

sb_ver  = "sb11"
elems   = ["Si", "Ca", "Al", "Fe", "Mg", "Na"]
n_el    = length(elems)

ss = JSON3.read("stx11_solution.json", Vector{ModelJSON}) 

n_ss = length(ss)
for i = 1:n_ss
    println(ss[i].endmembers)

    em          = String.(keys(ss[i].endmembers))
    n_em        = length(em)

    fml         = ss[i].endmembers[em[1]][2]

    matches     = eachmatch(r"\[([^\[\]]+)\]", fml)
    contents    = [m.captures[1] for m in matches]

    n_sf        = length(contents)

    mul         = zeros(Float64, n_sf)
    site_cmp    = zeros(Float64, n_sf,n_el, n_em)

    for j=1:n_em
        fml         = ss[i].endmembers[em[j]][2]
        matches     = eachmatch(r"\[([^\[\]]+)\]", fml)
        contents    = [m.captures[1] for m in matches]
        for k=1:length(contents)
            matches = eachmatch(r"_(\d+)", contents[k])
            n_atoms = [parse(Int, m.captures[1]) for m in matches]
            matches = eachmatch(r"[A-Za-z]+", contents[k])
            elements= [m.match for m in matches]
            mul[k]  = sum(n_atoms)

            if length(n_atoms) == 1
                id              = findfirst(elems .== elements[1])
                # site_cmp[id + (k-1)*n_el, j] = Float64(n_atoms[1])
                site_cmp[k,id, j] = Float64(n_atoms[1])
            else
                for l=1:length(n_atoms)
                    id              = findfirst(elems .== elements[l])
                    site_cmp[k,id, j] = Float64(n_atoms[1])
                end
            end
        end
    end

end