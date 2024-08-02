struct ModelJSON
    name            :: String
    abbrev          :: String
    endmembers      :: Dict{String, Vector{Vector{Float64}}}
    margules        :: Dict{String, Float64}
    van_laar        :: Vector{Float64}
end

using JSON3

ss = JSON3.read("stx11_solution.json", Vector{ModelJSON}) 

n_ss = length(ss)
for i = 1:n_ss
    println(ss[i].name)
end