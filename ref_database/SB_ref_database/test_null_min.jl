# NR 12-09-24 Script to test nullspace minimization with generic solution models 
using JSON3
using DataFrames
using JLD2


include("test_min/function_list.jl")
include("test_min/function_stixrude.jl")


# load stixrude informations
@load "STIX11.jld2" data2
ss      = JSON3.read("stx11_solution.json", Vector{ModelJSON}) 
elems   = ["Si", "Ca", "Al", "Fe", "Mg", "Na"]
n_el    = length(elems)


n_ss    = length(ss)

i       = 7

n_ox    = length(data2[1,:oxides])
n_em    = length(ss[i].endmembers)
em_list = String.(keys(ss[i].endmembers))

v       = []
if ~isempty(ss[i].van_laar)
    v   = [ ss[i].van_laar[j] for j in  keys(ss[i].van_laar)]
end

W       = []
if ~isempty(ss[i].margules)
    W   = [ ss[i].margules[j] for j in  keys(ss[i].margules)]
end

em_comp = zeros(Float64, n_em, n_ox)
for j=1:n_em
    em          = ss[i].endmembers[em_list[j]][1]
    id          = findfirst(data2[!,:abbrev] .== em)
    em_comp[j,:]= [ (data2[id,:oxides][k]) for k in keys(data2[id,:oxides]) ]
end

println("v: $v")
println("W: $W")
println("em_comp: $em_comp")
