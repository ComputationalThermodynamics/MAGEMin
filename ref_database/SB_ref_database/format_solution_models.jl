# Julia script to transform the STIX11.jld2 file into the format of the SB_ref_database.h file
# This reads the data exported as 


using JSON3, DataFrames, JLD2, Symbolics, Combinatorics
include("functions_ss.jl")

sb      = 11
elems   = sb==24 ? ["Si", "Ca", "Al", "Mg", "Na", "O", "Cr", "Fe"] : ["Si", "Ca", "Al", "Fe", "Mg", "Na"]
n_el    = length(elems)

if sb == 11

    sb_ver  = "sb11"

    data    = read_data("stx11_data.json")
    ss      = JSON3.read("stx11_solution.json", Vector{ModelJSON})
elseif sb == 21

    sb_ver  = "sb21"

    data    = read_data("stx21_data.json")
    ss      = JSON3.read("stx21_solution.json", Vector{ModelJSON})
elseif sb == 24

    sb_ver  = "sb24"

    data    = read_data("stx24_data.json")
    ss      = JSON3.read("stx24_solution.json", Vector{ModelJSON})
end

sb_gss_init_function, sb_gss_function, sb_objective_functions, sb_SS_xeos_PC, SB_NLopt_opt_functions = generate_C_files(sb_ver,ss,data);


# print(sb_gss_init_function)
# print(sb_gss_function)
# print(sb_objective_functions)
# print(sb_SS_xeos_PC)
# print(SB_NLopt_opt_functions)