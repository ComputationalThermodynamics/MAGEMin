# Julia script to transform the STIX11.jld2 file into the format of the SB_ref_database.h file
# This reads the data exported as 


using JSON3, DataFrames, JLD2, Symbolics, Combinatorics
include("functions_ss.jl")

sb_ver  = "sb11"
elems   = ["Si", "Ca", "Al", "Fe", "Mg", "Na"]
n_el    = length(elems)

ss      = JSON3.read("stx11_solution.json", Vector{ModelJSON}) 

sb_gss_init_function, sb_gss_function, sb_objective_functions, sb_SS_xeos_PC, SB_NLopt_opt_functions = generate_C_files(sb_ver,ss)


# print(sb_gss_init_function)

# print(sb_gss_function)

print(sb_objective_functions)

# print(sb_SS_xeos_PC)

# print(SB_NLopt_opt_functions)
