# NR 12-09-24 Script to test nullspace minimization with generic solution models 
using JSON3
using DataFrames
using JLD2
using LinearAlgebra
using Symbolics
using NLopt

include("test_min/function_list.jl")

n_em, N, f_config, f_grad_config, f_mu_Gex = get_functions(8);

x  = rand(n_em) .+ 0.1; x ./= sum(x);

opt = NLopt.Opt(:LD_LBFGS, n_em)
NLopt.lower_bounds!(opt, zeros(n_em) .+ 1e-8)
NLopt.upper_bounds!(opt, ones(n_em)  .- 1e-8)
NLopt.ftol_rel!(opt, 1e-4)
NLopt.min_objective!(opt, (x, g) -> objective(x, g, N, f_config, f_grad_config, f_mu_Gex))

min_f, min_x, ret   = NLopt.optimize(opt, x)
num_evals           = NLopt.numevals(opt)

println(
    """
    objective value       : $min_f
    solution              : $min_x
    solution status       : $ret
    # function evaluation : $num_evals
    """
)
