# NR 12-09-24 Script to test nullspace minimization with generic solution models 
using JSON3
using DataFrames
using JLD2
using LinearAlgebra
using Symbolics
using NLopt
using BenchmarkTools

include("test_min/function_list.jl")

n_em, N, f_config, f_grad_config, f_mu_Gex = get_functions(10);

active = ones(Float64, n_em)

x   = rand(n_em) .+ 0.1; x ./= sum(x);
eps = 1e-8
opt = NLopt.Opt(:LD_LBFGS, n_em)
NLopt.lower_bounds!(opt, zeros(n_em) .+ eps)
NLopt.upper_bounds!(opt, ones(n_em)  .- eps)
NLopt.ftol_rel!(opt, 1e-8)
NLopt.min_objective!(opt, (x, g) -> objective(x, g, N, f_config, f_grad_config, f_mu_Gex, active))

min_f, x, ret   = NLopt.optimize(opt, x)
num_evals           = NLopt.numevals(opt)

act                 = findall(x .!= eps)

if length(act) < n_em
    inact        = findall(x .== eps)
    x[act]      .= x[act] ./ sum(x[act])
    x[inact]    .= 1e-8
    active[inact] .= 0.0
    N = nullspace(active')
end

NLopt.min_objective!(opt, (x, g) -> objective(x, g, N, f_config, f_grad_config, f_mu_Gex, active))

min_f, x, ret   = NLopt.optimize(opt, x)
num_evals           = NLopt.numevals(opt)

println(
    """
    objective value       : $min_f
    solution              : $x
    solution status       : $ret
    # function evaluation : $num_evals
    """
)

a = [1,1,0,1,1,1]
b = [1,1,1,1,1]
Na = nullspace(a')
Nb = nullspace(b')

a = [0,1,1,1]
b = [1,1,1]
Na = nullspace(a')
Nb = nullspace(b')
