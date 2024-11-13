# NR 12-09-24 Script to test nullspace minimization with generic solution models 
using JSON3
using DataFrames
using JLD2
using LinearAlgebra
using Symbolics
using NLopt
using BenchmarkTools

include("test_min/function_list.jl")

n_em, N, f_config, f_grad_config, f_mu_Gex, f_grad = get_functions(6);

active = ones(Float64, n_em)
gb     = [73.1715026079, 6.9138347588, 69.1393955643, 30.4219113823].*1e3
# gb    = zeros(n_em)
# gb = [-3295812.69,-2585019.80,-3193633.69,-3310918.59]
# x   = rand(n_em) .+ 0.1; x ./= sum(x);
x = [ 0.0001000000, 0.6998000000, 0.0001000000 ,0.3000000000]
# x = [0.00010, 0.00010, 0.99970, 0.0001]
# x = [0.00010, 0.99970, 0.00010, 0.0001]
eps = 1e-8
opt = NLopt.Opt(:LD_SLSQP, n_em)
NLopt.lower_bounds!(opt, zeros(n_em) .+ eps)
NLopt.upper_bounds!(opt, ones(n_em)  .- eps)
NLopt.ftol_rel!(opt, 1e-5)
NLopt.min_objective!(opt, (x, g) -> objective(x, g, N, f_config, f_grad_config, f_mu_Gex, f_grad, gb, active))
# Define the equality constraint function
function equality_constraint(x, g)
    if g !== nothing
        g .= 1.0
    end
    return sum(x) - 1.0
end

# Add the equality constraint to the optimizer
equality_constraint!(opt, equality_constraint, 1e-8)

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

# objective value       : -16879.898968165835
# solution              : [0.04730673326327671, 0.49721125109705694, 0.06862191381345994, 0.3868601018262064]
# solution status       : FTOL_REACHED
# function evaluation : 51

# objective value       : -16879.89896836978
# solution              : [0.047307589116837634, 0.4972114871702633, 0.06862332110342888, 0.3868576026094702]
# solution status       : FTOL_REACHED
# function evaluation : 31

# act                 = findall(x .!= eps)

# if length(act) < n_em
#     inact        = findall(x .== eps)
#     x[act]      .= x[act] ./ sum(x[act])
#     x[inact]    .= 1e-8
#     active[inact] .= 0.0
#     N = nullspace(active')
# end

# NLopt.min_objective!(opt, (x, g) -> objective(x, g, N, f_config, f_grad_config, f_mu_Gex,  G0, active))

# min_f, x, ret   = NLopt.optimize(opt, x)
# num_evals           = NLopt.numevals(opt)

# println(
#     """
#     objective value       : $min_f
#     solution              : $x
#     solution status       : $ret
#     # function evaluation : $num_evals
#     """
# )

# a = [1,1,0,1,1,1]
# b = [1,1,1,1,1]
# Na = nullspace(a')
# Nb = nullspace(b')

# a = [0,1,1,1]
# b = [1,1,1]
# Na = nullspace(a')
# Nb = nullspace(b')
