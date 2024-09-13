
using NLopt

function my_objective_fn(x::Vector, grad::Vector, pouet)

    if length(grad) > 0
        grad[1] = 2.0 * x[1]
        grad[2] = 2.0 * x[2]
    end
    return x[1]^2 + x[2]^2 + pouet
end


pouet = 5.0

opt = NLopt.Opt(:LD_LBFGS, 2)

NLopt.lower_bounds!(opt, [-10.0, -10.0])
NLopt.upper_bounds!(opt, [+10.0, +10.0])
NLopt.xtol_rel!(opt, 1e-4)
NLopt.min_objective!(opt, (x, g) -> my_objective_fn(x, g, 5.0))

min_f, min_x, ret = NLopt.optimize(opt, [1.234, 5.678])
num_evals = NLopt.numevals(opt)
println(
    """
    objective value       : $min_f
    solution              : $min_x
    solution status       : $ret
    # function evaluation : $num_evals
    """
)