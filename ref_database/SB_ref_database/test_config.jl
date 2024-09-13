using LinearAlgebra
using Symbolics

function logish(x, meps = 1.0e-7)
    """
    2nd order series expansion of log(x) about eps:
    log(eps) - sum_k=1^infty (f_eps)^k / k
    Prevents infinities at x=0
    """
    f_eps = 1.0 - x / meps
    if x < meps
        ln = log(meps) - f_eps - f_eps * f_eps / 2.0
    else
        ln = log(x)
    end

    return ln
end


R = 8.31446261815324
T = 1000.0


# Si Ca Al Fe Mg Na

C       = [     3/4 0;    #Mg s1
                0 3/4;    #Fe s1
                1/4 1/4;  #Al s1
                1/8 0;    #Mg s2
                0 1/8;    #Fe s2
                7/8 7/8]  #Al s2

M       = [4.0,4.0,4.0,8.0,8.0,8.0]
n_em    = size(C,2)

@variables xv[1:n_em]
X               = Symbolics.scalarize(xv)
 
Xo              = C*X
config          = R * T * (M' * Diagonal(Xo) * log.(Xo))
grad_config     = Symbolics.gradient(config, X)

f_config        = build_function(config,        X, expression = Val{false})
f_grad_config,  = build_function(grad_config,   X, expression = Val{false});


x = [0.1, 0.9]
f_config(x)
f_grad_config(x)

 

# grad_config = Symbolics.gradient(config, [x1, x2])
# f_grad_config, grad_config! = build_function(grad_config, X, expression = Val{false});
# f_grad_config([0.1, 0.9])
# # Print the symbolic expression for the gradient
# println(grad_config)

# x1_val = 0.1
# x2_val = 0.9
# config_evaluated = Symbolics.substitute(config, Dict(x1 => x1_val, x2 => x2_val))
# grad_config_evaluated = Symbolics.substitute(grad_config, Dict(x1 => x1_val, x2 => x2_val))

# # Print the evaluated gradient
# println("Evaluated gradient at x1=$x1_val, x2=$x2_val: ", grad_config_evaluated)


# config_evaluated = Symbolics.substitute(config, Dict(x1 => x1_val, x2 => x2_val))





n_ss = length(ss)
for i = 1:n_ss

    mul, site_cmp = retrieve_site_cmp(ss, i)

    n_sf = size(site_cmp)[1]
    n_ox = size(site_cmp)[2]
    n_em = size(site_cmp)[3]

    C = zeros(n_sf*n_ox, n_em)
    M = zeros(n_sf*n_ox)

    for k=1:length(mul)
        M[(k-1)*n_ox+1:k*n_ox] .= mul[k]
    end

    for k=1:n_sf
        C[(k-1)*n_ox+1:k*n_ox,:] .= site_cmp[k,:,:] ./mul[k]
    end

    X = zeros(n_em) .+ 1.0/n_em


    Xo      = C*X
    config  = R*T* (M'*Diagonal(Xo)*log.(Xo))


    # config_em = R*T .* X .* log.(X)


    println("$i $(ss[i].abbrev): $config")
    
end









#=
def G_config_default(self):
    if self.ns==0: return 0
    C = sym.Matrix(self.C)
    X_o = C*self.X
    M = sym.Matrix(self.m)
    G = self.nT*self.Rsym*self.T*M.T*sym.diag(*X_o)*X_o.applyfunc(sym.log)
    G_config = G[0].simplify()
    return G_config
=#