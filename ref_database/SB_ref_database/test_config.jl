using LinearAlgebra

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
T = 1000.



# SiO2:0 CaO:1 Al2O3:2 FeO:3 MgO:4 Na2O:5
# Si Ca Al Fe Mg Na
C       = [     3/4 0;    #Mg s1
                0 0;      #Fe s1
                0 3/4;    #Fe s1
                1/4 1/4;  #Al s1
                1/8 0;    #Mg s2
                0 1/8;    #Fe s2
                7/8 7/8]  #Al s2

M       = [4.0,4.0,4.0,4.0,8.0,8.0,8.0]

X       = [0.5,0.5]

Xo      = C*X
config  = R*T* (M'*Diagonal(Xo)*logish.(Xo))


tmp     = M.*Diagonal(Xo)*logish.(Xo)
cfg1    = R*T* (sum(tmp[[1,3,4,6]]) - tmp[3]/2  - tmp[6]/2)
cfg2    = R*T* (sum(tmp[[2,3,5,6]]) - tmp[3]/2  - tmp[6]/2)

cfg1+cfg2


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
    config  = R*T* (M'*Diagonal(Xo)*logish.(Xo))

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