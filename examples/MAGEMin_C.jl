# set of examples to run MAGEMin_C

# metapelite database
using MAGEMin_C

data        =   Initialize_MAGEMin("mp", verbose=true, solver=0);
test        =   0
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   4.0
T           =   500.0
out         =   point_wise_minimization(P,T, data);
# P           =   8.0
# T           =   1900.0
# out         =   point_wise_minimization(P,T, data);

if 1 == 0
    # println(out.PP_vec[1].Comp_apfu);
    println(out.SS_vec[1].Comp_apfu);

    println(out.bulk_S)
    println(out.bulk_S_wt)
    println(out.bulk_M)
    println(out.bulk_M_wt)
    println(out.bulk_F)
    println(out.bulk_F_wt)

    println(sum(out.bulk_S))
    println(sum(out.bulk_S_wt))
    println(sum(out.bulk_M))
    println(sum(out.bulk_M_wt))
    println(sum(out.bulk_F))
    println(sum(out.bulk_F_wt))

    println(out.frac_F + out.frac_M + out.frac_S)
    println(out.frac_F_wt + out.frac_M_wt + out.frac_S_wt)

    println(out.rho," ", out.rho_M," ", out.rho_S," ", out.rho_F)
    println( out.rho - (out.rho_M*out.frac_M_wt + out.rho_S*out.frac_S_wt + out.rho_F*out.frac_F_wt))
end

id          = 1
el          = out.elements
n_el        = length(el)
cmd2eval    = "Mg / (Mg + Fe)"

for i = 1:n_el
    if occursin(el[i], cmd2eval)
        cmd2eval = replace(cmd2eval, el[i] => "out.SS_vec[id].Comp_apfu[$i]")
    end
end


command  = Meta.parse(cmd2eval)
eval(command)