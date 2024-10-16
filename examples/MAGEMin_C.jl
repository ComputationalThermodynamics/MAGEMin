# set of examples to run MAGEMin_C

# metapelite database
using MAGEMin_C

data        =   Initialize_MAGEMin("mp", verbose=true, solver=0);
test        =   1
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   4.0
T           =   800.0
out         =   point_wise_minimization(P,T, data);


println(out.PP_vec[1].Comp_apfu);
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
println( out.rho - (out.rho_M*out.frac_M_wt + out.rho_S*out.frac_S_wt))