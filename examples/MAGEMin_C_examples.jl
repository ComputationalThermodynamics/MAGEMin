#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#

using MAGEMin_C
rm_list =   remove_phases(["fl","H2O"],"mpe")
data    = Initialize_MAGEMin("mpe", verbose=1, solver=0);
P,T     = 0.8, 603.0
Xoxides = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "MnO", "H2O", "CO2","S"];
X       = [0.5430992417692342, 0.09796194931925375, 0.005897681874450057, 0.030429284690742315, 0.04851244934858918, 0.021337981308402636, 0.011328750786070731, 0.005798239767617565, 0.005578931675087756, 0.0, 0.20, 0.03,0.005]
sys_in  = "mol"
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, rm_list=rm_list)
Finalize_MAGEMin(data)

using MAGEMin_C
data    = Initialize_MAGEMin("all", verbose=-1, solver=0);
P,T     = 0.8, 1203.0
Xoxides = ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "MnO", "H2O", "CO2","S"];
X       = [0.5430992417692342, 0.09796194931925375, 0.005897681874450057, 0.030429284690742315, 0.04851244934858918, 0.021337981308402636, 0.011328750786070731, 0.005798239767617565, 0.005578931675087756, 0.0, 0.20, 0.03,0.005]
sys_in  = "mol"
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
Finalize_MAGEMin(data)


using Printf
@printf("P = %.3f kbar, T = %.3f C\n", P, T)
for i=1:length(Xoxides)
    @printf(" %-12s %10.6f\n", Xoxides[i], X[i])
end
for i=1:length(out.PH_vec[:fl_DEW].emFrac)
    if out.PH_vec[:fl_DEW].emFrac[i] > 0.0
        @printf(" %-12s %12.10f\n", out.PH_vec[:fl_DEW].emNames[i], out.PH_vec[:fl_DEW].emFrac[i])
    end
end


using MAGEMin_C
data    = Initialize_MAGEMin("mp", verbose=-1);
P,T     = 6.0, 710.0
Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
sys_in  = "wt"
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, scp = 1, dT=1.5)
Finalize_MAGEMin(data)

using MAGEMin_C
data        =   Initialize_MAGEMin("ig", verbose=-1);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);
P           =   8.0
T           =   1800.0
out         =   point_wise_minimization(P,T, data);
Finalize_MAGEMin(data)

using MAGEMin_C
data        =   Initialize_MAGEMin("mb", verbose=-1);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);
P           =   4.0
T           =   600.0
out         =   point_wise_minimization(P,T, data);
Finalize_MAGEMin(data)

using MAGEMin_C
data        =   Initialize_MAGEMin("xMELTS", verbose=-1);
test        =   0 #rough basalt
data        =   use_predefined_bulk_rock(data, test);
P           =   8.0
T           =   800.0
out         =   point_wise_minimization(P,T, data);

using MAGEMin_C
data        =   Initialize_MAGEMin("pMELTS", verbose=-1);
test        =   0 #rough basalt
data        =   use_predefined_bulk_rock(data, test);
P           =   8.0
T           =   800.0
out         =   point_wise_minimization(P,T, data);

using MAGEMin_C
data        =   Initialize_MAGEMin("pMELTS", verbose=1);
test        =   0 #rough basalt
data        =   use_predefined_bulk_rock(data, test);
P           =   8.0
T           =   800.0
#out         =   point_wise_minimization(P,T, data);
rm_list =   remove_phases(["spn","rhm"],"pMELTS")
out     = single_point_minimization(P, T, data, rm_list=rm_list)

using MAGEMin_C
data        =   Initialize_MAGEMin("rMELTS", verbose=-1);
test        =   0 #rough basalt
data        =   use_predefined_bulk_rock(data, test);
P           =   1.0
T           =   500.0
out         =   point_wise_minimization(P,T, data)

rm_list =   remove_phases(["spn","rhm"],"pMELTS")
out     = single_point_minimization(P, T, data, rm_list=rm_list)
