using MAGEMin_C
data    = Initialize_MAGEMin("mp", verbose=-1);
P,T     = 6.0, 930.0
Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
sys_in  = "wt"
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
Finalize_MAGEMin(data)


Zr_liq_sat = MAGEMin_C.zirconium_saturation(out; model = "CB")