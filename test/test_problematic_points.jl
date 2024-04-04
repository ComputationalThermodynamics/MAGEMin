# This part provides a list of points that have proven to be problematic in the past...
using MAGEMin_C

# This gives problems on apple M2 Max (and macOS in general):
data    = Initialize_MAGEMin("ig", verbose=-1, solver=2);
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO";  "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"];

print("\toxides: $(Xoxides)\n")

P,T     = 8.42, 705.0
X       = [31.8;  8.28 ; 2.1;  0.1;  13.9;  0.0 ; 2.2;  4.2;  0.0;  0.0; 37.4];
sys_in  = "mol"  
print("\tbulk: $(X)\n")
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
print("\tsuccess!\n")

P,T     = 8.24, 810.0
X       = [51.2;  8.6;  0.0;  0.0;  0.1;  3.57;  4.9;  0.2;  0.0;  0.0;  31.5];
sys_in  = "mol"  
print("\tbulk: $(X)\n")  
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
print("\tsuccess!\n")

P,T     = 8.29, 635.0
X       = [ 51.4,  6.1,  0.4,  0.0 , 0.0 , 2.83,  3.0,  0.4 , 0.0,  0.0 , 35.9];
sys_in  = "mol"  
print("\tbulk: $(X)\n")  
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
print("\tsuccess!\n\n")

Finalize_MAGEMin(data)

