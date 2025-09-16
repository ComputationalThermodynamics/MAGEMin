using MAGEMin_C

sys_in = "wt"
Xoxides = ["SiO2"; "CaO"; "Al2O3";  "FeO"; "MgO"; "Na2O"]
MAGEMin_db = Initialize_MAGEMin("sb21", solver=0, verbose=false)
p = [215.]
t = [754.]
bulk =  [0.49921349253700387, 0.0732104764346424, 0.0357867270241625, 0.11340340153185195, 0.27532164424930744, 0.003064258223031849]

out = multi_point_minimization(p, t, MAGEMin_db, X=bulk, Xoxides=Xoxides, sys_in=sys_in)

out[1].bulk_wt .- bulk


bulk_mol = wt2mol(bulk, Xoxides)
bulk_wt = mol2wt(bulk_mol, Xoxides)



MAGEMin_bulk = [49.921349253700384, 7.32104764346424, 3.57867270241625, 11.340340153185196, 27.532164424930745, 0.3064258223031849]
MAGEMin_ox   = ["SiO2", "CaO", "Al2O3", "FeO", "MgO", "Na2O"]