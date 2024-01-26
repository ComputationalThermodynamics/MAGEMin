# 15/01/24
# simple example of how to compute a liquid line of descent using the igneous database using MAGEMin_C

using MAGEMin_C

# Starting composition [mol fraction], here we used an hydrous basalt; composition taken from Blatter et al., 2013 (01SB-872, Table 1), with added O and water saturated
oxides  = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"]
bulk    = [38.448328757254195, 7.718376151972274, 8.254653357127351, 9.95911842561036, 5.97899305676308, 0.24079752710315697, 2.2556006776515964, 0.7244006013202644, 0.7233140004182841, 0.0, 25.696417444779453];


# Define bulk-rock composition unit
sys_in  = "mol"

# Choose database
data    = Initialize_MAGEMin("ig", verbose=false);

## begin loop for LLD
    # Define PT
    P       = 10.0
    T       = 1200.0

    # compute point
    out     = single_point_minimization(P, T, data, X=bulk, Xoxides=oxides, sys_in=sys_in) 

    # assign new bulk
    bulk           .= out.bulk_M 
## end loop for LLD

