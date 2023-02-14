using MAGEMin_C

gv, z_b, DB, splx_data      = init_MAGEMin();

test        = 0;
sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
gv          = use_predefined_bulk_rock(gv, test)

# Call optimization routine for given P & T & bulk_rock
P           = 8.0
T           = 800.0
gv.verbose  = -1        # switch off any verbose
out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
