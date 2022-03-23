using MAGEMin_C


# Initialize database 
gv, DB = init_MAGEMin();

# Bulk rock composition for a particular test 
test      = 0;
bulk_rock = get_bulk_rock(gv, test)

# Call optimization routine for given P & T
P       = 8.
T       = 800.
gv.verbose = -1
gv, z_b, time = point_wise_minimization(P,T, bulk_rock, gv, DB);


sgleP =0
rank=0;
LibMAGEMin.PrintOutput(gv, rank, sgleP, DB, time, z_b);	

