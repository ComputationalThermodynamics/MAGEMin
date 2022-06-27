using MAGEMin_C

# Initialize database 
gv, DB = init_MAGEMin();

# Bulk rock composition for a particular test 
test      = 0;
bulk_rock = get_bulk_rock(gv, test)

# Call optimization routine for given P & T
P       = 8.0;
T       = 800.0;
gv.verbose = -1;
out     = point_wise_minimization(P,T, bulk_rock, gv, DB);


#sgleP =0
#rank=0;
#gv.verbose=0
#LibMAGEMin.PrintOutput(gv, rank, sgleP, DB, out.time_ms/1e3, z_b);	



