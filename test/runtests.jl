# this tests the julia interface to MAGEMin
using Test

using MAGEMin_C

# Initialize database 
gv, DB = init_MAGEMin();

# Bulk rock composition for a particular test 
test      = 0;
bulk_rock = get_bulk_rock(gv, test)

# Call optimization routine for given P & T
P           = 8.
T           = 800.
gv.verbose  = 0
gv, z_b, time = point_wise_minimization(P,T, bulk_rock, gv, DB);

# check result
@test gv.G_system_mu ≈ -797.7491824869334