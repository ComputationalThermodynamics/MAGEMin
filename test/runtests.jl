# this tests the julia interface to MAGEMin
using Test

using MAGEMin_C

# Initialize database 
gv, DB = init_MAGEMin();

# Bulk rock composition for a particular test 
test      = 0;
bulk_rock = get_bulk_rock(gv, test)

# Call optimization routine for given P & T & bulk_rock
P           = 8.
T           = 800.
gv.verbose  = -1    # switch off any verbose
out         = point_wise_minimization(P,T, bulk_rock, gv, DB)
@show out

# check result
@test out.G_system ≈ -797.7491824869334
@test out.ph == [ "opx", "spn", "ol", "cpx"]
@test all(abs.(out.ph_frac - [ 0.24226960158631541, 0.027991246529842587, 0.5880694152724345, 0.1416697366114075])  .< 1e-6)

# print more detailed info about this point:
print_info(out)

