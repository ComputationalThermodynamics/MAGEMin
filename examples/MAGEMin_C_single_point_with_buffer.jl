using MAGEMin_C

buffer = "cco"

db          = "ig"
gv, z_b, DB, splx_data  = init_MAGEMin(db);
sys_in      =   "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
test        =   0         #KLB1
gv          =   use_predefined_bulk_rock(gv, test, db);

gv.buffer   = pointer(buffer)
gv.buffer_n = -1.0
gv.limitCaOpx = 1
gv.CaOpxLim  = 0.25

gv.verbose=-1
P           =   8.0
T           =   800.0
out         =   point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in);