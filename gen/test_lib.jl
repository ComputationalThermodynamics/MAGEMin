using MAGEMin_C

gv, z_b, DB, splx_data      = init_MAGEMin();


# Call optimization routine for given P & T & test 0
test        = 0;
sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
gv          = use_predefined_bulk_rock(gv, test);

P           = 8.0;
T           = 800.0;
gv.verbose  = -1;        # switch off any verbose
out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
print_info(out)
finalize_MAGEMin(gv,DB)


# Call optimization routine for given P & T & bulk_rock
using MAGEMin_C
gv, z_b, DB, splx_data      = init_MAGEMin();

# Specify the bulk rock composition
bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
bulk_in    = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
sys_in     = "wt"

bulk_rock  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,sys_in);
gv         = define_bulk_rock(gv, bulk_rock);

P           = 10.0;
T           = 1100.0;
gv.verbose  = -1;        # switch off any verbose
out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
finalize_MAGEMin(gv,DB)