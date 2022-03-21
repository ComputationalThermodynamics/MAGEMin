# This tests
include("magemin_library.jl")


# Initialize 
gv = LibMAGEMin.global_variable_init();


gv.verbose=0
EM_database=0;
DB = LibMAGEMin.InitializeDatabases(gv, EM_database)

function Perform_Minimisation(P::Float64,T::Float64, bulk_rock::Vector{Float64}, gv::LibMAGEMin.global_variables, DB::LibMAGEMin.Database)
    
    LibMAGEMin.norm_array(bulk_rock, gv.len_ox);	                    # normalize bulk_rock
    
    input_data      =   LibMAGEMin.io_data();                           # zero (not used actually)
    z_b             =   LibMAGEMin.zeros_in_bulk(	bulk_rock, P, T);

    z_b.T           =   T + 273.15    # in K
    z_b.P           =   P
    
    Mode            = 0;
    gv.Mode         = Mode;
    gv.BR_norm      = 1.0; 								# reset bulk rock norm 			*/
    gv.global_ite   = 0;              					# reset global iteration 			*/
    gv.numPoint     = 1; 							    # the number of the current point */

    gv      = LibMAGEMin.reset_global_variables(gv,DB.PP_ref_db, DB.SS_ref_db, DB.cp)
    gv      = LibMAGEMin.reset_phases(gv, z_b, DB.PP_ref_db, DB.SS_ref_db, DB.cp)
    gv      = LibMAGEMin.ComputeEquilibrium_Point(EM_database, input_data, Mode, z_b,gv,	DB.PP_ref_db,DB.SS_ref_db,DB.cp);

    return gv, z_b
end


# 
P       = 8.
T       = 800.
test    = 0;
bulk_rock =  zeros(gv.len_ox)
LibMAGEMin.get_bulk(bulk_rock, test, gv.len_ox)

# Call optimization routine
gv, z_b = Perform_Minimisation(P,T, bulk_rock, gv, DB);


sgleP =0
time_taken = 1.;
rank=0;
LibMAGEMin.PrintOutput(gv, rank, sgleP, DB, time_taken, z_b);	

