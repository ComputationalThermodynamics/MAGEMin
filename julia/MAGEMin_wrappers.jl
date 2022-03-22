# The full functionality of MAGEMin is wrapped in ../gen/magemin_library.jl
# Yet, the routines here make it more convenient to use this from julia

export init_MAGEMin, point_wise_minimization, get_bulk_rock


"""
    gv, DB = init_MAGEMin(;EM_database=0)

Initializes MAGEMin (including setting global options) and loads the Database.
"""
function  init_MAGEMin(;EM_database=0)

    gv = LibMAGEMin.global_variable_init();
    DB = LibMAGEMin.InitializeDatabases(gv, EM_database)

    return gv, DB
end


"""
    bulk_rock = get_bulk_rock(gv, test=0)

Returns the pre-defined bulk rock composition of a given test

"""
function get_bulk_rock(gv, test=0)

    bulk_rock =  zeros(gv.len_ox)
    LibMAGEMin.get_bulk(bulk_rock, test, gv.len_ox)

    return bulk_rock

end

"""
    point_wise_minimization(P::Float64,T::Float64, bulk_rock::Vector{Float64}, gv::LibMAGEMin.global_variables, DB::LibMAGEMin.Database)
    
Computes the stable assemblage at P[kbar], T[C] and for bulk rock composition bulk_rock
    
"""
function point_wise_minimization(P::Float64,T::Float64, bulk_rock::Vector{Float64}, gv, DB)
    
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

    EM_database = 0
    
    # Perform the point-wise minimization after resetting variables
    gv      = LibMAGEMin.reset_global_variables(gv,DB.PP_ref_db, DB.SS_ref_db, DB.cp)
    gv      = LibMAGEMin.reset_phases(gv, z_b, DB.PP_ref_db, DB.SS_ref_db, DB.cp)
    time = @elapsed  gv      = LibMAGEMin.ComputeEquilibrium_Point(EM_database, input_data, Mode, z_b,gv,	DB.PP_ref_db,DB.SS_ref_db,DB.cp);

    # Postprocessing (NOTE: we should switch off printing if gv.verbose=0)
    LibMAGEMin.ComputePostProcessing(0, z_b, gv, DB.PP_ref_db, DB.SS_ref_db, DB.cp)

    # Transform results to a more convenient julia struct
    # To be done

    return gv, z_b, time*1e3
end


