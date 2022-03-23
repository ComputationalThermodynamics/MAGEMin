# The full functionality of MAGEMin is wrapped in ../gen/magemin_library.jl
# Yet, the routines here make it more convenient to use this from julia

export init_MAGEMin, point_wise_minimization, get_bulk_rock, create_output


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
    gv      = LibMAGEMin.reset_gv(gv,z_b, DB.PP_ref_db, DB.SS_ref_db)
    LibMAGEMin.reset_cp(gv,z_b, DB.cp)
    LibMAGEMin.reset_SS(gv,z_b, DB.SS_ref_db)
    LibMAGEMin.reset_sp(gv, DB.sp)
    
    time = @elapsed  gv      = LibMAGEMin.ComputeEquilibrium_Point(EM_database, input_data, Mode, z_b,gv,	DB.PP_ref_db,DB.SS_ref_db,DB.cp);

    # Postprocessing (NOTE: we should switch off printing if gv.verbose=0)
    LibMAGEMin.ComputePostProcessing(0, z_b, gv, DB.PP_ref_db, DB.SS_ref_db, DB.cp)

    LibMAGEMin.fill_output_struct(	gv,	z_b, DB.PP_ref_db,DB.SS_ref_db,	DB.cp, DB.sp );

    # Transform results to a more convenient julia struct
    # To be done

    return gv, z_b, time*1e3
end

"""
    structure that holds the result of the pointwise minisation 

"""
struct output{len_ox,T}
    G_system::T             # G of system
    Gamma::Vector{T}        # Gamma
    P_kbar::T               # Pressure in kbar
    T_C::T                  # Temperature in Celcius
    iter::Int64             # number of iterations required
end

"""
    out = create_output(gv, z_b)

This extracts the output of a pointwise MAGEMin optimization and adds it into a julia structure
"""
function create_output(DB, gv)

    stb     =  unsafe_wrap(Vector{LibMAGEMin.stb_systems},DB.sp,1)[1]

    G_system = stb.G;
    Gamma    = unsafe_wrap(Vector{Cdouble},stb.gamma,gv.len_ox)
    P_kbar   = stb.P
    T_C      = stb.T-273.15 
    
    # Bulk rock info (total, melt, solid, fluid)
    bulk     = unsafe_wrap(Vector{Cdouble},stb.bulk,   gv.len_ox)
    bulk_M   = unsafe_wrap(Vector{Cdouble},stb.bulk_M, gv.len_ox)
    bulk_S   = unsafe_wrap(Vector{Cdouble},stb.bulk_S, gv.len_ox)
    bulk_F   = unsafe_wrap(Vector{Cdouble},stb.bulk_F, gv.len_ox)
    
    # Solid, melt, fluid fractions
    frac_M   = stb.frac_M      
    frac_S   = stb.frac_S
    frac_F   = stb.frac_F

    # Solid, melt, fluid densities
    rho_M   = stb.rho_M      
    rho_S   = stb.rho_S
    rho_F   = stb.rho_F

    # Numerics
    bulk_res_norm = stb.bulk_res_norm
    iter     =  gv.global_ite   

    # Stable assemblage info
    n_ph     =  stb.n_ph        # total # of stable phases
    n_PP     =  stb.n_PP        # number of pure phases
    n_SS     =  stb.n_SS        # number of solid solutions
    
    ph_frac  =  unsafe_wrap(Vector{Cdouble},stb.ph_frac,   n_ph)
    ph_type  =  unsafe_wrap(Vector{Cint},   stb.ph_type,   n_ph)
    ph_id    =  unsafe_wrap(Vector{Cint},   stb.ph_id  ,   n_ph)
    ph       =  unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, stb.ph, n_ph)) # stable phases
    
    # extract info about compositional variables of the solution models
    stb_SS_phase = unsafe_wrap(Vector{LibMAGEMin.stb_SS_phase},stb.SS,n_SS)
    stb_PP_phase = unsafe_wrap(Vector{LibMAGEMin.stb_PP_phase},stb.PP,n_PP)    

    # Names of oxides
    oxides   = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, stb.oxides, gv.len_ox))


   # return output{gv.len_ox, typeof(G_system)}(G_system, Gamma, P_kbar, T_C, iter)
end

