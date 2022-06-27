# The full functionality of MAGEMin is wrapped in ../gen/magemin_library.jl
# Yet, the routines here make it more convenient to use this from julia

import Base.show

export  init_MAGEMin, finalize_MAGEMin, point_wise_minimization, get_bulk_rock, create_output,
        print_info, create_gmin_struct


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
    finalize_MAGEMin(gv,DB)
Cleans up the memory 
"""
function  finalize_MAGEMin(gv,DB)
    LibMAGEMin.FreeDatabases(gv, DB)

    nothing
end


"""
    bulk_rock = get_bulk_rock(gv, test=0)

Returns the pre-defined bulk rock composition of a given test

"""
function get_bulk_rock(gv, test)
    bulk_rock   = zeros(gv.len_ox)
    LibMAGEMin.get_bulk(bulk_rock, test, gv.len_ox)

    return bulk_rock
end

"""
    point_wise_minimization(sys_in::String,P::Float64,T::Float64, bulk_rock::Vector{Float64}, gv::LibMAGEMin.global_variables, DB::LibMAGEMin.Database)
    
Computes the stable assemblage at P[kbar], T[C] and for bulk rock composition bulk_rock
    
"""
function point_wise_minimization(sys_in::String,P::Float64,T::Float64, bulk_rock::Vector{Float64}, gv, DB)
    
    input_data      =   LibMAGEMin.io_data();                           # zero (not used actually)
	z_b             =   LibMAGEMin.initialize_bulk_infos(P, T);

    z_b.T           =   T + 273.15    # in K
    z_b.P           =   P

    Mode            = 0;
    gv.Mode         = Mode;
    gv.BR_norm      = 1.0; 								# reset bulk rock norm 			*/
    gv.global_ite   = 0;              					# reset global iteration 			*/
    gv.numPoint     = 1; 							    # the number of the current point */

    EM_database     = 0
 
    # Declare LP structures
    splx_data       =   LibMAGEMin.simplex_data(); 

    LibMAGEMin.init_simplex_A( pointer_from_objref(splx_data), gv)
    LibMAGEMin.init_simplex_B_em( pointer_from_objref(splx_data), gv)

    LibMAGEMin.convert_system_comp(gv,sys_in,z_b,bulk_rock)
    LibMAGEMin.norm_array(bulk_rock, gv.len_ox)

    # Perform the point-wise minimization after resetting variables
    gv      = LibMAGEMin.reset_gv(gv,z_b, DB.PP_ref_db, DB.SS_ref_db)
    z_b     = LibMAGEMin.reset_z_b_bulk(	gv,	 bulk_rock,	z_b	   )	

    LibMAGEMin.reset_simplex_A(pointer_from_objref(splx_data), z_b, gv)
    LibMAGEMin.reset_simplex_B_em(pointer_from_objref(splx_data), gv)

    LibMAGEMin.reset_cp(gv,z_b, DB.cp)
    LibMAGEMin.reset_SS(gv,z_b, DB.SS_ref_db)
    LibMAGEMin.reset_sp(gv, DB.sp)

    time = @elapsed  gv      = LibMAGEMin.ComputeEquilibrium_Point(EM_database, input_data, Mode, z_b,gv, pointer_from_objref(splx_data),	DB.PP_ref_db,DB.SS_ref_db,DB.cp);

    # Postprocessing (NOTE: we should switch off printing if gv.verbose=0)
    gv = LibMAGEMin.ComputePostProcessing(0, z_b, gv, DB.PP_ref_db, DB.SS_ref_db, DB.cp)

    # Fill structure
    LibMAGEMin.fill_output_struct(	gv,	z_b, DB.PP_ref_db,DB.SS_ref_db,	DB.cp, DB.sp );

    # Print output to screen 
    LibMAGEMin.PrintOutput(gv, 0, 1, DB, time, z_b);	

    # Transform results to a more convenient julia struct
    out = create_gmin_struct(DB, gv, time);
    
    # LibMAGEMin.FreeDatabases(gv, DB);

    return out
end

point_wise_minimization(sys_in::String,P::Integer, T::Integer, bulk_rock::Vector{Float64}, gv, DB) = point_wise_minimization(String(sys_in),Float64(P),Float64(T), bulk_rock::Vector{Float64}, gv, DB)
point_wise_minimization(sys_in::String,P::Float64, T::Integer, bulk_rock::Vector{Float64}, gv, DB) = point_wise_minimization(String(sys_in),Float64(P),Float64(T), bulk_rock::Vector{Float64}, gv, DB)
point_wise_minimization(sys_in::String,P::Integer, T::Float64, bulk_rock::Vector{Float64}, gv, DB) = point_wise_minimization(String(sys_in),Float64(P),Float64(T), bulk_rock::Vector{Float64}, gv, DB)

"""
    structure that holds the result of the pointwise minisation 

"""
struct gmin_struct{T,I}
    G_system::T             # G of system
    Gamma::Vector{T}        # Gamma
    P_kbar::T               # Pressure in kbar
    T_C::T                  # Temperature in Celcius

    # bulk rock composition:
    bulk::Vector{T}   
    bulk_M::Vector{T}   
    bulk_S::Vector{T}   
    bulk_F::Vector{T}   

    bulk_wt::Vector{T}   
    bulk_M_wt::Vector{T}   
    bulk_S_wt::Vector{T}   
    bulk_F_wt::Vector{T}   
    
    # Fractions:
    # Solid, melt, fluid fractions
    frac_M::T   
    frac_S::T
    frac_F::T

    frac_M_wt::T   
    frac_S_wt::T
    frac_F_wt::T
 
    # Solid, melt, fluid densities
    rho::T
    rho_M::T      
    rho_S::T
    rho_F::T

    # Phase fractions and type:
    n_PP::Int64                 # number of pure phases
    n_SS::Int64                 # number of solid solutions
    
    ph_frac::Vector{T}          # phase fractions
    ph_frac_wt::Vector{T}          # phase fractions
    ph_type::Vector{I}      # type of phase (SS or PP)
    ph_id::Vector{I}        # id of phase
    ph::Vector{String}          # Name of phase
    
    SS_vec::Vector{LibMAGEMin.SS_data}
    PP_vec::Vector{LibMAGEMin.PP_data}
    
    oxides::Vector{String}

    # Seismic velocity info
    Vp::T               # P-wave velocity
    Vs::T               # S-wave velocity
    bulkMod::T          # Elastic bulk modulus
    shearMod::T         # Elastic shear modulus

    # Numerics:
    iter::I             # number of iterations required
    bulk_res_norm::T    # bulk residual norm
    time_ms::T          # computational time for this point
    status::I           # status of calculations
end

"""
    out = create_gmin_struct(gv, z_b)

This extracts the output of a pointwise MAGEMin optimization and adds it into a julia structure
"""
function create_gmin_struct(DB, gv, time)

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

    # Bulk rock info (total, melt, solid, fluid)
    bulk_wt     = unsafe_wrap(Vector{Cdouble},stb.bulk_wt,   gv.len_ox)
    bulk_M_wt   = unsafe_wrap(Vector{Cdouble},stb.bulk_M_wt, gv.len_ox)
    bulk_S_wt   = unsafe_wrap(Vector{Cdouble},stb.bulk_S_wt, gv.len_ox)
    bulk_F_wt   = unsafe_wrap(Vector{Cdouble},stb.bulk_F_wt, gv.len_ox)
    
    # Solid, melt, fluid fractions
    frac_M   = stb.frac_M      
    frac_S   = stb.frac_S
    frac_F   = stb.frac_F
    # Solid, melt, fluid fractions
    frac_M_wt   = stb.frac_M_wt     
    frac_S_wt   = stb.frac_S_wt
    frac_F_wt   = stb.frac_F_wt

    # Solid, melt, fluid densities
    rho     = stb.rho
    rho_M   = stb.rho_M      
    rho_S   = stb.rho_S
    rho_F   = stb.rho_F

    # Stable assemblage info
    n_ph     =  stb.n_ph        # total # of stable phases
    n_PP     =  stb.n_PP        # number of pure phases
    n_SS     =  stb.n_SS        # number of solid solutions
    
    ph_frac  =  unsafe_wrap(Vector{Cdouble},stb.ph_frac,   n_ph)
    ph_frac_wt  =  unsafe_wrap(Vector{Cdouble},stb.ph_frac_wt,   n_ph)
    ph_type  =  unsafe_wrap(Vector{Cint},   stb.ph_type,   n_ph)
    ph_id    =  unsafe_wrap(Vector{Cint},   stb.ph_id  ,   n_ph)
    ph       =  unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, stb.ph, n_ph)) # stable phases
    
    # extract info about compositional variables of the solution models:
    SS_vec  = convert.(LibMAGEMin.SS_data, unsafe_wrap(Vector{LibMAGEMin.stb_SS_phase},stb.SS,n_SS))
    
    # Info about the endmembers:
    PP_vec  = convert.(LibMAGEMin.PP_data, unsafe_wrap(Vector{LibMAGEMin.stb_PP_phase},stb.PP,n_PP))    

    # Names of oxides:
    oxides   = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, stb.oxides, gv.len_ox))
    
    # Numerics
    bulk_res_norm   =  gv.BR_norm
    iter            =  gv.global_ite   
    time_ms         =  time*1000.0

    # Store all in output struct 
    out = gmin_struct{Float64,Int64}( G_system, Gamma, P_kbar, T_C, 
                bulk, bulk_M, bulk_S, bulk_F,
                bulk_wt, bulk_M_wt, bulk_S_wt, bulk_F_wt,  
                frac_M, frac_S, frac_F, 
                frac_M_wt, frac_S_wt, frac_F_wt, 
                rho, rho_M, rho_S, rho_F,   
                n_PP, n_SS,
                ph_frac, ph_frac_wt, ph_type, ph_id, ph,
                SS_vec,  PP_vec, 
                oxides,  
                stb.Vp, stb.Vs, stb.bulkMod, stb.shearMod,
                iter, bulk_res_norm, time_ms, stb.status)
    
   return out
end


# Print brief info about pointwise calculation result 
function show(io::IO, g::gmin_struct)  
    println(io, "Pressure          : $(g.P_kbar)      [kbar]")  
    println(io, "Temperature       : $(round(g.T_C,digits=4))    [Celcius]")  
    
    println(io, "     Stable phase | Fraction ")  
    for i=1:length(g.ph)
        println(io, "   $(lpad(g.ph[i],14," "))   $( round(g.ph_frac[i], digits=5)) ")  
    end
    for i=1:length(g.ph)
        println(io, "   $(lpad(g.ph[i],14," "))   $( round(g.ph_frac_wt[i], digits=5)) ")  
    end
    println(io, "Gibbs free energy : $(round(g.G_system,digits=6))  ($(g.iter) iterations; $(round(g.time_ms,digits=2)) ms)")  
    if g.status>0
        println(io, "WARNING: calculation did not converge ----------------------------")  
    end
    
end

"""
    print_info(g::gmin_struct)

Prints a more extensive overview of the simulation results
"""
function print_info(g::gmin_struct)
   
    println("Stable phases @ {$(round(g.P_kbar,digits=4)), $(round(g.T_C,digits=4))}  kbar/°C :")
    for i=1:length(g.ph)
        print("$(lpad(g.ph[i],6," ")) ")  
    end
    print("  \n \n")  
   
    # ==
    println("Compositional variables (solution phase):")
    for i=1:g.n_SS
        print("$(lpad(g.ph[i],15," ")) ")  
        for j=1:length(g.SS_vec[i].compVariables)
            print("$(lpad(round(g.SS_vec[i].compVariables[j],digits=5),8," ")) ")  
        end
        print("\n")
    end
    print("\n")
    # ==

    # ==
    println("End-members fraction (solution phase):")
    for i=1:g.n_SS

        print("                ")
        for j=1:length(g.SS_vec[i].emNames)
            print("$(lpad(g.SS_vec[i].emNames[j],8," ")) ")  
        end
        print("\n")
        print("$(lpad(g.ph[i],15," ")) ")  
        for j=1:length(g.SS_vec[i].emFrac)
            print("$(lpad(round(g.SS_vec[i].emFrac[j],digits=5),8," ")) ")  
        end
        print("\n")
    end
    print("\n")
    # ==
      # ==
      println("Oxide compositions [mol%] (normalized):")
      print("                ")    
      for i=1:length(g.oxides)
          print("$(lpad(g.oxides[i],8," ")) ")  
      end
      print("\n")
  
      print("$(lpad("SYS",15)) ")  
      for i=1:length(g.oxides)
          print("$(lpad(round(g.bulk[i],digits=5),8," ")) ")  
      end
      print("\n")

      for i=1:g.n_SS
          print("$(lpad(g.ph[i],15," ")) ")  
          for j=1:length(g.oxides)
              print("$(lpad(round(g.SS_vec[i].Comp[j],digits=5),8," ")) ")  
          end
          print("\n")
      end

      for i=1:g.n_PP
        print("$(lpad(g.ph[i],15," ")) ")  
        for j=1:length(g.oxides)
            print("$(lpad(round(g.PP_vec[i].Comp[j],digits=5),8," ")) ")  
        end
        print("\n")
      end
      print("\n")

      # ==
    # ==

    # ==
      # ==
      println("Oxide compositions [wt%] (normalized):")
      print("                ")    
      for i=1:length(g.oxides)
          print("$(lpad(g.oxides[i],8," ")) ")  
      end
      print("\n")
  
      print("$(lpad("SYS",15)) ")  
      for i=1:length(g.oxides)
          print("$(lpad(round(g.bulk_wt[i],digits=5),8," ")) ")  
      end
      print("\n")

      for i=1:g.n_SS
        print("$(lpad(g.ph[i],15," ")) ")  
        for j=1:length(g.oxides)
            print("$(lpad(round(g.SS_vec[i].Comp_wt[j],digits=5),8," ")) ")  
        end
        print("\n")
    end

      print("\n")
      for i=1:g.n_PP
        print("$(lpad(g.ph[i],15," ")) ")  
        for j=1:length(g.oxides)
            print("$(lpad(round(g.PP_vec[i].Comp_wt[j],digits=5),8," ")) ")  
        end
        print("\n")
      end
      print("\n")

      # ==
    # ==
    


    println("Stable mineral assemblage:")
    println("          phase     mode        f           G        V       Cp  rho[kg/m3]  Thermal_Exp BulkMod[GPa] ShearMod[GPa]   Vp[km/s]   Vs[km/s]    ")
    for i=1:g.n_SS
        print("$(lpad(g.ph[i],15," ")) ")  
        print("$(lpad(round(g.ph_frac[i],digits=5),8," ")) ")  
        print("$(lpad(round(g.SS_vec[i].f,digits=5),8," ")) ")  
        print("$(lpad(round(g.SS_vec[i].G,digits=5),8," ")) ")  
        print("$(lpad(round(g.SS_vec[i].V,digits=5),8," ")) ")  
        print("$(lpad(round(g.SS_vec[i].cp,digits=5),8," ")) ")  
        print("$(lpad(round(g.SS_vec[i].rho,digits=5),11," "))   ")  
        print("$(lpad(round(g.SS_vec[i].alpha,digits=5),8," "))   ")  
        print("$(lpad(round(g.SS_vec[i].bulkMod,digits=5),12," ")) ")  
        print("$(lpad(round(g.SS_vec[i].shearMod,digits=5),13," ")) ")  
        print("$(lpad(round(g.SS_vec[i].Vp,digits=5),10," ")) ")  
        print("$(lpad(round(g.SS_vec[i].Vs,digits=5),10," ")) ")  
        
        print("\n")
    end
    for i=1:g.n_PP
        print("$(lpad(g.ph[i],15," ")) ")  
        print("$(lpad(round(g.ph_frac[i],digits=5),8," ")) ")  
        print("$(lpad(round(g.PP_vec[i].f,digits=5),8," ")) ")  
        print("$(lpad(round(g.PP_vec[i].G,digits=5),8," ")) ")  
        print("$(lpad(round(g.PP_vec[i].V,digits=5),8," ")) ")  
        print("$(lpad(round(g.PP_vec[i].cp,digits=5),8," ")) ")  
        print("$(lpad(round(g.PP_vec[i].rho,digits=5),11," "))   ")  
        print("$(lpad(round(g.PP_vec[i].alpha,digits=5),8," "))   ")  
        print("$(lpad(round(g.PP_vec[i].bulkMod,digits=5),12," ")) ")  
        print("$(lpad(round(g.PP_vec[i].shearMod,digits=5),13," ")) ")  
        print("$(lpad(round(g.PP_vec[i].Vp,digits=5),10," ")) ")  
        print("$(lpad(round(g.PP_vec[i].Vs,digits=5),10," ")) ")  
        
        print("\n")
    end
    
    print("$(lpad("SYS",15," ")) ")  
    print("$(lpad(round(sum(g.ph_frac),digits=5),4," ")) ")  
    print("$(lpad(round(g.G_system,digits=5),24," ")) ")  
    print("$(lpad(round(g.rho,digits=5),29," ")) ")  
    print("$(lpad(round(g.bulkMod,digits=5),25," ")) ")  
    print("$(lpad(round(g.shearMod,digits=5),13," ")) ")  
    print("$(lpad(round(g.Vp,digits=5),10," ")) ")  
    print("$(lpad(round(g.Vs,digits=5),10," ")) ")  
    
    print("\n")
    print("\n")
    # ==

    println("Gamma (chemical potential of oxides):")  
    for i=1:length(g.oxides)
        println("  $(lpad(g.oxides[i],6," "))  $(rpad(round(g.Gamma[i],digits=5),15," ")) ")  
    end
    print("\n")

    println("delta Gibbs energy (G-hyperplane distance):")  
    for i=1:g.n_SS
        println("  $(lpad(g.ph[i],6," "))  $(rpad(g.SS_vec[i].deltaG,15," ")) ")  
    end
    for i=1:g.n_PP
        println("  $(lpad(g.ph[i],6," "))  $(rpad(g.PP_vec[i].deltaG,15," ")) ")  
    end
    
    print("\n")

    println("mass residual :  $(g.bulk_res_norm)")  
    println("# iterations  :  $(g.iter)")  
    println("status        :  $(g.status)")  
    

    print("\n")

end



