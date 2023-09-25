# The full functionality of MAGEMin is wrapped in ../gen/magemin_library.jl
# Yet, the routines here make it more convenient to use this from julia
import Base.show
using Base.Threads: @threads 
using ProgressMeter

export  init_MAGEMin, finalize_MAGEMin, point_wise_minimization, convertBulk4MAGEMin, use_predefined_bulk_rock, define_bulk_rock, create_output,
        print_info, create_gmin_struct,
        single_point_minimization,
        multi_point_minimization, MAGEMin_Data,
        Initialize_MAGEMin, Finalize_MAGEMin
    

"""
Holds the MAGEMin databases & required structures for every thread
"""
struct MAGEMin_Data{TypeGV, TypeZB, TypeDB, TypeSplxData}
    db :: String
    gv :: TypeGV
    z_b :: TypeZB
    DB  :: TypeDB
    splx_data :: TypeSplxData
end

"""
    Dat = Initialize_MAGEMin(db = "ig"; verbose::Union{Bool, Int64} = true)

Initializes MAGEMin on one or more threads, for the database `db`. You can surpress all output with `verbose=false`. `verbose=true` will give a brief summary of the result, whereas `verbose=1` will give more details about the computations.
"""
function Initialize_MAGEMin(db = "ig"; verbose::Union{Int64,Bool} = 0)
    gv, z_b, DB, splx_data = init_MAGEMin(db);

    nt = Threads.nthreads()
    list_gv = Vector{typeof(gv)}(undef, nt)
    list_z_b = Vector{typeof(z_b)}(undef, nt)
    list_DB = Vector{typeof(DB)}(undef, nt)
    list_splx_data = Vector{typeof(splx_data)}(undef, nt)

    if isa(verbose,Bool)
        if verbose
            verbose=0
        else
            verbose=-1
        end
    end

    for id in 1:nt
        gv, z_b, DB, splx_data = init_MAGEMin(db)
        gv.verbose = verbose
        list_gv[id] = gv
        list_z_b[id] = z_b
        list_DB[id] = DB
        list_splx_data[id] = splx_data
    end

    return MAGEMin_Data(db, list_gv, list_z_b, list_DB, list_splx_data)
end


"""
    Finalize_MAGEMin(dat::MAGEMin_Data)
Finalizes MAGEMin and clears variables
"""
function Finalize_MAGEMin(dat::MAGEMin_Data)
    for id in 1:Threads.nthreads()
        gv = dat.gv[id]
        DB = dat.DB[id]
        LibMAGEMin.FreeDatabases(gv, DB)

        # These are indeed not freed yet (same with C-code), which should be added for completion
        # They are rather small structs compared to the others
        z_b = dat.z_b[id]
        splx_data = dat.splx_data[id]

     end
     return nothing
end


# Left for backwards compatibility
function  init_MAGEMin(db="ig")

    z_b         = LibMAGEMin.bulk_infos()
    gv          = LibMAGEMin.global_variables()
    splx_data   = LibMAGEMin.simplex_data(); 
    DB          = LibMAGEMin.Database()
    gv          = LibMAGEMin.global_variable_alloc( pointer_from_objref(z_b))
    
    if db == "ig"
        gv.EM_database = 2
    elseif db == "mp"
        gv.EM_database = 0
    elseif db == "um"
        gv.EM_database = 4
    else
        print("Database not implemented...\n")
    end

    gv          = LibMAGEMin.global_variable_init(gv, pointer_from_objref(z_b))
    DB          = LibMAGEMin.InitializeDatabases(gv, gv.EM_database)

    LibMAGEMin.init_simplex_A( pointer_from_objref(splx_data), gv)
    LibMAGEMin.init_simplex_B_em( pointer_from_objref(splx_data), gv)

    return gv, z_b, DB, splx_data
end

# left here for backwards compatibility
function finalize_MAGEMin(gv,DB)
    LibMAGEMin.FreeDatabases(gv,DB)
    return nothing
end

# wrapper for single point minimization
function single_point_minimization(     P::Float64,
                                        T::Float64,
                                        MAGEMin_db::MAGEMin_Data;  
                                        test::Int64 = 0, # if using a build-in test case
                                        X::Union{Nothing, Vector{_T}, Vector{Vector{_T}}} = nothing,
                                        Xoxides     = Vector{String},
                                        sys_in      = "mol",
                                        progressbar = true        # show a progress bar or not?
                                        ) where _T <: Float64

    P = [P];
    T = [T];
    if ~isnothing(X)
        X = [X]
    end


    Out_PT     =   multi_point_minimization(P,
                                            T,
                                            MAGEMin_db,
                                            test=test,
                                            X=X,
                                            Xoxides=Xoxides,
                                            sys_in=sys_in,
                                            progressbar=progressbar);

    return Out_PT[1]
end


"""
    Out_PT = multi_point_minimization(P:Vector{_T}, T::Vector, MAGEMin_db::MAGEMin_Data; sys_in="mol", test=0, X::Union{Nothing, Vector, Vector{Vector}}=nothing, progressbar=true)

Perform (parallel) MAGEMin calculations for a range of points as a function of pressure `P`, temperature `T` and/or composition `X`. The database `MAGEMin_db` must be initialised before calling the routine.
The bulk-rock composition can either be set to be one of the pre-defined build-in test cases, or can be specified specifically by passing `X`, `Xodides` and `sys_in` (that specifies whether the input is in "mol" or "wt").

Below a few examples:

Example 1 - build-in test vs. pressure and temperature
===
```julia
julia> data = Initialize_MAGEMin("ig", verbose=false);
julia> n = 10
julia> P = rand(8:40.0,n)
julia> T = rand(800:1500.0,n)
julia> out = multi_point_minimization(P, T, data, test=0)
julia> Finalize_MAGEMin(data)
```

Example 2 - Specify constant bulk rock composition for all points:
===
```julia
julia> data = Initialize_MAGEMin("ig", verbose=false);
julia> n = 10
julia> P = fill(10.0,n)
julia> T = fill(1100.0,n)
julia> Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
julia> X = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
julia> sys_in = "wt"    
julia> out = multi_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
julia> Finalize_MAGEMin(data)
```

Example 3 - Different bulk rock composition for different points
===
```julia
julia> data = Initialize_MAGEMin("ig", verbose=false);
julia> P = [10.0, 20.0]
julia> T = [1100.0, 1200]
julia> Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
julia> X1 = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
julia> X2 = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
julia> X = [X1,X2]
julia> sys_in = "wt"    
julia> out = multi_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
julia> Finalize_MAGEMin(data)
```

Activating multithreading on julia
===

To take advantage of multithreading, you need to start julia from the terminal with:
```julia
\$julia -t auto
```
which will automatically use all threads on your machine. Alternatively, use `julia -t 4` to start it on 4 threads.
If you are interested to see what you can do on your machine type:
```
julia> versioninfo()
``` 

"""
function multi_point_minimization(  P::Vector{Float64},
                                    T::Vector{Float64},
                                    MAGEMin_db::MAGEMin_Data;  
                                    test        = 0, # if using a build-in test case
                                    X::Union{Nothing, Vector{_T}, Vector{Vector{_T}}}=nothing,
                                    Xoxides     = Vector{String},
                                    sys_in      = "mol",
                                    progressbar = true        # show a progress bar or not?
                                    ) where _T <: Float64

    # Set the compositional info 
    CompositionType::Int64 = 0;
    if isnothing(X)
        # Use one of the build-in tests 
        # Create thread-local data
        for i in 1:Threads.nthreads()
            MAGEMin_db.gv[i] = use_predefined_bulk_rock(MAGEMin_db.gv[i], test, MAGEMin_db.db)
        end

        CompositionType = 0;    # build-in tests
    else
        if isa(X,Vector{Float64})
            # same bulk rock composition for the full diagram
            @assert length(X) == length(Xoxides)

            # Set the bulk rock composition for all points
            for i in 1:Threads.nthreads()
                MAGEMin_db.gv[i] = define_bulk_rock(MAGEMin_db.gv[i], X, Xoxides, sys_in, MAGEMin_db.db);
            end
            CompositionType = 1;    # specified bulk composition for all points
        else
            @assert length(X) == length(P)
            CompositionType = 2;    # different bulk rock composition for every point
        end

    end

    # initialize vectors
    Out_PT = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef, length(P))

    # Currently, there seem to be some type instabilities or something else so that
    # some compilation happens in the threaded loop below. This interferes badly
    # in some weird way with (libsc, p4est, t8code) - in particular on Linux where
    # we get segfaults. To avoid this, we force serial compilation by calling MAGEMin
    # once before the loop.
    let id = 1
        gv          = MAGEMin_db.gv[id]
        z_b         = MAGEMin_db.z_b[id]
        DB          = MAGEMin_db.DB[id]
        splx_data   = MAGEMin_db.splx_data[id]
        point_wise_minimization(P[1], T[1], gv, z_b, DB, splx_data, sys_in)
    end

    # main loop
    if progressbar
        progr = Progress(length(P), desc="Computing $(length(P)) points...") # progress meter
    end
    @threads :static for i in eachindex(P)
        # Get thread-local buffers. As of Julia v1.9, a dynamic scheduling of
        # the threads is the default setting. To avoid task migration and the
        # resulting concurrency issues, we restrict the loop to static scheduling.
        id          = Threads.threadid()
        gv          = MAGEMin_db.gv[id]
        z_b         = MAGEMin_db.z_b[id]
        DB          = MAGEMin_db.DB[id]
        splx_data   = MAGEMin_db.splx_data[id]

        if CompositionType==2
            # different bulk-rock composition for every point - specify it here
            gv = define_bulk_rock(gv, X[i], Xoxides, sys_in, MAGEMin_db.db);
        end

        # compute a new point using a ccall
        out         = point_wise_minimization(P[i], T[i], gv, z_b, DB, splx_data)
        Out_PT[i]   = deepcopy(out)  

        if progressbar
            next!(progr)
        end
    end
    if progressbar
        finish!(progr)
    end
    return Out_PT
end



"""
    bulk_rock = use_predefined_bulk_rock(gv, test=-1, db="ig")

Returns the pre-defined bulk rock composition of a given test
"""
function use_predefined_bulk_rock(gv, test=0, db="ig")

    if db == "ig"
        gv.test = test
        gv = LibMAGEMin.get_bulk_igneous(gv)
        LibMAGEMin.norm_array(gv.bulk_rock, gv.len_ox)
    elseif db == "mp"
        gv.test = test
        gv = LibMAGEMin.get_bulk_metapelite(gv)
        LibMAGEMin.norm_array(gv.bulk_rock, gv.len_ox)
    elseif db == "um"
        gv.test = test
        gv = LibMAGEMin.get_bulk_ultramafic(gv)
        LibMAGEMin.norm_array(gv.bulk_rock, gv.len_ox)
    else
        print("Database not implemented...\n")
    end

    return gv
end

"""
    data = use_predefined_bulk_rock(data::MAGEMin_Data, test=0)
Returns the pre-defined bulk rock composition of a given test
"""
function use_predefined_bulk_rock(data::MAGEMin_Data, test=0)
    nt = Threads.nthreads()
    for id in 1:nt
        data.gv[id] =  use_predefined_bulk_rock(data.gv[id], test, data.db)
    end
    return data
end


function define_bulk_rock(gv, bulk_in, bulk_in_ox, sys_in,db)

    bulk_rock       = convertBulk4MAGEMin(bulk_in,bulk_in_ox,sys_in,db)     # conversion changes the system unit to mol
    gv.bulk_rock    = pointer(bulk_rock)                                    # copy the bulk-rock

    LibMAGEMin.norm_array(gv.bulk_rock, gv.len_ox)

    return gv
end


function normalize(vector::Vector{Float64})
    return vector ./ sum(vector)
end


"""
convertBulk4MAGEMin( bulk_in, bulk_in_ox, sys_in)
    
receives bulk-rock composition in [mol,wt] fraction and associated oxide list and sends back bulk-rock composition converted for MAGEMin use
    
"""
function convertBulk4MAGEMin(bulk_in::Vector{Float64},bulk_in_ox::Vector{String},sys_in::String,db::String);

	ref_ox          = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "MnO"; "H2O"];
	ref_MolarMass   = [60.08; 101.96; 56.08; 40.30; 71.85; 79.85; 94.2; 61.98; 79.88; 16.0; 151.99; 70.937; 18.015];      #Molar mass of oxides

    if db == "ig"
	    MAGEMin_ox      = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"];
    elseif db == "mp"
        MAGEMin_ox      = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "H2O"];
    else
        print("Database not implemented...\n")
    end
    
	MAGEMin_bulk    = zeros(11);
    bulk            = zeros(11);
	# convert to mol, if system unit = wt
	if sys_in == "wt"
		for i=1:length(bulk_in_ox)
            id = findall(ref_ox .== bulk_in_ox[i]);
			bulk[i] = bulk_in[i]/ref_MolarMass[id[1]];
		end
	end
	bulk = normalize(bulk); 

	for i=1:length(MAGEMin_ox)
        id = findall(bulk_in_ox .== MAGEMin_ox[i]);
		if isempty(id) == 0
			MAGEMin_bulk[i] = bulk[id[1]];
		end
	end
    idFe2O3 = findall(bulk_in_ox .== "Fe2O3");

    if isempty(idFe2O3) == 0
        idFeO = findall(MAGEMin_ox .== "FeO");
        MAGEMin_bulk[idFeO[1]] += bulk[idFe2O3[1]]*2.0;

        idO = findall(MAGEMin_ox .== "O");
        MAGEMin_bulk[idO[1]] += bulk[idFe2O3[1]];
    end

    MAGEMin_bulk .= normalize(MAGEMin_bulk);

    idNonH2O = findall(MAGEMin_ox .!= "H2O");

    id0 = findall(MAGEMin_bulk[idNonH2O] .== 0.0)
    if isempty(id0) == 0
        MAGEMin_bulk[id0] .= 1e-4;
    end

    MAGEMin_bulk .= normalize(MAGEMin_bulk)*100.0;

    return MAGEMin_bulk;
end





"""
    point_wise_minimization(P::Float64,T::Float64, gv, z_b, DB, splx_data, sys_in::String="mol")
    
Computes the stable assemblage at `P` [kbar], `T` [C] and for a given bulk rock composition
    

# Example 1

This is an example of how to use it for a predefined bulk rock composition:
```julia
julia> db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
julia> gv, z_b, DB, splx_data      = init_MAGEMin(db);
julia> test        = 0;
julia> sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
julia> gv          = use_predefined_bulk_rock(gv, test, db)
julia> P           = 8.0;
julia> T           = 800.0;
julia> gv.verbose  = -1;        # switch off any verbose
julia> out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
Pressure          : 8.0      [kbar]
Temperature       : 800.0    [Celcius]
     Stable phase | Fraction (mol 1 atom basis) 
              opx   0.24229 
               ol   0.58808 
              cpx   0.14165 
              spn   0.02798 
     Stable phase | Fraction (wt fraction) 
              opx   0.23908 
               ol   0.58673 
              cpx   0.14583 
              spn   0.02846 
Gibbs free energy : -797.749183  (26 iterations; 94.95 ms)
Oxygen fugacity          : 9.645393319147175e-12
```

# Example 2
And here a case in which you specify your own bulk rock composition. 
We convert that in the correct format, using the `convertBulk4MAGEMin` function. 
```julia
julia> db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
julia> gv, z_b, DB, splx_data      = init_MAGEMin(db);
julia> bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
julia> bulk_in    = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
julia> sys_in     = "wt"
julia> gv         = define_bulk_rock(gv, bulk_in, bulk_in_ox, sys_in, db);
julia> P,T         = 10.0, 1100.0;
julia> gv.verbose  = -1;        # switch off any verbose
julia> out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
Pressure          : 10.0      [kbar]
Temperature       : 1100.0    [Celcius]
     Stable phase | Fraction (mol 1 atom basis) 
             pl4T   0.01114 
              liq   0.74789 
              cpx   0.21862 
              opx   0.02154 
     Stable phase | Fraction (wt fraction) 
             pl4T   0.01168 
              liq   0.72576 
              cpx   0.23872 
              opx   0.02277 
Gibbs free energy : -907.27887  (47 iterations; 187.11 ms)
Oxygen fugacity          : 0.02411835177808492
julia> finalize_MAGEMin(gv,DB)
```

"""
function point_wise_minimization(P::Float64,T::Float64, gv, z_b, DB, splx_data)
    
    input_data      =   LibMAGEMin.io_data();                           # zero (not used actually)

    z_b.T           =   T + 273.15    # in K
    z_b.P           =   P

    gv.numPoint     = 1; 							    # the number of the current point */

    # Perform the point-wise minimization after resetting variables
    gv      = LibMAGEMin.reset_gv(gv,z_b, DB.PP_ref_db, DB.SS_ref_db)
    z_b     = LibMAGEMin.reset_z_b_bulk(	gv,	z_b	   )	

    LibMAGEMin.reset_simplex_A(pointer_from_objref(splx_data), z_b, gv)
    LibMAGEMin.reset_simplex_B_em(pointer_from_objref(splx_data), gv)

    LibMAGEMin.reset_cp(gv,z_b, DB.cp)
    LibMAGEMin.reset_SS(gv,z_b, DB.SS_ref_db)
    LibMAGEMin.reset_sp(gv, DB.sp)

    time = @elapsed  gv      = LibMAGEMin.ComputeEquilibrium_Point(gv.EM_database, input_data, z_b, gv, pointer_from_objref(splx_data),	DB.PP_ref_db,DB.SS_ref_db,DB.cp);

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

point_wise_minimization(P::Number,T::Number, gv, z_b, DB, splx_data) = point_wise_minimization(Float64(P),Float64(T), gv, z_b, DB, splx_data)


point_wise_minimization(P::Number,T::Number, gv::MAGEMin_C.LibMAGEMin.global_variables, z_b::MAGEMin_C.LibMAGEMin.bulk_infos, DB::MAGEMin_C.LibMAGEMin.Database, splx_data::MAGEMin_C.LibMAGEMin.simplex_datas, sys_in::String) = point_wise_minimization(P,T, gv, z_b, DB, splx_data)

"""
    out = point_wise_minimization(P::Number,T::Number, data::MAGEMin_Data)

Performs a point-wise optimization for a given pressure `P` and temperature `T` foir the data specified in the MAGEMin database `MAGEMin_Data` (where also compoition is specified)
"""
point_wise_minimization(P::Number,T::Number, data::MAGEMin_Data) = point_wise_minimization(P,T, data.gv[1], data.z_b[1], data.DB[1], data.splx_data[1]) 



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

    # Oxygen fugacity
    fO2::T

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
    Vp_S::T               # P-wave velocity of solid aggregate
    Vs_S::T               # S-wave velocity of solid aggregate
    bulkMod::T          # Elastic bulk modulus
    shearMod::T         # Elastic shear modulus
    bulkModulus_M::T          # Elastic bulk modulus
    bulkModulus_S::T          # Elastic bulk modulus
    shearModulus_S::T         # Elastic shear modulus

    # thermodynamic properties
    entropy::T          # entropy
    enthalpy::T         # enthalpy

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

    stb      = unsafe_load(DB.sp)

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

    # Oxygen fugacity
    fO2     = stb.fO2

    # thermodynamic properties
    entropy = stb.entropy
    enthalpy= stb.enthalpy

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
                fO2, 
                n_PP, n_SS,
                ph_frac, ph_frac_wt, ph_type, ph_id, ph,
                SS_vec,  PP_vec, 
                oxides,  
                stb.Vp, stb.Vs, stb.Vp_S, stb.Vs_S, stb.bulkMod, stb.shearMod, stb.bulkModulus_M,  stb.bulkModulus_S, stb.shearModulus_S,
                entropy, enthalpy,
                iter, bulk_res_norm, time_ms, stb.status)
    
   return out
end


# Print brief info about pointwise calculation result 
function show(io::IO, g::gmin_struct)  
    println(io, "Pressure          : $(g.P_kbar)      [kbar]")  
    println(io, "Temperature       : $(round(g.T_C,digits=4))    [Celcius]")  
    
    println(io, "     Stable phase | Fraction (mol 1 atom basis) ")  
    for i=1:length(g.ph)
        println(io, "   $(lpad(g.ph[i],14," "))   $( round(g.ph_frac[i], digits=5)) ")  
    end
    println(io, "     Stable phase | Fraction (wt fraction) ")  
    for i=1:length(g.ph)
        println(io, "   $(lpad(g.ph[i],14," "))   $( round(g.ph_frac_wt[i], digits=5)) ")  
    end
    println(io, "Gibbs free energy : $(round(g.G_system,digits=6))  ($(g.iter) iterations; $(round(g.time_ms,digits=2)) ms)")  
    if g.status>0
        println(io, "WARNING: calculation did not converge ----------------------------")  
    end
    println(io, "Oxygen fugacity          : $(g.fO2)")  

    
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
    println("End-members fraction [mol% 1 atom basis](solution phase):")
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
    println("End-members fraction [wt%](solution phase):")
    for i=1:g.n_SS

        print("                ")
        for j=1:length(g.SS_vec[i].emNames)
            print("$(lpad(g.SS_vec[i].emNames[j],8," ")) ")  
        end
        print("\n")
        print("$(lpad(g.ph[i],15," ")) ")  
        for j=1:length(g.SS_vec[i].emFrac_wt)
            print("$(lpad(round(g.SS_vec[i].emFrac_wt[j],digits=5),8," ")) ")  
        end
        print("\n")
    end
    print("\n")
    # ==
      # ==
      println("Oxide compositions [mol% 1 atom basis] (normalized):")
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
    println("          phase  mode[mol1at] mode[wt]        f           G        V       Cp  rho[kg/m3]  Thermal_Exp Entropy[J/K] Enthalpy[J] BulkMod[GPa] ShearMod[GPa]   Vp[km/s]   Vs[km/s]")
    for i=1:g.n_SS
        print("$(lpad(g.ph[i],15," ")) ")  
        print("$(lpad(round(g.ph_frac[i],digits=5),13," ")) ")  
        print("$(lpad(round(g.ph_frac_wt[i],digits=5),8," ")) ")  
        print("$(lpad(round(g.SS_vec[i].f,digits=5),8," ")) ")  
        print("$(lpad(round(g.SS_vec[i].G,digits=5),8," ")) ")  
        print("$(lpad(round(g.SS_vec[i].V,digits=5),8," ")) ")  
        print("$(lpad(round(g.SS_vec[i].cp,digits=5),8," ")) ")  
        print("$(lpad(round(g.SS_vec[i].rho,digits=5),11," "))   ")  
        print("$(lpad(round(g.SS_vec[i].alpha,digits=5),8," "))   ")  
        print("$(lpad(round(g.SS_vec[i].entropy,digits=5),8," "))   ")  
        print("$(lpad(round(g.SS_vec[i].enthalpy,digits=5),8," "))   ")  
        print("$(lpad(round(g.SS_vec[i].bulkMod,digits=5),12," ")) ")  
        print("$(lpad(round(g.SS_vec[i].shearMod,digits=5),13," ")) ")  
        print("$(lpad(round(g.SS_vec[i].Vp,digits=5),10," ")) ")  
        print("$(lpad(round(g.SS_vec[i].Vs,digits=5),10," ")) ")  
        
        print("\n")
    end
    for i=1:g.n_PP
        print("$(lpad(g.ph[i],15," ")) ")  
        print("$(lpad(round(g.ph_frac[i],digits=5),13," ")) ")  
        print("$(lpad(round(g.ph_frac_wt[i],digits=5),8," ")) ") 
        print("$(lpad(round(g.PP_vec[i].f,digits=5),8," ")) ")  
        print("$(lpad(round(g.PP_vec[i].G,digits=5),8," ")) ")  
        print("$(lpad(round(g.PP_vec[i].V,digits=5),8," ")) ")  
        print("$(lpad(round(g.PP_vec[i].cp,digits=5),8," ")) ")  
        print("$(lpad(round(g.PP_vec[i].rho,digits=5),11," "))   ")  
        print("$(lpad(round(g.PP_vec[i].alpha,digits=5),8," "))   ")  
        print("$(lpad(round(g.PP_vec[i].entropy,digits=5),8," "))   ")  
        print("$(lpad(round(g.PP_vec[i].enthalpy,digits=5),8," "))   ")  
        print("$(lpad(round(g.PP_vec[i].bulkMod,digits=5),12," ")) ")  
        print("$(lpad(round(g.PP_vec[i].shearMod,digits=5),13," ")) ")  
        print("$(lpad(round(g.PP_vec[i].Vp,digits=5),10," ")) ")  
        print("$(lpad(round(g.PP_vec[i].Vs,digits=5),10," ")) ")  
        
        print("\n")
    end
    
    print("$(lpad("SYS",15," ")) ")  
    print("$(lpad(round(sum(g.ph_frac),digits=5),13," ")) ")  
    print("$(lpad(round(sum(g.ph_frac_wt),digits=5),8," ")) ")  
    print("$(lpad(round(g.G_system,digits=5),20," ")) ")  
    print("$(lpad(round(g.rho,digits=5),29," ")) ") 
    print("$(lpad(round(g.entropy,digits=5),21," ")) ")  
    print("$(lpad(round(g.enthalpy,digits=5),13," ")) ")  
    print("$(lpad(round(g.bulkMod,digits=5),13," ")) ")  
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
