# this tests the julia interface to MAGEMin
using Test

# Load MAGEMin (needs to be loaded from main directory to pick up correct
# library in case it is locally compiled). This is handled by the logic in
# runtests.jl
using MAGEMin_C         # load MAGEMin (needs to be loaded from main directory to pick up correct library in case it is locally compiled)


# generic test for sb11 database
data        =   Initialize_MAGEMin("sb11", verbose=true);
test        =   1         #KLB1
data        =   use_predefined_bulk_rock(data, test);
P           =   80.0
T           =   800.0
out         =   point_wise_minimization(P,T, data);

@test sort(out.ph) == sort(["gtmj", "hpcpx", "ol" ,"cpx"])
Finalize_MAGEMin(data)

# generic test for thermocalc database
data        =   Initialize_MAGEMin("ig", verbose=true);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);
P           =   8.0
T           =   800.0
out         =   point_wise_minimization(P,T, data);


@show out

@test out.G_system ≈ -797.7397668057848
@test sort(out.ph) == sort(["spn", "cpx",  "opx", "ol"])
@test abs(out.s_cp[1] - 1208.466551730128) < 2.0

# print more detailed info about this point:
print_info(out)
Finalize_MAGEMin(data)


# Initialize database  - new way
data        =   Initialize_MAGEMin("mtl", verbose=true);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   180.0
T           =   1400.0
out         =   point_wise_minimization(P,T, data);

@test sort(out.ph) == sort(["g", "ring", "wad"])
Finalize_MAGEMin(data)



@testset "test light output calculation" begin
    # Without a buffer at 1100.0 C
    data    = Initialize_MAGEMin("ig", verbose=false);
    P,T     = 10.0, 600.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [78.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_hT  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, light=true)

    @test sort(out_hT.ph_id_db) == sort([1, 1, 3, 4, 6, 10, 12])
    @test out_hT.ph_type == [1, 1, 1, 1, 1, 0, 0]

    Finalize_MAGEMin(data)
end

# Tests from L. Candioti - ETH - Oct 2024
@testset "test mass conservation" begin

    # Without a buffer at 1100.0 C
    data    = Initialize_MAGEMin("ig", verbose=false);
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_hT  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
    Δρ_hT = abs( out_hT.rho - ((out_hT.frac_M_wt * out_hT.rho_M + out_hT.frac_S_wt * out_hT.rho_S )) )
    @test Δρ_hT < 1e-10
    Finalize_MAGEMin(data)


    # Without a buffer at 800.0 C
    data    = Initialize_MAGEMin("ig", verbose=false);
    P,T     = 10.0, 800.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_lT  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
    Δρ_lT = abs( out_lT.rho - ((out_lT.frac_M_wt * out_lT.rho_M + out_lT.frac_S_wt * out_lT.rho_S )) )
    @test Δρ_lT < 1e-10
    Finalize_MAGEMin(data)

    # With a buffer at 1100.0 C
    data    = Initialize_MAGEMin("ig", buffer = "nno", verbose=false);
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_BhT = single_point_minimization(P, T, data, X=X, B=0.0, Xoxides=Xoxides, sys_in=sys_in)
    Δρ_BhT = abs( out_BhT.rho - ((out_BhT.frac_M_wt * out_BhT.rho_M + out_BhT.frac_S_wt * out_BhT.rho_S )) )
    @test Δρ_BhT < 1e-10
    Finalize_MAGEMin(data)
    
    # With a buffer at 800.0 C
    data    = Initialize_MAGEMin("ig", verbose=false);
    P,T     = 10.0, 800.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_BlT = single_point_minimization(P, T, data, X=X, B=0.0, Xoxides=Xoxides, sys_in=sys_in)
    Δρ_BlT = abs( out_BlT.rho - ((out_BlT.frac_M_wt * out_BlT.rho_M + out_BlT.frac_S_wt * out_BlT.rho_S )) )
    @test Δρ_BlT < 1e-10
    Finalize_MAGEMin(data)
end

@testset "test activity buffers" begin
    # Initialize database  - new way
    data        =   Initialize_MAGEMin("mp", verbose=true, buffer="aH2O");
    test        =   0        
    data        =   use_predefined_bulk_rock(data, test);

    # Call optimization routine for given P & T & bulk_rock
    P           =   8.0
    T           =   400.0
    out         =   point_wise_minimization(P,T, data, buffer_n=0.6);
    @test sort(out.ph) == sort(["chl", "sp", "mu", "mu", "fsp", "ep", "q", "ru", "aH2O"])
    Finalize_MAGEMin(data)

    # Initialize database  - new way
    data        =   Initialize_MAGEMin("mp", verbose=true, buffer="aTiO2");
    test        =   0        
    data        =   use_predefined_bulk_rock(data, test);

    # Call optimization routine for given P & T & bulk_rock
    P           =   8.0
    T           =   400.0
    out         =   point_wise_minimization(P,T, data, buffer_n=0.6);
    @test sort(out.ph) == sort(["H2O", "aTiO2", "chl", "ep", "fsp", "ilm", "mu", "mu", "q"])
    Finalize_MAGEMin(data)

    # Initialize database  - new way
    data        =   Initialize_MAGEMin("ig", verbose=true, buffer="aTiO2");
    test        =   0        
    data        =   use_predefined_bulk_rock(data, test);

    # Call optimization routine for given P & T & bulk_rock
    P           =   8.0
    T           =   1200.0
    out         =   point_wise_minimization(P,T, data, buffer_n=0.1);
    @test sort(out.ph) == sort(["aTiO2", "cpx", "fsp", "liq", "ol", "opx"])

    Finalize_MAGEMin(data)


    data    = Initialize_MAGEMin("ig", verbose=false, buffer="qfm");
    
    # One bulk rock for all points
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
    @test sort(out.ph) == sort(["cpx", "liq", "opx", "qfm"])

    Finalize_MAGEMin(data)

end

@testset "test zr saturation" begin
    data    = Initialize_MAGEMin("mp", verbose=false);

    P,T     = 6.0, 930.0
    Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
    X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
    sys_in  = "wt"
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

    # use compo from experiment from Boehnke et al., 2013
    compo1 = [54.2, 0.5, 16.9, 4.1, 0.0, 2, 7.6, 2.3, 0.8, 0, 0]
    bulk_melt = convertBulk4MAGEMin(compo1, Xoxides, "wt", "mp")[1]

    out.bulk_M .= bulk_melt
    zr_sat_B = MAGEMin_C.zirconium_saturation(out, model="B")
    zr_sat_WH = MAGEMin_C.zirconium_saturation(out, model="WH")

    @test zr_sat_B ≈ 1403.8755429428836 rtol=1e-5
    @test zr_sat_WH ≈ 1059.5976323423222 rtol=1e-5

    # test crisp and berry 2022, use compo from their example in the calculator from their paper
    P,T     = 20.0, 750.0
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

    bulk_melt = [61.26, 0, 13.12, 1.33, 0, 0.45, 2.51, 2.26, 2.15, 15.00, 0]
    # convert bulk_melt from wt to mol of oxides
    bulk_melt = convertBulk4MAGEMin(bulk_melt, Xoxides, "wt", "mp")[1]
    out.bulk_M .= bulk_melt

    zr_sat = MAGEMin_C.zirconium_saturation(out, model="CB")

    @test zr_sat ≈ 65.83158859091596 rtol=1e-5
end

@testset "test normalization" begin

    # Initialize database  - new way
    data        =   Initialize_MAGEMin("ig", verbose=true);
    test        =   5         #KLB1
    data        =   use_predefined_bulk_rock(data, test);
    P           =   8.0
    T           =   800.0
    out         =   point_wise_minimization(P,T, data);
    @test  out.frac_M    + out.frac_S    + out.frac_F        ≈ 1.0
    @test  out.frac_M_wt + out.frac_S_wt + out.frac_F_wt     ≈ 1.0
    @test  sum(out.bulk_M)                                   ≈ 1.0
    @test  sum(out.bulk_F)                                   ≈ 1.0
    @test  sum(out.bulk_S)                                   ≈ 1.0
    @test  sum(out.bulk_M_wt)                                ≈ 1.0
    @test  sum(out.bulk_F_wt)                                ≈ 1.0
    @test  sum(out.bulk_S_wt)                                ≈ 1.0

    test        =   0         #KLB1
    data        =   use_predefined_bulk_rock(data, test);

    P           =   8.0
    T           =   1500.0
    out         =   point_wise_minimization(P,T, data);
    @test  out.frac_M    + out.frac_S    + out.frac_F        ≈ 1.0
    @test  out.frac_M_wt + out.frac_S_wt + out.frac_F_wt     ≈ 1.0
    @test  sum(out.bulk_M)                                   ≈ 1.0
    @test  sum(out.bulk_S)                                   ≈ 1.0
    @test  sum(out.bulk_M_wt)                                ≈ 1.0
    @test  sum(out.bulk_S_wt)                                ≈ 1.0

    P           =   8.0
    T           =   800.0
    out         =   point_wise_minimization(P,T, data);
    @test  out.frac_M    + out.frac_S    + out.frac_F        ≈ 1.0
    @test  out.frac_M_wt + out.frac_S_wt + out.frac_F_wt     ≈ 1.0
    @test  sum(out.bulk_S)                                   ≈ 1.0
    @test  sum(out.bulk_S_wt)                                ≈ 1.0
  

    P           =   8.0
    T           =   1900.0
    out         =   point_wise_minimization(P,T, data);
    @test  out.frac_M    + out.frac_S    + out.frac_F        ≈ 1.0
    @test  out.frac_M_wt + out.frac_S_wt + out.frac_F_wt     ≈ 1.0
    @test  sum(out.bulk_M)                                   ≈ 1.0
    @test  sum(out.bulk_M_wt)                                ≈ 1.0
    
    Finalize_MAGEMin(data)
end



# previous way we defined this (left here for backwards compatibility)
db          = "ig"
gv, z_b, DB, splx_data  = init_MAGEMin(db);
sys_in      =   "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
test        =   0         #KLB1
gv          =   use_predefined_bulk_rock(gv, test, db);
gv.verbose=-1
P           =   8.0
T           =   800.0
out         =   point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in);
@test out.G_system ≈ -797.7397668057848
@test abs(out.s_cp[1] - 1208.466551730128) < 2.0
@test sort(out.ph) == sort(["spn", "cpx",  "opx", "ol"])
finalize_MAGEMin(gv,DB,z_b)

@testset "pointwise tests  " begin
    n       =   100;
    P       =   fill(8.0,n)
    T       =   fill(800.0,n)
    db      =   "ig" 
    data    =   Initialize_MAGEMin(db, verbose=false);
    out     =   multi_point_minimization(P, T, data, test=0);
    @test out[end].G_system ≈ -797.7397668057848
    @test sort(out[end].ph) == sort(["spn", "cpx",  "opx", "ol"])

    Finalize_MAGEMin(data)
end

@testset "specify bulk rock" begin
    

    data    = Initialize_MAGEMin("ig", verbose=false);
    
    # One bulk rock for all points
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

    @test abs(out.G_system + 916.8283889543869)/abs(916.8283889543869) < 2e-4


    # different bulk rock per point
    P       = [10.0, 10.0]
    T       = [1100.0, 1100.0]
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X1      = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X2      = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X       = [X1,X2]
    sys_in  = "wt"    
    out     = multi_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
    
    @test out[1].G_system ≈ -916.8283889543869 rtol=2e-4
    @test out[2].G_system ≈ -912.5920719174167 rtol=2e-4

    Finalize_MAGEMin(data)

    data    = Initialize_MAGEMin("um", verbose=-1, solver=0);
    # One bulk rock for all points
    P,T     = 10.0, 600.0
    Xoxides = ["SiO2", "Al2O3", "MgO", "FeO", "O", "H2O", "S"];
    X       = [20.044,0.6256,29.24,3.149,0.0,46.755,0.0]
    sys_in  = "mol"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)


end

@testset "PT adaptive refinement" begin
    data        = Initialize_MAGEMin("mp", verbose=-1, solver=0);

    init_sub    =  1
    ref_lvl     =  2
    Prange      = (1.0,10.0)
    Trange      = (400.0,800.0)
    Xoxides     = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","MnO","H2O"]
    X           = [70.999,12.805,0.771,3.978,6.342,2.7895,1.481,0.758,0.72933,0.075,30.0]
    sys_in      = "mol"    
    out         = AMR_minimization(init_sub, ref_lvl, Prange, Trange, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
    @test length(out) == 81
    @test sort(out[66].ph) == sort(["cd", "bi", "liq", "fsp", "sp", "ilm", "H2O"])
    Finalize_MAGEMin(data)
end


@testset "remove solution phase" begin

    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0);

    rm_list =   remove_phases(["liq","ilmm"],"mp")

    # One bulk rock for all points
    P,T     = 10.713125, 1177.34375
    Xoxides = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","MnO","H2O"]
    X       = [70.999,12.805,0.771,3.978,6.342,2.7895,1.481,0.758,0.72933,0.075,30.0]
    sys_in  = "mol"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in,rm_list=rm_list)
    @test sort(out.ph) == sort(["sp", "g", "fsp", "ilm", "q", "sill", "H2O"])
end

@testset "view array PT" begin

    data    = Initialize_MAGEMin("ig", verbose=false);

    # different bulk rock per point
    P       = [10.0, 10.0, 0]
    T       = [1100.0, 1100.0, 0]
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X1      = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X2      = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X       = [X1,X2]
    sys_in  = "wt"
    P_view = @view P[1:2]
    T_view = @view T[1:2]
    out     = multi_point_minimization(P_view, T_view, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

    Finalize_MAGEMin(data)
end

@testset "convert bulk rock" begin

    bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    bulk_in    = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    bulk_rock,ox  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"wt","ig");

    @test bulk_rock ≈ [46.12597764761598, 8.52489397284109, 11.805554573333653, 14.383528505497756, 6.471419573392541, 0.35839516987780934, 1.7264383468329216, 0.4871154401452383, 0.5876614012114892, 0.0, 9.529015369251512]

    bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"];
    bulk_in    = [69.64; 13.76; 1.77; 1.73; 4.32; 0.4; 2.61; 2.41; 0.80; 0.07; 0.0];
    bulk_rock,ox  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"wt","mp");

    @test bulk_rock ≈ [76.57038397179574, 8.914984523583415, 2.0849576977131403, 2.835783318610597, 4.30275071755529, 1.8302970975627948, 2.568605789798099, 0.6615823604771729, 0.16546809116073818, 0.06518643174302832, 0.0]
end


@testset "test Seismic velocities & modulus" begin
    # Call optimization routine for given P & T & bulk_rock
    data         = Initialize_MAGEMin("ig", verbose=false);
    test        = 0;
    data         = use_predefined_bulk_rock(data, test)
    P           = 8.0
    T           = 1200.0
    out         = point_wise_minimization(P,T, data)

    tol = 1.5e-2;
    @test abs(out.bulkMod - 95.36502708314566         )  < tol
    @test abs(out.shearMod - 29.908982525130096       )  < tol
    @test abs(out.Vs - 3.0562725380891456             )  < tol
    @test abs(out.Vp - 6.499047842613104              )  < tol
    @test abs(out.Vs_S -4.30859321743808              )  < tol
    @test abs(out.Vp_S - 7.392063221476717            )  < tol
    @test abs(out.bulkModulus_M - 27.239361903015595  )  < tol
    @test abs(out.bulkModulus_S - 95.74438231792864   )  < tol
    @test abs(out.shearModulus_S - 59.4633264822083   )  < tol

    Finalize_MAGEMin(data)
end

# Stores data of tests
mutable struct outP{ _T  } 
    P           ::  _T
    T           ::  _T 
    test        ::  Int64

    G           ::  _T
    ph          ::  Vector{String}
    ph_frac     ::  Vector{Float64}
end

print_error_msg(i,out) = println("ERROR for point $i with test=$(out.test); P=$(out.P); T=$(out.T); stable phases=$(out.ph), fractions=$(out.ph_frac)")


# Automatic testing of all points
function TestPoints(list, data::MAGEMin_Data)

    # Compute all points
    P = [ l.P for l in list]
    T = [ l.T for l in list]
    test = [ l.test for l in list]
    out_vec = multi_point_minimization(P, T, data, test = test[1]);

    # Check if the points this fit
    for (i,out) in enumerate(out_vec)
        VerifyPoint(out, list[i], i)
    end
    return nothing
end

# This checks whether a single point agrees with precomputed values & prints a message if not
function VerifyPoint(out, list, i)

     # We need to sort the phases (sometimes they are ordered differently)
     ind_sol = sortperm(list.ph)
     ind_out = sortperm(out.ph)
     
     result1 = @test out.G_system  ≈ list.G     rtol=1e-3
     result2 = @test out.ph[ind_out]        == list.ph[ind_sol]
     result3 = @test sort(out.ph_frac) ≈ sort(list.ph_frac) atol=5e-2       # ok, this is really large (needs fixing for test6!)
     
     # print more info about the point if one of the tests above fails
     if isa(result1,Test.Fail) || isa(result2,Test.Fail) || isa(result3,Test.Fail)
         print_error_msg(i,list)
     end
     
     return nothing
end

# load reference for built-in tests
println("Testing points from the reference diagrams:")
@testset verbose = true "Total tests" begin

    # Igneous database
    println("  Starting KLB-1 peridotite tests")
    db  = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
    data = Initialize_MAGEMin(db, verbose=false);
   
    gv.verbose=-1;
    @testset "IG-DB - KLB-1 peridotite" begin
        include("test_diagram_test0.jl")
        TestPoints(list, data)
    end
    Finalize_MAGEMin(data)

    println("  Starting RE-46 icelandic basalt tests")
    db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
    data = Initialize_MAGEMin(db, verbose=false);
    gv.verbose=-1;
    @testset "IG-DB - RE-46 icelandic basalt" begin
        include("test_diagram_test1.jl")
        TestPoints(list, data)
    end
   
    println("  Starting Wet MORB tests")
    db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
    data = Initialize_MAGEMin(db, verbose=false);
    @testset "IG-DB - Wet MORB" begin
        include("test_diagram_test6.jl")
        TestPoints(list, data)
    end
    Finalize_MAGEMin(data)

    # Metapelite database
    println("  Starting WM Pelite tests")
    db  = "mp"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
    data = Initialize_MAGEMin(db, verbose=false);
   
    gv.verbose=-1;
    @testset "MP-DB - WM Pelite" begin
        include("test_diagram_test0_mp.jl")
        TestPoints(list, data)
    end
    Finalize_MAGEMin(data)

    println("  Starting Gt-Migmatite tests")
    db  = "mp"  # database: mp, metapelite (White et al 2014b)
    data = Initialize_MAGEMin(db, verbose=false);
   
    gv.verbose=-1;
    @testset "MP-DB - Gt-Migmatite" begin
        include("test_diagram_test4_mp.jl")
        TestPoints(list, data)
    end
    Finalize_MAGEMin(data)

    # Metabasite database
    println("  Starting SQA Amphibole tests")
    db  = "mb"  # database: ig, igneous (Holland et al., 2018)
    data = Initialize_MAGEMin(db, verbose=false, mbCpx = 1);
   
    gv.verbose=-1;
    @testset "MB-DB - SQA Amphibole" begin
        include("test_diagram_test0_mb.jl")
        TestPoints(list, data)
    end
    Finalize_MAGEMin(data)

    # Igneous alkaline dry database
    println("  Starting Syenite tests")
    db  = "igad"  # database: igad, Igneous alkaline dry database (Weller et al., 2024)
    data = Initialize_MAGEMin(db, verbose=false);
   
    gv.verbose=-1;
    @testset "IGAD-DB - Syenite" begin
        include("test_diagram_test0_igad.jl")
        TestPoints(list, data)
    end
    Finalize_MAGEMin(data)


end

# a few tests that gave problems in the past
println("Testing problematic points:")
@testset verbose = true "Problematic points" begin
    include("test_problematic_points.jl")
end