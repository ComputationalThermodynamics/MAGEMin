#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
# this tests the julia interface to MAGEMin
using Test
using MAGEMin_C

function norm(vec :: Vector{Float64})
    return sqrt(sum(vec.^2))
end

#= A convenient testy test

using MAGEMin_C
data    = Initialize_MAGEMin("mp", verbose=-1);
P,T     = 6.0, 710.0
Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
sys_in  = "wt"
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
Finalize_MAGEMin(data)

using MAGEMin_C
data        =   Initialize_MAGEMin("ig", verbose=-1);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);
P           =   8.0
T           =   1800.0
out         =   point_wise_minimization(P,T, data);
Finalize_MAGEMin(data)

=#

data        =   Initialize_MAGEMin("sb21", verbose=-1);
test        =   1         #KLB1
data        =   use_predefined_bulk_rock(data, test);
P           =   80.0
T           =   800.0
out         =   point_wise_minimization(P,T, data);
@test sort(out.ph) == sort(["gtmj", "hpcpx", "ol" ,"cpx"])
Finalize_MAGEMin(data)

# generic test for thermocalc database
data        =   Initialize_MAGEMin("ig", verbose=-1);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);
P           =   8.0
T           =   800.0
out         =   point_wise_minimization(P,T, data);
Finalize_MAGEMin(data)

@test out.G_system ≈ -797.7873859283119
@test sort(out.ph) == sort(["spl", "cpx",  "opx", "ol"])
@test abs(out.s_cp[1] - 1208.466551730128) < 2.0

@testset "test external routines" begin
    ox              = ["SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "H2O"]
    mol_percents    = [62.38, 0.41, 11.79, 0.03, 0.02, 4.80, 9.73, 3.41, 0.59, 0.05, 6.80]
    T_C             = 1000.0
    viscosity       = compute_melt_viscosity_G08(ox, mol_percents, T_C)
    @test (viscosity) ≈ 4751.168588718496
end

@testset "test light output calculation" begin
    # Without a buffer at 1100.0 C
    data    = Initialize_MAGEMin("ig", verbose=-1);
    P,T     = 10.0, 600.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [78.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_hT  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, light=true);

    @test sort(out_hT.ph_id_db) == sort([1, 1, 3, 4, 6, 10, 12])
    @test out_hT.ph_type == [1, 1, 1, 1, 1, 0, 0]

    Finalize_MAGEMin(data)
end

# Tests from L. Candioti - ETH - Oct 2024
@testset "test mass conservation" begin

    # Without a buffer at 1100.0 C
    data    = Initialize_MAGEMin("ig", verbose=-1);
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_hT  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    Δρ_hT   = abs( out_hT.rho - ((out_hT.frac_M_wt * out_hT.rho_M + out_hT.frac_S_wt * out_hT.rho_S )) )
    @test Δρ_hT < 1e-10
    Finalize_MAGEMin(data)

    # Without a buffer at 800.0 C
    data    = Initialize_MAGEMin("ig", verbose=-1);
    P,T     = 10.0, 800.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_lT  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    Δρ_lT = abs( out_lT.rho - ((out_lT.frac_M_wt * out_lT.rho_M + out_lT.frac_S_wt * out_lT.rho_S )) )
    @test Δρ_lT < 1e-10
    Finalize_MAGEMin(data)

    # With a buffer at 1100.0 C
    data    = Initialize_MAGEMin("ig", buffer = "nno", verbose=-1);
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_BhT = single_point_minimization(P, T, data, X=X, B=0.0, Xoxides=Xoxides, sys_in=sys_in);
    Δρ_BhT  = abs( out_BhT.rho - ((out_BhT.frac_M_wt * out_BhT.rho_M + out_BhT.frac_S_wt * out_BhT.rho_S )) )
    @test Δρ_BhT < 1e-10
    Finalize_MAGEMin(data)
    
    # With a buffer at 800.0 C
    data    = Initialize_MAGEMin("ig", verbose=-1);
    P,T     = 10.0, 800.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_BlT = single_point_minimization(P, T, data, X=X, B=0.0, Xoxides=Xoxides, sys_in=sys_in);
    Δρ_BlT  = abs( out_BlT.rho - ((out_BlT.frac_M_wt * out_BlT.rho_M + out_BlT.frac_S_wt * out_BlT.rho_S )) )
    @test Δρ_BlT < 1e-10
    Finalize_MAGEMin(data)

end

@testset "test activity buffers" begin
    data        =   Initialize_MAGEMin("mp", verbose=-1, buffer="aH2O");
    test        =   0        
    data        =   use_predefined_bulk_rock(data, test);
    P           =   8.0
    T           =   400.0
    out         =   point_wise_minimization(P,T, data, buffer_n=0.6);
    @test sort(out.ph) == sort(["chl", "sp", "mu", "mu", "fsp", "ep", "q", "ru", "aH2O"])
    Finalize_MAGEMin(data)

    data        =   Initialize_MAGEMin("mp", verbose=-1, buffer="aTiO2");
    test        =   0        
    data        =   use_predefined_bulk_rock(data, test);
    P           =   8.0
    T           =   400.0
    out         =   point_wise_minimization(P,T, data, buffer_n=0.6);
    @test sort(out.ph) == sort(["H2O", "aTiO2", "chl", "ep", "fsp", "ilm", "mu", "mu", "q"])
    Finalize_MAGEMin(data)

    data        =   Initialize_MAGEMin("ig", verbose=-1, buffer="aTiO2");
    test        =   0        
    data        =   use_predefined_bulk_rock(data, test);
    P           =   8.0
    T           =   1200.0
    out         =   point_wise_minimization(P,T, data, buffer_n=0.1);
    @test sort(out.ph) == sort(["aTiO2", "cpx", "fsp", "liq", "ol", "opx", "spl"])

    Finalize_MAGEMin(data)

    data    = Initialize_MAGEMin("ig", verbose=-1, buffer="qfm");
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    @test sort(out.ph) == sort(["cpx", "liq", "opx", "qfm"])

    Finalize_MAGEMin(data)

    data    = Initialize_MAGEMin("ig", verbose=-1, buffer="iw");
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    @test sort(out.ph) == sort(["opx", "liq", "cpx", "iw"])

    Finalize_MAGEMin(data)

    data        =   Initialize_MAGEMin("mp", verbose=-1, buffer="iw");
    test        =   0        
    data        =   use_predefined_bulk_rock(data, test);
    P           =   8.0
    T           =   400.0
    out         =   point_wise_minimization(P,T, data, buffer_n=-5.0);
    @test sort(out.ph) == sort(["chl", "fsp", "mu", "mu", "q", "ru", "sph", "H2O", "iw"])
    Finalize_MAGEMin(data)
end

@testset "test sum frac_vol" begin
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0);
    P,T     = 10.713125, 1177.34375
    Xoxides = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","MnO","H2O"]
    X       = [70.999,12.805,0.771,3.978,6.342,2.7895,1.481,0.758,0.72933,0.075,30.0]
    sys_in  = "mol"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    @test sum(out.frac_M_vol + out.frac_F_vol + out.frac_S_vol) ≈ 1.0
    Finalize_MAGEMin(data)

    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0);
    P,T     = 5.713125, 477.34375
    Xoxides = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","MnO","H2O"]
    X       = [70.999,12.805,0.771,3.978,6.342,2.7895,1.481,0.758,0.72933,0.075,30.0]
    sys_in  = "mol"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    Finalize_MAGEMin(data)
    @test sum(out.frac_M_vol + out.frac_F_vol + out.frac_S_vol) ≈ 1.0
end

@testset "test zr saturation" begin
    data    = Initialize_MAGEMin("mp", verbose=-1);

    P,T     = 6.0, 930.0
    Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
    X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
    sys_in  = "wt"
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);

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
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);

    bulk_melt = [61.26, 0, 13.12, 1.33, 0, 0.45, 2.51, 2.26, 2.15, 15.00, 0]
    # convert bulk_melt from wt to mol of oxides
    bulk_melt = convertBulk4MAGEMin(bulk_melt, Xoxides, "wt", "mp")[1]
    out.bulk_M .= bulk_melt

    zr_sat = MAGEMin_C.zirconium_saturation(out, model="CB")

    @test zr_sat ≈ 65.83158859091596 rtol=1e-5
end

@testset "test normalization" begin
    data        =   Initialize_MAGEMin("ig", verbose=-1);
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
@test out.G_system ≈ -797.7873859283119
@test abs(out.s_cp[1] - 1208.466551730128) < 2.0
@test sort(out.ph) == sort(["spl", "cpx",  "opx", "ol"])
finalize_MAGEMin(gv,DB,z_b)

@testset "pointwise tests  " begin
    n       =   100;
    P       =   fill(8.0,n)
    T       =   fill(800.0,n)
    db      =   "ig" 
    data    =   Initialize_MAGEMin(db, verbose=-1);
    out     =   multi_point_minimization(P, T, data, test=0);
    @test out[end].G_system ≈ -797.7873859283119
    @test sort(out[end].ph) == sort(["spl", "cpx",  "opx", "ol"])

    Finalize_MAGEMin(data)
end

@testset "specify bulk rock" begin
    data    = Initialize_MAGEMin("ig", verbose=-1);
    
    # One bulk rock for all points
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);

    @test abs(out.G_system + 916.8283889543869)/abs(916.8283889543869) < 2e-4

    # different bulk rock per point
    P       = [10.0, 10.0]
    T       = [1100.0, 1100.0]
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X1      = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X2      = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X       = [X1,X2]
    sys_in  = "wt"    
    out     = multi_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    
    @test out[1].G_system ≈ -916.8283889543869 rtol=2e-4
    @test out[2].G_system ≈ -912.5920719174167 rtol=2e-4

    Finalize_MAGEMin(data)

    data    = Initialize_MAGEMin("um", verbose=-1, solver=0);
    # One bulk rock for all points
    P,T     = 10.0, 600.0
    Xoxides = ["SiO2", "Al2O3", "MgO", "FeO", "O", "H2O", "S"];
    X       = [20.044,0.6256,29.24,3.149,0.0,46.755,0.0]
    sys_in  = "mol"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
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
    out         = AMR_minimization(init_sub, ref_lvl, Prange, Trange, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    @test length(out) == 81
    @test sort(out[66].ph) == sort(["cd", "bi", "liq", "fsp", "sp", "ilm", "H2O"])
    Finalize_MAGEMin(data)
end

@testset "Trace-element partitioning model" begin
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0);
    P,T     = 6.0, 699.0
    Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
    X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.2];
    sys_in  = "wt"
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, name_solvus=true);
    Finalize_MAGEMin(data)


    # create database on the fly
    el      = ["Li","Zr"]
    ph      = ["q","afs","pl","bi","opx","cd","mu","amp","fl","cpx","g","zrn"]
    KDs     = ["0.17" "0.01";"0.14 * T_C/1000.0 + [:bi].compVariables[1]" "0.01";"0.33 + 0.01*P_kbar" "0.01";"1.67 * P_kbar / 10.0 + T_C/1000.0" "0.01";"0.2" "0.01";"125" "0.01";"0.82" "0.01";"0.2" "0.01";"0.65" "0.01";"0.26" "0.01";"0.01" "0.01";"0.01" "0.0"] 
    C0      = [100.0,400.0] #starting concentration of elements in ppm (ug/g)
    dtb     = "mp"

    KDs_database = create_custom_KDs_database(el, ph, KDs)

    out_TE = TE_prediction(out, C0, KDs_database, dtb; ZrSat_model = "CB");

    @test out_TE.Cliq[1] ≈ 154.51891525771387 rtol=1e-3
    @test out_TE.Cliq[2] ≈ 3107.727391290317 rtol=1e-3
end

@testset "remove solution phase" begin

    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0);
    rm_list =   remove_phases(["liq","ilm"],"mp")
    P,T     = 10.713125, 1177.34375
    Xoxides = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","MnO","H2O"]
    X       = [70.999,12.805,0.771,3.978,6.342,2.7895,1.481,0.758,0.72933,0.075,30.0]
    sys_in  = "mol"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in,rm_list=rm_list);
    @test sort(out.ph) == sort(["fsp", "g", "ilmm", "sp", "q", "sill", "H2O"])
    Finalize_MAGEMin(data)

    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0);
    rm_list =   remove_phases(["liq","ilmm","sill"],"mp")
    P,T     = 10.713125, 1177.34375
    Xoxides = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","MnO","H2O"]
    X       = [70.999,12.805,0.771,3.978,6.342,2.7895,1.481,0.758,0.72933,0.075,30.0]
    sys_in  = "mol"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in,rm_list=rm_list);
    @test sort(out.ph) == sort(["fsp", "cd", "sa", "ilm", "sp", "q", "H2O"])
    Finalize_MAGEMin(data)
end

@testset "view array PT" begin

    data    = Initialize_MAGEMin("ig", verbose=-1);

    # different bulk rock per point
    P       = [10.0, 10.0, 0]
    T       = [1100.0, 1100.0, 0]
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X1      = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X2      = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X       = [X1, X2] # only use first two points
    sys_in  = "wt"
    P_view  = @view P[1:2]
    T_view  = @view T[1:2]
    out     = multi_point_minimization(P_view, T_view, data, X=X, Xoxides=Xoxides, sys_in=sys_in);

    # test with a view of the bulk rock
    index_shufle   = [2,1,3,4,5,6,7,8,9,10,11]
    Xoxides_shufle = ["Al2O3"; "SiO2"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"]
    X1_view        = @view X1[index_shufle]

    # just run it to be sure it is not erroring
    out     = single_point_minimization(P[3], T[3], data, X=X1_view, Xoxides=Xoxides_shufle, sys_in=sys_in);
    mol2wt(X1_view, Xoxides_shufle) # convert to mol
    wt2mol(X1_view, Xoxides_shufle) # convert to mol

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
    data         = Initialize_MAGEMin("ig", verbose=-1);
    test        = 0;
    data         = use_predefined_bulk_rock(data, test)
    P           = 8.0
    T           = 1200.0
    out         = point_wise_minimization(P,T, data)

    tol = 1.5e-2;
    @test abs(out.bulkMod - 94.62309357990975          )  < tol
    @test abs(out.shearMod - 29.843843046045578        )  < tol
    @test abs(out.Vs - 3.0500442437065094              )  < tol
    @test abs(out.Vp - 6.472952899434848               )  < tol
    @test abs(out.Vs_S -4.303123606906489              )  < tol
    @test abs(out.Vp_S - 7.3759048706307055            )  < tol
    @test abs(out.bulkModulus_M - 27.774175695339732   )  < tol
    @test abs(out.bulkModulus_S - 95.39738730456645    )  < tol
    @test abs(out.shearModulus_S - 59.44716946888283   )  < tol

    Finalize_MAGEMin(data)
end


@testset "test Mantle HP13" begin

    data        =   Initialize_MAGEMin("mtl", verbose=-1);
    test        =   0         #KLB1
    data        =   use_predefined_bulk_rock(data, test);

    # Call optimization routine for given P & T & bulk_rock
    P           =   180.0
    T           =   1400.0
    out         =   point_wise_minimization(P,T, data);

    @test sort(out.ph) == sort(["g", "ring", "wad"])
    Finalize_MAGEMin(data)
end

@testset "test ume" begin

    data        =   Initialize_MAGEMin("ume", verbose=-1);
    test        =   0
    data        =   use_predefined_bulk_rock(data, test);
    P           =   20.0
    T           =   400.0
    out         =   point_wise_minimization(P,T, data);
    Finalize_MAGEMin(data)
    @test sort(out.ph) == sort(["amp", "atg", "chl", "fl", "hem", "pyr", "spi", "ta"])
end


# test from Philip Hartmeier
@testset "test apfu" begin

    data        =   Initialize_MAGEMin("mp", verbose=-1);
    T           = 580.0
    P           = 4.5
    X           = [64.13, 0.91, 19.63, 6.85, 0.08, 2.41, 0.65, 1.38, 3.95, 40.0]
    Xoxides     = ["SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "H2O"]
    sys_in      = "wt"
    out         =   single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    Finalize_MAGEMin(data)

    id_bi = findfirst( out.ph .== "bi" )
    @test sum(abs.(out.SS_vec[id_bi].Comp_apfu .- [2.7139545235947877, 1.572090952810424, 0.0, 1.1115442718289477, 1.5026608485480333, 1.0, 0.0, 0.08905581689654532, 12.0, 0.010693586321261687, 1.8218883662069094])) .< 1e-3
end


@testset "Text initial guess" begin

    MAGEMin_data    = Initialize_MAGEMin("ig", verbose=false, solver=0);

    Xoxides         = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"];
    X1              = [70.999,  12.805, 0.771,  3.978,  6.342,  2.7895, 1.481,  0.758,  0.72933,    0.1,    3.0];
    X2              = [70.999,  12.805, 0.771,  3.978,  6.342,  2.7895, 1.481,  0.758,  0.72933,    0.1,    9.0];
    X3              = [70.999,  12.805, 0.771,  3.978,  6.342,  2.7895, 1.481,  0.758,  0.72933,    0.1,    15.0];
    X4              = [70.999,  12.805, 0.771,  3.978,  6.342,  2.7895, 1.481,  0.758,  0.72933,    0.1,    21.0];
    sys_in          = "mol";

    P, T            = 19.0, 1350.0;

    Pvec,Tvec       = [19.0,19.0,19.5,19.5], [1325.0,1350.0,1325.0,1350.0]

    Xvec            = [X1,X2,X3,X4] # here the composition can also be slightly varied. how much I am not quite sure yet

    Out_XY          = Vector{MAGEMin_C.gmin_struct}(undef,length(Pvec))
    Out_XY_ig       = Vector{MAGEMin_C.gmin_struct}(undef,length(Pvec))
    Out_XY          = multi_point_minimization( Pvec, Tvec, MAGEMin_data;
                                                X=Xvec, Xoxides=Xoxides, sys_in=sys_in, 
                                                name_solvus=true); 

                            
    tmp             = [Out_XY[i].mSS_vec for i=1:length(Pvec)]
    Gig             = vcat(tmp...)                  

    Out_ig          = single_point_minimization(    19.25, 1337.5,
                                                    MAGEMin_data;
                                                    X           = sum(Xvec)./4.0,
                                                    Xoxides     = Xoxides,
                                                    sys_in      = sys_in, 
                                                    name_solvus = true,
                                                    iguess      = true,
                                                    G           = [Gig]);
                                            
    Finalize_MAGEMin(MAGEMin_data)

    @test sort(Out_ig.ph) == sort(["liq", "spl"])
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


@testset "Metastability function" begin
    data    = Initialize_MAGEMin("mp", verbose=-1; solver=0);
    P,T     = 6.0, 630.0
    Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
    X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
    sys_in  = "wt"

    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
    Pmeta, Tmeta       = 6.0, 500.0
    out2    = point_wise_metastability(out, Pmeta, Tmeta, data)

    Finalize_MAGEMin(data)

    @test abs(out.G_system + 806.7071168433587) < 1e-6
    @test abs(out2.G_system + 791.460287) < 1e-6
end

@testset verbose = true "Test Ws override" begin

    #= First we create a structure to store the data in memory =#
    dtb     = 0             # metapelite
    ss_id   = 3             # biotite
    n_Ws    = 21            # number of Margules parameters
    Ws      = [0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0]
    new_Ws      =  Vector{MAGEMin_C.W_data{Float64,Int64}}(undef, 1)
    new_Ws[1]   = MAGEMin_C.W_data(dtb, ss_id, n_Ws, Ws)   

    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0);
    P,T     = 4.0,650.0
    Xoxides = ["SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","MnO","H2O"]
    X       = [70.999,12.805,0.771,3.978,6.342,2.7895,1.481,0.758,0.72933,0.075,30.0]
    sys_in  = "mol"    
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in ,W=new_Ws)
    # out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
    Finalize_MAGEMin(data)
    @test norm(out.ph_frac) - 0.45682499466457954 < 0.01
end

#=
# The following part is not really a test yet, but more an example of how to use initial guesses

using MAGEMin_C

MAGEMin_data    = Initialize_MAGEMin("mp", verbose=-1);

Xoxides         = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
X1              = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
sys_in          = "wt"


Pvec,Tvec       = [6.0,6.1,6.0,6.1], [710.0,710.0,720.0,720.0]
Xvec            = [X1,X1,X1,X1] # here the composition can also be slightly varied. how much I am not quite sure yet

Out_XY          = Vector{out_struct}(undef,length(Pvec))
Out_XY          = multi_point_minimization( Pvec, Tvec, MAGEMin_data;
                                            X=Xvec, Xoxides=Xoxides, sys_in=sys_in, 
                                            name_solvus=true); 

# retrieve theinitial guesses
# The way it works is by retrieving the mSS_vec structure from the 4 previous minimizations and concatenating them into a single vector
tmp             = [Out_XY[i].mSS_vec for i=1:length(Pvec)]
Gig             = vcat(tmp...)                  

# note that below we use slightly different P,T conditions and that Gig is passed as an initial guess within square brackets
# A similar approach can be used in the case of multi_point_minimization but then a vector of P,T, iguess and G must be passed
# Note that G is a vector of vectors in that case => Gig = Vector{Vector{LibMAGEMin.mSS_data}}(undef,np);
Out_ig          = single_point_minimization(    6.05, 715.0, MAGEMin_data;
                                                X=X1, Xoxides=Xoxides, sys_in=sys_in, 
                                                name_solvus=true,
                                                iguess=true,G=[Gig]); 


Finalize_MAGEMin(MAGEMin_data)


=#