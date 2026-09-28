#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
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

@testset verbose=true "single point minimization" begin
    data        =   Initialize_MAGEMin("sb24", verbose=-1);
    test        =   1         #Pyrolite
    data        =   use_predefined_bulk_rock(data, test);
    P           =   60.0
    T           =   1000.0
    out         =   point_wise_minimization(P,T, data);
    @test sort(out.ph) == sort(["opx", "cpx", "gtmj", "ol"])
    Finalize_MAGEMin(data)

    data        =   Initialize_MAGEMin("ig", verbose=-1);
    test        =   0         #KLB1
    data        =   use_predefined_bulk_rock(data, test);
    P           =   8.0
    T           =   800.0
    out         =   point_wise_minimization(P,T, data);
    Finalize_MAGEMin(data)

    @test out.G_system ≈ -797.7873865220898
    @test sort(out.ph) == sort(["spl", "cpx",  "opx", "ol"])
    @test abs(out.s_cp[1] - 1208.466551730128) < 2.0
end

@testset verbose=true "test external routines" begin
    ox              = ["SiO2", "TiO2", "Al2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "P2O5", "H2O"]
    mol_percents    = [62.38, 0.41, 11.79, 0.03, 0.02, 4.80, 9.73, 3.41, 0.59, 0.05, 6.80]
    T_C             = 1000.0
    viscosity       = compute_melt_viscosity_G08(ox, mol_percents, T_C)
    @test (viscosity) ≈ 4751.168588718496
end

@testset verbose=true "test light output calculation" begin
    # Without a buffer at 1100.0 C
    data    = Initialize_MAGEMin("ig", verbose=-1);
    P,T     = 10.0, 600.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [78.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_hT  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, light=true);

    @test sort(out_hT.ph_name) == sort(["amp", "bi", "cpx", "ep", "fl", "q", "sph"])

    Finalize_MAGEMin(data)
end


@testset verbose=true "test light output calculation with initial guess" begin
    # Without a buffer at 1100.0 C
    data    = Initialize_MAGEMin("ig", verbose=-1);
    P,T     = 10.0, 600.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [78.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_hT  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, light=true, light_ig = true);

    @test sort(out_hT.ph_name) == sort(["amp", "bi", "cpx", "ep", "fl", "q", "sph"])

    Finalize_MAGEMin(data)
end


@testset verbose=true "test global TC database" begin

    using MAGEMin_C
    ss_list = ["liq_W14", "fsp_H22", "bi_W14", "g_W14", "ep_H11", "ma_W14", "mu_W14", "opx_W14", "sa_W14", "cd_W14", "st_W14", "chl_W14", "ctd_W14", "sp_W02", "ilm_W00", "DEW_S24"]
    pp_list = ["q", "crst", "trd", "coe", "stv", "law", "ky", "sill", "and", "ru", "sph","prl"]

    data    = Initialize_MAGEMin("all", verbose=false, solver=0);
    P, T    = 10.0, 400.0;
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "Cr2O3"; "H2O"; "CO2"; "S"];
    X       = [0.62212, 0.1122, 0.0, 0.03486, 0.05557, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.17525, 0.0, 0.0];
    sys_in  = "mol";
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, ss_list=ss_list, pp_list=pp_list)
    @test sort(out.ph) == ["DEW_S24", "chl_W14", "ctd_W14", "prl", "q"]

    ss_list = ["liq_W14", "fsp_H22", "bi_W14", "g_W14", "ep_H11", "ma_W14", "mu_W14", "opx_W14", "sa_W14", "cd_W14", "st_W14",  "ctd_W14", "sp_W02", "ilm_W00", "DEW_S24"]
    pp_list = ["q", "crst", "trd", "coe", "stv", "law", "ky", "sill", "and", "ru", "sph"]
    P, T    = 10.0, 400.0;
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "Cr2O3"; "H2O"; "CO2"; "S"];
    X       = [0.62212, 0.1122, 0.0, 0.03486, 0.05557, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.17525, 0.0, 0.0];
    sys_in  = "mol";
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, ss_list=ss_list, pp_list=pp_list)
    @test sort(out.ph) == ["DEW_S24", "ctd_W14", "ky", "q"]

    Finalize_MAGEMin(data)
end

# Tests from L. Candioti - ETH - Oct 2024
@testset verbose=true "test mass conservation" begin

    # Without a buffer at 1100.0 C
    data    = Initialize_MAGEMin("ig", verbose=-1);
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_hT  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    Δρ_hT   = abs( out_hT.rho - (1.0 / (out_hT.frac_M_wt / out_hT.rho_M + out_hT.frac_S_wt / out_hT.rho_S)) )
    @test Δρ_hT < 1e-10
    Finalize_MAGEMin(data)

    # Without a buffer at 800.0 C
    data    = Initialize_MAGEMin("ig", verbose=-1);
    P,T     = 10.0, 800.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_lT  = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);
    Δρ_lT = abs( out_lT.rho - (1.0 / (out_lT.frac_M_wt / out_lT.rho_M + out_lT.frac_S_wt / out_lT.rho_S)) )
    @test Δρ_lT < 1e-10
    Finalize_MAGEMin(data)

    # With a buffer at 1100.0 C
    data    = Initialize_MAGEMin("ig", buffer = "nno", verbose=-1);
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_BhT = single_point_minimization(P, T, data, X=X, B=0.0, Xoxides=Xoxides, sys_in=sys_in);
    Δρ_BhT  = abs( out_BhT.rho - (1.0 / (out_BhT.frac_M_wt / out_BhT.rho_M + out_BhT.frac_S_wt / out_BhT.rho_S)) )
    @test Δρ_BhT < 1e-10
    @test out_BhT.system_oxygen ≈ 0.004867165479051421    rtol=1e-3
    Finalize_MAGEMin(data)
    
    # With a buffer at 800.0 C
    data    = Initialize_MAGEMin("ig", verbose=-1);
    P,T     = 10.0, 800.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out_BlT = single_point_minimization(P, T, data, X=X, B=0.0, Xoxides=Xoxides, sys_in=sys_in);
    Δρ_BlT  = abs( out_BlT.rho - (1.0 / (out_BlT.frac_M_wt / out_BlT.rho_M + out_BlT.frac_S_wt / out_BlT.rho_S)) )
    @test Δρ_BlT < 1e-10
    Finalize_MAGEMin(data)

end


@testset verbose=true "test PH_vec keys" begin
    using MAGEMin_C
    data    = Initialize_MAGEMin("mp", verbose=-1);
    P,T     = 6.0, 710.0
    Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
    X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
    sys_in  = "wt"
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, scp = 1, dT=1.5)
    @test haskey(out.PH_vec, :ru)   == false
    @test haskey(out.PH_vec, :q)    == true
    @test haskey(out.PH_vec, :opx)  == true
    
    Finalize_MAGEMin(data)
end


@testset verbose=true "test activity buffers" begin
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
    out         =   point_wise_minimization(P,T, data, buffer_n=0.01);
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

@testset verbose=true "test native mu-mu chemical potential fixing" begin
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "H2O"];
    X       = [38.494; 1.977; 2.907; 50.264; 5.435; 0.017; 0.204; 0.077; 0.086; 0.469; 0.077];
    P, T    = 10.0, 1100.0
    sys_in  = "mol"

    data_base   = Initialize_MAGEMin("ig", verbose=-1);
    out_base    = single_point_minimization(P, T, data_base; X=X, Xoxides=Xoxides, sys_in=sys_in);
    Gamma_MgO   = out_base.Gamma[4]    # 1-based; MgO is Xoxides[4]
    Gamma_K2O   = out_base.Gamma[6]    # K2O is Xoxides[6]
    Finalize_MAGEMin(data_base)

    data = Initialize_MAGEMin("ig", verbose=-1, mu_fix_idx=["MgO"]);
    out  = single_point_minimization(P, T, data; X=X, Xoxides=Xoxides, sys_in=sys_in, mu_fix_val=[Gamma_MgO]);
    @test out.Gamma[4] ≈ Gamma_MgO atol=1e-3
    Finalize_MAGEMin(data)

    data     = Initialize_MAGEMin("ig", verbose=-1, mu_fix_idx=["MgO"]);
    X_oversat = copy(X); X_oversat[4] = 150.0
    out      = single_point_minimization(P, T, data; X=X_oversat, Xoxides=Xoxides, sys_in=sys_in, mu_fix_val=[-900.0]);
    @test out.Gamma[4] ≈ -900.0 atol=1e-6
    @test "mMgO" in out.ph
    @test isfinite(out.rho)
    @test isfinite(out.s_cp[1])
    Finalize_MAGEMin(data)

    # 3 oxides simultaneously, by name, all oversaturated: exercises both the
    # arbitrary-N-oxides support and oxide-name resolution together
    data = Initialize_MAGEMin("ig", verbose=-1, mu_fix_idx=["MgO", "K2O", "TiO2"]);
    X3   = [38.494; 1.977; 2.907; 70.0; 5.435; 0.1; 0.204; 0.2; 0.086; 0.469; 0.077];
    out  = single_point_minimization(P, T, data; X=X3, Xoxides=Xoxides, sys_in=sys_in,
                                      mu_fix_val=[Gamma_MgO, Gamma_K2O, -1077.756911]);
    @test out.Gamma[4] ≈ Gamma_MgO       atol=1e-3
    @test out.Gamma[6] ≈ Gamma_K2O       atol=1e-3
    @test out.Gamma[8] ≈ -1077.756911    atol=1e-3
    Finalize_MAGEMin(data)

    @test_throws ErrorException Initialize_MAGEMin("ig", verbose=-1, mu_fix_idx=["NotAnOxide"])
    data     = Initialize_MAGEMin("ig", verbose=-1, mu_fix_idx=["MgO"]);
    n        = 6
    Pvec     = fill(P, n)
    Tvec     = fill(T, n)
    targets  = collect(range(-950.0, -700.0, length=n))
    mu_grid  = [[t] for t in targets]
    out_grid = multi_point_minimization(Pvec, Tvec, data; X=X_oversat, Xoxides=Xoxides, sys_in=sys_in,
                                         mu_fix_val=mu_grid, progressbar=false);
    for i in 1:n
        @test out_grid[i].Gamma[4] ≈ targets[i] atol=1e-6
    end
    Finalize_MAGEMin(data)
end

@testset "test seismic corrections" begin

    using MAGEMin_C
    data    = Initialize_MAGEMin("mp", verbose=-1);
    P,T     = 6.0, 710.0
    Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
    X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0];
    sys_in  = "wt"
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, seismic_cor=true, aspect_ratio=0.3, seismic_water=0);

    @test out.Vp_cor ≈ 5.342675762785988    rtol=1e-4
    @test out.Vs_cor ≈ 3.415349818810448    rtol=1e-4

    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, seismic_cor=true, aspect_ratio=0.6, seismic_water=0);

    @test out.Vp_cor ≈ 5.35301103758362     rtol=1e-4
    @test out.Vs_cor ≈ 3.423002024620632    rtol=1e-4

    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, seismic_cor=true, aspect_ratio=0.1, seismic_water=0);

    @test out.Vp_cor ≈ 5.31879891895499     rtol=1e-4
    @test out.Vs_cor ≈ 3.39393959747249     rtol=1e-4

    Finalize_MAGEMin(data)
end


@testset verbose=true "test sum frac_vol" begin
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

@testset verbose=true "test zr saturation" begin
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
    zr_sat_B    = MAGEMin_C.zirconium_saturation(out, model="B")
    zr_sat_WH   = MAGEMin_C.zirconium_saturation(out, model="WH")

    @test zr_sat_B  ≈ 1403.8755429428836 rtol=1e-5
    @test zr_sat_WH ≈ 1059.5976323423222 rtol=1e-5

    # test crisp and berry 2022, use compo from their example in the calculator from their paper
    P,T     = 20.0, 750.0
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in);

    bulk_melt    = [61.26, 0, 13.12, 1.33, 0, 0.45, 2.51, 2.26, 2.15, 15.00, 0]
    bulk_melt   = convertBulk4MAGEMin(bulk_melt, Xoxides, "wt", "mp")[1]
    out.bulk_M .= bulk_melt

    zr_sat  = MAGEMin_C.zirconium_saturation(out, model="CB")

    @test zr_sat ≈ 65.83158859091596 rtol=1e-5
end

@testset verbose=true "test normalization" begin
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
@test out.G_system ≈ -797.7873865220898
@test abs(out.s_cp[1] - 1208.466551730128) < 2.0
@test sort(out.ph) == sort(["spl", "cpx",  "opx", "ol"])
finalize_MAGEMin(gv,DB,z_b,splx_data)

@testset verbose=true "pointwise tests  " begin
    n       =   100;
    P       =   fill(8.0,n)
    T       =   fill(800.0,n)
    db      =   "ig"
    data    =   Initialize_MAGEMin(db, verbose=-1);
    out     =   multi_point_minimization(P, T, data, test=0);
    @test out[end].G_system ≈ -797.7873865220898
    @test sort(out[end].ph) == sort(["spl", "cpx",  "opx", "ol"])

    Finalize_MAGEMin(data)
end

@testset verbose=true "calibration mode" begin
    # gv.calibration (default off): after the normal solve, additionally locally
    # minimizes every structurally-feasible-but-not-stable solution phase and appends
    # non-duplicate results to out.mSS_vec, tagged info="calib". See
    # calibration_output_struct in dump_function.c.
    data = Initialize_MAGEMin("ig", verbose=-1);

    # off by default: byte-identical to the existing "pointwise tests" reference values,
    # and no "calib"-tagged entries at all -- the new code path must be fully inert
    out_off = single_point_minimization(8.0, 800.0, data; test=0)
    @test out_off.G_system ≈ -797.7873865220898
    @test sort(out_off.ph) == sort(["spl", "cpx", "opx", "ol"])
    @test !any(m -> m.info == "calib", out_off.mSS_vec)

    # explicit calibration=false must match the implicit default above
    out_false = single_point_minimization(8.0, 800.0, data; test=0, calibration=false)
    @test out_false.G_system ≈ out_off.G_system
    @test !any(m -> m.info == "calib", out_false.mSS_vec)

    # on: same stable assemblage/energy (calibration must not perturb the real solve),
    # plus new "calib" entries for phases that are structurally feasible but not stable
    out_on = single_point_minimization(8.0, 800.0, data; test=0, calibration=true)
    @test out_on.G_system ≈ out_off.G_system
    @test sort(out_on.ph) == sort(out_off.ph)

    calib_entries = filter(m -> m.info == "calib", out_on.mSS_vec)
    @test length(calib_entries) > 0

    # dedup rule: no "calib" entry duplicates an already-stable phase (a distinct local
    # minimum of an already-stable MODEL is legitimate -- see plan doc -- but none of
    # KLB-1's own stable phases should reappear verbatim here)
    @test !any(m -> m.ph_name in out_on.ph, calib_entries)

    # every reported driving force is finite and small in magnitude -- not on the order
    # of gam_tot itself, which is what the double-subtraction bug (see plan doc) produced
    # before it was fixed (899 for a phase whose real answer was 4.4)
    @test all(m -> isfinite(m.deltaG) && abs(m.deltaG) < 30.0, calib_entries)

    Finalize_MAGEMin(data)
end

@testset verbose=true "specify bulk rock" begin
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

@testset verbose=true "PT adaptive refinement" begin
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

@testset verbose=true "Trace-element partitioning model" begin
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0);
    P,T     = 6.0, 699.0
    Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
    X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.2];
    sys_in  = "wt"
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, name_solvus=true);
    Finalize_MAGEMin(data)


    # create database on the fly
    el      = ["Li","Zr"]
    ph      = ["q","afs","pl","bi","opx","cd","mu","amp","fl","cpx","g","zrc"]
    KDs     = ["0.17" "0.01";"0.14 * T_C/1000.0 + [:bi].compVariables[1]" "0.01";"0.33 + 0.01*P_kbar" "0.01";"1.67 * P_kbar / 10.0 + T_C/1000.0" "0.01";"0.2" "0.01";"125" "0.01";"0.82" "0.01";"0.2" "0.01";"0.65" "0.01";"0.26" "0.01";"0.01" "0.01";"0.01" "0.0"] 
    C0      = [100.0,400.0] #starting concentration of elements in ppm (ug/g)
    dtb     = "mp"

    KDs_database = create_custom_KDs_database(el, ph, KDs)

    out_TE = TE_prediction(out, C0, KDs_database, dtb; ZrSat_model = "CB");

    @test out_TE.Cliq[1] ≈ 189.11851095903208 rtol=1e-3
    @test out_TE.Cliq[2] ≈ 47.86020212957779  rtol=1e-3
end


@testset verbose=true "Trace-element partitioning + saturation models" begin
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0)
    P, T    = 6.0, 699.0
    Xoxides = ["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","H2O","O"]
    X       = [58.509, 1.022, 14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.2]
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in="wt", name_solvus=true)

    # Li with real KDs across major phases; Zr/P2O5/S/CO2 saturation-controlled (KDs=0)
    # saturation phases (zrc, fapt, sulf, fl) are auto-added by SaturationConfig
    el  = ["Li","Zr","P2O5","S","CO2"]
    ph  = ["q","afs","pl","bi","opx","cd","mu","amp","fl","cpx","g"]
    KDs = ["0.17"   "0.01" "0.0" "0.0" "0.0";
           "0.14 * T_C/1000.0 + [:bi].compVariables[1]" "0.01" "0.0" "0.0" "0.0";
           "0.33 + 0.01*P_kbar" "0.01" "0.0" "0.0" "0.0";
           "1.67 * P_kbar / 10.0 + T_C/1000.0" "0.01" "0.0" "0.0" "0.0";
           "0.2"  "0.01" "0.0" "0.0" "0.0";
           "125"  "0.01" "0.0" "0.0" "0.0";
           "0.82" "0.01" "0.0" "0.0" "0.0";
           "0.2"  "0.01" "0.0" "0.0" "0.0";
           "0.65" "0.01" "0.0" "0.0" "0.0";
           "0.26" "0.01" "0.0" "0.0" "0.0";
           "0.01" "0.01" "0.0" "0.0" "0.0"]
    C0      = [100.0, 400.0, 1000.0, 1000.0, 500.0]
    KDs_dtb = create_custom_KDs_database(el, ph, KDs)
    sat     = SaturationConfig(Zr="CB", P2O5="HWBea92", S="Liu07", CO2="SY26")
    out_TE  = TE_prediction(out, C0, KDs_dtb, "mp"; sat=sat)

    tol = 1e-3
    @test out_TE.Cliq[1] ≈ 189.83559381921782    rtol=tol   # Li
    @test out_TE.Cliq[2] ≈ 47.86020212957779     rtol=tol   # Zr at saturation
    @test out_TE.Cliq[3] ≈ 133.18203710723262    rtol=tol   # P2O5 at saturation
    @test out_TE.Cliq[4] ≈ 16.185494729785756    rtol=tol   # S at saturation
    @test out_TE.Cliq[5] ≈ 1500.2696597838953    rtol=tol   # CO2 capped at saturation
    @test out_TE.zrc_wt  ≈ 0.0007395168468941766  rtol=tol
    @test out_TE.fapt_wt ≈ 0.0023231379042603197   rtol=tol
    @test out_TE.sulf_wt ≈ 0.0027220048406193923   rtol=tol
    @test out_TE.Sat_CO2_liq ≈ 1500.2696597838953  rtol=tol
    @test out_TE.fl_CO2_wt   ≈ 0.0003371727572556881 rtol=tol

    Finalize_MAGEMin(data)
end


@testset verbose=true "CO lattice-strain TE database (Cornet 2017)" begin
    # KLB-1 peridotite (predefined test=0), ig database.
    # At P=10.01 kbar, T=1300°C the stable assembly is liq + cpx
    # (verified in test_diagram_test0.jl), so D_cpx is exercised.
    data    = Initialize_MAGEMin("ig", verbose=-1, solver=0)
    data    = use_predefined_bulk_rock(data, 0)
    P, T    = 10.01, 1300.0
    out     = single_point_minimization(P, T, data, name_solvus=true)
    Finalize_MAGEMin(data)

    # build the CO KDs database — returns custom_KDs_database directly
    KDs_database = get_TE_database("CO")

    # structural checks on the database
    @test KDs_database isa custom_KDs_database
    @test length(KDs_database.element_name) == 28
    @test KDs_database.element_name == TE_names
    @test length(KDs_database.phase_name)   == 7   # cpx gt opx pl ol hb amp

    # run partitioning — flat 10 ppm for all 28 elements
    C0     = fill(10.0, 28)
    dtb    = "ig"
    out_TE = TE_prediction(out, C0, KDs_database, dtb)

    # output structure
    @test length(out_TE.elements) == 28
    @test out_TE.elements         == TE_names
    @test length(out_TE.Cliq)     == 28

    # melt is present and carries all elements
    @test out_TE.liq_wt_norm  > 0.0
    @test all(isfinite, out_TE.Cliq)
    @test all(>(0.0),   out_TE.Cliq)

    # mass balance: bulk_D is finite and positive
    @test isfinite(out_TE.bulk_D)
    @test out_TE.bulk_D > 0.0
end


@testset verbose=true "Saturation models" begin
    # using MAGEMin_C
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0);
    P,T     = 6.0, 699.0
    Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
    X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.2];
    X_mol, Xoxides  = convertBulk4MAGEMin(X, Xoxides,"wt","mp"); sys_in   = "mol"
    X_mol ./= sum(X_mol)                                                    # normalize to 1.0

    el      = ["Zr","P2O5","S"]
    ph      = ["zrc","fapt","sulf"]
    KDs     = [ "0.0" "0.0" "0.0";
                "0.0" "0.0" "0.0";
                "0.0" "0.0" "0.0"]                                          # phase crystallized from saturation models have 0.0 KDs

    C0      = [400.0, 1000, 1000.0]                                        # starting concentration of elements in ppm (ug/g)
    dtb     = "mp"

    KDs_dtb = create_custom_KDs_database(el, ph, KDs)

    out      = Vector{out_struct}(undef,1)
    out_TE   = Vector{out_TE_struct}(undef,1)

    X       = copy(X_mol)
    tol     = 1e-6
    res     = 1.0
    n0      = 0.0
    ite     = 0
    while res > tol && ite < 32
        out[1]     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, name_solvus=true);
        out_TE[1]  = TE_prediction(out[1] , C0, KDs_dtb, dtb; 
                                ZrSat_model     = "CB",
                                P2O5Sat_model   = "HWBea92",
                                SSat_model      = "Liu07",
                                norm_TE         = false);

        X       =  X_mol .- out_TE[1].bulk_cor_mol

        res     = abs(n0 - vec_norm(out_TE[1].bulk_cor_mol))
        n0      = vec_norm(out_TE[1].bulk_cor_mol)
        println("   Iteration $ite: residual = $res")
        ite    += 1
        if ite == 32
            @warn "Saturation model did not converge in 32 iterations, residual is $res"
        end
    end

    @test out_TE[1].zrc_wt  ≈ 0.0007885214282942672    rtol=1e-3
    @test out_TE[1].sulf_wt ≈ 0.002722767470774445      rtol=1e-3
    @test out_TE[1].fapt_wt ≈ 0.0023191756689226023      rtol=1e-3

    Finalize_MAGEMin(data)
end


@testset verbose=true "Saturation models — solve_with_saturation + CO2 (SY26)" begin
    # Same metapelite as above.  Uses the new SaturationConfig / solve_with_saturation API.
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0)
    P, T    = 6.0, 699.0
    Xoxides = ["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","H2O","O"]
    X       = [58.509, 1.022, 14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.2]
    X_mol, Xoxides = convertBulk4MAGEMin(X, Xoxides, "wt", "mp"); sys_in = "mol"
    X_mol ./= sum(X_mol)

    el      = ["Zr", "P2O5", "S", "CO2"]
    C0      = [400.0, 1000.0, 1000.0, 500.0]   # CO2 well below saturation → fl_CO2_wt = 0

    # phases are added automatically by SaturationConfig via _augment_KDs_for_saturation
    KDs_dtb = create_custom_KDs_database(el)

    sat = SaturationConfig(Zr="CB", P2O5="HWBea92", S="Liu07", CO2="SY26")

    out, out_TE, converged, n_iter = solve_with_saturation(P, T, data, X_mol, Xoxides, C0, KDs_dtb, "mp";
                                                            sat=sat, sys_in=sys_in)

    @test converged
    # Zr / S / P2O5 must be unchanged (CO2 undersaturated → zero fl weight, sum_wt unchanged)
    @test out_TE.zrc_wt  ≈ 0.0007882595805073063        rtol=1e-3
    @test out_TE.sulf_wt ≈ 0.002722767470774445         rtol=1e-3
    @test out_TE.fapt_wt ≈ 0.0023191756689226023        rtol=1e-3
    # CO2 saturation was computed and 50 ppm is well below the limit
    @test !isnan(out_TE.Sat_CO2_liq)
    @test out_TE.fl_CO2_wt ≈ 0.0003324204260014847      rtol=1e-2

    Finalize_MAGEMin(data)
end


@testset verbose=true "test monazite saturation" begin
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0)
    P, T    = 6.0, 800.0
    Xoxides = ["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","H2O","O"]
    X       = [58.509, 1.022, 14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.2]
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in="wt")

    Sat_LREE_1  = MAGEMin_C.monazite_saturation(out, 1.0)
    Sat_LREE_05 = MAGEMin_C.monazite_saturation(out, 0.5)

    @test Sat_LREE_1  ≈ 440.4193371817389 rtol=1e-3
    @test Sat_LREE_05  < Sat_LREE_1

    Cliq_LREE = [800.0, 1600.0, 180.0, 700.0, 140.0]
    mnz_wt, P2O5_wt, LREE_wt = MAGEMin_C.adjust_bulk_4_monazite(Cliq_LREE, Sat_LREE_1, out.frac_M_wt)

    @test mnz_wt   > 0.0
    @test P2O5_wt  > 0.0
    @test sum(LREE_wt) ≈ (mnz_wt - P2O5_wt) atol=1e-8

    Finalize_MAGEMin(data)
end


@testset verbose=true "Saturation models — solve_with_saturation + monazite (Stepanov12) vs. apatite competition" begin
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0)
    P, T    = 6.0, 730.0
    Xoxides = ["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","H2O","O"]
    X       = [58.509, 1.022, 14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 3.0, 0.2]
    X_mol, Xoxides = convertBulk4MAGEMin(X, Xoxides, "wt", "mp"); sys_in = "mol"
    X_mol ./= sum(X_mol)

    el      = ["La","Ce","Pr","Nd","Sm","Eu","Gd","Y","Th","U","P2O5"]
    C0      = [400.0, 800.0, 90.0, 350.0, 70.0, 16.0, 40.0, 200.0, 85.0, 20.0, 1200.0]
    KDs_dtb = create_custom_KDs_database(el)

    sat = SaturationConfig(P2O5="HWBea92", Mnz="Stepanov12")

    out, out_TE, converged, n_iter = solve_with_saturation(P, T, data, X_mol, Xoxides, C0, KDs_dtb, "mp";
                                                            sat=sat, sys_in=sys_in)

    @test converged
    @test sum(out_TE.Cliq[1:5]) ≈ out_TE.Sat_LREE_liq   rtol=1e-6
    @test out_TE.fapt_wt ≈ 0.0027272216144428353 rtol=1e-3
    @test out_TE.mnz_wt  ≈ 0.0001369575950076943 rtol=1e-2
    @test out_TE.X_mnz_LREE ≈ 0.8573455544311409 rtol=1e-2
    @test out_TE.X_mnz_LREE < 1.0

    Finalize_MAGEMin(data)
end


@testset verbose=true "test monazite saturation — Montel93 and Maimaiti19" begin
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0)
    P, T    = 4.0, 750.0
    Xoxides = ["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","H2O","O"]
    X       = [70.0, 0.3, 15.0, 2.5, 0.05, 1.0, 1.5, 3.5, 4.0, 3.0, 0.1]
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in="wt")

    M = MAGEMin_C._montel_maimaiti_M(out)
    @test M ≈ 1.024404401329306 rtol=1e-3
    @test M > 0.0

    Sat_Montel_1  = MAGEMin_C.monazite_saturation(out, 1.0; model="Montel93")
    Sat_Montel_07 = MAGEMin_C.monazite_saturation(out, 0.7; model="Montel93")
    @test Sat_Montel_1  ≈ 118.31340706544371 rtol=1e-3
    @test Sat_Montel_07 < Sat_Montel_1   # less pure REE-phosphate character -> lower solubility

    Sat_Maim_1  = MAGEMin_C.monazite_saturation(out, 1.0; model="Maimaiti19")
    Sat_Maim_07 = MAGEMin_C.monazite_saturation(out, 0.7; model="Maimaiti19")
    @test Sat_Maim_1  ≈ 252.07971509133452 rtol=1e-3
    @test Sat_Maim_07 < Sat_Maim_1

    @test MAGEMin_C.monazite_saturation(out; model="not_a_model") == -1

    Finalize_MAGEMin(data)
end


@testset verbose=true "Saturation models — solve_with_saturation + monazite (Montel93) vs. apatite competition" begin
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0)
    P, T    = 6.0, 730.0
    Xoxides = ["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","H2O","O"]
    X       = [58.509, 1.022, 14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 3.0, 0.2]
    X_mol, Xoxides = convertBulk4MAGEMin(X, Xoxides, "wt", "mp"); sys_in = "mol"
    X_mol ./= sum(X_mol)

    el      = ["La","Ce","Pr","Nd","Sm","Eu","Gd","Y","Th","U","P2O5"]
    C0      = [400.0, 800.0, 90.0, 350.0, 70.0, 16.0, 40.0, 200.0, 85.0, 20.0, 1200.0]
    KDs_dtb = create_custom_KDs_database(el)

    sat = SaturationConfig(P2O5="HWBea92", Mnz="Montel93")

    out, out_TE, converged, n_iter = solve_with_saturation(P, T, data, X_mol, Xoxides, C0, KDs_dtb, "mp";
                                                            sat=sat, sys_in=sys_in)

    @test converged
    @test sum(out_TE.Cliq[1:5]) ≈ out_TE.Sat_LREE_liq   rtol=1e-6
    @test out_TE.Sat_LREE_liq ≈ 150.08389000223016 rtol=1e-2
    # the P2O5 joint clamp binds identically regardless of which saturation-front model is used
    @test out_TE.fapt_wt ≈ 0.0027272216144428353 rtol=1e-3
    @test out_TE.mnz_wt  ≈ 0.0001369575950076943 rtol=1e-2
    @test out_TE.X_mnz_LREE ≈ 0.8573455544311409 rtol=1e-2

    Finalize_MAGEMin(data)
end


@testset verbose=true "Saturation models — solve_with_saturation + monazite (Maimaiti19) vs. apatite competition" begin
    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0)
    P, T    = 6.0, 730.0
    Xoxides = ["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","H2O","O"]
    X       = [58.509, 1.022, 14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 3.0, 0.2]
    X_mol, Xoxides = convertBulk4MAGEMin(X, Xoxides, "wt", "mp"); sys_in = "mol"
    X_mol ./= sum(X_mol)

    el      = ["La","Ce","Pr","Nd","Sm","Eu","Gd","Y","Th","U","P2O5"]
    C0      = [400.0, 800.0, 90.0, 350.0, 70.0, 16.0, 40.0, 200.0, 85.0, 20.0, 1200.0]
    KDs_dtb = create_custom_KDs_database(el)

    sat = SaturationConfig(P2O5="HWBea92", Mnz="Maimaiti19")

    out, out_TE, converged, n_iter = solve_with_saturation(P, T, data, X_mol, Xoxides, C0, KDs_dtb, "mp";
                                                            sat=sat, sys_in=sys_in)

    @test converged
    @test sum(out_TE.Cliq[1:5]) ≈ out_TE.Sat_LREE_liq   rtol=1e-6
    @test out_TE.Sat_LREE_liq ≈ 296.1406721714649 rtol=1e-2
    @test out_TE.fapt_wt ≈ 0.0027272216144428353 rtol=1e-3
    @test out_TE.mnz_wt  ≈ 0.0001369575950076943 rtol=1e-2
    @test out_TE.X_mnz_LREE ≈ 0.8573455544311409 rtol=1e-2

    Finalize_MAGEMin(data)
end


@testset verbose=true "Yak25 accessory-phase TE database — zircon/apatite/monazite competition" begin
    KDs_dtb = get_TE_database("Yak25")

    @test KDs_dtb.element_name == ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Y", "Th", "U", "Zr", "P2O5"]
    @test "zrc" in KDs_dtb.phase_name
    @test "fapt" in KDs_dtb.phase_name
    @test "mnz" in KDs_dtb.phase_name

    dummy_out = single_point_minimization(6.0, 800.0, Initialize_MAGEMin("mp", verbose=-1, solver=0),
                                            X=[58.509, 1.022, 14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.2],
                                            Xoxides=["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","H2O","O"],
                                            sys_in="wt")

    id_zrc = findfirst(KDs_dtb.phase_name .== "zrc")
    id_mnz = findfirst(KDs_dtb.phase_name .== "mnz")
    id_fapt = findfirst(KDs_dtb.phase_name .== "fapt")
    id_Th = findfirst(KDs_dtb.element_name .== "Th")
    id_Zr = findfirst(KDs_dtb.element_name .== "Zr")
    id_P2O5 = findfirst(KDs_dtb.element_name .== "P2O5")

    @test KDs_dtb.KDs_expr[id_zrc,  id_Th](dummy_out) ≈ 62.0
    @test KDs_dtb.KDs_expr[id_mnz,  id_Th](dummy_out) ≈ 5000.0
    @test KDs_dtb.KDs_expr[id_fapt, id_Th](dummy_out) ≈ 23.0
    @test all(KDs_dtb.KDs_expr[i, id_Zr](dummy_out)   == 0.0 for i in eachindex(KDs_dtb.phase_name))
    @test all(KDs_dtb.KDs_expr[i, id_P2O5](dummy_out) == 0.0 for i in eachindex(KDs_dtb.phase_name))

    P, T    = 6.0, 900.0
    Xoxides = ["SiO2","TiO2","Al2O3","FeO","MnO","MgO","CaO","Na2O","K2O","H2O","O"]
    X       = [45.0, 1.0, 25.0, 12.0, 0.2, 10.0, 1.0, 0.5, 1.0, 3.0, 0.2]
    X_mol, Xoxides = convertBulk4MAGEMin(X, Xoxides, "wt", "mp"); sys_in = "mol"
    X_mol ./= sum(X_mol)

    el   = KDs_dtb.element_name
    vals = Dict("La"=>40.0, "Ce"=>80.0, "Pr"=>9.0, "Nd"=>35.0, "Sm"=>7.0, "Eu"=>1.6, "Gd"=>4.0,
                "Y"=>20.0, "Th"=>8.5, "U"=>2.0, "Zr"=>150.0, "P2O5"=>1200.0)
    C0   = [vals[e] for e in el]

    sat = SaturationConfig(Zr="CB", P2O5="HWBea92", Mnz="Stepanov12")
    data = Initialize_MAGEMin("mp", verbose=-1, solver=0)
    out, out_TE, converged, n_iter = solve_with_saturation(P, T, data, X_mol, Xoxides, C0, KDs_dtb, "mp";
                                                            sat=sat, sys_in=sys_in)

    @test converged
    mnz_idx  = findfirst(out_TE.ph_TE .== "mnz")
    zrc_idx  = findfirst(out_TE.ph_TE .== "zrc")
    fapt_idx = findfirst(out_TE.ph_TE .== "fapt")
    @test !isnothing(mnz_idx) && !isnothing(zrc_idx) && !isnothing(fapt_idx)
    @test out_TE.Cmin[zrc_idx,  id_Th] > 0.0
    @test out_TE.Cmin[fapt_idx, id_Th] > 0.0
    @test out_TE.Cmin[mnz_idx,  id_Th] > 0.0

    Finalize_MAGEMin(data)
end


@testset verbose=true "remove solution phase" begin

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
    @test sort(out.ph) == sort(["H2O", "cd", "cor", "fsp", "g", "ilm", "q", "sp"])
    Finalize_MAGEMin(data)
end

@testset verbose=true "view array PT" begin

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

@testset verbose=true "convert bulk rock" begin

    bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    bulk_in    = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    bulk_rock,ox  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"wt","ig");

    @test bulk_rock ≈ [46.12597764761598, 8.52489397284109, 11.805554573333653, 14.383528505497756, 6.471419573392541, 0.35839516987780934, 1.7264383468329216, 0.4871154401452383, 0.5876614012114892, 0.0, 9.529015369251512]

    bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"];
    bulk_in    = [69.64; 13.76; 1.77; 1.73; 4.32; 0.4; 2.61; 2.41; 0.80; 0.07; 0.0];
    bulk_rock,ox  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"wt","mp");

    @test bulk_rock ≈ [76.57038397179574, 8.914984523583415, 2.0849576977131403, 2.835783318610597, 4.30275071755529, 1.8302970975627948, 2.568605789798099, 0.6615823604771729, 0.16546809116073818, 0.06518643174302832, 0.0]
end

@testset verbose=true "FeO + Oᵉˣᵗʳᵃ -> Feᵀᵒᵗᵃˡ + Oᵀᵒᵗᵃˡ" begin

    bulk_in_ox  =  ["SiO2"; "CaO"; "Al2O3"; "MgO"; "Na2O"; "FeO"; "Cr2O3"; "O"]    
    bulk_in     = [38.83, 2.94, 2.03, 50.03, 0.11, 5.69, 0.19, 0.17] # Pyrolite
    bulk_mod, bulk_ox  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"mol","sb24");
    @test bulk_mod ≈ [36.742997728993195, 2.7819833459500383, 1.9208932626797883, 47.34102952308858, 0.10408781226343683, 5.545041635124907, 0.17978803936411814, 5.384178652535959]
    @test bulk_ox == ["SiO2"; "CaO"; "Al2O3"; "MgO"; "Na2O"; "O"; "Cr2O3"; "Fe"]

    bulk_in_ox  =  ["SiO2"; "CaO"; "Al2O3"; "MgO"; "Na2O"; "FeO"; "Cr2O3"; "O"]    
    bulk_in     = [45.13, 3.19, 4.02, 39.00, 0.13, 7.91, 0.56, 0.05] # Pyrolite in wt%
    bulk_mod, bulk_ox  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"wt","sb24");
    @test bulk_mod ≈ [36.744275604478226, 2.7825114393856785, 1.9286370551405936, 47.33842901832809, 0.10259969288150061, 5.538090389905078, 0.1802300870140965, 5.385226712866728]
    @test bulk_ox == ["SiO2"; "CaO"; "Al2O3"; "MgO"; "Na2O"; "O"; "Cr2O3"; "Fe"]

    bulk_in_ox  =  ["SiO2"; "CaO"; "Al2O3"; "MgO"; "Na2O"; "FeO"; "Cr2O3"; "Fe2O3"]    
    bulk_in     = [38.83, 2.94, 2.03, 50.03, 0.11, 5.35, 0.19, 0.17] # Pyrolite
    bulk_mod, bulk_ox  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"mol","sb24");
    @test bulk_mod ≈ [36.74299772899319, 2.7819833459500383, 1.9208932626797879, 47.34102952308858, 0.1040878122634368, 5.545041635124905, 0.17978803936411814, 5.384178652535957]
    @test bulk_ox == ["SiO2"; "CaO"; "Al2O3"; "MgO"; "Na2O"; "O"; "Cr2O3"; "Fe"]

    bulk_in_ox  =  ["SiO2"; "CaO"; "Al2O3"; "MgO"; "Na2O"; "Fe"; "Cr2O3"; "O"]    
    bulk_in     = [38.83, 2.94, 2.03, 50.03, 0.11, 5.69, 0.19, 5.86] # Pyrolite
    bulk_mod, bulk_ox  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"mol","sb24");
    @test bulk_mod ≈ [36.74299772899319, 2.7819833459500383, 1.9208932626797879, 47.34102952308857, 0.1040878122634368, 5.545041635124906, 0.1797880393641181, 5.384178652535958]
    @test bulk_ox == ["SiO2"; "CaO"; "Al2O3"; "MgO"; "Na2O"; "O"; "Cr2O3"; "Fe"]

end


@testset verbose=true "test Seismic velocities & modulus - VRH" begin
    # Call optimization routine for given P & T & bulk_rock
    data         = Initialize_MAGEMin("ig", verbose=-1; seismicScheme="VRH", seismicWeightFactor=0.5);
    test        = 0;
    data         = use_predefined_bulk_rock(data, test)
    P           = 8.0
    T           = 1200.0
    out         = point_wise_minimization(P,T, data)
    tol         = 1.5e-2;

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


@testset verbose=true "test Seismic velocities & modulus - HS" begin
    # Call optimization routine for given P & T & bulk_rock
    data         = Initialize_MAGEMin("ig", verbose=-1; seismicScheme="HS", seismicWeightFactor=0.95);
    test        = 0;
    data         = use_predefined_bulk_rock(data, test)
    P           = 8.0
    T           = 1200.0
    out         = point_wise_minimization(P,T, data)
    tol         = 1.5e-2;

    @test abs(out.bulkMod - 94.9207227314658           )  < tol
    @test abs(out.shearMod - 56.23048827772405         )  < tol
    @test abs(out.Vs - 4.189889782462432               )  < tol
    @test abs(out.Vp - 7.282937716875037               )  < tol
    @test abs(out.Vs_S - 4.305197788922828             )  < tol
    @test abs(out.Vp_S - 7.373574860960485             )  < tol
    @test abs(out.bulkModulus_M - 27.774116805966923   )  < tol
    @test abs(out.bulkModulus_S - 95.4968207451378     )  < tol
    @test abs(out.shearModulus_S - 59.68335800904911   )  < tol

    Finalize_MAGEMin(data)
end


@testset verbose=true "test Mantle HP13" begin

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

@testset verbose=true "test ume" begin

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
@testset verbose=true "test apfu" begin

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


@testset verbose=true "Text initial guess" begin

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


@testset verbose=true "Metastability function" begin
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
    @test abs(out2.G_system + 791.4602868457285) < 1e-6
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

@testset verbose = true "Test gbase override" begin

    #= Additive shift on the reference Gibbs energy of one endmember of "liq" =#
    dtb         = 2             # igneous
    ss_id       = 9             # liq
    n_Gs        = 1
    em_ids      = [1]
    dG          = reshape([-50.0, 0.0, 0.0], 1, 3)   # constant -50 kJ shift, no T/P dependence
    new_gbase   = Vector{MAGEMin_C.gbase_data{Float64,Int64}}(undef, 1)
    new_gbase[1]= MAGEMin_C.gbase_data(dtb, ss_id, n_Gs, em_ids, dG)

    data        = Initialize_MAGEMin("ig", verbose=-1);
    P,T         = 10.0, 1100.0
    Xoxides     = ["SiO2","Al2O3","CaO","MgO","FeO","Fe2O3","K2O","Na2O","TiO2","Cr2O3","H2O"]
    X           = [48.43,15.19,11.57,10.13,6.65,1.64,0.59,1.87,0.68,0.0,3.0]
    sys_in      = "wt"

    out_base    = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
    out_shift   = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, gbase=new_gbase)
    Finalize_MAGEMin(data)

    # lowering one endmember's reference G must lower (or keep equal) the system's total G
    @test out_shift.G_system < out_base.G_system
    @test abs(out_shift.G_system + 918.2105250940679) < 1e-4
end

@testset verbose = true "exclude_DEW_species" begin

    #= resolution is purely db_infos-based - no Initialize_MAGEMin/minimization needed =#
    excl                    = exclude_DEW_species("all", "H4SiO4")
    @test length(excl)      == 1
    @test excl[1].dtb       == 8    # "all"
    @test excl[1].ss_ids    == 44   # DEW_S24's position in "all"'s solution-phase list
    @test excl[1].n_Gs      == 1
    @test excl[1].em_ids    == [38]
    @test excl[1].dG        == reshape([1.0e6, 0.0, 0.0], 1, 3)

    excl2                   = exclude_DEW_species("all", ["HCOOH", "H4SiO4"])
    @test length(excl2)     == 1
    @test excl2[1].n_Gs     == 2
    @test excl2[1].em_ids   == [43, 38]

    #= unknown species name / database without a DEW phase both error =#
    @test_throws ErrorException exclude_DEW_species("all", "NOT_A_SPECIES")
    @test_throws ErrorException exclude_DEW_species("sb11", "HCOOH")

    ss_list = ["liq_W14", "fsp_H22", "bi_W14", "g_W14", "ep_H11", "ma_W14", "mu_W14", "opx_W14",
               "sa_W14", "cd_W14", "st_W14", "chl_W14", "ctd_W14", "sp_W02", "ilm_W00", "DEW_S24"]
    pp_list = ["q", "crst", "trd", "coe", "stv", "law", "ky", "sill", "and", "ru", "sph", "prl"]

    data    = Initialize_MAGEMin("all", verbose=false, solver=0)
    P, T    = 10.0, 400.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "Cr2O3"; "H2O"; "CO2"; "S"]

    X1      = [0.62212, 0.1122, 0.0, 0.03486, 0.05557, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.17525, 0.0, 0.0]

    out_base    = single_point_minimization(P, T, data, X=X1, Xoxides=Xoxides, sys_in="mol", ss_list=ss_list, pp_list=pp_list)
    idx_b       = findfirst(n -> startswith(n, "DEW"), out_base.ph)
    @test !isnothing(idx_b)
    ssb         = out_base.SS_vec[idx_b]
    i_h4        = findfirst(==("H4SiO4"), ssb.emNames)
    @test i_h4  == 38               # matches excl[1].em_ids above - db_infos and the runtime index agree
    @test ssb.emFrac[i_h4] > 1e-5  # meaningfully present in the baseline
    @test abs(out_base.G_system + 865.7668984405251) < 1e-4

    out_shift   = single_point_minimization(P, T, data, X=X1, Xoxides=Xoxides, sys_in="mol", ss_list=ss_list, pp_list=pp_list, gbase=excl)
    idx_s       = findfirst(n -> startswith(n, "DEW"), out_shift.ph)
    @test !isnothing(idx_s)
    sss         = out_shift.SS_vec[idx_s]
    @test sss.emFrac[findfirst(==("H4SiO4"), sss.emNames)] == 0.0
    @test out_shift.G_system > out_base.G_system
    @test abs(out_shift.G_system + 865.7666083401803) < 1e-4

    X2          = [0.61, 0.1122, 0.0, 0.03486, 0.05557, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.17525, 0.012, 0.0]

    out_base2   = single_point_minimization(P, T, data, X=X2, Xoxides=Xoxides, sys_in="mol", ss_list=ss_list, pp_list=pp_list)
    out_shift2  = single_point_minimization(P, T, data, X=X2, Xoxides=Xoxides, sys_in="mol", ss_list=ss_list, pp_list=pp_list, gbase=excl2)
    idx_s2      = findfirst(n -> startswith(n, "DEW"), out_shift2.ph)
    @test !isnothing(idx_s2)
    sss2        = out_shift2.SS_vec[idx_s2]
    @test sss2.emFrac[findfirst(==("HCOOH"),  sss2.emNames)] == 0.0
    @test sss2.emFrac[findfirst(==("H4SiO4"), sss2.emNames)] == 0.0
    @test out_shift2.G_system > out_base2.G_system

    Finalize_MAGEMin(data)
end

@testset verbose=true "test matrix (2D grid) input for multi_point_minimization" begin
    data    = Initialize_MAGEMin("ig", verbose=-1);
    Xoxides = ["SiO2","Al2O3","CaO","MgO","FeO","Fe2O3","K2O","Na2O","TiO2","Cr2O3","H2O"]

    # 2×2 P-T grid
    P_grid  = [8.0  9.0; 10.0 11.0]
    T_grid  = [800.0 850.0; 900.0 950.0]

    # per-node bulk composition matrix: 4 rows (one per grid node, column-major order), 11 columns
    X1      = [48.43, 15.19, 11.57, 10.13, 6.65, 1.64, 0.59, 1.87, 0.68, 0.0, 3.0]
    X2      = [50.0,  14.0,  10.0,  12.0,  7.0,  1.5,  0.6,  2.0,  0.7,  0.0, 2.0]
    X3      = [47.0,  16.0,  12.0,   9.0,  6.0,  2.0,  0.5,  1.5,  0.6,  0.0, 4.0]
    X4      = [52.0,  13.0,   9.0,  11.0,  8.0,  1.2,  0.7,  2.2,  0.8,  0.0, 1.5]
    X_grid  = [X1'; X2'; X3'; X4']   # (4, 11)

    out     = multi_point_minimization(P_grid, T_grid, data, X=X_grid, Xoxides=Xoxides, sys_in="wt", progressbar=false)

    @test size(out) == (2, 2)
    @test out[1,1].P_kbar ≈ 8.0
    @test out[1,1].T_C    ≈ 800.0
    @test out[2,2].P_kbar ≈ 11.0
    @test out[2,2].T_C    ≈ 950.0

    Finalize_MAGEMin(data)
end

@testset verbose=true "test SY26 volatile saturation (Sun & Yao 2026)" begin
    # Low-MgO, low-K2O rhyolite matching Table S3 of Sun & Yao (2026).
    # MgO≈0 and K2O≈2 mol% are critical: b3=-1019 for X_MgO and b4=-452 for X_K2O²
    # dominate the molecular-CO2 term; typical rhyolites with K2O>3% yield far too little CO2.
    data    = Initialize_MAGEMin("ig", verbose=-1, solver=0)
    P, T    = 2.5, 900.0
    Xoxides = ["SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O"]
    X       = [78.0, 0.30, 12.0, 2.0, 0.01, 1.0, 4.0, 2.0, 5.0]
    sys_in  = "mol"
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
    P_H2O_auto, P_CO2_auto, S_CO2_auto = MAGEMin_C.CO2_from_dissolved_H2O(out)


    @test out.frac_M > 0.0

    P_bar = P * 1000.0

    # pure H2O fluid: only H2O partial pressure provided
    S_H2O_pure, _ = MAGEMin_C.volatile_saturation_SY26(out; P_H2O = P_bar)
    @test !isnan(S_H2O_pure)
    @test S_H2O_pure > 0.0

    # mixed H2O–CO2 fluid: P_H2O = 0.7·P, P_CO2 = 0.3·P
    S_H2O_mix, S_CO2_mix = MAGEMin_C.volatile_saturation_SY26(out; P_H2O = 0.7*P_bar, P_CO2 = 0.3*P_bar)
    @test !isnan(S_H2O_mix) && S_H2O_mix > 0.0
    @test !isnan(S_CO2_mix) && S_CO2_mix > 0.0
    @test S_H2O_mix < S_H2O_pure   # lower P_H2O → less H2O solubility

    # CO2_from_dissolved_H2O: use half the pure-H2O saturation so P_CO2 > 0
    S_H2O_in              = 0.5 * S_H2O_pure
    P_H2O_inv, P_CO2_inv, S_CO2_inv = MAGEMin_C.CO2_from_dissolved_H2O(out, S_H2O_in)
    @test !isnan(P_H2O_inv) && P_H2O_inv > 0.0
    @test !isnan(P_CO2_inv) && P_CO2_inv > 0.0
    @test !isnan(S_CO2_inv) && S_CO2_inv > 0.0
    @test P_H2O_inv + P_CO2_inv ≈ P_bar  rtol=1e-4

    # round-trip: re-evaluate S_H2O at the inverted P_H2O — must recover S_H2O_in
    S_H2O_rt, _ = MAGEMin_C.volatile_saturation_SY26(out; P_H2O = P_H2O_inv)
    @test S_H2O_rt ≈ S_H2O_in  rtol=1e-4

    # no-arg overload: reads dissolved H2O from melt directly
    P_H2O_auto, P_CO2_auto, S_CO2_auto = MAGEMin_C.CO2_from_dissolved_H2O(out)
    H2O_idx      = findfirst(==("H2O"), out.oxides)
    S_H2O_melt   = out.SS_vec[out.SS_syms[:liq]].Comp_wt[H2O_idx] * 100.0
    @test !isnan(P_H2O_auto) && P_H2O_auto > 0.0
    @test P_H2O_auto + P_CO2_auto ≈ P_bar  rtol=1e-4
    # re-evaluate S_H2O at auto P_H2O — must match the melt dissolved H2O
    S_H2O_auto, _ = MAGEMin_C.volatile_saturation_SY26(out; P_H2O = P_H2O_auto)
    @test S_H2O_auto ≈ S_H2O_melt  rtol=1e-4

    Finalize_MAGEMin(data)
end

@testset verbose=true "test p2x_convert/pc_convert/lm_convert (endmember fractions -> phase Gibbs energy / local minimization)" begin
    gv, z_b, DB, splx_data      = init_MAGEMin("all")
    gv                          = use_predefined_bulk_rock(gv, 0, "all")
    gv, z_b, DB, splx_data      = pwm_init(2.0, 700.0, gv, z_b, DB, splx_data)

    p                           = Dict("ab"=>0.2, "an"=>0.2, "san"=>0.6)
    SS_ref_db                   = p2x_convert(gv, DB, "fsp_H22", p)

    em_names                    = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, SS_ref_db.EM_list, SS_ref_db.n_em))
    cv_names                    = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, SS_ref_db.CV_list, SS_ref_db.n_xeos))

    @test em_names == ["ab", "an", "san"]
    @test cv_names == ["ca", "k"]
    @test unsafe_wrap(Vector{Float64}, SS_ref_db.p, SS_ref_db.n_em) ≈ [0.2, 0.2, 0.6]
    @test unsafe_wrap(Vector{Float64}, SS_ref_db.xeos, SS_ref_db.n_xeos) ≈ [0.2, 0.6] atol=1e-6

    SS_ref_db                   = pc_convert(gv, z_b, DB, "fsp_H22", SS_ref_db)

    @test SS_ref_db.sf_ok == 1
    @test isfinite(SS_ref_db.df)
    @test SS_ref_db.df < 0.0   # sanity: molar Gibbs energy of a stable silicate is negative

    p_bad                       = Dict("ab"=>-0.4, "an"=>0.7, "san"=>0.7)
    SS_ref_bad                  = p2x_convert(gv, DB, "fsp_H22", p_bad)
    SS_ref_bad                  = pc_convert(gv, z_b, DB, "fsp_H22", SS_ref_bad)

    @test SS_ref_bad.sf_ok == 0
    @test_throws ErrorException p2x_convert(gv, DB, "fsp_H22", Dict("ab"=>0.5, "an"=>0.5))               # missing endmember
    @test_throws ErrorException p2x_convert(gv, DB, "fsp_H22", Dict("ab"=>0.2,"an"=>0.2,"san"=>0.5,"xx"=>0.1)) # unknown endmember
    @test_throws ErrorException p2x_convert(gv, DB, "not_a_phase", p)                                     # unknown phase

    gamma0                      = [-973.680237,-1768.944009,-820.744391,-693.029850,-366.575047,-949.528471,-875.814737,-1025.018224,-240.991436,-507.666146,-1306.489776,-385.336776,-592.598340,-88.853592]
    xeos0                       = [0.2, 0.6, 0.1]
    SS_ref_min                  = lm_convert(gv, z_b, DB, "fsp_H22op", gamma0, xeos0)
    @test SS_ref_min.status     == 3
    @test SS_ref_min.sf_ok      == 1
    @test isfinite(SS_ref_min.df)

    xeos_min                    = copy(unsafe_wrap(Vector{Float64}, SS_ref_min.xeos, SS_ref_min.n_xeos))
    @test xeos_min              ≈ [0.3125135893680029, 0.031900521359415714, 0.07754632896992686] atol=1e-4
    @test !isapprox(xeos_min, xeos0, atol=1e-3)

    # re-running from the converged point under the same Gamma is a fixed point
    SS_ref_min2                 = lm_convert(gv, z_b, DB, "fsp_H22op", gamma0, xeos_min)
    @test SS_ref_min2.status    == 3
    @test unsafe_wrap(Vector{Float64}, SS_ref_min2.xeos, SS_ref_min2.n_xeos) ≈ xeos_min atol=1e-6

    @test_throws ErrorException lm_convert(gv, z_b, DB, "fsp_H22op", zeros(gv.len_ox - 1), xeos0) # wrong-length gamma
    @test_throws ErrorException lm_convert(gv, z_b, DB, "fsp_H22op", gamma0, [0.2])                # wrong-length xeos
    @test_throws ErrorException lm_convert(gv, z_b, DB, "not_a_phase", gamma0, xeos0)             # unknown phase

    finalize_MAGEMin(gv, DB, z_b, splx_data)
end



@testset verbose=true "filter_DEW_species — drop chemically infeasible DEW_S24 species" begin
    data    = Initialize_MAGEMin("all", verbose=false, solver=0)
    ss_list = ["liq_W14", "fsp_H22", "bi_W14", "g_W14", "ep_H11", "ma_W14", "mu_W14", "opx_W14",
               "sa_W14", "cd_W14", "st_W14", "chl_W14", "ctd_W14", "sp_W02", "ilm_W00", "DEW_S24"]
    pp_list = ["q", "crst", "trd", "coe", "stv", "law", "ky", "sill", "and", "ru", "sph", "prl"]

    P, T    = 10.0, 400.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "O"; "MnO"; "Cr2O3"; "H2O"; "CO2"; "S"]
    X       = [0.62212, 0.1122, 0.0, 0.03486, 0.05557, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.17525, 0.0, 0.0]

    out_unfiltered = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in="mol", ss_list=ss_list, pp_list=pp_list)
    out_filtered   = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in="mol", ss_list=ss_list, pp_list=pp_list, filter_DEW_species=true)

    idx_u = findfirst(==("DEW_S24"), out_unfiltered.ph)
    idx_f = findfirst(==("DEW_S24"), out_filtered.ph)
    @test !isnothing(idx_u) && !isnothing(idx_f)

    ss_u = out_unfiltered.SS_vec[idx_u]
    ss_f = out_filtered.SS_vec[idx_f]

    # unfiltered: full DEW_S24 species list regardless of bulk feasibility
    @test length(ss_u.emNames) == 107

    # filtered: only the Al/Fe/Mg/Si/H/O-bearing species survive, since the bulk has
    # no CaO/K2O/Na2O/MnO/Cr2O3/CO2/S — matches the exact species Ca/K/Na/Mn/Cr/C/S-free chemistry predicts
    expected_species = ["Al(OH)3", "Al(OH)4-", "Al(OH)Si(OH)-", "Al+3", "Fe(H3SiO4)+", "Fe(OH)+",
                         "Fe(OH)2", "Fe(OH)3-", "Fe+2", "H+", "H3SiO4-", "H4SiO4", "H6Si2O7",
                         "H8Si3O10", "Mg(H3SiO4)+", "Mg(OH)+", "Mg(OH)2", "Mg+2", "OH-", "H2O"]
    @test sort(ss_f.emNames) == sort(expected_species)
    @test all(ss_f.emFrac .> 0.0)

    # all per-endmember fields filtered in lockstep
    n_kept = length(ss_f.emNames)
    @test length(ss_f.molality)    == n_kept
    @test length(ss_f.activity)    == n_kept
    @test length(ss_f.emFrac_wt)   == n_kept
    @test length(ss_f.emChemPot)   == n_kept
    @test length(ss_f.emComp)      == n_kept
    @test length(ss_f.emComp_wt)   == n_kept
    @test length(ss_f.emComp_apfu) == n_kept

    # phase-level fields (not per-endmember) must be untouched
    @test ss_f.Comp == ss_u.Comp
    @test ss_f.compVariables == ss_u.compVariables

    # other stable solution phases in the same run must be completely unaffected
    n_SS = length(out_unfiltered.SS_vec)
    for i in 1:n_SS
        if out_unfiltered.ph[i] != "DEW_S24"
            @test length(out_unfiltered.SS_vec[i].emNames) == length(out_filtered.SS_vec[i].emNames)
        end
    end

    Finalize_MAGEMin(data)
end


