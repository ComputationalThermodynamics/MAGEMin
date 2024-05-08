# tests for TE partitioning and Zircons saturation calculations

using MAGEMin_C

include("julia/TE_partitioning.jl")
include("julia/Zircon_saturation.jl")

using MAGEMin_C

data    = Initialize_MAGEMin("ig", verbose=0, solver=0);
P,T     = 10.0, 1500.0
Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 0.0];
sys_in  = "wt"    
out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

X       = [48.4; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 0.0];
out2     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, data_in = out)




"""
    structure that holds the result of the trace element predictive model
"""
struct tepm_struct{T}
    te          :: Vector{String}                   # Name of the trace elements
    ph          :: Union{Vector{String}, Nothing}   # Name of the phases bearing trace elements

    C0          :: Union{Vector{T}, Nothing}        # starting TE composition
    Cliq        :: Union{Vector{T}, Nothing}        # partitioned trace element composition for the liquid
    Cmin        :: Union{Matrix{T}, Nothing}        # partinioned trace element composition for the minerals

    te_pm       :: String                           # predictive model used to compute trace elements partitioning
 
    zr_sat_pm   :: String                           # used predictive model to computate zircon saturation
    zr_liq_sat  :: Union{T, Nothing}                # zircon saturation in ptr_comp_pc
    zr_wt_pc    :: Union{T, Nothing}                # zircon wt crystallized from melt
end



using MAGEMin_C
# First we provide a bulk
TE_db       = get_TE_database("TE_OL_felsic")
components  = ["SiO2","TiO2","Al2O3","O","FeO","MnO","MgO","CaO","K2O","Na2O","H2O","Li","Be","B","Sc","V","Cr","Ni","Cu","Zn","Rb","Sr","Y","Zr","Nb","Cs","Ba","La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","Pb","Th","U"]
Ton	        = [69.20111492,0.385480479,15.36802994,0.08,2.835169639,0.052720082,1.173738821,3.375444597,1.767078881,4.503466421,2,0.,0.,0.,5.524762752,39.87762745,49.26129576,15.04133712,15.15293785,63.74940359,64.5373756,484.5753262,7.459369222,150.2402264,5.430265866,2.474580356,499.0093971,29.06202556,53.78324684,5.797620808,20.54095673,3.317738791,0.961548784,2.418869313,0.315058021,1.593968618,0.296925128,0.799927928,0.128677572,0.713354049,0.133287066,3.95645511,0.722386615,15.17205921,5.69152854,1.064824497]
Bas         = [50.60777553,0.953497243,13.70435413,0.19,11.28130762,0.202560796,8.496024312,9.502380068,0.700881685,2.07434927,4,29.14258603,0.434482763,29.69200003,38.23663423,257.4346716,529.333066,208.2057375,88.87615683,91.7592182,16.60777308,163.4533209,20.74016207,66.90677472,3.808354064,1.529226981,122.8449739,6.938172601,16.04827796,2.253943183,10.18276823,3.3471043,0.915941652,3.28230146,1.417695298,3.851230952,0.914966282,2.20425,0.343734976,2.136202593,0.323405135,1.841502082,0.330971265,5.452969044,1.074692058,0.290233271]

C0_TE_idx   = [findfirst(isequal(x), components) for x in TE_db.element_name]


C0_TE       = Ton[C0_TE_idx]

# Then we call MAGENin
db          = "ig"
sys_in      = "wt"    
data        = Initialize_MAGEMin(db, verbose=false);

P       = [1.0, 2.0,1.0, 2.0]
T       = [800.0, 800.0,800.0, 800.0]
Xoxides     = ["SiO2","TiO2","Al2O3","O","FeO","MgO","CaO","K2O","Na2O","H2O"];
X1           = [69.20111492,0.385480479,15.36802994,0.08,2.835169639,1.173738821,3.375444597,1.767078881,4.503466421,2];
X2           = [69.20111492,0.385480479,15.36802994,0.08,2.835169639,1.173738821,3.375444597,1.767078881,4.503466421,2];
X       = [X1,X2,X1,X2]
XTE     = [C0_TE,C0_TE,C0_TE,C0_TE]

out, out_te     = multi_point_minimization(    P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, scp=0,
                                                tepm    = 1,
                                                te_db   = "OL",
                                                zr_sat  = "CB",
                                                te_X    = XTE )




    # # here we compute trace element partitioning and zircon saturation
    # if (tepm == 1)
    #     if (out.frac_M > 0.0 && out.frac_S > 0.0)
    #         Cliq, Cmin, ph_TE, ph_wt_norm, liq_wt_norm, Cliq_Zr, te_names  = compute_TE_partitioning(   te_X,
    #                                                                                                     out,
    #                                                                                                     dtb;
    #                                                                                                     TE_db = te_db)

    #         # Then we compute zirconium saturation
    #         Sat_zr_liq  = zirconium_saturation( out; 
    #                                             model = zr_sat)   

    #         if Cliq_Zr > Sat_zr_liq
    #             zircon_wt, SiO2_wt, O_wt  = adjust_bulk_4_zircon(Cliq_Zr, Sat_zr_liq)
    #             SiO2_id     = findall(out.oxides .== "SiO2")[1]

    #             bulk_act    = copy(out.bulk_wt)
    #             bulk_act[SiO2_id]    = out.bulk_wt[SiO2_id] - SiO2_wt 
    #             bulk_act  ./= sum(bulk_act)
    #             gv          = define_bulk_rock(gv, bulk_act, out.oxides, "wt", dtb);
    #             mSS_vec     = deepcopy(out.mSS_vec)
    #             out_cor     = point_wise_minimization_with_guess(mSS_vec, P, T, gv, z_b, DB, splx_data)

    #             Cliq, Cmin, ph_TE, ph_wt_norm, liq_wt_norm, Cliq_Zr, te_names = compute_TE_partitioning(    te_X,
    #                                                                                                         out_cor,
    #                                                                                                         dtb;
    #                                                                                                         TE_db = te_db)

    #             # Then we compute zirconium saturation
    #             Sat_zr_liq  = zirconium_saturation( out; 
    #                                                 model = zr_sat)     

    #             zircon_wt, SiO2_wt, O_wt  = adjust_bulk_4_zircon(Cliq_Zr, Sat_zr_liq)
    #         else
    #             zircon_wt = 0.0;
    #         end

    #     elseif out.frac_M == 1.0
    #         TE_dtb      =  get_TE_database("TE_OL_felsic")

    #         Sat_zr_liq  = zirconium_saturation( out; 
    #                                             model = zr_sat)     

    #         zircon_wt, SiO2_wt, O_wt  = adjust_bulk_4_zircon(Cliq_Zr, Sat_zr_liq)

    #         te_names    = TE_dtb.element_name
    #         ph_TE       = nothing
    #         Cliq        = te_X 
    #         Cmin        = nothing
    #         te_db       = te_db
    #         zr_sat      = zr_sat
    #     else 
    #         TE_dtb      =  get_TE_database("TE_OL_felsic")

    #         te_names    = TE_dtb.element_name
    #         ph_TE       = nothing
    #         Cliq        = nothing 
    #         Cmin        = nothing
    #         te_db       = te_db
    #         zr_sat      = zr_sat
    #         Sat_zr_liq  = nothing
    #         zircon_wt   = 0.0
    #     end

    #     out_te = tepm_struct{Float64}(  te_names, ph_TE, te_X, Cliq, Cmin,
    #                                     te_db, zr_sat, 
    #                                     Sat_zr_liq, zircon_wt)

    #     return out, out_te
    # else
