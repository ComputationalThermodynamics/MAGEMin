"""
Lattice strain models for trace element partitioning (Brice 1975).
Mineral compositions (wt% oxides) are converted to crystallographic site
fractions following Thermocalc a-x solution-model conventions defined in
igG23_H18w_MAGEMin_descriptions.txt (Green, Powell, Holland, Riel, Weller):
  cpx  → cpx_G23  (Green et al., in prep; after Holland et al. 2018)
  gt   → g_G23    (Green et al., in prep; after Holland et al. 2018)
  opx  → opx_G23  (Green et al., in prep; after Holland et al. 2018)
  ol   → ol_H18   (Holland et al. 2018)
  fsp  → fsp_H21  (Holland et al. 2021)
  amph → hb_G16   (Green et al. 2016)

Ported from TEPM v02.02 (J. Cornet, 2017).

Output element order (28 elements, same for every mineral):
  LILE 1+:   Cs  Rb  K
  LILE 2+:   Ba  Sr
  REE+Sc:    La Ce Pr Nd Sm Eu Gd Tb Dy Y Ho Er Tm Yb Lu Sc
  HFSE 4+:   Ti  Hf  Zr
  Actinides: U   Th
  HFSE 5+:   Ta  Nb
"""

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
const R_gas = 8.3144598        # J mol⁻¹ K⁻¹
# N_avo in GPa·Å³ units: 6.022140857e23 × 1e-21 (GPa→Pa) × (Å³→m³).
# This matches the MATLAB constant_input.m convention (N_avo = 602.2140857).
const fpe   = -4π * 602.2140857 / R_gas   # Brice prefactor ≈ -909.8

# ---------------------------------------------------------------------------
# Ionic radii (Å) — Shannon 1976
# ---------------------------------------------------------------------------

# 8-fold coordination
const ri8_REE_Sc = [1.160, 1.143, 1.126, 1.090, 1.079, 1.066, 1.053, 1.040,
                    1.027, 1.019, 1.015, 1.004, 0.994, 0.985, 0.977, 0.870]
                 #   La    Ce     Pr     Nd     Sm     Eu     Gd     Tb
                 #   Dy    Y      Ho     Er     Tm     Yb     Lu     Sc
const ri8_LILE1  = [1.74, 1.61, 1.51]          # Cs  Rb  K
const ri8_LILE2  = [1.42, 1.26]                 # Ba  Sr
const ri8_HFSE4  = [0.7399, 0.83, 0.84, 1.00, 1.05]  # Ti Hf Zr U Th
const ri8_HFSE5  = [0.7295, 0.74]              # Ta  Nb

# 6-fold coordination
const ri6_REE_Sc = [1.032, 1.010, 0.990, 0.983, 0.958, 0.947, 0.938, 0.923,
                    0.912, 0.900, 0.901, 0.890, 0.880, 0.868, 0.861, 0.745]
const ri6_LILE1  = [1.67, 1.52, 1.38]
const ri6_LILE2  = [1.35, 1.18]
const ri6_HFSE4  = [0.605, 0.71, 0.72, 0.89, 0.94]   # Ti Hf Zr U Th
const ri6_HFSE5  = [0.6295, 0.64]

# ---------------------------------------------------------------------------
# Molar masses (g mol⁻¹)
# Column order used throughout: SiO2 TiO2 Al2O3 FeO MnO MgO CaO Na2O K2O H2O
# ---------------------------------------------------------------------------
const MM = (SiO2=60.084, TiO2=79.979, Al2O3=101.961, FeO=71.846, MnO=70.937,
            MgO=40.304, CaO=56.077, Na2O=61.979, K2O=94.196, H2O=18.01528)

# ---------------------------------------------------------------------------
# Brice (1975) lattice-strain partition coefficient
# ---------------------------------------------------------------------------
@inline function D_brice(D0, E, r0, ri, T)
    D0 * exp(fpe * E * (r0/2 * (r0^2 - ri^2) - (r0^3 - ri^3)/3) / T)
end

D_brice_vec(D0, E, r0, ri_arr, T) =
    D0 .* exp.(fpe .* E .* (r0/2 .* (r0^2 .- ri_arr.^2) .- (r0^3 .- ri_arr.^3)/3) ./ T)

# "Mixed" Brice variant used in MATLAB TE_*.m for LILE and some HFSE:
# elastic radius r_e modulates the quadratic term; r_ref is the fixed reference
# cation whose D0 is the peak (the formula returns D0 exactly at ri = r_ref).
@inline function D_brice_mx(D0, E, r_e, r_ref, ri, T)
    D0 * exp(fpe * E * (r_e/2*(r_ref^2 - ri^2) - (r_ref^3 - ri^3)/3) / T)
end
D_brice_mx_vec(D0, E, r_e, r_ref, ri_arr, T) =
    D0 .* exp.(fpe .* E .* (r_e/2 .* (r_ref^2 .- ri_arr.^2) .- (r_ref^3 .- ri_arr.^3)/3) ./ T)

# ===========================================================================
# Section 1 — Structural formula and site fractions
#
# Each function accepts a NamedTuple of wt% oxides and returns site occupancies
# following the TC a-x conventions in igG23_H18w_MAGEMin_descriptions.txt.
#
# TC composition variables (lowercase) → site fractions (x prefix):
#   y  = 2·xAlT   (Al on tetrahedral site × 2, i.e. total Al_T count per f.u.)
#   x  = Fe/(Fe+Mg) on M sites (overall iron number)
#   c  = xCaM2 or xCaM1 depending on mineral
#   n  = xNaM2,  k = xKM2,  o = xFeM2+xMgM2  (for pyroxenes)
# ===========================================================================

"""
    cpx_sites(wt)

cpx_G23 — Green et al. (in prep), after Holland et al. (2018).
Sites: T*(Si Al), M1(Mg Fe Al Fe³⁺ Cr Ti), M2(Ca Na K Mg Fe).
Normalised to 6 oxygens. TC composition variables:
  y = 2·xAlT (= Al_T count per f.u. = XAl4)
  x = Fe/(Fe+Mg), o = xMgM2+xFeM2, n = xNaM2, k = xKM2
  Q = order variable → set to 0 (disordered approximation).
"""
function cpx_sites(wt)
    n_SiO2  = wt.SiO2  / MM.SiO2
    n_TiO2  = wt.TiO2  / MM.TiO2
    n_Al2O3 = wt.Al2O3 / MM.Al2O3
    n_FeO   = wt.FeO   / MM.FeO
    n_MnO   = wt.MnO   / MM.MnO
    n_MgO   = wt.MgO   / MM.MgO
    n_CaO   = wt.CaO   / MM.CaO
    n_Na2O  = wt.Na2O  / MM.Na2O
    n_K2O   = wt.K2O   / MM.K2O

    O_sum = 2n_SiO2 + 2n_TiO2 + 3n_Al2O3 + n_FeO + n_MnO + n_MgO +
            n_CaO + n_Na2O + n_K2O
    f6 = 6.0 / O_sum   # normalisation factor to 6 O per formula unit

    Si = f6 * n_SiO2
    Ti = f6 * n_TiO2
    Al = f6 * 2n_Al2O3
    Fe = f6 * n_FeO
    Mg = f6 * n_MgO
    Ca = f6 * n_CaO
    Na = f6 * 2n_Na2O
    K  = f6 * 2n_K2O

    # TC: y = 2·xAlT → total Al on T per formula unit
    y  = max(2.0 - Si, 0.0)   # Al_T = XAl4
    # TC: n = xNaM2, k = xKM2
    n  = Na
    k  = K
    # TC: o = xMgM2 + xFeM2 = M2 remainder after large cations
    o  = max(1.0 - Ca - n - k, 0.0)
    # TC: x = Fe/(Fe+Mg) overall (order Q = 0)
    x  = (Fe + Mg) > 0 ? Fe / (Fe + Mg) : 0.0

    # --- M2 site fractions (TC cpx_G23, Q=0) ---
    xCaM2 = Ca
    xNaM2 = n
    xKM2  = k
    xMgM2 = o * (1.0 - x)
    xFeM2 = o * x

    # --- T site fractions ---
    xAlT = y / 2.0   # fraction of T sites occupied by Al
    xSiT = 1.0 - xAlT

    # --- M1 site fractions: excess Al after T-site filling ---
    xAlM1 = max(Al - y, 0.0)
    xTiM1 = Ti
    xMgM1 = Mg - xMgM2
    xFeM1 = Fe - xFeM2

    XWo = Ca / (Ca + Mg + Fe)
    XEn = Mg / (Ca + Mg + Fe)
    XFs = Fe / (Ca + Mg + Fe)

    return (xSiT=xSiT, xAlT=xAlT, Al_T=y,
            xAlM1=xAlM1, xTiM1=xTiM1, xMgM1=xMgM1, xFeM1=xFeM1,
            xCaM2=xCaM2, xNaM2=xNaM2, xKM2=xKM2, xMgM2=xMgM2, xFeM2=xFeM2,
            XWo=XWo, XEn=XEn, XFs=XFs)
end

"""
    gt_sites(wt)

g_G23 — Green et al. (in prep), after Holland et al. (2018).
Sites: M1(Mg Fe Ca) dodecahedral, M2(Al Cr Fe³⁺ Mg Ti) octahedral.
Normalised to 12 oxygens. TC composition variables:
  c = xCaM1 (= XGr), x = xFeM1/(xFeM1+xMgM1) = Fe/(Fe+Mg).
Note: TC labels the dodecahedral site M1 and octahedral site M2;
xAlM2 ≈ 1 for typical igneous garnets.
"""
function gt_sites(wt)
    n_SiO2  = wt.SiO2  / MM.SiO2
    n_TiO2  = wt.TiO2  / MM.TiO2
    n_Al2O3 = wt.Al2O3 / MM.Al2O3
    n_FeO   = wt.FeO   / MM.FeO
    n_MnO   = wt.MnO   / MM.MnO
    n_MgO   = wt.MgO   / MM.MgO
    n_CaO   = wt.CaO   / MM.CaO

    O_sum = 2n_SiO2 + 2n_TiO2 + 3n_Al2O3 + n_FeO + n_MnO + n_MgO + n_CaO
    f12 = 12.0 / O_sum

    Si  = f12 * n_SiO2
    Ti  = f12 * n_TiO2
    Al  = f12 * 2n_Al2O3
    Fe  = f12 * n_FeO
    Mn  = f12 * n_MnO
    Mg  = f12 * n_MgO
    Ca  = f12 * n_CaO

    tot = Ca + Mg + Fe + Mn

    # TC: c = xCaM1, x = Fe/(Fe+Mg) within non-Ca fraction
    c = Ca / tot   # xCaM1 = XGr
    x = (Mg + Fe) > 0 ? Fe / (Mg + Fe) : 0.0

    # TC site fractions (M1 = dodecahedral X site):
    xCaM1 = c                    # XGr
    xMgM1 = (1.0 - c) * (1.0 - x)  # XPy (normalised)
    xFeM1 = (1.0 - c) * x          # XAlm (normalised)

    XGr  = xCaM1
    XPy  = xMgM1
    XAlm = xFeM1
    XSpss = Mn / tot

    # TC M2 (octahedral Y site): xAlM2 = 1 - cr - f - 2t
    # From structural formula: Al on Z site fills 3 positions (Si_Z + Al_Z = 3)
    Si_Z  = min(Si, 3.0)
    Al_Z  = 3.0 - Si_Z
    Al_Y  = max(Al - Al_Z, 0.0)   # Al remaining for Y (M2) site
    xAlM2 = Al_Y / 2.0            # fraction of 2 Y-site positions

    return (xCaM1=xCaM1, xMgM1=xMgM1, xFeM1=xFeM1, xAlM2=xAlM2, xTiM2=Ti/2,
            XGr=XGr, XPy=XPy, XAlm=XAlm, XSpss=XSpss)
end

"""
    opx_sites(wt)

opx_G23 — Green et al. (in prep), after Holland et al. (2018).
Sites: T*(Si Al), M1(Mg Fe Al Fe³⁺ Cr Ti), M2(Ca Na Mg Fe).
Normalised to 6 oxygens. TC composition variables:
  y = 2·xAlT, c = xCaM2, j = xNaM2, x = Fe/(Fe+Mg), Q=0.
"""
function opx_sites(wt)
    n_SiO2  = wt.SiO2  / MM.SiO2
    n_TiO2  = wt.TiO2  / MM.TiO2
    n_Al2O3 = wt.Al2O3 / MM.Al2O3
    n_FeO   = wt.FeO   / MM.FeO
    n_MnO   = wt.MnO   / MM.MnO
    n_MgO   = wt.MgO   / MM.MgO
    n_CaO   = wt.CaO   / MM.CaO
    n_Na2O  = wt.Na2O  / MM.Na2O

    O_sum = 2n_SiO2 + 2n_TiO2 + 3n_Al2O3 + n_FeO + n_MnO + n_MgO +
            n_CaO + n_Na2O
    f6 = 6.0 / O_sum

    Si = f6 * n_SiO2
    Ti = f6 * n_TiO2
    Al = f6 * 2n_Al2O3
    Fe = f6 * n_FeO
    Mg = f6 * n_MgO
    Ca = f6 * n_CaO
    Na = f6 * 2n_Na2O

    y  = max(2.0 - Si, 0.0)   # Al_T = XAl4
    c  = Ca                    # xCaM2
    j  = Na                    # xNaM2
    x  = (Fe + Mg) > 0 ? Fe / (Fe + Mg) : 0.0

    o     = max(1.0 - c - j, 0.0)
    xMgM2 = o * (1.0 - x)
    xFeM2 = o * x
    xMgM1 = Mg - xMgM2
    xFeM1 = Fe - xFeM2

    # M1 Al: excess Al after T-site filling (consistent with cpx_G23)
    xAlM1 = max(Al - y, 0.0)

    XWo = Ca / (Ca + Mg + Fe)
    XEn = Mg / (Ca + Mg + Fe)
    XFs = Fe / (Ca + Mg + Fe)

    return (xAlT=y/2, Al_T=y,
            xAlM1=xAlM1, xTiM1=Ti,
            xCaM2=Ca, xNaM2=Na, xMgM2=xMgM2, xFeM2=xFeM2,
            xMgM1=xMgM1, xFeM1=xFeM1,
            XWo=XWo, XEn=XEn, XFs=XFs)
end

"""
    fsp_sites(wt)

fsp_H21 — Holland et al. (2021), ternary feldspar.
Sites: A(Na Ca K), TB*(Si Al with 1/4 entropy contribution).
Normalised to 8 oxygens. TC variables: ca = xCaA = XAn, k = xKA = XOr.
  xAlTB = (1 + ca)/4,  xSiTB = (3 - ca)/4
"""
function fsp_sites(wt)
    n_SiO2  = wt.SiO2  / MM.SiO2
    n_Al2O3 = wt.Al2O3 / MM.Al2O3
    n_CaO   = wt.CaO   / MM.CaO
    n_Na2O  = wt.Na2O  / MM.Na2O
    n_K2O   = wt.K2O   / MM.K2O

    O_sum = 2n_SiO2 + 3n_Al2O3 + n_CaO + n_Na2O + n_K2O
    f8 = 8.0 / O_sum

    Ca = f8 * n_CaO
    Na = f8 * 2n_Na2O
    K  = f8 * 2n_K2O

    tot_A = Ca + Na + K
    XAn = Ca / tot_A
    XAb = Na / tot_A
    XOr = K  / tot_A

    # TC fsp_H21: xAlTB = 1/4 + 1/4*ca = (1+XAn)/4
    xAlTB = (1.0 + XAn) / 4.0
    xSiTB = (3.0 - XAn) / 4.0

    return (xCaA=XAn, xNaA=XAb, xKA=XOr,
            xAlTB=xAlTB, xSiTB=xSiTB,
            XAn=XAn, XAb=XAb, XOr=XOr)
end

"""
    ol_sites(wt)

ol_H18 — Holland et al. (2018).
Sites: M1(Mg Fe), M2(Mg Fe Ca).
Normalised to 4 oxygens. TC variables: x = Fe/(Fe+Mg), c = xCaM2, Q=0.
  xMgM1 = 1-x, xFeM1 = x, xMgM2 = 1-c-x, xFeM2 = x*(1-c), xCaM2 = c
"""
function ol_sites(wt)
    n_SiO2 = wt.SiO2 / MM.SiO2
    n_FeO  = wt.FeO  / MM.FeO
    n_MnO  = wt.MnO  / MM.MnO
    n_MgO  = wt.MgO  / MM.MgO
    n_CaO  = wt.CaO  / MM.CaO

    O_sum = 2n_SiO2 + n_FeO + n_MnO + n_MgO + n_CaO
    f4 = 4.0 / O_sum

    Fe = f4 * n_FeO
    Mn = f4 * n_MnO
    Mg = f4 * n_MgO
    Ca = f4 * n_CaO

    # TC ol_H18: x = Fe/(Fe+Mg), c = xCaM2
    x = (Fe + Mg) > 0 ? Fe / (Fe + Mg) : 0.0
    c = Ca   # xCaM2

    # Site fractions with Q=0 (disordered)
    xMgM1 = 1.0 - x
    xFeM1 = x
    xMgM2 = 1.0 - c - x
    xFeM2 = x * (1.0 - c)
    xCaM2 = c

    X_Fo = 1.0 - x   # = Mg/(Mg+Fe), ignoring Mn

    return (xMgM1=xMgM1, xFeM1=xFeM1,
            xMgM2=xMgM2, xFeM2=xFeM2, xCaM2=xCaM2,
            X_Fo=X_Fo, Mn=Mn)
end

"""
    amph_sites(wt)

hb_G16 — Green et al. (2016).
Sites: A(v Na K), M13(Mg Fe), M2(Mg Fe Al Fe³⁺ Ti), M4(Ca Mg Fe Na), T1*(Si Al).
Normalised to 23 oxygens (8 T1 positions per formula unit). TC variables:
  y = xAlM2, z = xNaM4, c = xCaM4, a = total A-site occupancy.
  xAlT1 ≈ y/2 + a/4  (simplified: ignoring small f, t, z contributions)
"""
function amph_sites(wt)
    n_SiO2  = wt.SiO2  / MM.SiO2
    n_TiO2  = wt.TiO2  / MM.TiO2
    n_Al2O3 = wt.Al2O3 / MM.Al2O3
    n_FeO   = wt.FeO   / MM.FeO
    n_MgO   = wt.MgO   / MM.MgO
    n_CaO   = wt.CaO   / MM.CaO
    n_Na2O  = wt.Na2O  / MM.Na2O
    n_K2O   = wt.K2O   / MM.K2O

    O_sum = 2n_SiO2 + 2n_TiO2 + 3n_Al2O3 + n_FeO + n_MgO +
            n_CaO + n_Na2O + n_K2O
    f23 = 23.0 / O_sum   # 23 O per formula unit (excluding OH)

    Si = f23 * n_SiO2
    Ti = f23 * n_TiO2
    Al = f23 * 2n_Al2O3
    Fe = f23 * n_FeO
    Mg = f23 * n_MgO
    Ca = f23 * n_CaO
    Na = f23 * 2n_Na2O
    K  = f23 * 2n_K2O

    # M4 site (2 positions): Ca preferentially, then Na
    xCaM4 = min(Ca, 2.0) / 2.0          # c = xCaM4
    xNaM4 = min(Na, 2.0 - 2*xCaM4) / 2.0  # z = xNaM4
    Na_A  = Na - 2*xNaM4               # Na on A site
    K_A   = K                            # K always on A site

    # A site occupancy: a = xNaA + xKA
    a_tot = min(Na_A + K_A, 1.0)
    xNaA  = Na_A / max(Na_A + K_A, 1e-12) * a_tot
    xKA   = K_A  / max(Na_A + K_A, 1e-12) * a_tot
    xvA   = 1.0 - a_tot

    # T1 site (8 positions): Si fills first, Al completes
    # TC hb_G16: xAlT1 = 1/2*f + 1/2*t + 1/2*y - 1/2*z + 1/4*a
    # Simplified (ignoring f, t ≈ 0):
    # First compute Al_T from structural formula
    Al_T8 = max(8.0 - Si, 0.0)          # total Al on T (8 positions)
    xAlT1 = Al_T8 / 8.0                  # fraction

    # y = xAlM2: Al remaining after T1 assignment
    Al_oct = max(Al - Al_T8, 0.0)
    xAlM2  = Al_oct / 2.0               # normalised to M2 site

    return (xCaM4=xCaM4, xNaM4=xNaM4, xNaA=xNaA, xKA=xKA, xvA=xvA,
            xAlT1=xAlT1, Al_T=Al_T8,
            xAlM2=xAlM2, xTiM2=Ti/2,
            Mg=Mg, Fe=Fe)
end

# ===========================================================================
# Section 2 — Partition coefficient models
#
# Each function returns a 28-element Vector{Float64} with D values in the
# fixed output order: Cs Rb K | Ba Sr | REE(La→Lu) Sc | Ti Hf Zr | U Th | Ta Nb
# ===========================================================================

"""
    D_cpx(T, P, melt, min_wt)

Clinopyroxene/melt partition coefficients after Blundy & Wood (1994),
Sun & Liang (2012), Corgne et al. (2012), and Hill et al. (2011).
Site fractions from cpx_G23 (Green et al., in prep).
"""
function D_cpx(T, P, melt, min_wt)
    s    = cpx_sites(min_wt)
    XAl4 = s.Al_T          # = y in TC notation
    XAl6 = s.xAlM1
    XMg_M2 = s.xMgM2
    XMg_M1 = s.xMgM1

    H2O_melt = melt.H2O / MM.H2O

    # Melt mole fractions for z0 (MATLAB TE_cpx.m: oxide moles, 4×Ti not 4×Si)
    n_Ti    = melt.TiO2 / MM.TiO2
    n_Fe    = melt.FeO  / MM.FeO
    n_Mg    = melt.MgO  / MM.MgO
    n_Ca    = melt.CaO  / MM.CaO
    n_Mn    = melt.MnO  / MM.MnO
    n_Na_ox = melt.Na2O / MM.Na2O
    n_K_ox  = melt.K2O  / MM.K2O
    z0 = (H2O_melt + n_K_ox + n_Na_ox + 2n_Ca + 2n_Mg + 2n_Mn + 2n_Fe + 4n_Ti) / 15.0

    DCa = min_wt.CaO / melt.CaO
    DNa = exp((10367 + 2100P - 165P^2) / T - 10.27 + 0.36P - 0.018P^2)

    # D_Ti (Hill et al. 2011, first model)
    # Uses per-formula-unit Si count (XSi = 2 - Al_T), not T-site fraction
    fpe_Ti = exp((46837 - 4257P - 9967P^2) / (R_gas*T))
    XNa = s.xNaM2;  XK = s.xKM2
    XSi = 2.0 - XAl4   # Si per formula unit (= MATLAB XSi = 2 - XAl4)
    Y = ((1-(XK+XNa))*XAl4^2 + 2*(XK+XNa)*XAl4*XSi) +
        ((XK+XNa)*XSi^2 + 2*(1-(XK+XNa))*XAl6*XSi*XAl4) *
        exp(-16500 / (R_gas*T)) +
        ((1-(XNa+XK))*XSi^2) * exp(-4*16500 / (R_gas*T))
    D_Ti = Y * fpe_Ti

    # REE 3+ on M2 — Sun & Liang (2012)
    ro_M2_3 = 1.066 - 0.104*XAl6 - 0.212*XMg_M2
    E_3     = (2.27*ro_M2_3 - 2.0) * 1e3
    Do_REE  = exp(-7.14 + 71900/(R_gas*T) + 4.37*XAl4 + 1.98*XMg_M2 - 0.91*H2O_melt)
    Di_REE  = D_brice_vec(Do_REE, E_3, ro_M2_3, ri8_REE_Sc[1:15], T)

    # Sc — Hill et al. (2011); uses per-f.u. XSi; -28.000 matches MATLAB
    fpe_Sc = exp((255.646 - 0.149T + 4.233P^2 - 2.280*(3-z0)^2) / (R_gas*T))
    L = ((XK+XNa)*XSi^2 + 2*(1-(XK+XNa))*XAl6*XSi*XAl4) +
        ((XAl4^2*(1-(XK+XNa))) + 2*(XK+XNa)*XAl4*XSi +
         (1-(XK+XNa))*XSi) * exp(-28.000 / (R_gas*T))
    D_Sc = L * fpe_Sc
    Di_REE_Sc = [Di_REE; D_Sc]

    # LILE 1+ on M2 — mixed Brice: elastic r=ro_1, fixed ref = Na (1.18 Å)
    ro_1 = ro_M2_3 + 0.20;  E_1 = E_3 / 3.0
    Di_LILE1 = D_brice_mx_vec(DNa, E_1, ro_1, 1.18, ri8_LILE1, T)

    # LILE 2+ on M2 — r_ref = Ca
    ro_2 = ro_M2_3 + 0.06;  E_2 = 2.0/3.0 * E_3
    Di_LILE2 = DCa .* exp.(fpe .* E_2 .*
        (ro_2/2 .* (1.12^2 .- ri8_LILE2.^2) .- (1.12^3 .- ri8_LILE2.^3)/3) ./ T)

    # HFSE 4+ on M1 — Corgne et al. (2012)
    # Mixed Brice: elastic r=ro_M1_4, fixed ref = Ti⁶ᶠ (0.605 Å)
    ro_M1_4 = 0.64 - 0.008P + 0.071*XAl6 + 0.004
    E_4     = 10473 - 5.09T - 201.54P + 14633*XAl4 + 331
    D_Hf    = D_brice_mx(D_Ti, E_4, ro_M1_4, ri6_HFSE4[1], ri6_HFSE4[2], T)
    D_Zr    = D_brice_mx(D_Ti, E_4, ro_M1_4, ri6_HFSE4[1], ri6_HFSE4[3], T)
    Di_HFSE4 = [D_Ti, D_Hf, D_Zr]

    # HFSE 5+
    D_Ta = exp(-2.127 + 3.769*XAl4)
    D_Nb = 0.003 + 0.292*D_Ta
    Di_HFSE5 = [D_Ta, D_Nb]

    # Actinides on M2 — Blundy & Wood (2003)
    E_4_M2   = 4.0/3.0 * E_3
    gamma_Mg = exp(7500 / (R_gas*T))
    gamma_Th = exp(fpe * E_4_M2 *
        (ro_M2_3/2*(ri8_HFSE4[5]-ro_M2_3)^2 + (ri8_HFSE4[5]-ro_M2_3)^3/3) / T)
    D_Th = n_Mg / (gamma_Mg * gamma_Th * XMg_M1) *
           exp((214.79 - 0.757T + 16.42P - 1.5P^2) / (R_gas*T))
    D_U  = D_Th * exp(-910.17 * E_4_M2 *
           (ro_M2_3/2*(ri8_HFSE4[5]^2-ri8_HFSE4[4]^2) -
            (ri8_HFSE4[5]^3-ri8_HFSE4[4]^3)/3) / T)
    Di_Act = [D_U, D_Th]

    return vcat(Di_LILE1, Di_LILE2, Di_REE_Sc, Di_HFSE4, Di_Act, Di_HFSE5)
end

"""
    D_gt(T, P, melt, min_wt)

Garnet/melt partition coefficients after van Westrenen & Draper (2007),
and Sun & Liang (2013).
Site fractions from g_G23 (Green et al., in prep).
"""
function D_gt(T, P, melt, min_wt)
    s   = gt_sites(min_wt)
    XCa = s.XGr
    XMg = s.XPy
    Melt_poly = 0.008458 * melt.SiO2 + 0.1734

    n_Fe = melt.FeO / MM.FeO
    n_Mg = melt.MgO / MM.MgO
    n_Si = melt.SiO2 / MM.SiO2

    # D_Mg — van Westrenen & Draper (2007); 19.000 matches MATLAB constant
    DMg = exp((25.820 - 0.1415T + 5.418P) / (3*R_gas*T)) /
          exp((19.000 * XCa^2) / (R_gas*T))

    # REE 3+ on M1 (dodecahedral)
    ro_3 = 0.780 + 0.155*XCa
    E_3  = (-1.62 + 2.29*ro_3) * 1e3
    Do   = exp(-2.75 + (91700 - 91.35P*(38-P)) / (R_gas*T) - 1.42*XCa)
    Di_REE    = D_brice_vec(Do, E_3, ro_3, ri8_REE_Sc[1:15], T)
    D_Sc      = 5.79   # Adam & Green 2006
    Di_REE_Sc = [Di_REE; D_Sc]

    # LILE 1+
    D_Cs = 0.001
    D_Rb = exp(0.47029 - 4.644*XMg + 0.656*XMg^2)
    D_K  = exp(4.72 - 1.4P - 50*XCa + 8.09P*XCa)
    Di_LILE1 = [D_Cs, D_Rb, D_K]

    # LILE 2+ — mixed Brice: elastic r=ro_2, fixed ref = Mg²⁺ (0.89 Å)
    ro_2 = ro_3 + 0.12;  E_2 = 2.0/3.0 * E_3
    Di_LILE2 = D_brice_mx_vec(DMg, E_2, ro_2, 0.89, ri8_LILE2, T)

    # Actinides — Salters & Longhi (1999) / Elkins et al. (2008)
    # Computed before HFSE4 because D_Th is used in the XCa > 0.19 HFSE branch
    D_Th = (n_Fe + n_Mg)^4 * n_Si^2 *
           ((11.46 - 24200/T + 8.6*(1-(n_Fe+n_Mg))^2 - 2.08*(1-XCa)^2)^2)
    D_U  = 3.6*D_Th + 0.003
    Di_Act = [D_U, D_Th]

    # HFSE 4+ — mixed Brice; X-site uses Ti⁸ᶠ (0.7399 Å), Y-site uses Ti⁶ᶠ (0.605 Å)
    D_Ti  = 0.0037 * exp(7.386 * Melt_poly)
    ro_X4 = ro_3 - 0.1;  ro_Y4 = 0.67
    E_X4  = 1325.0
    E_Y4  = 15870 * exp(-0.1086*(XCa*100)) + 409.6 * exp(0.01191*(XCa*100))

    if XCa > 0.19 && XCa < 0.4
        # X-site (8-fold) uses D_Th as D0 and ri8 Ti reference
        D_Hf = D_brice_mx(D_Th, E_X4, ro_X4, ri8_HFSE4[1], ri8_HFSE4[2], T) +
               D_brice_mx(D_Ti,  E_Y4, ro_Y4, ri6_HFSE4[1], ri6_HFSE4[2], T)
        D_Zr = D_brice_mx(D_Th, E_X4, ro_X4, ri8_HFSE4[1], ri8_HFSE4[3], T) +
               D_brice_mx(D_Ti,  E_Y4, ro_Y4, ri6_HFSE4[1], ri6_HFSE4[3], T)
    else
        D_Hf = D_brice_mx(D_Ti, E_Y4, ro_Y4, ri6_HFSE4[1], ri6_HFSE4[2], T)
        D_Zr = D_brice_mx(D_Ti, E_Y4, ro_Y4, ri6_HFSE4[1], ri6_HFSE4[3], T)
    end
    Di_HFSE4 = [D_Ti, D_Hf, D_Zr]

    # HFSE 5+
    Di_HFSE5 = [0.0146, 0.0290]   # Ta, Nb — fixed values (Adam & Green 2006)

    return vcat(Di_LILE1, Di_LILE2, Di_REE_Sc, Di_HFSE4, Di_Act, Di_HFSE5)
end

"""
    D_opx(T, P, melt, min_wt)

Orthopyroxene/melt partition coefficients after Bédard (2007),
Frei et al. (2009), and Wood & Blundy (2013).
Site fractions from opx_G23 (Green et al., in prep).
"""
function D_opx(T, P, melt, min_wt)
    s       = opx_sites(min_wt)
    XAl4    = s.Al_T
    XCa     = s.xCaM2
    XMg_M2  = s.xMgM2
    XMg_M1  = s.xMgM1
    XFe_M1  = s.xFeM1

    Mg_n_opx = min_wt.MgO / (min_wt.MgO + min_wt.FeO)
    Mg_n     = melt.MgO   / (melt.MgO   + melt.FeO)

    DMg = 1.0   # Wood & Blundy (2013): DMg ≈ 1 between basalt and opx

    # REE 3+ on M2
    ro_3 = 0.69 + 0.43*XCa + 0.23*XMg_M2
    E_3  = (-1.37 + 1.85*ro_3 - 0.53*XCa) * 1e3
    Do   = exp(-5.37 + 38700/(R_gas*T) + 3.56*XCa + 3.54*XAl4)
    Di_REE    = D_brice_vec(Do, E_3, ro_3, ri8_REE_Sc[1:15], T)
    D_Sc      = exp(5.04 - 4.37*Mg_n_opx - 0.41*log(min_wt.FeO) -
                    0.318*log(melt.MgO))
    Di_REE_Sc = [Di_REE; D_Sc]

    # LILE 1+
    D_Cs = exp(23.6 - 1.46P - 0.26*melt.SiO2 - 32.49*Mg_n)
    if D_Cs > 1; D_Cs = exp(2.901 - 17.32*Mg_n); end
    D_Rb = exp(6.995148 - 0.1015*melt.SiO2 - 0.0571*melt.Al2O3 +
               0.4567*log(melt.MgO) - 16.0222*Mg_n)
    D_K  = exp(-13.04737 + 0.9458P + 0.4307*s.XWo + 0.13*melt.SiO2 -
               0.9374*log(melt.FeO) + 1.3376*log(melt.MgO))
    Di_LILE1 = [D_Cs, D_Rb, D_K]

    # LILE 2+ on M2 — mixed Brice: elastic r=ro_2, fixed ref = Mg²⁺ (0.89 Å)
    ro_2 = ro_3 + 0.12;  E_2 = 2.0/3.0 * E_3
    Di_LILE2 = D_brice_mx_vec(DMg, E_2, ro_2, 0.89, ri8_LILE2, T)

    # HFSE 4+ on M1 — Frei et al. (2009)
    ro_4 = 0.618 + 0.032*XCa + 0.030*XMg_M1
    E_4  = 2203.0
    Do_4 = exp(-4.825 + 31780/(R_gas*T) + 4.17*XAl4 + 8.55*XCa*XMg_M2 - 2.62*XFe_M1)
    Di_HFSE4 = D_brice_vec(Do_4, E_4, ro_4, ri6_HFSE4[1:3], T)

    # HFSE 5+
    D_Ta = exp(-1.68 + 0.62*XAl4)
    D_Nb = D_Ta * (1.17 - 3.16*XAl4)
    if D_Nb < 0; D_Nb = D_Ta * (1.17 - 3.16 * 0.369); end
    Di_HFSE5 = [D_Ta, D_Nb]

    # Actinides
    D_U  = exp(3.967 - 11.668*Mg_n_opx)
    D_Th = exp(-6.713332 + 19.48971*XAl4 + 2.11*log(melt.FeO) - 2.55*log(melt.MgO))
    Di_Act = [D_U, D_Th]

    return vcat(Di_LILE1, Di_LILE2, Di_REE_Sc, Di_HFSE4, Di_Act, Di_HFSE5)
end

"""
    D_pl(T, P, melt, min_wt)

Plagioclase/melt partition coefficients after Dohmen & Blundy (2014).
Site fractions from fsp_H21 (Holland et al. 2021).
"""
function D_pl(T, _, melt, min_wt)
    s   = fsp_sites(min_wt)
    XAn = s.XAn

    DCa = exp((-19107 + 467*melt.SiO2) / (R_gas*T))
    DNa = exp((13621 - 18990*XAn)      / (R_gas*T))

    # REE 3+ — Dohmen & Blundy (2014)
    # Mixed Brice: elastic r=ro_3, fixed ref = La (ri8_REE_Sc[1] = 1.160 Å)
    ro_3 = 1.331 - 0.068*XAn
    E_3  = 152 - 31*XAn
    DLa  = (DCa^2 / DNa) * exp(4400/(R_gas*T) - 30.8/R_gas)
    Di_REE    = D_brice_mx_vec(DLa, E_3, ro_3, ri8_REE_Sc[1], ri8_REE_Sc[1:15], T)
    D_Sc      = exp((-23000 - 7550*melt.MgO) / (R_gas*T))
    Di_REE_Sc = [Di_REE; D_Sc]

    # LILE 1+ on A site — mixed Brice: elastic r=ro_1, fixed ref = Na (1.18 Å)
    ro_1 = 1.24 - 0.017*XAn
    E_1  = 49.05 + 17.16*XAn
    Di_LILE1 = D_brice_mx_vec(DNa, E_1, ro_1, 1.18, ri8_LILE1, T)

    # LILE 2+ on A site — Dohmen & Blundy (2014) r0 uses fsp_H21 T-correction
    Tr   = 1565.0
    ro_2 = 1.2895 + 1.3e-4*(T-Tr) + XAn*(-0.0952 - 4e-5*(T-Tr))
    E_2  = 120.03824 - 0.3686*XAn
    Di_LILE2 = DCa .* exp.(fpe .* E_2 .*
        (ro_2/2 .* (1.12^2 .- ri8_LILE2.^2) .- (1.12^3 .- ri8_LILE2.^3)/3) ./ T)

    # HFSE 4+
    D_Ti = exp((-15400 - 28900*XAn) / (R_gas*T))
    D_Zr = exp((-187270 + 2260*melt.SiO2) / (R_gas*T))
    D_Hf = exp((-126860 + 1399*melt.SiO2) / (R_gas*T))
    Di_HFSE4 = [D_Ti, D_Hf, D_Zr]

    # HFSE 5+
    D_Ta = exp((-271600 + 3157*melt.SiO2 + 8926*melt.MgO + 26501*XAn) / (R_gas*T))
    D_Nb = exp((-417350 + 5300*melt.SiO2 + 10340*melt.MgO + 40880*XAn) / (R_gas*T))
    Di_HFSE5 = [D_Ta, D_Nb]

    # Actinides
    D_U  = exp((-7450  - 48400*XAn) / (R_gas*T))
    D_Th = exp((-11646 - 31369*XAn) / (R_gas*T))
    Di_Act = [D_U, D_Th]

    return vcat(Di_LILE1, Di_LILE2, Di_REE_Sc, Di_HFSE4, Di_Act, Di_HFSE5)
end

"""
    D_ol(T, P, melt, min_wt)

Olivine/melt partition coefficients after Yao et al. (2012) and Bédard (2005).
Site fractions from ol_H18 (Holland et al. 2018).
"""
function D_ol(T, P, melt, min_wt)
    s    = ol_sites(min_wt)
    X_Fo = s.X_Fo

    DMg = min_wt.MgO / melt.MgO

    # REE 3+ on M site (8-fold, Yao et al. 2012)
    # MATLAB uses mc(1,3) = Al2O3_wt/MM_Al2O3 (moles), not wt%
    ro_3 = 0.72;  E_3 = 426.0
    Do   = exp(-0.45 - 0.11P + 1.54*(melt.Al2O3/MM.Al2O3) - 0.0194*X_Fo)
    Di_REE_Sc = D_brice_vec(Do, E_3, ro_3, ri8_REE_Sc, T)

    # LILE 1+
    if melt.MgO > 1.11
        val = exp(0.055 - 3.05*log(melt.MgO))
        Di_LILE1 = [val, val, val]
    else
        Di_LILE1 = [0.035, 0.035, 0.035]
    end

    # LILE 2+ on M site
    ro_2 = 0.89;  E_2 = 240.0
    Di_LILE2 = D_brice_vec(DMg, E_2, ro_2, ri8_LILE2, T)

    # HFSE 4+ — Bédard (2005); Mg_n is melt wt% Mg# × 100 (MATLAB TE_ol.m convention)
    Mg_n_melt = 100.0 * melt.MgO / (melt.MgO + melt.FeO)
    if melt.SiO2 > 57.1
        D_Ti = exp(-9.72 + 0.11*melt.SiO2)
    elseif melt.MgO < 0.493
        D_Ti = 0.29 - 0.52*melt.MgO
    elseif Mg_n_melt < 8.35
        D_Ti = 0.32 - 0.034*Mg_n_melt
    else
        D_Ti = 0.0302
    end
    D_Hf = melt.MgO > 5.2 ? exp(-2.98 - 0.76*log(melt.MgO)) : 0.04
    D_Zr = melt.MgO > 5.2 ? exp(-0.25 - 1.88*log(melt.MgO)) : 0.036
    Di_HFSE4 = [D_Ti, D_Hf, D_Zr]

    # HFSE 5+
    D_Ta = melt.MgO > 3.79 ? exp(-1.5*log(melt.MgO))         : 0.126
    D_Nb = melt.MgO > 4.4  ? exp(0.28 - 3.29*log(melt.MgO))  : 0.01025
    Di_HFSE5 = [D_Ta, D_Nb]

    # Actinides
    D_U  = melt.MgO > 5.91 ? exp(7.08 - 5.67*log(melt.MgO)) : 0.048
    D_Th = melt.MgO > 4.9  ? exp(3.8  - 4.22*log(melt.MgO)) : 0.0542
    Di_Act = [D_U, D_Th]

    return vcat(Di_LILE1, Di_LILE2, Di_REE_Sc, Di_HFSE4, Di_Act, Di_HFSE5)
end

"""
    D_amph(T, P, melt, min_wt)

Amphibole/melt partition coefficients after Tiepolo et al. (2007) and
Dalpe & Baker (2000).
Site fractions from hb_G16 (Green et al. 2016).
"""
function D_amph(T, _, melt, min_wt)
    s         = amph_sites(min_wt)
    Melt_poly = 0.008458*melt.SiO2 + 0.1734
    DCa       = min_wt.CaO / melt.CaO

    # REE 3+ — empirical regression to melt SiO2 (Tiepolo et al. 2007)
    a_REE = [0.0135, 0.0194, 0.0266, 0.0388, 0.0417, 0.0445, 0.0466, 0.0482,
             0.0492, 0.0495, 0.0496, 0.0495, 0.0491, 0.0485, 0.0478]
    b_REE = [0.0478, 0.0578, 0.0596, 0.0607, 0.0607, 0.0607, 0.0606, 0.0605,
             0.0603, 0.0601, 0.0600, 0.0597, 0.0593, 0.0589, 0.0585]
    Di_REE    = a_REE .* exp.(b_REE .* melt.SiO2)
    D_Sc      = 0.0169 * exp(9.33*Melt_poly)
    Di_REE_Sc = [Di_REE; D_Sc]

    # LILE 1+ on A site — lattice strain, Dalpe & Baker (2000)
    # hb_G16: A site (12-fold), ro from TC mean values
    r_CN12  = [1.88, 1.72, 1.64]   # Cs Rb K in 12-fold coordination
    ro_1    = 1.52;  E_1 = 86.75;  Do_1 = 2.99
    Di_LILE1 = D_brice_vec(Do_1, E_1, ro_1, r_CN12, T)

    # LILE 2+ on A site — Ba, Sr
    r_Ba_12 = 1.61
    D_Ba    = D_brice(2.6, 315.0, 1.52, r_Ba_12, T)
    D_Sr    = exp(0.251*log(DCa) - 0.998)
    Di_LILE2 = [D_Ba, D_Sr]

    # HFSE 4+ on M2/T sites — Tiepolo et al. (2007)
    if melt.SiO2 < 65
        D_Ti = 0.0273 * exp(6.9531*Melt_poly)
        D_Hf = 0.0035 * exp(8.894 *Melt_poly)
        D_Zr = 0.001493 * exp(9.352*Melt_poly)
    else
        E_4  = 1393.2;  ro_4 = 0.65;  Do_4 = 1.27
        D_Ti = D_brice(Do_4, E_4, ro_4, ri6_HFSE4[1], T)
        D_Hf = D_brice(Do_4, E_4, ro_4, ri6_HFSE4[2], T)
        D_Zr = D_brice(Do_4, E_4, ro_4, ri6_HFSE4[3], T)
    end
    Di_HFSE4 = [D_Ti, D_Hf, D_Zr]

    # HFSE 5+ on M2/T sites
    if melt.SiO2 < 56
        D_Nb = 6.5e-6 * exp(18.46*Melt_poly)
        D_Ta = D_Nb   * (3.48e-5 * exp(16.51*Melt_poly))
    else
        E_5  = 1985.0;  ro_5 = 0.73;  Do_5 = 2.10
        D_Nb = D_brice(Do_5, E_5, ro_5, ri6_HFSE5[2], T)
        D_Ta = D_brice(Do_5, E_5, ro_5, ri6_HFSE5[1], T)
    end
    Di_HFSE5 = [D_Ta, D_Nb]

    # Actinides
    D_U  = 9e-6 * exp(12.15*Melt_poly)
    D_Th = 1e-5 * exp(12.24*Melt_poly)
    Di_Act = [D_U, D_Th]

    # use s to suppress unused-binding warning (xAlT1 documents Al-T saturation)
    _ = s.xAlT1

    return vcat(Di_LILE1, Di_LILE2, Di_REE_Sc, Di_HFSE4, Di_Act, Di_HFSE5)
end

# ---------------------------------------------------------------------------
# Trace element names in output order
# ---------------------------------------------------------------------------
const TE_names = ["Cs", "Rb", "K",
                  "Ba", "Sr",
                  "La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb",
                  "Dy", "Y",  "Ho", "Er", "Tm", "Yb", "Lu", "Sc",
                  "Ti", "Hf", "Zr",
                  "U",  "Th",
                  "Ta", "Nb"]

# ---------------------------------------------------------------------------
# Composition constructor
# ---------------------------------------------------------------------------

"""
    make_composition(names, wt; Fe2O3=0.0) -> NamedTuple

Build the oxide wt% NamedTuple expected by every `D_*` and `*_sites`
function from parallel vectors of oxide names and wt% values.

Names recognised (case-sensitive):
  "SiO2"  "TiO2"  "Al2O3"  "FeO"  "MnO"  "MgO"  "CaO"  "Na2O"  "K2O"  "H2O"

Any oxide not listed defaults to 0.0.

If your analysis reports total iron as Fe₂O₃, pass `Fe2O3=<wt%>` and it is
converted to FeO-equivalent (×0.8998) and added to any existing "FeO" entry.

# Example
```julia
oxides = ["SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O"]
cpx  = make_composition(oxides, [50.45, 1.45, 6.30, 9.38, 13.53, 18.10, 0.70])
melt = make_composition(
    ["SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O","H2O"],
    [67.67,  1.27,  15.30,  2.79,  1.39,  4.17,  3.13,  1.94, 13.8])

D = D_cpx(1273.0, 1.0, melt, cpx)
```
"""
function make_composition(names::Vector{String}, wt::Vector{Float64}; Fe2O3::Real=0.0)
    g(k) = (i = findfirst(==(k), names); i === nothing ? 0.0 : wt[i])
    feo_total = g("FeO") + Fe2O3 * 0.8998   # Fe2O3 → FeO equiv: 2×71.846/159.692
    (SiO2  = g("SiO2"),
     TiO2  = g("TiO2"),
     Al2O3 = g("Al2O3"),
     FeO   = feo_total,
     MnO   = g("MnO"),
     MgO   = g("MgO"),
     CaO   = g("CaO"),
     Na2O  = g("Na2O"),
     K2O   = g("K2O"),
     H2O   = g("H2O"))
end
