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


"""
    zirconium_saturation(out; model="WH")

    Compute the zirconium saturation concentration in the melt phase [ppm].

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output (must contain a melt phase).
    model : String, optional
        Saturation model (default: "WH"). Valid options:
        - "WH"  — Watson & Harrison (1983)
        - "B"   — Boehnke et al. (2013)
        - "CB"  — Crisp & Berry (2022)

    Returns
    -------
    C_zr_liq : Float64
        Zr saturation concentration in the melt [ppm], or -1 if no melt is present.
"""
function zirconium_saturation(  out     :: MAGEMin_C.gmin_struct{Float64, Int64};
                                model   :: String = "WH"    )

    if out.frac_M > 0.0
        if model == "WH" || model == "B"
            ref_ox      = ("SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O", "Na2O", "TiO2", "O", "Cr2O3", "MnO", "H2O", "S")
            n_cation    = (1.0, 2, 1, 1, 1, 2, 2, 1, 1, 2, 1, 2, 1)
            cation_name = ["Si", "Al", "Ca", "Mg", "Fe", "K", "Na", "Ti", "O", "Cr", "Mn", "H", "S"]

            # Calculate cation_idx using a comprehension with the tuple ref_ox
            cation_idx  = [findfirst(isequal(x), ref_ox) for x in out.oxides]
            bulk_M_dry  = anhydrous_renormalization(out.bulk_M,out.oxides)

            cation      = zeros(length(cation_idx))
            # # Perform in-place multiplication to fill the cation array
            @inbounds for i in axes(cation, 1)
                cation[i] = bulk_M_dry[i] * n_cation[cation_idx[i]]
            end

            sum_cation  = sum(cation);

            Na          = findfirst(==("Na"), cation_name[cation_idx]);
            K           = findfirst(==("K"), cation_name[cation_idx]);
            Ca          = findfirst(==("Ca"), cation_name[cation_idx]);
            Al          = findfirst(==("Al"), cation_name[cation_idx]);
            Si          = findfirst(==("Si"), cation_name[cation_idx]);

            M           = (cation[Na] + cation[K]+ 2.0*cation[Ca]) / (cation[Al] * cation[Si]) * sum_cation;
            C_zr_zrc    = 497644;

            # Watson & Harrison (1983)
            if model == "WH"
                C_zr_liq    = C_zr_zrc/exp((-3.80 - 0.85*(M - 1.0)) + (12900/(out.T_C+273.15)))
            end
            # Boehnke et al., 2013
            if model == "B"
                C_zr_liq    = C_zr_zrc/exp((-1.48 - 1.16*(M - 1.0)) + (10108/(out.T_C+273.15)))
            end
        # Crisp and Berry 2022
        elseif model == "CB"
            opt_basi_oxides = ("SiO2","TiO2","Al2O3","MgO","MnO","FeO","Fe2O3","CaO","Na2O","K2O","P2O5","H2O")
            optical_basicity= (0.48,0.75,0.60,0.78,0.96,1.00,0.77,1.00,1.15,1.40,0.33,0.40)
            n_oxygen        = (2.0,2,3,1,1,1,3,1,1,1,5,1)

            commonOxide   =  intersect(Set(out.oxides),Set(opt_basi_oxides))
            idOx_MM       = [findfirst(isequal(x), out.oxides) for x in commonOxide]
            idOx_OB       = [findfirst(isequal(x), opt_basi_oxides) for x in commonOxide]

            liqCompNorm   = out.bulk_M[idOx_MM] ./ sum(out.bulk_M[idOx_MM])
            oxListDry     = findall(commonOxide .!= "H2O")

            # Compute optical basicity directly without intermediate arrays
            opt_basicity_numerator = 0.0
            opt_basicity_denominator = 0.0

            @inbounds for i in oxListDry
                idx = idOx_OB[i]
                ob = optical_basicity[idx]
                no = n_oxygen[idx]
                comp = liqCompNorm[i]
                opt_basicity_numerator += comp * ob * no
                opt_basicity_denominator += comp * no
            end
            opt_basicity = opt_basicity_numerator / opt_basicity_denominator

            # # Directly access the water content without creating a temporary array
            xH2O_index = findfirst(==("H2O"), String.(commonOxide))
            if !isnothing(xH2O_index)
                xH2O = liqCompNorm[xH2O_index]
            else
                xH2O = 0.0
            end

            # # Compute C_zr_liq
            C_zr_liq = 10^(0.96 - 5790.0 / (out.T_C + 273.15) - 1.28 * (out.P_kbar * 0.1) + 12.39 * opt_basicity + 0.83 * xH2O + 2.06 * (out.P_kbar * 0.1) * opt_basicity)
        else
            print("Model $model for zirconium saturation is invalid\n")
        end
    else
        print("Cannot compute zirconium saturation in liquid if melt is not predicted!\n")
        C_zr_liq = -1
    end

    return C_zr_liq
end

"""
    adjust_bulk_4_zircon(zr_liq, sat_liq, liq_wt)

    Compute the weight fractions of zircon, SiO₂, and O returned to the bulk when the melt exceeds Zr saturation.

    Parameters
    ----------
    zr_liq : Float64
        Zr concentration in the melt [ppm].
    sat_liq : Float64
        Zr saturation concentration in the melt [ppm].
    liq_wt : Float64
        Melt weight fraction.

    Returns
    -------
    zircon_wt : Float64
        Weight fraction of precipitated zircon (scaled by `liq_wt`).
    SiO2_wt : Float64
        SiO₂ weight returned to the bulk (scaled by `liq_wt`).
    O2_wt : Float64
        O weight returned to the bulk (scaled by `liq_wt`).
"""
function adjust_bulk_4_zircon(  zr_liq  ::  Float64,
                                sat_liq ::  Float64,
                                liq_wt  ::  Float64 )

    SiO2_wt         = 0.0
    O2_wt           = 0.0
    zircon_wt       = 0.0
    zircon_excess   = (zr_liq - sat_liq)/1e6

    zircon_wt   = zircon_excess*0.497644
    SiO2_wt     = (zircon_wt *0.327765)
    O2_wt       = (zircon_wt *0.174570)


    return zircon_wt*liq_wt, SiO2_wt*liq_wt, O2_wt*liq_wt
end

"""
    phosphate_saturation(out; model="Klein26")

    Compute the P₂O₅ saturation concentration in the melt phase [ppm].

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output (must contain a melt phase with a "liq" solution phase).
    model : String, optional
        Saturation model (default: "Klein26"). Valid options:
        - "Klein26"  — Klein et al. (2026)
        - "HWBea92"  — Harrison & Watson (1984) + correction from Bea et al. (1992); for hydrous systems
        - "Tollari06" — Tollari et al. (2006); calibrated for dry systems

    Returns
    -------
    C_P2O5_liq : Float64
        P₂O₅ saturation concentration in the melt [ppm], or -1 if the model is unrecognized.
"""
function phosphate_saturation(      out     :: MAGEMin_C.gmin_struct{Float64, Int64};
                                    model   :: String = "Klein26"    )

    if model == "Tollari06"
        ox_list         = ["SiO2","TiO2","Al2O3","FeO","MgO","CaO","Na2O","K2O"]
        T_K             = out.T_C + 273.15

        M_liq_dry       = [out.SS_vec[out.SS_syms[:liq]].Comp[i] * get_molar_mass(out.oxides[i]) for i in eachindex(out.oxides) if out.oxides[i] in ox_list]
        bulk_dry        = [out.SS_vec[out.SS_syms[:liq]].Comp[i] for i in eachindex(out.oxides) if out.oxides[i] in ox_list]
        ox_dry          = [out.oxides[i] for i in eachindex(out.oxides) if out.oxides[i] in ox_list]
        bulk_dry_norm   = bulk_dry./sum(bulk_dry)

        MCaO            = bulk_dry_norm[findfirst(ox_dry .== "CaO" )] * 100.0
        MSiO2           = bulk_dry_norm[findfirst(ox_dry .== "SiO2")] * 100.0

        # C_P_liq_molar = exp( T_K * ( (-0.8579)/(139.0 - MSiO2) + 0.0165)  - 3.3333 * log(MCaO) )  # abstract equation
        C_P_liq_molar   = exp( 2/3 * ( T_K*((-1.2868)/(139.0 - MSiO2) + 0.0247 ) - 5.0*log(MCaO)) )   # equation p1532
        M               = get_molar_mass("P2O5") * C_P_liq_molar/100    # g/mol for P2O5

        M_factor        = M/(sum(M_liq_dry) + M)                        # molar fraction of P2O5 in melt

        C_P2O5_liq      = C_P_liq_molar/M_factor * 1e4                             # ppm   

    elseif model == "HWBea92"

        SiO2_wt     = out.SS_vec[out.SS_syms[:liq]].Comp_wt[findfirst(isequal("SiO2"), out.oxides)]
        Al2O3_mol   = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("Al2O3"), out.oxides)]
        CaO_mol     = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("CaO"), out.oxides)]
        Na2O_mol    = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("Na2O"), out.oxides)]
        K2O_mol     = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("K2O"), out.oxides)]
        T_K         = out.T_C + 273.15

        AoCNK       = Al2O3_mol / (2.0*CaO_mol + Na2O_mol + K2O_mol)    # A over CNK

        C_P2O5_WH   = 52.5525567/ exp( (8400.0 + 2.64e4(SiO2_wt - 0.5))/T_K - (3.1 + 12.4(SiO2_wt - 0.5)) )  # wt% as in Harrison and Watson 1984
        C_P2O5_liq  =  maximum([C_P2O5_WH * 1e4, (C_P2O5_WH * (AoCNK -1.0) * 6429.0/out.T_C) * 1e4])   # ppm Harrison and Watson 1984 + Correction from Bea et al. 1992
    elseif model == "Klein26"
        SiO2_wt     = out.SS_vec[out.SS_syms[:liq]].Comp_wt[findfirst(isequal("SiO2"), out.oxides)]
        Al2O3_mol   = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("Al2O3"), out.oxides)]
        CaO_mol     = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("CaO"), out.oxides)]
        Na2O_mol    = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("Na2O"), out.oxides)]
        K2O_mol     = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("K2O"), out.oxides)]
        T_K         = out.T_C + 273.15

        ASI         = Al2O3_mol / (2.0*CaO_mol + Na2O_mol + K2O_mol)    # A over CNK
        log10_P2O5  = (1e4 * (-1.47 + 1.28 * SiO2_wt))/(T_K) + 12.79 + 1.06*ASI - 14.06 * SiO2_wt
        C_P2O5_liq  = exp10(log10_P2O5) * 1e4
    else
        print("Cannot compute phosphate saturation in liquid if melt is not predicted!\n")
        C_P2O5_liq = -1
    end
    
    return C_P2O5_liq
end


"""
    adjust_bulk_4_fapatite(P2O5_liq, sat_liq, liq_wt)

    Compute the weight fractions of fluorapatite and CaO returned to the bulk when the melt exceeds P₂O₅ saturation.

    Parameters
    ----------
    P2O5_liq : Float64
        P₂O₅ concentration in the melt [ppm].
    sat_liq : Float64
        P₂O₅ saturation concentration in the melt [ppm].
    liq_wt : Float64
        Melt weight fraction.

    Returns
    -------
    fapt_wt : Float64
        Weight fraction of precipitated fluorapatite (scaled by `liq_wt`).
    CaO_wt : Float64
        CaO weight returned to the bulk (scaled by `liq_wt`).
"""
function adjust_bulk_4_fapatite(    P2O5_liq   ::  Float64,
                                    sat_liq    ::  Float64,
                                    liq_wt     ::  Float64 )

    P2O5        = (P2O5_liq - sat_liq)/1e6
    fapt_wt     =  P2O5 * 2.37
    CaO_wt      =  P2O5 * 1.316

    return fapt_wt*liq_wt, CaO_wt*liq_wt
end

#=
    Routine to compute Sulfur saturation and adjust bulk-rock composition when sulfide crystallizes 
    based on different empirical models.
    Models implemented:
    - Liu et al. (2007) - Sulfur concentration at sulfide saturation (SCSS) in magmatic silicate melts.
    - O'Neill (2021) - Thermodynamic controls on sulfide saturation in silicate melts.
    - Bockrath et al. (2004) - Sulfur solubility in mafic melts: Implications for sulfur degassing of basaltic volcanoes.
    NR 7/11/2025
=#
"""
    logish(x)

    Safe natural logarithm that returns 0.0 for non-positive inputs instead of `-Inf` or `NaN`.

    Parameters
    ----------
    x : Float64
        Input value.

    Returns
    -------
    Float64
        `log(x)` if `x > 0`, otherwise `0.0`.
"""
function logish(x::Float64)
    if x <= 0.0
        return 0.0
    else
        return log(x)
    end
end

"""
    sulfur_saturation(out; model="1000ppm")

    Compute the sulfur concentration at sulfide saturation (SCSS) in the melt phase [ppm].

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output (must contain a melt phase with a "liq" solution phase).
    model : String, optional
        Saturation model (default: "1000ppm"). Valid options:
        - `"<N>ppm"`   — Fixed saturation value of N ppm (e.g., "1000ppm", "500ppm")
        - "Liu07"      — Liu et al. (2007); SCSS model, may be unreliable for anhydrous conditions
        - "Oneill21"   — O'Neill (2021); thermodynamic SCSS model using fO₂ from the minimization

    Returns
    -------
    C_s_liq : Float64
        Sulfur saturation concentration in the melt [ppm].
"""
function sulfur_saturation(     out     :: MAGEMin_C.gmin_struct{Float64, Int64};
                                model   :: String = "1000ppm"    )

    model_type = model[end-2:end]
    if model_type == "ppm"
        model = "fixed"
    end

    if model == "fixed"
        ppm_value   = parse(Float64, model[1:end-3])
        C_s_liq     = ppm_value
    elseif model == "Liu07"
        #=
            Liu et al., 2007
            Sulfur concentration at sulfide saturation (SCSS) in magmatic silicate melts
            NR: does not seem to work for anhydrous conditions
        =#
        X_ox        = out.SS_vec[out.SS_syms[:liq]].Comp
        oxides      = out.oxides
        oxide_data  = Dict(
            "SiO2"  => ("Si",   1),
            "Al2O3" => ("Al",   2),
            "CaO"   => ("Ca",   1),
            "MgO"   => ("Mg",   1),
            "FeO"   => ("Fe2",  1),
            "K2O"   => ("K",    2),
            "Na2O"  => ("Na",   2),
            "TiO2"  => ("Ti",   1),
            "Cr2O3" => ("Cr",   2),
            "H2O"   => ("H",    2),
            "O"     => ("Fe3",  1)
        )
        # compute unnormalized cation amounts
        cation_moles = Dict{String, Float64}()

        for (oxide, X) in zip(oxides, X_ox)
            if haskey(oxide_data, oxide)
                cation, ncat = oxide_data[oxide]
                cation_moles[cation] = get(cation_moles, cation, 0.0) + ncat * X
            end
        end
        # normalize to total cations
        total_cat = sum(values(cation_moles))
        cat_mol_F = Dict(k => v / total_cat for (k, v) in cation_moles)

        # Calculate MFM compositional parameter (cation mole fraction)
        MFM     = (cat_mol_F["Na"] + cat_mol_F["K"] + 2 * (cat_mol_F["Ca"] + cat_mol_F["Mg"] + cat_mol_F["Fe2"])) / (cat_mol_F["Si"] * (cat_mol_F["Al"] + cat_mol_F["Fe3"]))
        P_bar   = out.P_kbar * 1000.0
        T_K     = out.T_C + 273.15

        X_H2O_melt = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("H2O"), out.oxides)]
        X_FeO_melt = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("FeO"), out.oxides)]

        # Function to calculate ln(S in ppm)_SCSS
        ln_S_SCSS = 11.35251 -
                    4454.6 / T_K -
                    0.03190 * (P_bar / T_K) +
                    0.71006 * log(MFM) -
                    1.98063 * (MFM * X_H2O_melt) +
                    0.21867 * logish(X_H2O_melt) +
                    0.36192 * logish(X_FeO_melt)

        C_s_liq     = exp(ln_S_SCSS)
    elseif model == "Oneill21"
        #=
            O'Neill, H.S.C. (2021) The thermodynamic controls on sulfide saturation in silicate melts with application to ocean floor basalts.
        =#
        X_ox        = out.SS_vec[out.SS_syms[:liq]].Comp
        oxides      = out.oxides
        oxide_data  = Dict(
            "SiO2"  => ("Si",   1),
            "Al2O3" => ("Al",   2),
            "CaO"   => ("Ca",   1),
            "MgO"   => ("Mg",   1),
            "FeO"   => ("Fe2",  1),
            "K2O"   => ("K",    2),
            "Na2O"  => ("Na",   2),
            "TiO2"  => ("Ti",   1),
            "Cr2O3" => ("Cr",   2),
            "H2O"   => ("H",    2),
            "O"     => ("Fe3",  1)
        )
        # compute unnormalized cation amounts
        cation_moles = Dict{String, Float64}()

        for (oxide, X) in zip(oxides, X_ox)
            if haskey(oxide_data, oxide)
                cation, ncat = oxide_data[oxide]
                cation_moles[cation] = get(cation_moles, cation, 0.0) + ncat * X
            end
        end
        # normalize to total cations
        total_cat = sum(values(cation_moles))
        cat_mol_F = Dict(k => v / total_cat for (k, v) in cation_moles)
        P_bar       = out.P_kbar * 1000.0
        T_K         = out.T_C   + 273.15

        # compute lnCs2 capacity from O'Neill 2021
        lnCs2        = -23590.0 / T_K + 8.77 + (1673.0 / T_K) * (   6.7   * (cat_mol_F["Na"]
                                                                       + cat_mol_F["K"])
                                                                + 1.8  * cat_mol_F["Al"]
                                                                + 4.9  * cat_mol_F["Mg"]
                                                                + 8.1  * cat_mol_F["Ca"] 
                                                                + 5    * cat_mol_F["Ti"]
                                                                + 8.9  * cat_mol_F["Fe2"]
                                                                - 22.2 * cat_mol_F["Fe2"] * cat_mol_F["Ti"]
                                                                + 7.2  * cat_mol_F["Fe2"] * cat_mol_F["Si"])
                                                                - 2.06 * erf(-7.2 * cat_mol_F["Fe2"])

        #= Liu et al., 2007 relationship to get logfS2, extended after Bockrath et al. 2004 =#
        R           = 8.314
        ΔVr         = 0.904 #J/bar
        ΔFMQ        = out.dQFM
        logfO2      = out.fO2
        X_FeO_melt  = out.SS_vec[out.SS_syms[:liq]].Comp[findfirst(isequal("FeO"), out.oxides)]
        logfS2      = 6.7 - 12800 / T_K - 2 * log10(X_FeO_melt) + ΔFMQ + (ΔVr * (P_bar - 1)) / (2.303 * R * T_K)

        C_s_liq      = exp(lnCs2) * (exp10(logfS2)/exp10(logfO2))^(0.5)

    else
        println("Sulfur saturation model not recognized, using 1000ppm as default.")
        C_s_liq     = 1000.0
    end

    return C_s_liq
end

"""
    adjust_bulk_4_sulfide(S_liq, sat_liq, liq_wt)

    Compute the weight fractions of sulfide, FeO, and O returned to the bulk when the melt exceeds sulfur saturation.

    Parameters
    ----------
    S_liq : Float64
        Sulfur concentration in the melt [ppm].
    sat_liq : Float64
        Sulfur saturation concentration in the melt [ppm].
    liq_wt : Float64
        Melt weight fraction.

    Returns
    -------
    sulfide_wt : Float64
        Weight fraction of precipitated sulfide (scaled by `liq_wt`).
    FeO_wt : Float64
        FeO weight returned to the bulk (scaled by `liq_wt`).
    O_wt : Float64
        O weight correction (negative, scaled by `liq_wt`).
"""
function adjust_bulk_4_sulfide( S_liq  ::  Float64,
                                sat_liq ::  Float64,
                                liq_wt  ::  Float64 )


    sulfur_excess   = (S_liq - sat_liq)/1e6

    sulfide_wt      = sulfur_excess * 2.742
    FeO_wt          = sulfur_excess * 1.742
    O_wt            = - FeO_wt * 0.2227

    return sulfide_wt*liq_wt, FeO_wt*liq_wt, O_wt*liq_wt
end

#=
    H2O–CO2 solubility model: Sun & Yao (2026), Earth Planet. Sci. Lett.
    "M2Fluid": unified model calibrated on 2090 experiments,
    1 bar – 6 GPa, 660 – 1924 °C, carbonatite to rhyolite.

    Eq. 7 — H2O (wt%):
        S_H2O = S_H2O(m) + S_OH-
        ln S_H2O(m) = a0 ln P_H2O + a1 ln T + a2 (P²/T) + a3 (√P · X_Na2O/T) + a4 (X_Na2O/T)
        ln S_OH-    = a5 ln P_H2O + a6 ln T + (a7 + a8 X_Al2O3 + a9 X_Na2O)(P/T)

    Eq. 8 — CO2 (ppm):
        S_CO2 = S_CO2(m) + S_CO3²⁻
        ln S_CO2(m)  = b0 ln P_CO2 + b1 ln P + b2 X_Al2O3 √P + b3 X_MgO² + b4 X_K2O²
        ln S_CO3²⁻   = b0 ln P_CO2 + b5 ln T + b6 ln P + Ω + F + Q · S_H2O
        Ω = (b7 X_SiO2 + b8 X_CaO)(P/T) + (b9 X_CaO + b10 X_Na2O) √P
        F = b11 X_SiO2 + b12 X_Al2O3 + b13 X_FeOt + b14 X_MgO + b15 X_K2O
          + b16 ln(1 − (X_SiO2 + X_TiO2 + X_Al2O3 + X_P2O5))
        Q = b17 X_SiO2 + b18 X_Al2O3 + b19 X_Na2O + b20 X_K2O

    Units: T in K, P in bar, X = oxide molar fractions on volatile-free basis.
    NR 05/2026
=#

# Table 1 coefficients
const _SY26_a = (
     0.807805,      # a0
    -0.700336,      # a1
     2.949559e-6,   # a2
    -74.546597,     # a3
     7146.275452,   # a4
     0.449387,      # a5
    -0.360194,      # a6
    -0.308259,      # a7
     2.984703,      # a8
    -1.693710,      # a9
)

const _SY26_b = (
     0.711812,      # b0
     0.285783,      # b1
    -0.044440,      # b2
 -1019.065935,      # b3
  -451.908354,      # b4
     1.282791,      # b5
     0.651135,      # b6
     0.040306,      # b7
     0.583632,      # b8
    -0.116595,      # b9
    -0.036834,      # b10
    -9.761378,      # b11
    -9.184586,      # b12
   -18.300993,      # b13
   -11.489536,      # b14
     9.856334,      # b15
     3.126584,      # b16
     0.584610,      # b17
    -2.027465,      # b18
     0.552689,      # b19
    -1.570268,      # b20
)

"""
    volatile_saturation_SY26(out; P_H2O=NaN, P_CO2=NaN)

Compute H2O and CO2 solubility in the melt using the M2Fluid model of Sun & Yao (2026).

If `P_H2O` or `P_CO2` are not provided and a fluid phase is present in `out`, partial
pressures are estimated from the fluid mole fractions assuming ideal mixing.  If neither
a value nor a fluid phase is available for a volatile species, the corresponding
saturation is returned as `NaN`.

Parameters
----------
out : MAGEMin_C.gmin_struct{Float64, Int64}
    MAGEMin minimization output (must contain a melt phase).
P_H2O : Float64, optional
    Partial pressure of H₂O [bar].  When `NaN` (default), derived from fluid composition
    if a fluid phase is present.
P_CO2 : Float64, optional
    Partial pressure of CO₂ [bar].  When `NaN` (default), derived from fluid composition
    if a fluid phase is present.

Returns
-------
S_H2O : Float64
    H₂O solubility in the melt [wt%], or `NaN` if P_H2O cannot be determined.
S_CO2 : Float64
    CO₂ solubility in the melt [ppm], or `NaN` if P_CO2 cannot be determined.
"""
function volatile_saturation_SY26(  out     :: MAGEMin_C.gmin_struct{Float64, Int64};
                                    P_H2O   :: Float64 = NaN,
                                    P_CO2   :: Float64 = NaN,
                                    P_total :: Float64 = NaN  )

    out.frac_M == 0.0 && return NaN, NaN

    T_K   = out.T_C + 273.15
    P_bar = isnan(P_total) ? out.P_kbar * 1000.0 : P_total
    oxides = out.oxides

    # --- volatile-free mole fractions from the melt ---
    volatile_ox = ("H2O", "CO2", "S", "O")
    liq_comp    = out.SS_vec[out.SS_syms[:liq]].Comp   # mole fractions (include volatiles)

    dry_mol  = [liq_comp[i] for i in eachindex(oxides) if !(oxides[i] in volatile_ox)]
    dry_oxs  = [oxides[i]   for i in eachindex(oxides) if !(oxides[i] in volatile_ox)]
    tot_dry  = sum(dry_mol)
    Xd(name) = (i = findfirst(==(name), dry_oxs); isnothing(i) ? 0.0 : dry_mol[i] / tot_dry)

    # Fe³⁺ is tracked as "O" (16 g/mol) in MAGEMin: 1 mol O ↔ 2 mol Fe as FeO-equiv.
    # "O" is excluded from tot_dry, so add its contribution explicitly.
    idx_O   = findfirst(==("O"), oxides)
    n_O_raw = isnothing(idx_O) ? 0.0 : liq_comp[idx_O]

    X_SiO2  = Xd("SiO2")
    X_TiO2  = Xd("TiO2")
    X_Al2O3 = Xd("Al2O3")
    X_FeOt  = Xd("FeO") + 2.0 * n_O_raw / tot_dry
    X_MgO   = Xd("MgO")
    X_CaO   = Xd("CaO")
    X_Na2O  = Xd("Na2O")
    X_K2O   = Xd("K2O")
    X_P2O5  = Xd("P2O5")

    # --- derive P_H2O / P_CO2 from fluid phase when not supplied ---
    if (isnan(P_H2O) || isnan(P_CO2)) && out.frac_F > 0.0
        fl_mol   = out.bulk_F
        tot_fl   = sum(fl_mol)
        idx_H2O  = findfirst(==("H2O"), oxides)
        idx_CO2  = findfirst(==("CO2"), oxides)

        if isnan(P_H2O) && !isnothing(idx_H2O)
            P_H2O = (fl_mol[idx_H2O] / tot_fl) * P_bar
        end
        if isnan(P_CO2) && !isnothing(idx_CO2)
            P_CO2 = (fl_mol[idx_CO2] / tot_fl) * P_bar
        end
    end

    a = _SY26_a
    b = _SY26_b

    # --- H2O solubility (wt%) ---
    S_H2O = NaN
    if !isnan(P_H2O) && P_H2O > 0.0
        ln_Sm  = a[1]*log(P_H2O) + a[2]*log(T_K) + a[3]*(P_bar^2/T_K) +
                 a[4]*(sqrt(P_bar)*X_Na2O/T_K) + a[5]*(X_Na2O/T_K)
        ln_SOH = a[6]*log(P_H2O) + a[7]*log(T_K) + (a[8] + a[9]*X_Al2O3 + a[10]*X_Na2O)*(P_bar/T_K)
        S_H2O  = exp(ln_Sm) + exp(ln_SOH)
    end

    # --- CO2 solubility (ppm) ---
    S_CO2 = NaN
    if !isnan(P_CO2) && P_CO2 > 0.0
        S_H2O_eff = isnan(S_H2O) ? 0.0 : S_H2O   # coupling term; 0 if H2O unavailable

        inner = max(1.0 - (X_SiO2 + X_TiO2 + X_Al2O3 + X_P2O5), 1e-12)

        ln_Sm_co2 = b[1]*log(P_CO2) + b[2]*log(P_bar) + b[3]*X_Al2O3*sqrt(P_bar) +
                    b[4]*X_MgO^2 + b[5]*X_K2O^2

        Ω  = (b[8]*X_SiO2 + b[9]*X_CaO)*(P_bar/T_K) + (b[10]*X_CaO + b[11]*X_Na2O)*sqrt(P_bar)
        F  = b[12]*X_SiO2 + b[13]*X_Al2O3 + b[14]*X_FeOt + b[15]*X_MgO + b[16]*X_K2O +
             b[17]*log(inner)
        Q  = b[18]*X_SiO2 + b[19]*X_Al2O3 + b[20]*X_Na2O + b[21]*X_K2O

        ln_Sco3 = b[1]*log(P_CO2) + b[6]*log(T_K) + b[7]*log(P_bar) + Ω + F + Q*S_H2O_eff

        S_CO2 = exp(ln_Sm_co2) + exp(ln_Sco3)
    end

    return S_H2O, S_CO2
end

"""
    CO2_from_dissolved_H2O(out, S_H2O_wt; tol=1e-6)

Given a known dissolved H₂O content in the melt, compute the CO₂ saturation
concentration using the M2Fluid model of Sun & Yao (2026).

Assumes a binary H₂O–CO₂ fluid so that P_CO₂ = P − P_H₂O.  P_H₂O is found by
numerically inverting Eq. (7) via bisection on [0, P].

Parameters
----------
out : MAGEMin_C.gmin_struct{Float64, Int64}
    MAGEMin minimization output (must contain a melt phase).
S_H2O_wt : Float64
    Total dissolved H₂O content in the melt [wt%].
tol : Float64, optional
    Convergence tolerance on P_H₂O [bar] (default 1e-6).

Returns
-------
P_H2O : Float64
    H₂O partial pressure [bar].
P_CO2 : Float64
    CO₂ partial pressure [bar]  (= P − P_H₂O).
S_CO2 : Float64
    CO₂ saturation in the melt [ppm], or NaN if P_CO₂ ≤ 0.
"""
function CO2_from_dissolved_H2O(    out         :: MAGEMin_C.gmin_struct{Float64, Int64},
                                    S_H2O_wt    :: Float64;
                                    tol         :: Float64 = 1e-6   )

    out.frac_M == 0.0 && return NaN, NaN, NaN

    P_bar = out.P_kbar * 1000.0

    # upper bound: pure H2O fluid at total pressure
    S_max, _ = volatile_saturation_SY26(out; P_H2O = P_bar, P_CO2 = 0.0)

    if isnan(S_max) || S_H2O_wt >= S_max
        # melt is at or above pure-H2O saturation — no CO2 in fluid
        return P_bar, 0.0, NaN
    end

    # bisect on P_H2O in (0, P_bar] to match S_H2O_wt
    lo = 0.0
    hi = P_bar
    while (hi - lo) > tol
        mid = 0.5*(lo + hi)
        S, _ = volatile_saturation_SY26(out; P_H2O = mid, P_CO2 = 0.0)
        (isnan(S) || S < S_H2O_wt) ? lo = mid : hi = mid
    end
    P_H2O = 0.5*(lo + hi)

    P_CO2 = max(P_bar - P_H2O, 0.0)

    _, S_CO2 = volatile_saturation_SY26(out; P_H2O = P_H2O, P_CO2 = P_CO2)

    return P_H2O, P_CO2, S_CO2
end

"""
    CO2_from_dissolved_H2O(out)

Convenience overload: reads dissolved H₂O directly from the melt phase of `out`
(in wt%) and forwards to `CO2_from_dissolved_H2O(out, S_H2O_wt)`.

Returns `(NaN, NaN, NaN)` if no melt is present or the melt contains no H₂O.
"""
function CO2_from_dissolved_H2O(out :: MAGEMin_C.gmin_struct{Float64, Int64})
    liq_idx = get(out.SS_syms, :liq, 0)
    liq_idx == 0 && return NaN, NaN, NaN

    H2O_idx = findfirst(==("H2O"), out.oxides)
    isnothing(H2O_idx) && return NaN, NaN, NaN

    S_H2O_wt = out.SS_vec[liq_idx].Comp_wt[H2O_idx] * 100.0
    S_H2O_wt <= 0.0 && return NaN, NaN, NaN

    return CO2_from_dissolved_H2O(out, S_H2O_wt)
end

"""
    co2_saturation(out; model="SY26")

    Compute the CO₂ saturation concentration in the melt phase [ppm].

    Reads dissolved H₂O from the melt phase of `out`, inverts the SY26 H₂O
    equation to find P_H₂O, then evaluates CO₂ solubility at
    P_CO₂ = P_total − P_H₂O (binary H₂O–CO₂ fluid assumption).

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output (must contain a melt phase with dissolved H₂O).
    model : String, optional
        Saturation model (default: "SY26"). Currently only "SY26" is supported.

    Returns
    -------
    S_CO2 : Float64
        CO₂ saturation concentration in the melt [ppm], or NaN if H₂O is
        unavailable or the model is unrecognised.
"""
function co2_saturation(    out     :: MAGEMin_C.gmin_struct{Float64, Int64};
                            model   :: String = "SY26"  )
    if model == "SY26"
        _, _, S_CO2 = CO2_from_dissolved_H2O(out)
        return S_CO2
    else
        print("CO2 saturation model $model not recognized.\n")
        return NaN
    end
end

"""
    adjust_bulk_4_fluid(CO2_liq, sat_liq, liq_wt)

    Compute the weight fractions of CO₂ fluid and the oxide correction when the
    melt exceeds CO₂ saturation.

    Unlike mineral saturation phases (zircon, sulfide, apatite), the excess CO₂
    is not converted to a stoichiometrically distinct solid — it degasses into a
    CO₂-bearing fluid.  The same CO₂ weight is therefore both the fluid weight
    and the amount returned to the bulk CO₂ oxide budget for the next iteration.

    Parameters
    ----------
    CO2_liq : Float64
        CO₂ concentration in the melt [ppm].
    sat_liq : Float64
        CO₂ saturation concentration in the melt [ppm].
    liq_wt : Float64
        Melt weight fraction.

    Returns
    -------
    fluid_wt : Float64
        Weight fraction of CO₂ fluid formed (scaled by `liq_wt`).
    CO2_wt : Float64
        CO₂ weight returned to the bulk oxide budget (scaled by `liq_wt`);
        equals `fluid_wt` for a pure CO₂ fluid.
"""
function adjust_bulk_4_fluid(   CO2_liq ::  Float64,
                                sat_liq ::  Float64,
                                liq_wt  ::  Float64 )

    CO2_excess  = (CO2_liq - sat_liq) / 1e6   # ppm → dimensionless weight fraction
    fluid_wt    = CO2_excess
    CO2_wt      = CO2_excess

    return fluid_wt * liq_wt, CO2_wt * liq_wt
end
