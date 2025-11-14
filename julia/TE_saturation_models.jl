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
            C_zr_zrn    = 497644;

            # Watson & Harrison (1983)
            if model == "WH"
                C_zr_liq    = C_zr_zrn/exp((-3.80 - 0.85*(M - 1.0)) + (12900/(out.T_C+273.15)))
            end
            # Boehnke et al., 2013
            if model == "B"
                C_zr_liq    = C_zr_zrn/exp((-1.48 - 1.16*(M - 1.0)) + (10108/(out.T_C+273.15)))
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

function phosphate_saturation(      out     :: MAGEMin_C.gmin_struct{Float64, Int64};
                                    model   :: String = "Tollari06"    )

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

        AoCNK       = Al2O3_mol / (CaO_mol + Na2O_mol + K2O_mol)    # A over CNK

        C_P2O5_WH   = 52.5525567/ exp( (8400.0 + 2.64e4(SiO2_wt - 0.5))/T_K - (3.1 + 12.4(SiO2_wt - 0.5)) )  # wt% as in Harrison and Watson 1984

        C_P2O5_liq  =  (C_P2O5_WH * (AoCNK -1.0) * 6429.0/out.T_C) * 1e4  # ppm Harrison and Watson 1984 + Correction from Bea et al. 1992
    else
        print("Cannot compute phosphate saturation in liquid if melt is not predicted!\n")
        C_P2O5_liq = -1
    end
    
    return C_P2O5_liq
end


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
function logish(x::Float64)
    if x <= 0.0
        return 0.0
    else
        return log(x)
    end
end

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

function adjust_bulk_4_sulfide( S_liq  ::  Float64,
                                sat_liq ::  Float64,
                                liq_wt  ::  Float64 )

    Fe_wt           = 0.0
    sulfur_excess   = (S_liq - sat_liq)/1e6

    sulfide_wt      = sulfur_excess * 2.742
    Fe_wt           = sulfur_excess * 1.742


    return sulfide_wt*liq_wt, Fe_wt*liq_wt
end
