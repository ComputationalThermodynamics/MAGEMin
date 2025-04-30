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
# Routine to compute zircon saturation and adjust bulk-rock composition when zircon crystallizes
# NR 12/04/2024, benchmarked HD 25/06/2024

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
            xH2O = liqCompNorm[xH2O_index]

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
                                sat_liq ::  Float64 )

    SiO2_wt         = 0.0
    O2_wt           = 0.0
    zircon_wt       = 0.0
    zircon_excess   = (zr_liq - sat_liq)/1e4

    zircon_wt   = zircon_excess*0.497644
    SiO2_wt     = (zircon_wt *0.327765)
    O2_wt       = (zircon_wt *0.174570)


    return zircon_wt ,SiO2_wt ,O2_wt
end
