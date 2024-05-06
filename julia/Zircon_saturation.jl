# Routine to compute zircon saturation and adjust bulk-rock composition when zircon crystallizes
# NR 12/04/2023

function zirconium_saturation(  out     :: gmin_struct{Float64, Int64}; 
                                model   :: String = "WH"    )

    if out.frac_M > 0.0                            
        if model == "WH" || model == "B"
            ref_ox      = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "O"; "Cr2O3"; "MnO"; "H2O"; "S"];
            ratio_cation= [1/3.0,2/5,1/2,1/2,1/2,2/5,2/3,2/3,1/3,1,2/5,1/2,2/3,1]
            
            cation_name = ["Si"; "Al"; "Ca"; "Mg"; "Fe"; "Fe"; "K"; "Na"; "Ti"; "O"; "Cr"; "Mn"; "H"; "S"]
            
            cation_idx  = [findfirst(isequal(x), ref_ox) for x in out.oxides];
            cation      = out.bulk_M_wt .*ratio_cation[cation_idx];
            cation      = cation ./ sum(cation);

            Na          = findall(cation_name[cation_idx] .== "Na")[1];
            K           = findall(cation_name[cation_idx] .== "K")[1];
            Ca          = findall(cation_name[cation_idx] .== "Ca")[1];
            Al          = findall(cation_name[cation_idx] .== "Al")[1];
            Si          = findall(cation_name[cation_idx] .== "Si")[1];

            M           = (cation[Na] + cation[K]+ 2.0*cation[Ca])/(cation[Al]*cation[Si]);
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
            opt_basi_oxides = ["SiO2","TiO2","Al2O3","MgO","MnO","FeO","Fe2O3","CaO","Na2O","K2O","P2O5","H2O"]
            optical_basicity= [0.48,0.75,0.60,0.78,0.96,1.00,0.77,1.00,1.15,1.40,0.33,0.40]
            n_oxygen        = [2.0,2,3,1,1,1,3,1,1,1,5,1]
            
            commonOxide   =  intersect(out.oxides,opt_basi_oxides)
            idOx_MM       = [findfirst(isequal(x), out.oxides) for x in commonOxide]
            idOx_OB       = [findfirst(isequal(x), opt_basi_oxides) for x in commonOxide]

            liqCompNorm   = out.bulk[idOx_MM] ./ sum(out.bulk[idOx_MM])
            oxListDry     = findall(commonOxide .!= "H2O")
            opt_basicity  = sum(liqCompNorm[oxListDry] .* optical_basicity[idOx_OB[oxListDry]].* n_oxygen[idOx_OB[oxListDry]]) / sum(liqCompNorm[oxListDry] .* n_oxygen[idOx_OB[oxListDry]])
            xH2O          = liqCompNorm[findall(commonOxide .== "H2O")[1]]

            C_zr_liq   = exp(0.96 - 5790.0/(out.T_C+273.15) - 1.28*(out.P_kbar/10.0) + 12.39*opt_basicity + 0.83*xH2O + 2.06*(out.P_kbar/10.0)*opt_basicity)
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
