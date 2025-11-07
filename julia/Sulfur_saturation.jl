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
# Routine to compute zircon saturation and adjust bulk-rock composition when sulfide crystallizes
# NR 7/11/2025

function sulfur_saturation(    out     :: MAGEMin_C.gmin_struct{Float64, Int64};
                                model   :: String = "1000ppm"    )

    model_type = model[end-2:end]
    if model_type == "ppm"
        ppm_value   = parse(Float64, model[1:end-3])
        C_s_liq     = ppm_value
    else 
        println("Sulfur saturation model not recognized, using 1000ppm as default.")
        C_s_liq     = 1000.0
    end

    return C_s_liq
end

function adjust_bulk_4_sulfide( S_liq  ::  Float64,
                                sat_liq ::  Float64 )

    Fe_wt       = 0.0
    sulfur_excess   = (S_liq - sat_liq)/1e4

    sulfide_wt  = sulfur_excess*2.742
    Fe_wt       = sulfur_excess * 1.742


    return sulfide_wt, Fe_wt
end
