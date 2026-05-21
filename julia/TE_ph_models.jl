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
    bi_Li_CB_model(y, f, T)

    Compute the Li partition coefficient between biotite and melt, after Beard (2025).

    Parameters
    ----------
    y : Float64
        Biotite Al-content compositional variable (Al on T site).
    f : Float64
        Biotite Fe/(Fe+Mg) compositional variable.
    T : Float64
        Temperature [°C].

    Returns
    -------
    KD_Li : Float64
        Li partition coefficient (biotite/melt).
"""
function bi_Li_CB_model(    y :: Float64,
                            f :: Float64,
                            T :: Float64 )

    bi_xSiT     = 0.5-f/2.0 - y/2.0
    bi_SiT      = (2.0 * bi_xSiT)+2.0 # need cations, not fraction
    ln_D_Li     = (-22.359) + (5.592*bi_SiT) + (10000*0.67*1/(T+273.15))
    KD_Li       = exp(ln_D_Li)

    return KD_Li 
end

"""
    bi_Li_IL_model(T)

    Compute the Li partition coefficient between biotite and melt using a linear temperature model (Imai & Liang, unpublished).

    Parameters
    ----------
    T : Float64
        Temperature [°C].

    Returns
    -------
    KD_Li : Float64
        Li partition coefficient (biotite/melt).
"""
function bi_Li_IL_model(   T :: Float64 )
    return -0.0076*(T)+6.5775
end

"""
    cd_Li_IL_model(T)

    Compute the Li partition coefficient between cordierite and melt using a linear temperature model (Imai & Liang, unpublished).

    Parameters
    ----------
    T : Float64
        Temperature [°C].

    Returns
    -------
    KD_Li : Float64
        Li partition coefficient (cordierite/melt).
"""
function cd_Li_IL_model(   T :: Float64 )
    return -0.0021*(T)+1.933
end