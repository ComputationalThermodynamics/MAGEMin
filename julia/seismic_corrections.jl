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
    wave_melt_correction(gv, P_kbar, solid_Vp, solid_Vs, aspectRatio=0.3)

    Melt-fraction correction for P- and S-wave velocities of the solid aggregate.

    Uses the reduction formulation of Clark et al. (2017), based on the
    equilibrium geometry model for the solid skeleton of Takei (1997).

    If `gv.melt_fraction == 0.0`, a fixed near-surface porosity correction
    (aspect ratio 0.25, Poisson ratio 0.25) is applied instead, using `P_kbar`
    to estimate depth.

    Parameters
    ----------
    gv : LibMAGEMin.global_variables
        Global variables structure (melt/solid fractions, densities and moduli).
    P_kbar : Float64
        Pressure [kbar], used for the depth/porosity estimate when `melt_fraction == 0`.
    solid_Vp : Float64
        P-wave velocity of the solid aggregate [km/s] (uncorrected, or anelastically corrected).
    solid_Vs : Float64
        S-wave velocity of the solid aggregate [km/s] (uncorrected, or anelastically corrected).
    aspectRatio : Float64, optional
        Coefficient defining the geometry of the solid framework (contiguity):
        0.0 (layered melt) < 0.1 (grain boundary melt) < 1.0 (melt in separated
        bubble pockets) (default: 0.3).

    Returns
    -------
    (Vp_cor, Vs_cor) : Tuple{Float64, Float64}
        Melt-corrected P- and S-wave velocities of the solid aggregate [km/s].

    Reference
    ---------
    Clark A.N. and Lesher C.E. (2017)
"""
function wave_melt_correction(  gv, 
                                P_kbar      :: Float64,
                                solid_Vp    :: Float64,
                                solid_Vs    :: Float64,
                                frac_M_vol  :: Float64,
                                frac_S_vol  :: Float64,
                                frac_F_vol  :: Float64;
                                aspectRatio :: Float64 = 0.25,
                                shallow_correction  :: Bool = false,
                                fluid_as_melt       :: Bool = false)

    Vp_cor = solid_Vp
    Vs_cor = solid_Vs

    if frac_M_vol > 0.0

        if fluid_as_melt
            melt_fraction = frac_M_vol + frac_F_vol
            solid_fraction = frac_S_vol
        else 
            melt_fraction = frac_M_vol / ( 1.0 - frac_F_vol)
            solid_fraction = frac_S_vol  / ( 1.0 - frac_F_vol)
        end

        poisson = 0.25

        aij = ( (0.318, 6.780, 57.560,  0.182),
                (0.164, 4.290, 26.658,  0.464),
                (1.549, 4.814,  8.777, -0.290) )

        bij = ( (-0.3238,  0.2341),
                (-0.1819,  0.5103) )

        a = ntuple(i -> aij[i][1]*exp( aij[i][2]*(poisson-0.25) + aij[i][3]*(poisson-0.25)^3.0 ) + aij[i][4], 3)
        b = ntuple(i -> bij[i][1]*poisson + bij[i][2], 2)

        nk  = a[1]*aspectRatio + a[2]*(1.0 - aspectRatio) + a[3]*aspectRatio*(1.0 - aspectRatio)*(0.5 - aspectRatio)
        nmu = b[1]*aspectRatio + b[2]*(1.0 - aspectRatio)

        ksk_k   = aspectRatio^nk
        musk_mu = aspectRatio^nmu

        ksk  = ksk_k  *gv.solid_bulkModulus
        musk = musk_mu*gv.solid_shearModulus

        kb = (1.0 - melt_fraction)*ksk
        mu = (1.0 - melt_fraction)*musk

        LambdaK = gv.solid_bulkModulus /kb
        LambdaG = gv.solid_shearModulus/mu

        beta  = gv.solid_bulkModulus /gv.melt_bulkModulus
        gamma = gv.solid_shearModulus/gv.solid_bulkModulus

        deltaVp = (((((beta-1.0)*LambdaK)/((beta-1.0) + LambdaK) + 4.0/3.0*gamma*LambdaG)/(1.0 + 4.0/3.0*gamma)) - (1.0 - gv.melt_density/gv.solid_density))*(melt_fraction/2.0)
        deltaVs = (LambdaG - (1.0 - gv.melt_density/gv.solid_density))*(melt_fraction/2.0)

        Vp_cor = solid_Vp - deltaVp*solid_Vp
        Vs_cor = solid_Vs - deltaVs*solid_Vs

        if Vp_cor < 0.0; Vp_cor = 0.0; end
        if Vs_cor < 0.0; Vs_cor = 0.0; end
    end

    if shallow_correction
        # aspect ratio and Poisson ratio fixed to fit Vs = 2.0 km/s at the surface
        aspectRatio = 0.25
        poisson     = 0.25

        depth_m            = P_kbar*1e5/(2600.0*9.81)
        porosity_fraction  = 0.474/(1.0 + 0.071*depth_m)^5.989

        porosity_density    = 1000.0

        aij = ( (0.318, 6.780, 57.560,  0.182),
                (0.164, 4.290, 26.658,  0.464),
                (1.549, 4.814,  8.777, -0.290) )

        bij = ( (-0.3238,  0.2341),
                (-0.1819,  0.5103) )

        a = ntuple(i -> aij[i][1]*exp( aij[i][2]*(poisson-0.25) + aij[i][3]*(poisson-0.25)^3.0 ) + aij[i][4], 3)
        b = ntuple(i -> bij[i][1]*poisson + bij[i][2], 2)

        nk  = a[1]*aspectRatio + a[2]*(1.0 - aspectRatio) + a[3]*aspectRatio*(1.0 - aspectRatio)*(0.5 - aspectRatio)
        nmu = b[1]*aspectRatio + b[2]*(1.0 - aspectRatio)

        ksk_k   = aspectRatio^nk
        musk_mu = aspectRatio^nmu

        ksk  = ksk_k  *gv.solid_bulkModulus
        musk = musk_mu*gv.solid_shearModulus

        mu = (1.0 - porosity_fraction)*musk

        LambdaG = gv.solid_shearModulus/mu

        deltaVs = (LambdaG - (1.0 - porosity_density/gv.solid_density))*(porosity_fraction/2.0)

        Vs_cor = Vs_cor - deltaVs*Vs_cor
        if Vs_cor < 0.0; Vs_cor = 0.0; end
    end

    return Vp_cor, Vs_cor
end


"""
    anelastic_correction(water, Vs0, P_kbar, T_K)

    Anelastic (attenuation) correction for the S-wave velocity, following the
    model of Behn et al. (2009), with frequency/grain-size terms from
    Cobden et al. (2018).

    Parameters
    ----------
    water : Int
        Water content mode: `0` = dry mantle, `1` = damp mantle, `2` = wet mantle (saturated).
    Vs0 : Float64
        Unrelaxed (elastic) S-wave velocity [km/s].
    P_kbar : Float64
        Pressure [kbar].
    T_K : Float64
        Temperature [K].

    Returns
    -------
    Vs_anel : Float64
        Anelastically corrected S-wave velocity [km/s].

    Reference
    ---------
    Behn M.D., Hirth G., Elsenbeck J.R. (2009); Cobden L. et al. (2018)
"""
function anelastic_correction(water::Int, Vs0::Float64, P_kbar::Float64, T_K::Float64)

    kbar2pa = 100.0e3
    Pref    = P_kbar*kbar2pa            # Pa
    R       = 8.31446261815324

    # values based on fitting experimental constraints (Behn et al., 2009)
    alpha  = 0.27
    B0     = 1.28e8                     # m/s
    dref   = 1.24e-5                    # m
    COHref = 50.0/1e6                   # 50 H/1e6 Si
    Gref   = 1.09
    Eref   = 505.0e3                    # J/mol
    Vref   = 1.2e-5                     # m3/mol

    G = 1.00
    E = 420.0e3                         # J/mol (activation energy)
    V = 1.2e-5                          # m3/mol (activation volume)

    omega = 3.0                         # Hz (frequency for Toba)
    d     = 1e-2                        # m  (grain size)

    if water == 0
        COH, rH = 50.0/1e6,   0.0       # dry mantle
    elseif water == 1
        COH, rH = 1000.0/1e6, 1.0       # damp mantle
    elseif water == 2
        COH, rH = 3000.0/1e6, 2.0       # wet mantle (saturated)
    else
        @warn "anelastic_correction: water mode $water is not implemented, defaulting to dry mantle"
        COH, rH = 50.0/1e6,   0.0
    end

    B    = B0*dref^(G-Gref)*(COH/COHref)^rH * exp(((E + Pref*V) - (Eref + Pref*Vref))/(R*T_K))
    Qinv = (B*d^(-G)*(1.0/omega) * exp(-(E + Pref*V)/(R*T_K)))^alpha

    Vs_anel = Vs0*(1.0 - Qinv/(2.0*tan(pi*alpha/2.0)))

    return Vs_anel
end
