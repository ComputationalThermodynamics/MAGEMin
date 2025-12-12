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
# The full functionality of MAGEMin is wrapped in ../gen/magemin_library.jl
# Yet, the routines here make it more convenient to use this from julia

# This file, gathers external function to compute other parameters such as melt viscosity

#=
- Giordano D, Russell JK, & Dingwell DB (2008). Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science Letters, 271, 123-134. (https://dx.doi.org/10.1016/j.epsl.2008.03.038)
=#
"""
    compute_melt_viscosity_G08(oxides, M_mol, T_C; A = -4.55)

    Takes as input arguments:
        oxides :: Vector{String}    -> oxide list of the melt composition
        M_mol  :: Vector{Float64}   -> melt composition in mol
        T_C    :: Float64           -> temperature in °C

    returns melt viscosity in Pa.s

    Formulation after Giordano et al., 2008
"""
function compute_melt_viscosity_G08(oxides, M_mol, T_C; A = -4.55)

    G08_ox_list = ["SiO2", "Al2O3", "TiO2", "FeO", "CaO", "MgO", "MnO", "Na2O", "K2O", "P2O5", "H2O", "F2O-1"]
    n_G08_ox    = length(G08_ox_list)
    G_mol       = zeros(Float64,n_G08_ox)

    for i=1:n_G08_ox
        if G08_ox_list[i] in oxides
            idx = findfirst(==(G08_ox_list[i]), oxides)
            G_mol[i] = M_mol[idx]
        end
    end
    G_mol = G_mol ./ (sum(G_mol)) .* 100.0            #normalize to 100

    b, b1,   = [159.6, -173.3, 72.13, 75.69, -38.9, -84.08, 141.54], [-2.43, -0.91, 17.62]
    c, c1    = [2.75, 15.72, 8.32, 10.2, -12.29, -99.54], [0.30]
    V        = G_mol[11] + G_mol[12]
    FM       = G_mol[4] + G_mol[6] + G_mol[7] 
    TA       = G_mol[3] + G_mol[2]
    NK       = G_mol[8] + G_mol[9]

    M        = [G_mol[1] + G_mol[3], G_mol[2], G_mol[4] + G_mol[7] + G_mol[10], G_mol[6], G_mol[5], G_mol[8] + V, V + log(1.0 + G_mol[11])]
    M1       = [(G_mol[1] + G_mol[3])*FM, (G_mol[1] + TA + G_mol[10]) * (NK + G_mol[11]), G_mol[2] * NK ]

    N        = [G_mol[1], TA, FM, G_mol[5], NK, log(1.0 + V)]
    N1       = [(G_mol[2] + FM + G_mol[5] - G_mol[10]) * (NK + V)]

    B        = b'* M + b1'*(M1)
    C        = c'* N + c1'*(N1)

    TK       = T_C + 273.12

    eta = 10^( A + B/(TK - C) )
    eta > 1e14 ? eta = 1e14 : eta
    
    return eta
end

