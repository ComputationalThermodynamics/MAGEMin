#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
using Test
using MAGEMin_C
import MAGEMin_jll

const GH_EOS_LIB = isfile("libMAGEMin.dylib") ? abspath("libMAGEMin.dylib") : MAGEMin_jll.libMAGEMin

function gh_mix_G(x_h2o, T, P)
    dGdx = Ref{Cdouble}(0.0)
    g = ccall((:GH_pitzer_sterner_mix_G, GH_EOS_LIB), Cdouble,
              (Cdouble, Cdouble, Cdouble, Ref{Cdouble}), x_h2o, T, P, dGdx)
    return g, dGdx[]
end

function gh_G(is_H2O, T, P)
    return ccall((:GH_pitzer_sterner_G, GH_EOS_LIB), Cdouble,
                 (Cint, Cdouble, Cdouble), is_H2O ? 1 : 0, T, P)
end

@testset "GH fluid EOS regression" begin
    T = 800.0
    P = 1000.0

    g_h2o_mix, _ = gh_mix_G(1.0, T, P)
    @test isapprox(g_h2o_mix, gh_G(1, T, P); atol=1e-10, rtol=0)

    g_co2_mix, _ = gh_mix_G(0.0, T, P)
    @test isapprox(g_co2_mix, gh_G(0, T, P); atol=1e-10, rtol=0)

    x = 0.25
    g_center, dg_center = gh_mix_G(x, T, P)
    dx = 1e-4
    g_plus, _  = gh_mix_G(x + dx, T, P)
    g_minus, _ = gh_mix_G(x - dx, T, P)
    @test isapprox(dg_center, (g_plus - g_minus) / (2.0 * dx); atol=1e-2, rtol=1e-4)
end
