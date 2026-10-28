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

@testset verbose=true "ρ-H sections" begin
    @test MAGEMin_C.SECTION_XY[:rh] == (:rho, :H)

    data    = Initialize_MAGEMin("mp", verbose=-1, solver=0)
    Pbox, Tbox = (2.0, 10.0), (550.0, 750.0)
    ref     = single_point_minimization(6.0, 650.0, data, test=0)
    M       = sum(ref.bulk .* [get_molar_mass(o) for o in ref.oxides])
    @test ref.rho ≈ 1000.0 * M / ref.V  rtol=1e-6

    gv, z_b, DB, splx = data.gv[1], data.z_b[1], data.DB[1], data.splx_data[1]
    gv      = use_predefined_bulk_rock(gv, 0, data.db)
    o_rh, i_rh = section_point(:rh, ref.rho, ref.enthalpy[1], gv, z_b, DB, splx; Pbox = Pbox, Tbox = Tbox, P0 = 5.0, T0 = 620.0)
    o_vh, i_vh = section_point(:vh, ref.V,   ref.enthalpy[1], gv, z_b, DB, splx; Pbox = Pbox, Tbox = Tbox, P0 = 5.0, T0 = 620.0)
    @test i_rh.status == 0
    @test o_rh.rho          ≈ ref.rho           rtol=1e-6
    @test o_rh.enthalpy[1]  ≈ ref.enthalpy[1]   rtol=1e-6
    @test o_rh.P_kbar       ≈ o_vh.P_kbar       rtol=1e-4
    @test o_rh.T_C          ≈ o_vh.T_C          rtol=1e-4

    xr, yr, grid = section_ranges(:rh, data; Pbox = Pbox, Tbox = Tbox, n = 5, test = 0)
    xv, _, _     = section_ranges(:vh, data; Pbox = Pbox, Tbox = Tbox, n = 5, test = 0)
    @test xr[1] ≈ 1000.0 * M / xv[2]  rtol=1e-6
    @test xr[2] ≈ 1000.0 * M / xv[1]  rtol=1e-6

    ρ       = [xr[1] + f*(xr[2] - xr[1]) for f in (0.3, 0.5, 0.7)]
    H       = [yr[1] + f*(yr[2] - yr[1]) for f in (0.5, 0.5, 0.5)]
    outs, infos = multi_point_section(:rh, ρ, H, data; Pbox = Pbox, Tbox = Tbox, test = 0, grid = grid, progressbar = false)
    for (o, inf, r, h) in zip(outs, infos, ρ, H)
        inf.status == 0 || continue
        @test o.rho         ≈ r  rtol=1e-5
        @test o.enthalpy[1] ≈ h  rtol=1e-5
    end
    @test count(i -> i.status == 0, infos) >= 1
    Finalize_MAGEMin(data)
end
