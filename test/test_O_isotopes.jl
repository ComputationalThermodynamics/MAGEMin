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

const O18M = MAGEMin_C

@testset verbose=true "DBOxygen database" begin
    db = get_O18_database()
    @test length(db.phase) == 154
    @test db.version == "DBOxygen 2.0.3"
    @test all(==(0.0), db.C)
    @test O18_fractionation(db, "Quartz", 600.0) == 0.0
    @test O18_fractionation(db, "Water", 500.0) ≈ 4.11e6 / 773.15^2 - 3.33e3 / 773.15  rtol=1e-12
    @test O18_fractionation(db, "Water", 500.0) ≈ 2.57  atol=5e-3
    @test O18_fractionation(db, "Forsterite", 1200.0) ≈ O18_fractionation(db, "Forsterite", 900.0) * (1173.15 / 1473.15)^2  rtol=1e-12
    @test O18_fractionation(db, "Forsterite", 1200.0; T_max=Inf) ≈ 1.50e6 / 1473.15^2 + 2.84e3 / 1473.15  rtol=1e-12
    @test_throws ErrorException get_O18_database("none")

    mapping = get_O18_mapping()
    for (key, parts) in mapping.map
        @test sum(last, parts) ≈ 1.0  atol=1e-3
        @test all(p -> haskey(db.id, first(p)), parts)
    end
end

@testset verbose=true "DBOxygen mapping coverage (TC and BR databases)" begin
    mapping     = get_O18_mapping().map
    not_phases  = ["O2", "qfm", "mw", "qif", "nno", "hm", "iw", "cco", "aH2O", "aO2", "aMgO", "aFeO", "aAl2O3", "aTiO2"]
    no_oxygen   = ["pyr", "gph", "po:trov", "po:trot", "fl:H2"]
    no_analogue = ["stv", "mpm", "fper:per", "fper:wu", "fsp:op", "flc:CO2", "fl:CO2",
                   "sa:spr4", "sa:spr5", "sa:fspm", "sa:spro", "sa:ospr",
                   "carp:mcar", "carp:fcar", "car:mcar", "car:fcar", "scp:me", "scp:coma", "scp:mizz"]
    allowed     = Set(vcat(not_phases, no_oxygen, no_analogue))
    for d in MAGEMin_C.DB_INFOS_GENERATED
        d.db_name in ("mp", "mpe", "mb", "mbe", "ig", "igad", "igd", "um", "ume", "po") || continue
        missing_ss = ["$(s.ss_name):$e" for s in d.data_ss if !(s.ss_name in ("liq", "DEW"))
                        for e in s.ss_em[2:end] if !(haskey(mapping, (s.ss_name, e)) || haskey(mapping, ("*", e)))]
        missing_pp = [p for p in d.data_pp if !haskey(mapping, (p, ""))]
        @test isempty(setdiff(vcat(missing_ss, missing_pp), allowed))
    end
end

@testset verbose=true "Vho et al. (2019) Table 3 pelite" begin
    db      = get_O18_database()
    fsp(ab, an, sa) = [("Albite", ab), ("Anorthite", an), ("K-Feldspar", sa)]
    stages  = [
        (800.0, [("Quartz", 1.0)], [
            (fsp(0.26, 0.02, 0.72), 1.0),
            (fsp(0.64, 0.30, 0.06), 1.1),
            ([("Grossular", 0.06), ("Pyrope", 0.29), ("Almandine", 0.63), ("Spessartine", 0.02)], 3.1),
            ([("Biotite", 1.0)], 3.3),
            ([("Sillimanite", 1.0)], 1.7),
            ([("Rutile", 1.0)], 4.7)]),
        (550.0, [("Quartz", 1.0)], [
            ([("Grossular", 0.18), ("Pyrope", 0.14), ("Almandine", 0.66), ("Spessartine", 0.02)], 4.4),
            ([("Phengite", 1.0)], 2.8),
            ([("Paragonite", 1.0)], 2.8),
            ([("Glaucophane", 1.0)], 3.4),
            ([("Rutile", 1.0)], 6.6)]) ]
    for (T, qz, phases) in stages
        for (parts, Δ_ref) in phases
            Δ = O18M._O18_delta(O18M._O18_mix(parts, db)..., T + 273.15)
            @test Δ ≈ Δ_ref  atol=0.15
        end
    end
end

@testset verbose=true "oxygen isotope mass balance" begin
    Δ       = [0.0, 1.2, 3.5, NaN, 5.0]
    nO      = [2.0, 8.0, 12.0, 4.0, 2.0]
    δ_sys, δ_qz, δ_ph = O18_mass_balance(Δ, nO, 10.0)
    use     = .!isnan.(Δ)
    @test δ_sys == 10.0
    @test sum(nO[use] .* δ_ph[use]) / sum(nO[use]) ≈ 10.0  rtol=1e-14
    @test δ_qz - δ_ph[3] ≈ 3.5
    δ_sys2, δ_qz2, _ = O18_mass_balance(Δ, nO, δ_ph[3]; anchor=3)
    @test δ_sys2 ≈ 10.0  rtol=1e-14
    @test δ_qz2 ≈ δ_qz  rtol=1e-14
end

@testset verbose=true "CIPW norm (Zhao & Zheng, 2003, Table 2)" begin
    ox      = ["SiO2", "TiO2", "Al2O3", "Fe2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O"]
    rocks   = [
        ([71.30, 0.31, 14.32, 1.21, 1.64, 0.05, 0.71, 1.84, 3.68, 4.07], 0.12, 0.05,
            Dict("Q" => 29.31, "C" => 0.94, "Or" => 24.23, "Ab" => 31.36, "An" => 8.09, "Hy" => 3.38, "Mt" => 1.77, "Il" => 0.59)),
        ([49.20, 1.84, 15.74, 3.79, 7.13, 0.20, 6.73, 9.47, 2.91, 1.10], 0.35, 0.11,
            Dict("Or" => 6.63, "Ab" => 24.97, "An" => 27.07, "Di" => 14.24, "Hy" => 15.55, "Ol" => 1.43, "Mt" => 5.58, "Il" => 3.55)),
        ([40.60, 2.66, 14.33, 5.48, 6.17, 0.26, 6.39, 11.89, 4.79, 3.46], 1.07, 0.60,
            Dict("Or" => 3.01, "An" => 7.61, "Lc" => 14.14, "Ne" => 22.60, "Di" => 33.38, "Ol" => 2.54, "Mt" => 8.19, "Il" => 5.21)),
        ([53.41, 1.22, 15.28, 6.66, 2.30, 0.46, 1.35, 1.79, 9.48, 4.65], 0.14, 0.03,
            Dict("Or" => 28.42, "Ab" => 19.44, "Ne" => 18.98, "Ac" => 19.92, "Ns" => 1.36, "Di" => 6.66, "Ol" => 2.46, "Il" => 2.40)),
        ([42.26, 0.63, 4.23, 3.61, 6.58, 0.41, 31.24, 5.05, 0.49, 0.34], 0.10, 0.30,
            Dict("Or" => 2.10, "Ab" => 4.33, "An" => 8.82, "Di" => 11.79, "Hy" => 16.70, "Ol" => 48.86, "Mt" => 5.51, "Il" => 1.27)) ]
    groups  = Dict("Hy" => ("En", "Fs"), "Ol" => ("Fo", "Fa"), "Di" => ("Di", "Hd"))
    mass(st) = sum(v * get_molar_mass(o) for (o, v) in st)
    for (wt, p2o5, co2, ref) in rocks
        mol     = wt ./ [get_molar_mass(o) for o in ox]
        mol[8] -= 10.0 / 3.0 * p2o5 / get_molar_mass("P2O5") + co2 / get_molar_mass("CO2")
        n, D    = cipw_norm(mol, ox)
        @test D == 0.0
        w       = Dict(m[1] => n[m[1]] * mass(m[2]) for m in O18M._CIPW_minerals)
        norm_wt(k) = haskey(groups, k) ? w[groups[k][1]] + w[groups[k][2]] : w[k]
        s_ref   = sum(values(ref))
        s_calc  = sum(norm_wt(k) for k in keys(ref))
        for (k, v) in ref
            @test 100.0 * norm_wt(k) / s_calc ≈ 100.0 * v / s_ref  atol=0.3
        end
        O_ox    = sum(mol[j] * O18M._O18_nu_O[ox[j]] for j in eachindex(ox))
        O_norm  = sum(n[m[1]] * sum(v * O18M._O18_nu_O[o] for (o, v) in m[2]) for m in O18M._CIPW_minerals)
        @test O_norm ≈ O_ox  rtol=1e-12
    end
end

@testset verbose=true "O18_prediction" begin
    for (db, P, T) in [("mp", 8.0, 650.0), ("mb", 12.0, 600.0), ("ig", 5.0, 1300.0), ("um", 20.0, 600.0), ("po", 20.0, 550.0)]
        data    = Initialize_MAGEMin(db, verbose=-1, solver=0)
        out     = single_point_minimization(P, T, data, test=0)
        for melt_model in ("norm", "ZZ03", "liq_em"), unmapped in (:exclude, :index)
            o18     = O18_prediction(out, 6.0; melt_model, unmapped)
            use     = o18.f_O .> 0.0
            @test sum(o18.f_O) ≈ 1.0  rtol=1e-12
            @test sum(o18.f_O[use] .* o18.δ_ph[use]) ≈ 6.0  rtol=1e-12
            @test all(isfinite, o18.δ_ph[use])
            @test 0.0 <= o18.f_O_unmapped < 0.05
            iq      = findfirst(==("q"), o18.ph)
            isnothing(iq) || @test o18.δ_ph[iq] ≈ o18.δ_qz  rtol=1e-12
            ia      = findfirst(use)
            o18b    = O18_prediction(out, o18.δ_ph[ia]; anchor=o18.ph[ia], melt_model, unmapped)
            @test o18b.δ_sys ≈ 6.0  rtol=1e-10
            @test o18b.δ_ph[use] ≈ o18.δ_ph[use]  rtol=1e-10
        end
        Finalize_MAGEMin(data)
    end

    data    = Initialize_MAGEMin("ig", verbose=-1, solver=0)
    out     = single_point_minimization(5.0, 1300.0, data, test=0)
    o_norm  = O18_prediction(out, 5.5)
    o_zz    = O18_prediction(out, 5.5; melt_model="ZZ03")
    @test Δ18O(o_norm, "liq", "ol") > 0.0
    @test Δ18O(o_norm, "opx", "ol") > 0.0
    @test abs(Δ18O(o_norm, "liq", "ol") - Δ18O(o_zz, "liq", "ol")) < 0.5
    @test 5.0 < o_norm.δ_ph[findfirst(==("ol"), o_norm.ph)] < 5.5
    @test any(f -> occursin("1/T²", f), o_norm.flags)
    @test isnan(Δ18O(o_norm, "liq", "q"))
    @test_throws ErrorException O18_prediction(out, 5.5; anchor="q")
    @test_throws ErrorException O18_prediction(out, 5.5; melt_model="none")
    @test_throws ErrorException O18_prediction(out, 5.5; unmapped=:none)
    Finalize_MAGEMin(data)
end
