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
# Oxygen isotope (δ18O) equilibrium fractionation between stable phases
# References: Vho, Lanari & Rubatto (2019) J. Petrol. doi:10.1093/petrology/egaa001 (DBOxygen 2.0.3)
#             Zhao & Zheng (2003) Chem. Geol. 193, 59-80. doi:10.1016/s0009-2541(02)00226-7
#             Bindeman (2008) Rev. Mineral. Geochem. 69, 445-478. doi:10.2138/rmg.2008.69.12

"""
    O18_database

    Oxygen isotope fractionation database, quartz referenced:
    `1000 ln α(Qz-phase) = A·10⁶/T² + B·10³/T + C` with T in K.

    Fields
    ------
    info : String
        Citation of the database.
    version : String
        Database version.
    phase : Vector{String}
        Phase names.
    group : Vector{String}
        Mineral group of each phase.
    A, sA, B, sB, C, sC : Vector{Float64}
        Fractionation parameters and their 2σ uncertainties.
    id : Dict{String,Int64}
        Phase name to row index.
"""
struct O18_database
    info    :: String
    version :: String
    phase   :: Vector{String}
    group   :: Vector{String}
    A       :: Vector{Float64}
    sA      :: Vector{Float64}
    B       :: Vector{Float64}
    sB      :: Vector{Float64}
    C       :: Vector{Float64}
    sC      :: Vector{Float64}
    id      :: Dict{String,Int64}
end

"""
    O18_mapping

    Mapping of MAGEMin phases and end-members onto `O18_database` phases.

    Fields
    ------
    map : Dict{Tuple{String,String},Vector{Tuple{String,Float64}}}
        Key `(phase, em)` returns the list of `(database phase, oxygen fraction)`. Pure phases use `em = ""`,
        `phase = "*"` applies to any solution phase.
"""
struct O18_mapping
    map :: Dict{Tuple{String,String},Vector{Tuple{String,Float64}}}
end

"""
    out_O18

    Result of `O18_prediction` for one MAGEMin point.

    Fields
    ------
    T_C : Float64
        Temperature [°C].
    δ_sys : Float64
        δ18O of the oxygen accounted for in the system [‰ VSMOW].
    δ_qz : Float64
        δ18O of a virtual quartz in equilibrium with the assemblage [‰].
    N_O : Float64
        Moles of oxygen accounted for per 100 g of system.
    ph : Vector{String}
        Phase names, same order as `out.ph`.
    f_O : Vector{Float64}
        Fraction of the accounted oxygen held by each phase.
    A, B, C : Vector{Float64}
        Quartz-phase fractionation parameters of each phase.
    δ_ph : Vector{Float64}
        δ18O of each phase [‰], NaN when the phase is not mapped.
    em : Vector{Vector{String}}
        End-member names of each phase (empty for pure phases).
    δ_em : Vector{Vector{Float64}}
        δ18O of each end-member [‰], NaN when not mapped.
    δ_liq, δ_fl, δ_sol : Float64
        Oxygen-weighted δ18O of melt, fluid and solid phases [‰], NaN when absent.
    f_O_unmapped : Float64
        Fraction of the total phase oxygen left out because it could not be mapped.
    flags : Vector{String}
        Warnings raised during the computation.
"""
struct out_O18
    T_C             :: Float64
    δ_sys           :: Float64
    δ_qz            :: Float64
    N_O             :: Float64
    ph              :: Vector{String}
    f_O             :: Vector{Float64}
    A               :: Vector{Float64}
    B               :: Vector{Float64}
    C               :: Vector{Float64}
    δ_ph            :: Vector{Float64}
    em              :: Vector{Vector{String}}
    δ_em            :: Vector{Vector{Float64}}
    δ_liq           :: Float64
    δ_fl            :: Float64
    δ_sol           :: Float64
    f_O_unmapped    :: Float64
    flags           :: Vector{String}
end

const _O18_nu_O = Dict( "SiO2" => 2.0, "TiO2" => 2.0, "Al2O3" => 3.0, "Cr2O3" => 3.0, "Fe2O3" => 3.0, "CO2" => 2.0,
                        "FeO" => 1.0, "MgO" => 1.0, "CaO" => 1.0, "MnO" => 1.0, "NiO" => 1.0, "Na2O" => 1.0,
                        "K2O" => 1.0, "H2O" => 1.0, "O" => 1.0, "P2O5" => 5.0, "S" => 0.0, "Fe" => 0.0 )

const _ZZ03_index = Dict(   "SiO2" => 1.0000, "TiO2" => 0.6322, "Al2O3" => 0.8697, "Fe2O3" => 0.4809, "FeO" => 0.4140,
                            "MnO" => 0.4066, "MgO" => 0.3555, "CaO" => 0.2374, "Na2O" => 0.0701, "K2O" => 0.0578 )

const _ZZ03_A_qz = 11.83

const _O18_fluid_names = ("fl", "flc", "H2O", "h2o", "DEW")

function _read_O18_database(path :: String)
    rows = Vector{Vector{SubString{String}}}()
    for line in eachline(path)
        startswith(line, '#') && continue
        isempty(strip(line)) && continue
        parts = split(line, ',')
        parts[1] == "phase" && continue
        push!(rows, strip.(parts))
    end
    num(i)  = [parse(Float64, r[i]) for r in rows]
    phase   = [String(r[1]) for r in rows]
    return O18_database("Vho, A., Lanari, P. & Rubatto, D. (2019). An internally-consistent database for oxygen isotope fractionation between minerals. Journal of Petrology 60, 2101-2129.",
                        "DBOxygen 2.0.3", phase, [String(r[2]) for r in rows],
                        num(3), num(4), num(5), num(6), num(7), num(8),
                        Dict(p => i for (i, p) in enumerate(phase)))
end

function _read_O18_mapping(path :: String)
    map = Dict{Tuple{String,String},Vector{Tuple{String,Float64}}}()
    for line in eachline(path)
        startswith(line, '#') && continue
        isempty(strip(line)) && continue
        parts = strip.(split(line, ','))
        parts[1] == "phase" && continue
        key = (String(parts[1]), String(parts[2]))
        push!(get!(map, key, Tuple{String,Float64}[]), (String(parts[3]), parse(Float64, parts[4])))
    end
    return O18_mapping(map)
end

const _O18_DB = let
    path = joinpath(@__DIR__, "O18_DBOxygen_2.0.3.csv")
    include_dependency(path)
    _read_O18_database(path)
end

const _O18_MAP = let
    path = joinpath(@__DIR__, "O18_mapping.csv")
    include_dependency(path)
    _read_O18_mapping(path)
end

"""
    get_O18_database(name="DBOxygen_2.0.3")

    Return the built-in oxygen isotope fractionation database.

    Parameters
    ----------
    name : String, optional
        Database identifier. Valid option: "DBOxygen_2.0.3" (Vho et al., 2019).

    Returns
    -------
    db : O18_database
        Quartz-referenced fractionation parameters for 153 minerals and water.
"""
function get_O18_database(name :: String = "DBOxygen_2.0.3")
    name == "DBOxygen_2.0.3" || error("Unknown oxygen isotope database \"$name\", valid option: \"DBOxygen_2.0.3\"")
    return _O18_DB
end

"""
    get_O18_mapping()

    Return the built-in mapping of MAGEMin phases and end-members onto DBOxygen phases.

    Returns
    -------
    mapping : O18_mapping
        Mapping read from `O18_mapping.csv`.
"""
get_O18_mapping() = _O18_MAP

"""
    O18_fractionation(db, phase, T_C; T_max=900.0)

    Fractionation between quartz and a database phase. Above `T_max` the fractionation is continued as `1/T²` from its
    value at `T_max`, the high-temperature limit of isotope exchange theory (Bottinga & Javoy, 1973).

    Parameters
    ----------
    db : O18_database
        Fractionation database.
    phase : String
        Database phase name (e.g. "Water", "Almandine").
    T_C : Float64
        Temperature [°C].
    T_max : Float64, optional
        Upper calibration temperature [°C] (default: 900, use `Inf` to extrapolate the polynomial).

    Returns
    -------
    Δ : Float64
        `1000 ln α(Qz-phase)` [‰].
"""
function O18_fractionation(db :: O18_database, phase :: String, T_C :: Float64; T_max :: Float64 = 900.0)
    i = db.id[phase]
    return _O18_delta(db.A[i], db.B[i], db.C[i], T_C + 273.15, T_max + 273.15)
end

function _O18_delta(A, B, C, TK, TK_max = Inf)
    TK <= TK_max && return A * 1e6 / TK^2 + B * 1e3 / TK + C
    return (A * 1e6 / TK_max^2 + B * 1e3 / TK_max + C) * (TK_max / TK)^2
end

function _O18_oxygen_mol(comp :: AbstractVector{Float64}, oxides :: Vector{String})
    return sum(comp[j] * get(_O18_nu_O, oxides[j], 0.0) for j in eachindex(oxides))
end

function _O18_oxygen_wt(comp_wt :: AbstractVector{Float64}, oxides :: Vector{String})
    return sum(comp_wt[j] * get(_O18_nu_O, oxides[j], 0.0) / get_molar_mass(oxides[j]) for j in eachindex(oxides))
end

function _O18_oxide_dict(comp :: AbstractVector{Float64}, oxides :: Vector{String})
    d = Dict{String,Float64}(ox => 0.0 for ox in ("SiO2", "TiO2", "Al2O3", "Cr2O3", "Fe2O3", "FeO", "MnO", "MgO", "CaO", "Na2O", "K2O", "H2O", "CO2"))
    for j in eachindex(oxides)
        ox = oxides[j]
        if ox == "O"
            continue
        elseif haskey(d, ox)
            d[ox] += comp[j]
        end
    end
    jO = findfirst(==("O"), oxides)
    if !isnothing(jO)
        fe3         = min(comp[jO], d["FeO"] / 2.0)
        d["Fe2O3"] += fe3
        d["FeO"]   -= 2.0 * fe3
    end
    return d
end

function _O18_mix(parts :: Vector{Tuple{String,Float64}}, db :: O18_database)
    A = B = C = 0.0
    for (p, w) in parts
        i  = db.id[p]
        A += w * db.A[i]
        B += w * db.B[i]
        C += w * db.C[i]
    end
    return A, B, C
end

"""
    ZZ03_coefficients(comp, oxides; db=get_O18_database())

    Quartz-phase fractionation parameters of a melt or rock from its oxide composition using the oxygen isotope
    index of Zhao & Zheng (2003): `1000 ln α(Qz-melt) = (1 - I)·11.83·10⁶/T²`, with `I` the oxygen-weighted sum of the
    oxide indices (their Table 4). H2O is assigned the DBOxygen water parameters, CO2 and oxides without index are ignored.

    Parameters
    ----------
    comp : Vector{Float64}
        Molar oxide composition.
    oxides : Vector{String}
        Oxide names matching `comp`.
    db : O18_database, optional
        Database used for the water parameters.

    Returns
    -------
    A, B, C : Float64
        Fractionation parameters.
    w_skip : Float64
        Fraction of the oxygen ignored.
"""
function ZZ03_coefficients( comp    :: AbstractVector{Float64},
                            oxides  :: Vector{String};
                            db      :: O18_database = get_O18_database())

    d       = _O18_oxide_dict(comp, oxides)
    O_dry   = 0.0
    I       = 0.0
    for (ox, idx) in _ZZ03_index
        o       = d[ox] * _O18_nu_O[ox]
        O_dry  += o
        I      += o * idx
    end
    O_h2o   = d["H2O"]
    O_skip  = _O18_oxygen_mol(comp, oxides) - O_dry - O_h2o
    O_tot   = O_dry + O_h2o
    O_tot > 0.0 || return 0.0, 0.0, 0.0, 1.0
    iw      = db.id["Water"]
    A_dry   = O_dry > 0.0 ? (1.0 - I / O_dry) * _ZZ03_A_qz : 0.0
    A       = (O_dry * A_dry + O_h2o * db.A[iw]) / O_tot
    B       = O_h2o * db.B[iw] / O_tot
    C       = O_h2o * db.C[iw] / O_tot
    return A, B, C, max(O_skip, 0.0) / (O_tot + max(O_skip, 0.0))
end

const _CIPW_minerals = (
    ("Q",   Dict("SiO2" => 1.0),                                        "Quartz"),
    ("C",   Dict("Al2O3" => 1.0),                                       "Corundum"),
    ("Or",  Dict("K2O" => 0.5, "Al2O3" => 0.5, "SiO2" => 3.0),          "K-Feldspar"),
    ("Ab",  Dict("Na2O" => 0.5, "Al2O3" => 0.5, "SiO2" => 3.0),         "Albite"),
    ("An",  Dict("CaO" => 1.0, "Al2O3" => 1.0, "SiO2" => 2.0),          "Anorthite"),
    ("Lc",  Dict("K2O" => 0.5, "Al2O3" => 0.5, "SiO2" => 2.0),          "Leucite"),
    ("Ne",  Dict("Na2O" => 0.5, "Al2O3" => 0.5, "SiO2" => 1.0),         "Nepheline"),
    ("Ac",  Dict("Na2O" => 0.5, "Fe2O3" => 0.5, "SiO2" => 2.0),         "Acmite"),
    ("Ns",  Dict("Na2O" => 1.0, "SiO2" => 1.0),                         ""),
    ("Di",  Dict("CaO" => 1.0, "MgO" => 1.0, "SiO2" => 2.0),            "Diopside"),
    ("Hd",  Dict("CaO" => 1.0, "FeO" => 1.0, "SiO2" => 2.0),            "Hedenbergite"),
    ("Wo",  Dict("CaO" => 1.0, "SiO2" => 1.0),                          "Wollastonite"),
    ("En",  Dict("MgO" => 1.0, "SiO2" => 1.0),                          "Enstatite"),
    ("Fs",  Dict("FeO" => 1.0, "SiO2" => 1.0),                          "Ferrosilite"),
    ("Fo",  Dict("MgO" => 2.0, "SiO2" => 1.0),                          "Forsterite"),
    ("Fa",  Dict("FeO" => 2.0, "SiO2" => 1.0),                          "Fayalite"),
    ("Mt",  Dict("FeO" => 1.0, "Fe2O3" => 1.0),                         "Magnetite"),
    ("Hm",  Dict("Fe2O3" => 1.0),                                       "Hematite"),
    ("Il",  Dict("FeO" => 1.0, "TiO2" => 1.0),                          "Ilmenite"),
    ("Ru",  Dict("TiO2" => 1.0),                                        "Rutile"),
    ("Chr", Dict("FeO" => 1.0, "Cr2O3" => 1.0),                         "Chromite"),
)

"""
    cipw_norm(comp, oxides)

    CIPW norm of an anhydrous oxide composition (MnO is added to FeO, Fe3+ is read from Fe2O3 or from the "O" component).
    Normative minerals are returned per formula unit, Mg and Fe pyroxene and olivine components are kept separate
    (Di/Hd, En/Fs, Fo/Fa). Silica undersaturation is resolved by converting hypersthene to olivine, then albite to
    nepheline and orthoclase to leucite.

    Parameters
    ----------
    comp : Vector{Float64}
        Molar oxide composition.
    oxides : Vector{String}
        Oxide names matching `comp`.

    Returns
    -------
    norm : Dict{String,Float64}
        Moles of each normative mineral (formula units: Q SiO2, C Al2O3, Or KAlSi3O8, Ab NaAlSi3O8, An CaAl2Si2O8,
        Lc KAlSi2O6, Ne NaAlSiO4, Ac NaFeSi2O6, Ns Na2SiO3, Di, Hd, Wo CaSiO3, En, Fs, Fo, Fa, Mt, Hm, Il, Ru, Chr FeCr2O4).
    si_deficit : Float64
        Moles of SiO2 that could not be supplied (0 for a complete norm).
"""
function cipw_norm( comp    :: AbstractVector{Float64},
                    oxides  :: Vector{String})

    d   = _O18_oxide_dict(comp, oxides)
    n   = Dict{String,Float64}(m[1] => 0.0 for m in _CIPW_minerals)
    si, ti, al, cr  = d["SiO2"], d["TiO2"], d["Al2O3"], d["Cr2O3"]
    fe3, fe, mg, ca = d["Fe2O3"], d["FeO"] + d["MnO"], d["MgO"], d["CaO"]
    na, k           = d["Na2O"], d["K2O"]

    x       = min(cr, fe);      n["Chr"] = x;   fe -= x
    x       = min(ti, fe);      n["Il"]  = x;   fe -= x;    ti -= x
    n["Ru"] = ti
    x       = min(k, al);       n["Or"]  = 2x;  al -= x
    x       = min(na, al);      n["Ab"]  = 2x;  al -= x;    na -= x
    x       = min(na, fe3);     n["Ac"]  = 2x;  fe3 -= x;   na -= x
    n["Ns"] = na
    x       = min(ca, al);      n["An"]  = x;   ca -= x;    al -= x
    n["C"]  = al
    x       = min(fe3, fe);     n["Mt"]  = x;   fe3 -= x;   fe -= x
    n["Hm"] = fe3
    m       = mg + fe
    x       = min(ca, m)
    xmg     = m > 0.0 ? mg / m : 0.0
    n["Di"] = x * xmg;          n["Hd"]  = x * (1.0 - xmg)
    ca     -= x;                mg      -= n["Di"];         fe -= n["Hd"]
    n["Wo"] = ca
    n["En"] = mg;               n["Fs"]  = fe

    si_unit = Dict("Or" => 3.0, "Ab" => 3.0, "An" => 2.0, "Ac" => 2.0, "Ns" => 1.0, "Di" => 2.0, "Hd" => 2.0, "Wo" => 1.0, "En" => 1.0, "Fs" => 1.0)
    need    = sum(n[s] * v for (s, v) in si_unit)
    if si >= need
        n["Q"] = si - need
        return n, 0.0
    end

    D   = need - si
    hy  = n["En"] + n["Fs"]
    xmg = hy > 0.0 ? n["En"] / hy : 0.0
    h   = min(hy, 2.0 * D)
    n["En"] -= xmg * h;         n["Fs"] -= (1.0 - xmg) * h
    n["Fo"]  = xmg * h / 2.0;   n["Fa"]  = (1.0 - xmg) * h / 2.0
    D       -= h / 2.0

    x        = min(n["Ab"], D / 2.0)
    n["Ab"] -= x;               n["Ne"] += x;               D -= 2.0 * x
    x        = min(n["Or"], D)
    n["Or"] -= x;               n["Lc"] += x;               D -= x

    return n, max(D, 0.0)
end

"""
    norm_coefficients(comp, oxides; db=get_O18_database())

    Quartz-phase fractionation parameters of a melt treated as a mixture of its CIPW normative minerals (Bindeman, 2008),
    each normative mineral weighted by its share of the oxygen. Dissolved H2O is assigned the DBOxygen water parameters.

    Parameters
    ----------
    comp : Vector{Float64}
        Molar oxide composition.
    oxides : Vector{String}
        Oxide names matching `comp`.
    db : O18_database, optional
        Fractionation database.

    Returns
    -------
    A, B, C : Float64
        Fractionation parameters.
    w_skip : Float64
        Fraction of the oxygen that could not be assigned (Na2SiO3, silica deficit, CO2).
"""
function norm_coefficients( comp    :: AbstractVector{Float64},
                            oxides  :: Vector{String};
                            db      :: O18_database = get_O18_database())

    n, _    = cipw_norm(comp, oxides)
    d       = _O18_oxide_dict(comp, oxides)
    A = B = C = 0.0
    O_map   = 0.0
    for (name, stoich, dbox) in _CIPW_minerals
        isempty(dbox) && continue
        o       = n[name] * sum(v * _O18_nu_O[ox] for (ox, v) in stoich)
        o > 0.0 || continue
        i       = db.id[dbox]
        A      += o * db.A[i]
        B      += o * db.B[i]
        C      += o * db.C[i]
        O_map  += o
    end
    O_h2o   = d["H2O"]
    iw      = db.id["Water"]
    A      += O_h2o * db.A[iw]
    B      += O_h2o * db.B[iw]
    C      += O_h2o * db.C[iw]
    O_map  += O_h2o
    O_tot   = _O18_oxygen_mol(comp, oxides)
    O_map > 0.0 || return 0.0, 0.0, 0.0, 1.0
    return A / O_map, B / O_map, C / O_map, max(O_tot - O_map, 0.0) / O_tot
end

function _O18_em_coefficients(  ph      :: String,
                                em      :: String,
                                emComp  :: AbstractVector{Float64},
                                oxides  :: Vector{String},
                                db      :: O18_database,
                                mapping :: O18_mapping,
                                unmapped:: Symbol)

    parts = get(mapping.map, (ph, em), nothing)
    isnothing(parts) && (parts = get(mapping.map, ("*", em), nothing))
    if !isnothing(parts)
        A, B, C = _O18_mix(parts, db)
        return A, B, C, true
    end
    if unmapped == :index
        A, B, C, _ = ZZ03_coefficients(emComp, oxides; db)
        return A, B, C, true
    end
    return NaN, NaN, NaN, false
end

function _O18_phase(out, i, db, mapping, melt_model, unmapped, TK, TK_max, flags)
    oxides  = out.oxides
    name    = out.ph[i]
    if out.ph_type[i] == 1
        ss      = out.SS_vec[out.ph_id[i] + 1]
        comp    = ss.Comp
        if name == "liq" && melt_model != "liq_em"
            if melt_model == "norm"
                A, B, C, w_skip = norm_coefficients(comp, oxides; db)
            elseif melt_model == "ZZ03"
                A, B, C, w_skip = ZZ03_coefficients(comp, oxides; db)
            else
                error("Unknown melt model \"$melt_model\", valid options: \"norm\", \"ZZ03\", \"liq_em\"")
            end
            w_skip > 0.01 && push!(flags, "liq: $(round(100w_skip, digits=1)) % of the melt oxygen not assigned")
            δem = fill(NaN, length(ss.emNames))
            return A, B, C, 1.0 - w_skip, ss.emNames, δem
        end
        n_em    = length(ss.emNames)
        v       = [ss.emFrac[k] * _O18_oxygen_mol(ss.emComp[k], oxides) for k in 1:n_em]
        v_tot   = sum(v)
        A = B = C = 0.0
        v_map   = 0.0
        Δem     = fill(NaN, n_em)
        missing = String[]
        for k in 1:n_em
            v[k] == 0.0 && continue
            Ak, Bk, Ck, ok = _O18_em_coefficients(name, ss.emNames[k], ss.emComp[k], oxides, db, mapping, unmapped)
            if ok
                A      += v[k] * Ak
                B      += v[k] * Bk
                C      += v[k] * Ck
                v_map  += v[k]
                Δem[k]  = _O18_delta(Ak, Bk, Ck, TK, TK_max)
            elseif abs(v[k] / v_tot) > 1e-3
                push!(missing, ss.emNames[k])
            end
        end
        isempty(missing) || push!(flags, "$name: end-member(s) $(join(missing, ", ")) not mapped")
        if v_tot <= 0.0 || v_map <= 0.0
            return NaN, NaN, NaN, 0.0, ss.emNames, Δem
        end
        minimum(v ./ v_tot) < -0.05 && push!(flags, "$name: negative end-member oxygen fraction")
        return A / v_map, B / v_map, C / v_map, v_map / v_tot, ss.emNames, Δem
    else
        pp      = out.PP_vec[out.ph_id[i] + 1]
        parts   = get(mapping.map, (name, ""), nothing)
        if !isnothing(parts)
            A, B, C = _O18_mix(parts, db)
            return A, B, C, 1.0, String[], Float64[]
        elseif unmapped == :index && _O18_oxygen_mol(pp.Comp, oxides) > 0.0
            A, B, C, _ = ZZ03_coefficients(pp.Comp, oxides; db)
            return A, B, C, 1.0, String[], Float64[]
        end
        return NaN, NaN, NaN, 0.0, String[], Float64[]
    end
end

"""
    O18_mass_balance(Δ, nO, δ; anchor=0)

    Closed-form oxygen isotope mass balance of an assemblage in isotopic equilibrium. With `δ_i = δ_qz - Δ_i` and
    `Σ nO_i δ_i = δ_sys Σ nO_i`, the virtual quartz value is `δ_qz = δ_sys + Σ f_i Δ_i`.

    Parameters
    ----------
    Δ : Vector{Float64}
        `1000 ln α(Qz-i)` of each phase [‰]; NaN entries are ignored.
    nO : Vector{Float64}
        Moles of oxygen held by each phase (any consistent unit).
    δ : Float64
        δ18O of the system [‰], or of phase `anchor` when `anchor > 0`.
    anchor : Int64, optional
        Index of the phase whose δ18O is `δ` (0: `δ` is the system value).

    Returns
    -------
    δ_sys : Float64
        δ18O of the system [‰].
    δ_qz : Float64
        δ18O of the virtual quartz [‰].
    δ_ph : Vector{Float64}
        δ18O of each phase [‰].
"""
function O18_mass_balance(  Δ       :: Vector{Float64},
                            nO      :: Vector{Float64},
                            δ       :: Float64;
                            anchor  :: Int64 = 0)

    use     = [!isnan(Δ[i]) && nO[i] > 0.0 for i in eachindex(Δ)]
    ΣfΔ     = sum(nO[use] .* Δ[use]) / sum(nO[use])
    if anchor == 0
        δ_sys   = δ
        δ_qz    = δ + ΣfΔ
    else
        δ_qz    = δ + Δ[anchor]
        δ_sys   = δ_qz - ΣfΔ
    end
    return δ_sys, δ_qz, δ_qz .- Δ
end

"""
    O18_prediction(out, δ; db=get_O18_database(), mapping=get_O18_mapping(), anchor="", melt_model="norm", unmapped=:exclude, T_max=900.0)

    Equilibrium oxygen isotope composition of every stable phase of a MAGEMin point (Vho et al., 2019, Eq. 7-11 and 20).
    Each phase is written as `δ_i = δ_qz - 1000 ln α(Qz-i)`; the virtual quartz value `δ_qz` follows in closed form from
    oxygen mass balance `Σ f_i δ_i = δ_sys`, with `f_i` the fraction of the oxygen held by phase `i`. Solid solutions combine
    their end-member parameters weighted by end-member oxygen fractions. The pressure effect is neglected (< 0.2 ‰ at 1 GPa).
    Above `T_max` every fractionation is continued as `1/T²` from its value at `T_max`.

    Parameters
    ----------
    out : MAGEMin_C.gmin_struct{Float64, Int64}
        MAGEMin minimization output.
    δ : Float64
        δ18O of the bulk system [‰ VSMOW], or of the `anchor` phase when `anchor` is given.
    db : O18_database, optional
        Fractionation database (default: DBOxygen 2.0.3).
    mapping : O18_mapping, optional
        Phase and end-member mapping onto `db`.
    anchor : String, optional
        Name of a stable phase whose δ18O is `δ` (inverse mode); the bulk δ18O is then an output.
    melt_model : String, optional
        "norm" (CIPW normative minerals, default), "ZZ03" (Zhao & Zheng, 2003 oxide index) or "liq_em"
        (melt model end-members mapped in `mapping`).
    unmapped : Symbol, optional
        `:exclude` (default) leaves unmapped phases and end-members out of the oxygen budget; `:index` assigns them
        parameters from their oxide index (Zhao & Zheng, 2003).
    T_max : Float64, optional
        Upper calibration temperature of `db` [°C] (default: 900, use `Inf` to extrapolate the polynomials).

    Returns
    -------
    o18 : out_O18
        δ18O of phases, end-members and subsystems with oxygen fractions and flags.
"""
function O18_prediction(out         :: MAGEMin_C.gmin_struct{Float64, Int64},
                        δ           :: Float64;
                        db          :: O18_database     = get_O18_database(),
                        mapping     :: O18_mapping      = get_O18_mapping(),
                        anchor      :: String           = "",
                        melt_model  :: String           = "norm",
                        unmapped    :: Symbol           = :exclude,
                        T_max       :: Float64          = 900.0)

    unmapped in (:exclude, :index) || error("Unknown unmapped option :$unmapped, valid options: :exclude, :index")

    TK      = out.T_C + 273.15
    TK_max  = T_max + 273.15
    flags   = String[]
    out.T_C > T_max && push!(flags, "T > $(T_max) °C: fractionations continued as 1/T²")
    out.T_C > 900.0 && isinf(T_max) && push!(flags, "T > 900 °C: DBOxygen polynomials extrapolated")
    out.T_C < 200.0 && push!(flags, "T < 200 °C: outside the recommended DBOxygen range")

    n_ph    = length(out.ph)
    A       = fill(NaN, n_ph)
    B       = fill(NaN, n_ph)
    C       = fill(NaN, n_ph)
    nO      = zeros(n_ph)
    nO_all  = zeros(n_ph)
    em      = Vector{Vector{String}}(undef, n_ph)
    Δem     = Vector{Vector{Float64}}(undef, n_ph)

    for i in 1:n_ph
        comp_wt     = out.ph_type[i] == 1 ? out.SS_vec[out.ph_id[i] + 1].Comp_wt : out.PP_vec[out.ph_id[i] + 1].Comp_wt
        nO_all[i]   = 100.0 * out.ph_frac_wt[i] * _O18_oxygen_wt(comp_wt, out.oxides)
        A[i], B[i], C[i], w_map, em[i], Δem[i] = _O18_phase(out, i, db, mapping, melt_model, unmapped, TK, TK_max, flags)
        nO[i]       = isnan(A[i]) ? 0.0 : nO_all[i] * w_map
        if isnan(A[i]) && nO_all[i] > 0.0
            push!(flags, "$(out.ph[i]): not mapped, excluded")
        end
    end

    N_O     = sum(nO)
    N_O > 0.0 || error("O18_prediction: no oxygen-bearing phase could be mapped")
    f_O     = nO ./ N_O
    Δ       = [isnan(A[i]) ? NaN : _O18_delta(A[i], B[i], C[i], TK, TK_max) for i in 1:n_ph]

    ia      = 0
    if !isempty(anchor)
        ia  = something(findfirst(==(anchor), out.ph), 0)
        ia == 0         && error("O18_prediction: anchor phase \"$anchor\" is not stable")
        isnan(Δ[ia])    && error("O18_prediction: anchor phase \"$anchor\" is not mapped")
    end
    δ_sys, δ_qz, δ_ph = O18_mass_balance(Δ, nO, δ; anchor = ia)
    δ_em    = [δ_qz .- d for d in Δem]

    function subsystem(sel)
        idx = [i for i in 1:n_ph if sel(out.ph[i]) && nO[i] > 0.0]
        isempty(idx) && return NaN
        return sum(nO[i] * δ_ph[i] for i in idx) / sum(nO[i] for i in idx)
    end
    is_liq(p)   = p == "liq"
    is_fl(p)    = p in _O18_fluid_names
    δ_liq   = subsystem(is_liq)
    δ_fl    = subsystem(is_fl)
    δ_sol   = subsystem(p -> !is_liq(p) && !is_fl(p))

    f_unmap = (sum(nO_all) - N_O) / sum(nO_all)

    return out_O18( out.T_C, δ_sys, δ_qz, N_O, copy(out.ph), f_O, A, B, C, δ_ph, em, δ_em,
                    δ_liq, δ_fl, δ_sol, f_unmap, unique(flags))
end

"""
    Δ18O(o18, i, j)

    Oxygen isotope fractionation between two phases of an `O18_prediction` result, independent of the bulk δ18O.

    Parameters
    ----------
    o18 : out_O18
        Result of `O18_prediction`.
    i, j : String
        Phase names (first occurrence is used).

    Returns
    -------
    Δ : Float64
        `δ18O_i - δ18O_j` [‰], NaN when a phase is absent or not mapped.
"""
function Δ18O(o18 :: out_O18, i :: String, j :: String)
    a = findfirst(==(i), o18.ph)
    b = findfirst(==(j), o18.ph)
    (isnothing(a) || isnothing(b)) && return NaN
    return o18.δ_ph[a] - o18.δ_ph[b]
end
