#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
#
#   Regenerates the phase and end-member tables of database.md from the library, so they
#   cannot drift from the code the way the hand-written ones did. The solvus abbreviation
#   list (problem.md) and the Yak25 Kd table (database.md) are maintained by hand.
#
#   Usage:  julia --project=. tools/generate_docs_tables.jl [docs/src dir]
#
#   The block is located by its heading and the existing ':::tabs' fences, so the
#   markdown carries no generator markers.

using MAGEMin_C
using MAGEMin_C: MAGEMIN_DATABASES, SOLVUS_FAMILY, retrieve_solution_phase_information,
                 get_Warr_name, ALWAYS_ACTIVE_PP

const DOCS = length(ARGS) >= 1 ? ARGS[1] :
             normpath(joinpath(@__DIR__, "..", "..", "MAGEMin_C.jl_2.3.7", "docs", "src"))

# MELTS databases are beta and deliberately kept out of the documentation
const SKIP   = ("xMELTS", "pMELTS", "rMELTS")
const TITLES = Dict(
    "mp"   => "Metapelite (mp)",            "mb"   => "Metabasite (mb)",
    "mbe"  => "Metabasite extended (mbe)",  "ig"   => "Igneous (ig)",
    "igaw" => "Igneous Alkaline Wet (igaw)",
    "igd"  => "Igneous (igd)",              "igad" => "Igneous Alkaline (igad)",
    "um"   => "Ultramafic (um)",            "ume"  => "Ultramafic extended (ume)",
    "mtl"  => "Mantle (mtl)",               "mpe"  => "Metapelite extended (mpe)",
    "po"   => "HP/LT Metapelite (po)",      "sb11" => "Mantle SB11",
    "sb21" => "Mantle SB21",                "sb24" => "Mantle SB24",
    "all"  => "Global (all)")
const ORDER = ["mp","mb","mbe","ig","igaw","igd","igad","um","ume","mtl","mpe","po",
               "sb11","sb21","sb24","all"]

"""
    NOTES

    Hand-written commentary that belongs to a database's tab, carried through
    regeneration. Only the tables themselves come from the library.
"""
const NOTES = Dict(
    "mbe"  => "Extends `mb` with two additional solution phases.",
    "ume"  => "Extends `um` with plagioclase, amphibole, augite, spinel, carbonated fluid and carbonates.",
    "mpe"  => "Extends `mp` with phases from Green et al. (2016), Evans & Frost (2021), and Diener et al. (2007).",
    "po"   => "Berman (1988) formalism (`br` research group), not the Holland-Powell-family datasets used elsewhere in this table.",
    "sb11" => "Stixrude & Lithgow-Bertelloni (2011). End-member names follow the SLB internal convention.",
    "sb21" => "Stixrude & Lithgow-Bertelloni (2021). Extends SB11 with new-aluminium-phase (`nal`) and an extra end-member in `mw`.",
    "sb24" => "Stixrude & Lithgow-Bertelloni (2024). Expanded solid solutions throughout; new iron polymorphs and high-pressure phases as pure phases.")

const FOOTNOTES = Dict(
    "all"  => "\\* The two `fsp_H22` rows share the same internal citation tag despite differing endmember\n" *
              "counts (a known, pre-existing data labelling quirk carried over from the source databases, not\n" *
              "introduced by `all`) - distinguish them by endmember count/list in practice.")

const BUFFERS    = ["qfm","mw","qif","nno","hm","iw","cco"]
const ACTIVITIES = ["aH2O","aO2","aMgO","aFeO","aAl2O3","aTiO2"]

warr(n) = (w = get_Warr_name(n); w == n ? n : w)
dotjoin(v) = join(v, " · ")

function phase_tabs()
    io = IOBuffer()
    for db in ORDER
        db in SKIP && continue
        inf = retrieve_solution_phase_information(db)
        println(io, "== ", TITLES[db], "\n")
        haskey(NOTES, db) && println(io, NOTES[db], "\n")
        println(io, "| Phase | Warr (2021) | Model | em | End-members |")
        println(io, "|---|---|---|:---:|---|")
        for s in inf.data_ss
            ems   = s.ss_em[2:end]
            model = isempty(s.ss_fName) ? "-" : s.ss_fName
            emtxt = startswith(s.ss_name, "DEW") || startswith(s.ss_fName, "DEW") ?
                    "*$(length(ems)) aqueous species — see the DEW section below*" : dotjoin(ems)
            println(io, "| `", s.ss_name, "` | ", warr(s.ss_name), " | ", model, " | ", s.n_em, " | ", emtxt, " |")
        end
        println(io)
        buf = [p for p in inf.data_pp if p in BUFFERS]
        act = [p for p in inf.data_pp if p in ACTIVITIES]
        pur = [p for p in inf.data_pp if !(p in BUFFERS) && !(p in ACTIVITIES)]
        isempty(pur) || println(io, "**Pure phases:** ", dotjoin(pur), "\n")
        if !isempty(buf) || !isempty(act)
            parts = String[]
            isempty(buf) || push!(parts, "**Buffers:** " * dotjoin(buf))
            isempty(act) || push!(parts, "**Activities:** " * dotjoin(act))
            println(io, join(parts, " &nbsp;&nbsp; "), "\n")
        end
        haskey(FOOTNOTES, db) && println(io, FOOTNOTES[db], "\n")
    end
    return rstrip(String(take!(io)))
end

const _FULL = Dict{String,String}()
function get_full_name(n)
    if isempty(_FULL)
        for line in eachline(normpath(joinpath(@__DIR__, "..", "julia", "MAGEMin_Warr2021_mapping.csv")))
            (startswith(line, '#') || isempty(strip(line))) && continue
            p = split(line, ','); length(p) < 4 && continue
            k = strip(p[1]); k == "magemin_name" && continue
            haskey(_FULL, k) || (_FULL[k] = String(strip(p[3])))
        end
    end
    return get(_FULL, n, "-")
end

"""
    replace_tabs_block!(path, heading, body)

    Replace the contents of the `:::tabs` block that follows `heading` in `path`.

    Anchored on the heading and the existing tab fences rather than on marker comments,
    so nothing has to be added to the markdown: the fences are already part of the page.
"""
function replace_tabs_block!(path, heading, body)
    lines = readlines(path)

    h = findfirst(l -> startswith(l, heading), lines)
    isnothing(h) && error("heading '$heading' not found in $path")
    o = findnext(l -> strip(l) == ":::tabs", lines, h)
    isnothing(o) && error("no ':::tabs' after '$heading' in $path")
    c = findnext(l -> strip(l) == ":::", lines, o + 1)
    isnothing(c) && error("unterminated ':::tabs' after '$heading' in $path")

    open(path, "w") do io
        for l in lines[1:o];      println(io, l); end
        println(io)
        println(io, body)
        println(io)
        for l in lines[c:end];    println(io, l); end
    end
    println("  rewrote the ", strip(heading), " tables in ", basename(path))
end

println(" Generating documentation tables from the library...")
replace_tabs_block!(joinpath(DOCS, "database.md"), "## Phase and End-member Listing", phase_tabs())
println(" done.")
