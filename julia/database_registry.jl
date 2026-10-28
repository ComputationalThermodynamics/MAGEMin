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

"""
    db_registry

    One row of [`MAGEMIN_DATABASES`], the single hand-maintained description of a
    MAGEMin thermodynamic database.

    Fields
    ------
    db_name : String
        Database acronym, as passed to `Initialize_MAGEMin`/`init_MAGEMin` and to the
        CLI `--db` flag (e.g. "mp", "ig", "sb21").
    research_group : String
        Research group the database belongs to, as expected by `gv.research_group`:
        "tc" (THERMOCALC), "sb" (Stixrude & Lithgow-Bertelloni), "gh" (Ghiorso/MELTS)
        or "br" (Berman).
    EM_database : Int64
        Numeric database code written into `gv.EM_database`. Unique within a research
        group only - the same code means different databases in different groups.
    db_info : String
        Human-readable database name and reference, used by `print_phase_info` and
        written into the metadata header of every CSV export.
    db_dataset : Int64
        End-member dataset the library falls back to when none is requested
        (`gv.EM_dataset == -1`). `-1` where the research group has no dataset choice.
    dataset_opt : NTuple{5,Int64}
        End-member datasets that may be requested for this database.
"""
struct db_registry
    db_name         :: String
    research_group  :: String
    EM_database     :: Int64
    db_info         :: String
    db_dataset      :: Int64
    dataset_opt     :: NTuple{5,Int64}
end

const GH_ENABLED   = false
const IGAW_ENABLED = false

"""
    MAGEMIN_DATABASES

    Registry of every thermodynamic database the library implements - the single place
    where a database acronym, its research group, its `gv.EM_database` code, its
    human-readable name and its dataset options are written down.

    Everything else about a database (oxide list, solution- and pure-phase lists,
    end-member / compositional-variable / site-fraction names) is *not* listed here:
    it is harvested straight from the C library by `gen/generate_db_infos.jl` into
    `julia/db_infos_generated.jl`, and `test/test_db_infos.jl` fails the build if the
    two ever disagree.

    The row order is significant: `db_infos_generated.jl` is emitted in this order.
"""
const MAGEMIN_DATABASES = filter(d -> (GH_ENABLED || d.research_group != "gh") && (IGAW_ENABLED || d.db_name != "igaw"), db_registry[
    db_registry("mp",     "tc",  0, "Metapelite (White et al., 2014)",                                                                              62, (62, 633, 634, 635, 636)),
    db_registry("mb",     "tc",  1, "Metabasite (Green et al., 2016)",                                                                              62, (62, 633, 634, 635, 636)),
    db_registry("mbe",    "tc", 11, "Metabasite extended (Green et al., 2016 with oamp from Diener et al., 2007 and ta from Rebay et al., 2022)",   62, (62, 633, 634, 635, 636)),
    db_registry("ig",     "tc",  2, "Igneous (Green et al., 2025, corrected after Holland et al., 2018)",                                          636, (62, 633, 634, 635, 636)),
    db_registry("igaw",   "tc", 10, "Igneous alkaline wet (Holland et al., 2026, in progress)",                636, (62, 633, 634, 635, 636)),
    db_registry("igad",   "tc",  3, "Igneous alkaline dry (Weller et al., 2024)",                                                                  636, (62, 633, 634, 635, 636)),
    db_registry("igd",    "tc", 22, "Igneous dry (Su et al., 2026, corrected after Tomlinson & Holland, 2021)",                                    634, (62, 633, 634, 635, 636)),
    db_registry("um",     "tc",  4, "Ultramafic (Evans & Frost., 2021)",                                                                           633, (62, 633, 634, 635, 636)),
    db_registry("ume",    "tc",  5, "Ultramafic extended (Evans & Frost., 2021 with pl, amp and aug from Green et al., 2016)",                     633, (62, 633, 634, 635, 636)),
    db_registry("mtl",    "tc",  6, "Mantle (Holland et al., 2013)",                                                                               633, (62, 633, 634, 635, 636)),
    db_registry("mpe",    "tc",  7, "Metapelite extended (White et al., 2014 with po from Evans & Frost., 2021, amp dio and aug from Green et al., 2016)", 62, (62, 633, 634, 635, 636)),
    db_registry("all",    "tc",  8, "Global solution dataset",                                                                                     636, (62, 633, 634, 635, 636)),
    db_registry("po",     "br",  0, "HP/LT (Pourteau et al., 2014)",                                                                                 1, ( 1,   1,   1,   1,   1)),
    db_registry("sb11",   "sb",  0, "Stixrude & Lithgow-Bertelloni (2011)",                                                                         -1, (-1,  -1,  -1,  -1,  -1)),
    db_registry("sb21",   "sb",  1, "Stixrude & Lithgow-Bertelloni (2021)",                                                                         -1, (-1,  -1,  -1,  -1,  -1)),
    db_registry("sb24",   "sb",  2, "Stixrude & Lithgow-Bertelloni (2024)",                                                                         -1, (-1,  -1,  -1,  -1,  -1)),
    db_registry("xMELTS", "gh",  0, "xMELTS dev",                                                                                                    1, ( 1,   1,   1,   1,   1)),
    db_registry("pMELTS", "gh",  2, "pMELTS 5.6.1",                                                                                                  1, ( 1,   1,   1,   1,   1)),
    db_registry("rMELTS", "gh",  1, "rMELTS 1.2.0",                                                                                                  1, ( 1,   1,   1,   1,   1)),
])

const _DB_REGISTRY_BY_NAME = Dict{String,db_registry}(d.db_name => d for d in MAGEMIN_DATABASES)
const _DB_REGISTRY_INDEX   = Dict{String,Int64}(d.db_name => i for (i,d) in enumerate(MAGEMIN_DATABASES))

"""
    get_db_list(; research_group=nothing)

    Acronyms of every implemented thermodynamic database, in registry order, optionally
    restricted to one research group ("tc", "sb", "gh" or "br").

    Examples
    --------
    ```julia
    get_db_list()                       # all enabled databases
    get_db_list(research_group = "sb")  # ["sb11", "sb21", "sb24"]
    ```
"""
function get_db_list(; research_group :: Union{Nothing,String} = nothing)
    isnothing(research_group) && return [d.db_name for d in MAGEMIN_DATABASES]
    return [d.db_name for d in MAGEMIN_DATABASES if d.research_group == research_group]
end

"""
    get_db(dtb)

    Registry entry ([`db_registry`](@ref)) of database `dtb`. Throws if the acronym is
    unknown - use [`is_db`](@ref) to test first.
"""
function get_db(dtb :: String)
    if !haskey(_DB_REGISTRY_BY_NAME, dtb)
        known = join(get_db_list(), ", ")
        error("get_db: unknown database \"$dtb\". Implemented databases: $known.")
    end
    return _DB_REGISTRY_BY_NAME[dtb]
end

"""
    is_db(dtb)

    Whether `dtb` is an implemented database acronym.
"""
is_db(dtb :: String) = haskey(_DB_REGISTRY_BY_NAME, dtb)

"""
    get_db_label(dtb)

    Human-readable name and reference of database `dtb` (`db_registry.db_info`).
"""
get_db_label(dtb :: String) = get_db(dtb).db_info

"""
    get_research_group(dtb)

    Research group ("tc", "sb", "gh" or "br") database `dtb` belongs to.
"""
get_research_group(dtb :: String) = get_db(dtb).research_group

"""
    has_dataset_choice(dtb)

    Whether the end-member dataset can be chosen for database `dtb`. Only the "tc"
    research group offers a choice (`available_TC_ds`); the others pin the dataset to
    the database itself, so a dataset selector should be hidden for them.
"""
has_dataset_choice(dtb :: String) = get_research_group(dtb) == "tc"
