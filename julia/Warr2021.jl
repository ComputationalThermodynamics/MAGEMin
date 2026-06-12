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

# Load Warr (2021) mapping at module init from the CSV file alongside this source.
# Reference: Warr L.N. (2021) IMA-CNMNC approved mineral symbols.
#            Mineralogical Magazine 85, 291-320. doi:10.1180/mgm.2021.43
const _warr_dict = let
    dict = Dict{String,String}()
    csv_path = joinpath(@__DIR__, "MAGEMin_Warr2021_mapping.csv")
    for line in eachline(csv_path)
        startswith(line, '#') && continue
        isempty(strip(line))  && continue
        parts = split(line, ',')
        length(parts) < 4    && continue
        name = strip(parts[1])
        name == "magemin_name" && continue   # header row
        sym  = strip(parts[4])
        # First occurrence wins — avoids collision from duplicate names
        haskey(dict, name) || (dict[name] = sym)
    end
    dict
end

"""
    get_Warr_name(name)

    Convert a MAGEMin mineral/endmember abbreviation to its IMA-CNMNC approved
    symbol following Warr (2021).

    Names with no official Warr equivalent (composite endmembers, model-specific
    components, buffer assemblages, activities) are returned as `name * "*"` so
    that downstream code always receives a plain `String` without special-casing.

    Parameters
    ----------
    name : String
        MAGEMin internal phase or endmember name (e.g. `"ab"`, `"phl"`, `"q4L"`).

    Returns
    -------
    symbol : String
        Warr (2021) symbol (e.g. `"Ab"`, `"Phl"`) or `name * "*"` when not listed.

    Reference
    ---------
    Warr L.N. (2021) IMA-CNMNC approved mineral symbols.
    Mineralogical Magazine 85, 291–320. doi:10.1180/mgm.2021.43
"""
function get_Warr_name(name::String)
    return get(_warr_dict, name, name)
end

"""
    get_Warr_names(names)

    Vectorised form of `get_Warr_name`: convert a vector of MAGEMin names to
    their Warr (2021) symbols in one call. Entries with no mapping are returned
    as `name * "*"`.

    Parameters
    ----------
    names : Vector{String}
        Vector of MAGEMin phase or endmember names.

    Returns
    -------
    symbols : Vector{String}
        Corresponding Warr (2021) symbols.
"""
function get_Warr_names(names::Vector{String})
    return get_Warr_name.(names)
end

"""
    get_mineral_name_Warr(db, ss, SS_vec)

    Return the Warr (2021) IMA-CNMNC symbol for the solvus-disambiguated name of
    a solution phase.  Internally calls `get_mineral_name` for solvus logic, then
    maps the result through `get_Warr_name`.

    Parameters
    ----------
    db     : String  — database identifier (e.g. `"ig"`, `"mp"`).
    ss     : String  — solution phase short name (e.g. `"fsp"`, `"amp"`).
    SS_vec : LibMAGEMin.SS_data — solution phase data structure.

    Returns
    -------
    symbol : String — Warr (2021) symbol, or the disambiguated name * "*"` if
             no official symbol is defined.
"""
function get_mineral_name_Warr(db, ss, SS_vec)
    return get_Warr_name(get_mineral_name(db, ss, SS_vec))
end
