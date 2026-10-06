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
    SOLVUS_FAMILY

    For each database and solution phase that straddles a solvus, every mineral name
    [`get_mineral_name`](@ref) can return for it, including the parent name where that is
    still reachable.

    This is the declared counterpart of `get_mineral_name`'s branches: the function picks
    *which* name applies at a given composition, this table states which names are
    *possible*. Code that needs the whole family up front - building a phase list, or
    matching a phase along a path where `name_solvus=true` has already renamed it - reads
    this instead of keeping its own copy.

    Keys are `(database, solution phase)`, the phase spelled as that database spells it, so
    the "all" database uses citation-tagged names. Only pairs the database actually has are
    listed. `test/test_db_infos.jl` sweeps `get_mineral_name` over random compositional
    variables and fails if it ever returns a name this table does not list.
"""
const SOLVUS_FAMILY = Dict{Tuple{String,String},Vector{String}}(
    # ig, igad, igd
    ("ig",   "spl")         => ["cm", "mgt", "spl", "usp"],
    ("ig",   "fsp")         => ["afs", "pl"],
    ("ig",   "mu")          => ["mu", "pat"],
    ("ig",   "amp")         => ["act", "amp", "cumm", "gl", "tr"],
    ("ig",   "ilm")         => ["hem", "ilm"],
    ("ig",   "cpx")         => ["Na-cpx", "cpx", "pig"],
    ("igad", "spl")         => ["cm", "mgt", "spl", "usp"],
    ("igad", "fsp")         => ["afs", "pl"],
    ("igad", "ilm")         => ["hem", "ilm"],
    ("igad", "nph")         => ["K-nph", "nph"],
    ("igad", "cpx")         => ["Na-cpx", "cpx", "pig"],
    ("igd",  "spl")         => ["cm", "mgt", "spl", "usp"],
    ("igd",  "fsp")         => ["afs", "pl"],
    ("igd",  "ilm")         => ["hem", "ilm"],
    ("igd",  "cpx")         => ["Na-cpx", "cpx", "pig"],
    # mp, mpe, mb, ume, mbe
    ("mp",   "sp")          => ["smt", "sp"],
    ("mp",   "fsp")         => ["afs", "pl"],
    ("mp",   "mu")          => ["mu", "pat"],
    ("mp",   "ilmm")        => ["hemm", "ilmm"],
    ("mp",   "ilm")         => ["hem", "ilm"],
    ("mpe",  "sp")          => ["smt", "sp"],
    ("mpe",  "fsp")         => ["afs", "pl"],
    ("mpe",  "mu")          => ["mu", "pat"],
    ("mpe",  "amp")         => ["act", "amp", "cumm", "gl", "tr"],
    ("mpe",  "ilmm")        => ["hemm", "ilmm"],
    ("mpe",  "ilm")         => ["hem", "ilm"],
    ("mpe",  "dio")         => ["dio", "jd", "omph"],
    ("mpe",  "occm")        => ["ank", "cc", "mag", "sid"],
    ("mpe",  "oamp")        => ["anth", "ged"],
    ("mb",   "sp")          => ["smt", "sp"],
    ("mb",   "spl")         => ["cm", "mgt", "spl"],
    ("mb",   "fsp")         => ["afs", "pl"],
    ("mb",   "mu")          => ["mu", "pat"],
    ("mb",   "amp")         => ["act", "amp", "cumm", "gl", "tr"],
    ("mb",   "ilmm")        => ["hemm", "ilmm"],
    ("mb",   "ilm")         => ["hem", "ilm"],
    ("mb",   "dio")         => ["dio", "jd", "omph"],
    ("ume",  "spl")         => ["cm", "mgt", "spl"],
    ("ume",  "amp")         => ["act", "amp", "cumm", "gl", "tr"],
    ("ume",  "occm")        => ["ank", "cc", "mag", "sid"],
    ("mbe",  "sp")          => ["smt", "sp"],
    ("mbe",  "spl")         => ["cm", "mgt", "spl"],
    ("mbe",  "fsp")         => ["afs", "pl"],
    ("mbe",  "mu")          => ["mu", "pat"],
    ("mbe",  "amp")         => ["act", "amp", "cumm", "gl", "tr"],
    ("mbe",  "ilmm")        => ["hemm", "ilmm"],
    ("mbe",  "ilm")         => ["hem", "ilm"],
    ("mbe",  "dio")         => ["dio", "jd", "omph"],
    ("mbe",  "oamp")        => ["anth", "ged"],
    # all
    ("all",  "fsp_H22")     => ["afs", "pl"],
    ("all",  "fsp_H22op")   => ["afs", "pl"],
    ("all",  "spl_T21")     => ["cm", "mgt", "spl"],
    ("all",  "sp_W02")      => ["smt", "sp"],
    ("all",  "ilm_T21")     => ["hem", "ilm"],
    ("all",  "ilm_W24")     => ["hem", "ilm"],
    ("all",  "ilm_W00")     => ["hem", "ilm"],
    ("all",  "ilmm_W14")    => ["hemm", "ilmm"],
    ("all",  "ilmmn_W14")   => ["hemm", "ilmm"],
    ("all",  "amp_G16")     => ["act", "amp", "cumm", "gl", "tr"],
    ("all",  "mu_W14")      => ["mu", "pat"],
    ("all",  "cpx_T21")     => ["Na-cpx", "cpx", "pig"],
    ("all",  "cpx_W24")     => ["Na-cpx", "cpx", "pig"],
    ("all",  "nph_W24")     => ["K-nph", "nph"],
    ("all",  "dio_G16")     => ["dio", "jd", "omph"],
    ("all",  "occm_F11")    => ["ank", "cc", "mag", "sid"],
    ("all",  "oamp_D07")    => ["anth", "ged"],
)

"""
    get_mineral_name(db, ss, SS_vec)

    Return a mineralogically meaningful name for a solution phase based on its compositional variables (solvus disambiguation).

    For solution phases that straddle a solvus (e.g., feldspar → plagioclase or alkali feldspar; spinel → spinel, magnetite, or ulvöspinel), the returned name reflects the dominant endmember rather than the generic solution phase label.

    Parameters
    ----------
    db : String
        Database identifier (e.g., "ig", "igad", "mp", "mpe", "mb", "ume", "mbe").
    ss : String
        Solution phase short name (e.g., "fsp", "spl", "amp", "ilm").
    SS_vec : LibMAGEMin.SS_data
        Solution phase data structure containing `compVariables`.

    Returns
    -------
    mineral_name : String
        Disambiguated mineral name (e.g., "pl", "afs", "mgt", "ilm", "hem").
"""
function get_mineral_name(db, ss, SS_vec)

    # @warn "Breaking changes in v2.1.3 by disambiguation of solvus names: 'spl': 'sp' > 'spl'; 'sp': 'mt' > 'smt'."
    mineral_name = ss
   
    if db == "ig" || db == "igad"  || db == "igd"
        x = SS_vec.compVariables
        if ss == "spl"
            if x[3] - 0.5 > 0.0;        mineral_name = "cm";
            elseif x[4] - 0.5 > 0.0;    mineral_name = "usp";
            elseif x[2] - 0.5 > 0.0;    mineral_name = "mgt";
            else                        mineral_name = "spl";    end
        elseif ss == "fsp"
            if x[2] - 0.5 > 0.0;       mineral_name = "afs";
            else                        mineral_name = "pl";    end
        elseif ss == "mu"
            if x[4] - 0.5 > 0.0;        mineral_name = "pat";
            else                        mineral_name = "mu";    end
        elseif ss == "amp"
            if x[3] - 0.5 > 0.0;        mineral_name = "gl";
            elseif -x[3] -x[4] + 0.2 > 0.0;   mineral_name = "act";
            else
                if x[6] < 0.1;          mineral_name = "cumm"; 
                elseif -1/2*x[4]+x[6]-x[7]-x[8]-x[2]+x[3]>0.5;      mineral_name = "tr";       
                else                    mineral_name = "amp";    end
            end  
        elseif ss == "ilm"
            if -x[1] + 0.5 > 0.0;       mineral_name = "hem";
            else                        mineral_name = "ilm";   end 
        elseif ss == "nph"
            if x[2] - 0.5 > 0.0;       mineral_name = "K-nph";
            else                        mineral_name = "nph";   end 
        elseif ss == "cpx"
            if x[3] - 0.6 > 0.0;        mineral_name = "pig";
            elseif x[4] - 0.5 > 0.0;    mineral_name = "Na-cpx";
            else                        mineral_name = "cpx";   end 
        end

    elseif db == "mp" || db == "mpe" || db == "mb" || db == "ume" || db == "mbe"
        x = SS_vec.compVariables
        if ss == "sp"
            if x[2] - 0.5 > 0.0;        mineral_name = "sp";
            else                        mineral_name = "smt";    end
        elseif ss == "spl" && db == "ume"
            if x[3] - 0.5 > 0.0;        mineral_name = "cm";
            elseif x[2] - 0.5 > 0.0;    mineral_name = "mgt";
            else                        mineral_name = "spl";    end
        elseif ss == "fsp"
            if x[2] - 0.5 > 0.0;       mineral_name = "afs";
            else                        mineral_name = "pl";    end
        elseif ss == "mu"
            if x[4] - 0.5 > 0.0;        mineral_name = "pat";
            else                        mineral_name = "mu";    end
        elseif ss == "amp"
            if x[3] - 0.5 > 0.0;        mineral_name = "gl";
            elseif -x[3]-x[4]+0.2>0.0;  mineral_name = "act";
            else
                if x[6] < 0.1;          mineral_name = "cumm"; 
                elseif -1/2*x[4]+x[6]-x[7]-x[8]-x[2]+x[3]>0.5;      mineral_name = "tr";     
                else                    mineral_name = "amp";    end
            end  
        elseif ss == "ilmm"
            if x[1] - 0.5 > 0.0;        mineral_name = "ilmm";
            else                        mineral_name = "hemm";   end 
        elseif ss == "ilm"
            if 1.0 - x[1] > 0.5;        mineral_name = "hem";
            else                        mineral_name = "ilm";   end 
        elseif ss == "dio"
            if x[2] > 0.0 && x[2] <= 0.3;       mineral_name = "dio";
            elseif x[2] > 0.3 && x[2] <= 0.7;   mineral_name = "omph";
            else                                mineral_name = "jd";   end 
        elseif ss == "occm"
            if x[2] > 0.5;              mineral_name = "sid";
            elseif x[3] > 0.5;          mineral_name = "ank";  
            elseif x[1] > 0.25 && x[3] < 0.01;         mineral_name = "mag";  
            else                        mineral_name = "cc";   end
        elseif ss == "oamp"
            if x[2] < 0.3;              mineral_name = "anth";  #compositional variable y
            else                        mineral_name = "ged";   end
        end

    elseif db == "all"
        x = SS_vec.compVariables
        if ss == "fsp_H22" || ss == "fsp_H22op"
            if x[2] - 0.5 > 0.0;       mineral_name = "afs";
            else                        mineral_name = "pl";    end
        elseif ss == "spl_T21"
            if x[3] - 0.5 > 0.0;        mineral_name = "cm";
            elseif x[2] - 0.5 > 0.0;    mineral_name = "mgt";
            else                        mineral_name = "spl";    end
        elseif ss == "sp_W02"
            if x[2] - 0.5 > 0.0;        mineral_name = "sp";
            else                        mineral_name = "smt";    end
        elseif ss == "ilm_T21"
            if -x[1] + 0.5 > 0.0;       mineral_name = "hem";
            else                        mineral_name = "ilm";   end
        elseif ss == "ilm_W24"
            if -x[1] + 0.5 > 0.0;       mineral_name = "hem";
            else                        mineral_name = "ilm";   end
        elseif ss == "ilm_W00"
            if 1.0 - x[1] > 0.5;        mineral_name = "hem";
            else                        mineral_name = "ilm";   end
        elseif ss == "ilmm_W14" || ss == "ilmmn_W14"
            if x[1] - 0.5 > 0.0;        mineral_name = "ilmm";
            else                        mineral_name = "hemm";   end
        elseif ss == "amp_G16"
            if x[3] - 0.5 > 0.0;        mineral_name = "gl";
            elseif -x[3]-x[4]+0.2>0.0;  mineral_name = "act";
            else
                if x[6] < 0.1;          mineral_name = "cumm";
                elseif -1/2*x[4]+x[6]-x[7]-x[8]-x[2]+x[3]>0.5;      mineral_name = "tr";
                else                    mineral_name = "amp";    end
            end
        elseif ss == "mu_W14"
            if x[4] - 0.5 > 0.0;        mineral_name = "pat";
            else                        mineral_name = "mu";    end
        elseif ss == "cpx_T21"
            if x[3] - 0.6 > 0.0;        mineral_name = "pig";
            elseif x[4] - 0.5 > 0.0;    mineral_name = "Na-cpx";
            else                        mineral_name = "cpx";   end
        elseif ss == "cpx_W24"
            if x[3] - 0.6 > 0.0;        mineral_name = "pig";
            elseif x[4] - 0.5 > 0.0;    mineral_name = "Na-cpx";
            else                        mineral_name = "cpx";   end
        elseif ss == "nph_W24"
            if x[2] - 0.5 > 0.0;       mineral_name = "K-nph";
            else                        mineral_name = "nph";   end
        elseif ss == "dio_G16"
            if x[2] > 0.0 && x[2] <= 0.3;       mineral_name = "dio";
            elseif x[2] > 0.3 && x[2] <= 0.7;   mineral_name = "omph";
            else                                mineral_name = "jd";   end
        elseif ss == "occm_F11"
            if x[2] > 0.5;              mineral_name = "sid";
            elseif x[3] > 0.5;          mineral_name = "ank";
            elseif x[1] > 0.25 && x[3] < 0.01;         mineral_name = "mag";
            else                        mineral_name = "cc";   end
        elseif ss == "oamp_D07"
            if x[2] < 0.3;              mineral_name = "anth";
            else                        mineral_name = "ged";   end
        elseif occursin("_", ss)
            mineral_name = split(ss, "_")[1]
        end

    end

    return mineral_name
end

"""
    get_ss_from_mineral(db, mrl, mbCpx, active_ss=String[])

    Return the solution phase name corresponding to a disambiguated mineral name (inverse of `get_mineral_name`).

    For `db == "all"`, some generic mineral names are ambiguous: several solution-phase models
    disambiguate to the same mineral name (e.g., `ilm_T21`, `ilm_W24`, and `ilm_W00` all disambiguate
    to "ilm"/"hem"; `cpx_T21` and `cpx_W24` both disambiguate to "cpx"/"pig"/"Na-cpx"). When `active_ss`
    is supplied (e.g., the list of solution phases actually available/used in the current database
    configuration or computation), the first matching candidate is returned; otherwise a fixed default
    is used for backward compatibility. For `db == "all"` mineral names with no dedicated solvus rule
    (e.g. "opx", "ol"), `active_ss` is searched for an entry whose prefix before "_" equals `mrl`, so the
    exact underlying solution-phase name (including its dataset tag) is recovered whenever available.

    Parameters
    ----------
    db : String
        Database identifier (e.g., "ig", "igad", "mp", "mpe", "mb", "ume", "mbe").
    mrl : String
        Disambiguated mineral name (e.g., "pl", "afs", "mgt", "hem", "omph").
    mbCpx : Int64
        Metabasite clinopyroxene model flag. Controls whether omphacite/diopside maps to "dio" (0) or "aug" (1).
    active_ss : AbstractVector{<:AbstractString}
        Optional list of solution phase names to disambiguate against when a mineral name maps to
        more than one candidate solution phase (only relevant for `db == "all"`). Defaults to empty.

    Returns
    -------
    ss : String
        Solution phase short name (e.g., "fsp", "spl", "amp", "ilm").
"""
function get_ss_from_mineral(db, mrl, mbCpx, active_ss::AbstractVector{<:AbstractString}=String[])

    ss = mrl
    pick(candidates) = candidates[something(findfirst(in(active_ss), candidates), firstindex(candidates))]

    if db =="ig" || db == "igad" || db == "igd"

        if mrl == "cm" || mrl == "mgt" || mrl == "usp" || mrl == "spl"
            ss = "spl"
        elseif mrl == "pat" || mrl == "mu"
            ss = "mu"
        elseif mrl == "afs" || mrl == "pl"
            ss = "fsp"
        elseif mrl == "gl" || mrl == "act" || mrl == "amp" || mrl == "cumm" || mrl == "tr"
            ss = "amp"
        elseif mrl == "hem" || mrl == "ilm"
            ss = "ilm"
        elseif mrl == "pig" || mrl == "Na-cpx"
            ss = "cpx"
        elseif mrl == "K-nph"
            ss = "nph"
        end

    elseif db == "mp" || db == "mpe" || db == "mb" || db == "ume" || db == "mbe"

        if mrl == "smt" || mrl == "sp"
            ss = "sp"
        elseif mrl == "cm" || mrl == "mgt" || mrl == "usp"
            ss = "spl"
        elseif mrl == "afs" || mrl == "pl"
            ss = "fsp"
        elseif mrl == "pat" || mrl == "mu"
            ss = "mu"
        elseif mrl == "gl" || mrl == "act" || mrl == "amp" || mrl == "cumm" || mrl == "tr"
            ss = "amp"
        elseif mrl == "hem" || mrl == "ilm"
            ss = "ilm"
        elseif mrl == "hemm" || mrl == "ilmm"
            ss = "ilmm"
        elseif mrl == "omph" || mrl == "dio" || mrl == "jd"
            if mbCpx == 0
                ss = "dio"
            else
                ss = "aug"
            end
        elseif mrl == "sid" || mrl == "mag" || mrl == "ank" || mrl == "cc"
            ss = "occm"
        elseif mrl == "anth" || mrl == "ged"
            ss = "oamp"
        end

    elseif db == "all"
        if mrl == "afs" || mrl == "pl"
            ss = "fsp_H22"
        elseif mrl == "cm" || mrl == "mgt" || mrl == "spl"
            ss = "spl_T21"
        elseif mrl == "sp" || mrl == "smt"
            ss = "sp_W02"
        elseif mrl == "hem" || mrl == "ilm"
            ss = pick(["ilm_W24", "ilm_T21", "ilm_W00"])
        elseif mrl == "hemm" || mrl == "ilmm"
            ss = pick(["ilmm_W14", "ilmmn_W14"])
        elseif mrl == "gl" || mrl == "act" || mrl == "amp" || mrl == "cumm" || mrl == "tr"
            ss = "amp_G16"
        elseif mrl == "pat" || mrl == "mu"
            ss = "mu_W14"
        elseif mrl == "pig" || mrl == "Na-cpx" || mrl == "cpx"
            ss = pick(["cpx_W24", "cpx_T21"])
        elseif mrl == "K-nph" || mrl == "nph"
            ss = "nph_W24"
        elseif mrl == "omph" || mrl == "dio" || mrl == "jd"
            ss = "dio_G16"
        elseif mrl == "sid" || mrl == "mag" || mrl == "ank" || mrl == "cc"
            ss = "occm_F11"
        elseif mrl == "anth" || mrl == "ged"
            ss = "oamp_D07"
        else
            match = findfirst(s -> split(s, "_")[1] == mrl, active_ss)
            if !isnothing(match)
                ss = active_ss[match]
            end
        end

    elseif occursin("_", mrl)
        ss = split(mrl, "_")[1]
    end

    return ss
end
