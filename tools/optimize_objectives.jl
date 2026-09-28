
#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# julia tools/optimize_objectives.jl src/TC_database/objective_functions.c -o src/TC_database/objective_functions.c

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#
const HELPERS_MARK = "static inline void mu_Gex_sym_n2("

const HELPERS = "\n" * """
static inline void mu_Gex_sym_n2(SS_ref *d, double *mu_Gex)
{
    int n = d->n_em, jmax = d->n_xeos, it = 0;
    double *p = d->p, *Wv = d->W, Q = 0.0, Wp[n];
    for (int i = 0; i < n; i++) Wp[i] = 0.0;
    for (int j = 0; j < jmax; j++) {
        for (int k = j + 1; k < n; k++) {
            double w = Wv[it++];
            Wp[j] += w * p[k];
            Wp[k] += w * p[j];
            Q += w * p[j] * p[k];
        }
    }
    for (int i = 0; i < n; i++) mu_Gex[i] = Wp[i] - Q;
}

static inline void mu_Gex_asym_n2(SS_ref *d, double *mu_Gex)
{
    int n = d->n_em, jmax = d->n_xeos, it = 0;
    double *phi = d->mat_phi, *v = d->v, *Wv = d->W, Q = 0.0, Bf[n];
    for (int i = 0; i < n; i++) Bf[i] = 0.0;
    for (int j = 0; j < jmax; j++) {
        for (int k = j + 1; k < n; k++) {
            double b = Wv[it++] * 2.0 / (v[j] + v[k]);
            Bf[j] += b * phi[k];
            Bf[k] += b * phi[j];
            Q += b * phi[j] * phi[k];
        }
    }
    for (int i = 0; i < n; i++) mu_Gex[i] = v[i] * (Bf[i] - Q);
}
"""

const NEM = raw"(?:d->)?n_em"
const LOOP_A = Regex(
    raw"[ \t]*for \(int i = 0; i < " * NEM * raw"; i\+\+\)\s*\{\s*" *
    raw"Gex = 0\.0;\s*int it\s*=\s*0;\s*" *
    raw"for \(int j = 0; j < d->n_xeos; j\+\+\)\s*\{\s*" *
    raw"tmp = \((?<t>[^;]*)\);\s*" *
    raw"for \(int k = j\s*\+\s*1; k < " * NEM * raw"; k\+\+\)\s*\{\s*" *
    raw"Gex -= (?<e>[^;]*);\s*it \+= 1;\s*\}\s*\}\s*" *
    raw"mu_Gex\[i\] = Gex;\s*\}")
const LOOP_B = Regex(
    raw"[ \t]*for \(int i = 0; i < " * NEM * raw"; i\+\+\)\s*\{\s*" *
    raw"mu_Gex\[i\] = 0\.0;\s*int it\s*=\s*0;\s*" *
    raw"for \(int j = 0; j < d->n_xeos; j\+\+\)\s*\{\s*" *
    raw"for \(int k = j\s*\+\s*1; k < " * NEM * raw"; k\+\+\)\s*\{\s*" *
    raw"mu_Gex\[i\] -= (?<e>[^;]*);\s*it \+= 1;\s*\}\s*\}\s*\}")

const SYM_A  = ("d->eye[i][j]-d->p[j]", "tmp*(d->eye[i][k]-d->p[k])*(d->W[it])")
const ASYM_A = ("d->eye[i][j]-d->mat_phi[j]", "tmp*(d->eye[i][k]-d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]))")
const SYM_B  = "(d->eye[i][j]-d->p[j])*(d->eye[i][k]-d->p[k])*(d->W[it])"
const ASYM_B = "(d->eye[i][j]-d->mat_phi[j])*(d->eye[i][k]-d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]))"

const POW       = r"cpow\(sf\[(\d+)\],\s*([-\d.]+)\)"
const SQRT      = r"csqrt\(sf\[(\d+)\]\)"
const SF_ASSIGN = r"^\s*sf\[\d+\]\s*="
const DECL      = r"^\s*double complex c[ps]_\w+ = "
const FUNC_HEAD = r"^double obj_\w+\("m
const CALL_SYM  = "    mu_Gex_sym_n2(d, mu_Gex);"
const CALL_ASYM = "    mu_Gex_asym_n2(d, mu_Gex);"

nows(s) = replace(s, r"\s+" => "")
tag(y) = replace(replace(y, "-" => "m"), "." => "p")

function mu_gex_pass(func::AbstractString)
    func = replace(func, LOOP_A => function (s)
        m = match(LOOP_A, s)
        key = (nows(m[:t]), nows(m[:e]))
        key == SYM_A  && return CALL_SYM
        key == ASYM_A && return CALL_ASYM
        return s
    end)
    return replace(func, LOOP_B => function (s)
        e = nows(match(LOOP_B, s)[:e])
        e == SYM_B  && return CALL_SYM
        e == ASYM_B && return CALL_ASYM
        return s
    end)
end

function cse_pass(func::AbstractString)
    lines = split(func, '\n')
    uses = [i for (i, l) in enumerate(lines) if (occursin(POW, l) || occursin(SQRT, l)) && !occursin(DECL, l)]
    if isempty(uses)
        return func, any(l -> occursin(DECL, l), lines) ? "hoisted" : "none"
    end
    first_use = uses[1]
    l1 = lines[first_use]
    if !startswith(l1, "    ") || startswith(l1, "     ")
        return func, "first use not at function level"
    end
    sf_lines = [i for (i, l) in enumerate(lines) if occursin(SF_ASSIGN, l)]
    if !isempty(sf_lines) && maximum(sf_lines) > first_use
        return func, "sf assigned after first use"
    end
    body = join(lines[first_use:end], "\n")
    pows = Tuple{String,String}[]
    sqs  = String[]
    for m in eachmatch(POW, body)
        k = (String(m[1]), String(m[2]))
        k in pows || push!(pows, k)
    end
    for m in eachmatch(SQRT, body)
        String(m[1]) in sqs || push!(sqs, String(m[1]))
    end
    decl = String[]
    for (k, y) in pows
        push!(decl, "    double complex cp_$(k)_$(tag(y)) = cpow(sf[$(k)], $(y));")
    end
    for k in sqs
        push!(decl, "    double complex cs_$(k) = csqrt(sf[$(k)]);")
    end
    body = replace(body, POW => s -> (m = match(POW, s); "cp_$(m[1])_$(tag(m[2]))"))
    body = replace(body, SQRT => s -> "cs_$(match(SQRT, s)[1])")
    return join(vcat(lines[1:first_use-1], decl, [body]), "\n"), "hoisted"
end

function classify_mu(func::AbstractString)
    body = first(split(func, "\n}\n"))
    (occursin("mu_Gex_sym_n2(", body) || occursin("mu_Gex_asym_n2(", body)) && return "O(n^2)"
    (!occursin("mu_Gex", body) && !occursin("Gex", body)) && return "no mu_Gex"
    occursin(r"mu_Gex\[\d+\]\s*=", body) && return "explicit mu_Gex[k] lines"
    return "loop not recognised"
end

mugex_done(func::AbstractString) = occursin("mu_Gex_sym_n2(", func) || occursin("mu_Gex_asym_n2(", func)
cse_done(func::AbstractString) = any(l -> occursin(DECL, l), split(func, '\n'))

function split_functions(src::AbstractString)
    starts = [m.offset for m in eachmatch(FUNC_HEAD, src)]
    isempty(starts) && return src, String[]
    head  = src[1:prevind(src, starts[1])]
    ends  = vcat([prevind(src, s) for s in starts[2:end]], [lastindex(src)])
    funcs = [src[starts[i]:ends[i]] for i in eachindex(starts)]
    return head, funcs
end

function optimize_objectives(src::AbstractString; cse::Bool = true, mugex::Bool = true)
    head, funcs = split_functions(src)
    stats = Dict{String,Dict{String,Int}}()
    out = String[]
    for f in funcs
        name = match(r"double obj_(\w+)\(", f)[1]
        db = String(first(split(name, "_")))
        c = get!(stats, db, Dict{String,Int}())
        c["objectives"] = get(c, "objectives", 0) + 1
        g = String(f)
        if mugex
            if mugex_done(g)
                c["mu_Gex already O(n^2)"] = get(c, "mu_Gex already O(n^2)", 0) + 1
            else
                g = mu_gex_pass(g)
                k = "mu_Gex " * classify_mu(g)
                k == "mu_Gex O(n^2)" && (k = "mu_Gex O(n^2) new")
                c[k] = get(c, k, 0) + 1
            end
        end
        if cse
            if cse_done(g)
                c["cpow/csqrt already hoisted"] = get(c, "cpow/csqrt already hoisted", 0) + 1
            else
                g, why = cse_pass(g)
                k = "cpow/csqrt " * (why == "hoisted" ? "hoisted new" : why)
                c[k] = get(c, k, 0) + 1
            end
        end
        push!(out, g)
    end
    res = head * join(out)
    if mugex && occursin("mu_Gex_", res) && !occursin(HELPERS_MARK, res)
        pos = match(r"^(double obj_|void px_)"m, res).offset
        res = res[1:prevind(res, pos)] * HELPERS * "\n" * res[pos:end]
    end
    return res, stats
end

function print_report(io::IO, stats)
    cols = sort(unique([k for c in values(stats) for k in keys(c) if k != "objectives"]))
    println(io, rpad("db", 8), " ", lpad("obj", 4), "  ", join(cols, "  "))
    for db in sort(collect(keys(stats)))
        c = stats[db]
        println(io, rpad(db, 8), " ", lpad(c["objectives"], 4), "  ", join([lpad(get(c, k, 0), length(k)) for k in cols], "  "))
    end
end

function default_save_dir(file::AbstractString)
    d = dirname(abspath(file))
    while true
        basename(d) == "src" && return joinpath(d, "saves")
        p = dirname(d)
        p == d && return joinpath(dirname(abspath(file)), "saves")
        d = p
    end
end

function save_original(file::AbstractString, src::AbstractString, save_dir::AbstractString)
    mkpath(save_dir)
    stem, ext = splitext(basename(file))
    stamp = Libc.strftime("%Y%m%d_%H%M%S", time())
    path = joinpath(save_dir, "$(stem)_$(stamp)$(ext)")
    n = 1
    while isfile(path)
        path = joinpath(save_dir, "$(stem)_$(stamp)_$(n)$(ext)")
        n += 1
    end
    write(path, src)
    return path
end

function main(args)
    file = nothing; output = nothing; cse = true; mugex = true; report = false; save = true; save_dir = nothing
    i = 1
    while i <= length(args)
        a = args[i]
        if a in ("-o", "--output")
            output = args[i+1]; i += 1
        elseif a == "--no-cse"
            cse = false
        elseif a == "--no-mugex"
            mugex = false
        elseif a == "--report"
            report = true
        elseif a == "--no-save"
            save = false
        elseif a == "--save-dir"
            save_dir = args[i+1]; i += 1
        else
            file = a
        end
        i += 1
    end
    file === nothing && error("usage: julia optimize_objectives.jl FILE [-o OUT] [--no-cse] [--no-mugex] [--report] [--no-save] [--save-dir DIR]")
    src = read(file, String)
    res, stats = optimize_objectives(src; cse, mugex)
    if output !== nothing
        if res == src
            println(stderr, "already optimized: nothing to do, file not written, nothing saved")
            abspath(output) == abspath(file) || write(output, res)
        else
            if save
                path = save_original(file, src, save_dir === nothing ? default_save_dir(file) : save_dir)
                println(stderr, "saved unmodified copy: ", path)
            end
            write(output, res)
        end
    elseif !report
        print(stdout, res)
    end
    print_report(stderr, stats)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
