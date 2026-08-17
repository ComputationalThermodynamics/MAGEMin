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
# Generate a pseudocompound (PC) composition grid for a "gh" (Ghiorso/MELTS)
# solution phase, in the same struct literal format already hand-written in
# src/GH_database/SS_xeos_PC_gh.c (e.g. "struct ss_pc gh_liq_pc_xeos[26] = {...}").
#
# generate_simplex/adjust_simplex below are ported verbatim from
# ref_database/SB_ref_database/functions_ss.jl (get_sb_SS_xeos_PC and its
# helpers) - same algorithm, unchanged. The only real difference is the
# function signature: gh's endmember counts are small fixed constants
# hardcoded across the codebase (gh_gss_init_function.c), not read from a
# JSON endmember database the way sb's phases are, so the entry point here
# takes (ph_name, step_size) directly instead of (sb_ver, ss).
#
# Usage (CLI):
#   julia tools/generate_pseudocompounds.jl <ph_name> <step_size>
#   julia tools/generate_pseudocompounds.jl liq 0.3
#
# Usage (REPL):
#   include("tools/generate_pseudocompounds.jl")
#   print(generate_gh_pseudocompounds("ol", 0.1))

# n_em per gh solution phase - keep in sync with gh_gss_init_function.c's
# n_em assignments (liq=13, ol=2, fsp=3, bi=2, g=3, hb=3, lc=3, mel=4, cum=2,
# spn=5, cpx=7, opx=7 (shares cpx's grid), fl=2).

const GH_N_EM = Dict(
    "liq" => 13,
    "ol"  => 2,
    "fsp" => 3,
    "bi"  => 2,
    "g"   => 3,
    "hb"  => 3,
    "lc"  => 3,
    "mel" => 4,
    "cum" => 2,
    "spn" => 5,
    "cpx" => 7,
    "opx" => 7,
    "fl"  => 2,
    "rhm" => 5,
    "nph" => 4,
    "kls" => 4,
    "liq_pmelts_dataset" => 12,
)
 
# Function to generate a simplex in n dimensions with given step size
# (verbatim from functions_ss.jl)
function generate_simplex(n, step)
    # Generate all possible points with the given step size
    points = collect(0:step:1)

    # Generate all combinations of points in n dimensions
    combinations = Iterators.product(ntuple(_ -> points, n)...)

    # Filter combinations to keep only those where the sum of coordinates is 1
    simplex = [collect(comb) for comb in combinations if sum(comb) ≈ 1.0]

    return simplex
end

# Function to adjust the simplex points
# (verbatim from functions_ss.jl)
function adjust_simplex(simplex, eps)
    adjusted_simplex = []
    for point in simplex
        adjusted_point = copy(point)
        for i in 1:length(point)
            if point[i] == 0.0
                adjusted_point[i] = eps
            elseif point[i] == 1.0
                adjusted_point[i] = 1.0 - eps
            end
        end
        # Adjust the remaining coordinates to ensure the sum is still 1
        sum_adjusted = sum(adjusted_point)
        if sum_adjusted ≈ 1.0
            push!(adjusted_simplex, adjusted_point)
        else
            # Find the index of the largest coordinate
            max_index = argmax(adjusted_point)
            adjusted_point[max_index] += 1.0 - sum_adjusted
            push!(adjusted_simplex, adjusted_point)
        end
    end
    return adjusted_simplex
end

"""
    generate_gh_pseudocompounds(ph_name::String, step_size::Float64; eps=1e-4)

Generate the pseudocompound composition grid for gh solution phase `ph_name`
at `step_size`, and return it as a ready-to-paste C struct literal matching
SS_xeos_PC_gh.c's existing convention, e.g.:

    struct ss_pc gh_ol_pc_xeos[9] = {
        {{0.1,0.9}},
        ...
    };
"""
function generate_gh_pseudocompounds(ph_name::String, step_size::Float64; eps::Float64 = 1e-2)
    haskey(GH_N_EM, ph_name) ||
        error("unknown gh phase \"$ph_name\" - add it to GH_N_EM first (known: $(join(sort(collect(keys(GH_N_EM))), ", ")))")
    n_em = GH_N_EM[ph_name]

    simplex          = generate_simplex(n_em, step_size)
    adjusted_simplex = adjust_simplex(simplex, eps)
    n_pc             = length(adjusted_simplex)

    tab = "    "
    out = "struct ss_pc gh_$(ph_name)_pc_xeos[$n_pc] = {\n"
    for i = 1:n_pc
        out *= tab * "{{"
        out *= join(round.(adjusted_simplex[i], digits = 6), ",")
        out *= "}},\n"
    end
    out *= "};\n"
    return out
end

# if abspath(PROGRAM_FILE) == @__FILE__
#     if length(ARGS) < 2
#         println("usage: julia generate_pseudocompounds.jl <ph_name> <step_size>")
#         println("known phases: ", join(sort(collect(keys(GH_N_EM))), ", "))
#         exit(1)
#     end
#     ph_name   = ARGS[1]
#     step_size = parse(Float64, ARGS[2])
#     print(generate_gh_pseudocompounds(ph_name, step_size))
# end

# out = generate_gh_pseudocompounds("spn", 1/10.0)

out = generate_gh_pseudocompounds("liq", 1/4.0;eps=1e-2)
write("gh_liq_v2_pc_0.25.txt", out)