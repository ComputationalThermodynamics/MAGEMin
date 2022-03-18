using Clang.Generators
using Clang.LibClang.Clang_jll
using Pkg
using Pkg.Artifacts
using MAGEMin_jll
using NLopt_jll
using MPICH_jll

cd(@__DIR__)

# Exclude a few things

const HASH_JEN = 0

# headers 
MAGEMin_toml = joinpath(dirname(pathof(MAGEMin_jll)), "..", "Artifacts.toml")
MAGEMin_dir = Pkg.Artifacts.ensure_artifact_installed("MAGEMin", MAGEMin_toml)

NLopt_toml = joinpath(dirname(pathof(NLopt_jll)), "..", "Artifacts.toml")
NLopt_dir = Pkg.Artifacts.ensure_artifact_installed("NLopt", NLopt_toml)

MPICH_toml = joinpath(dirname(pathof(MPICH_jll)), "..", "Artifacts.toml")       # not for windows
MPICH_dir = Pkg.Artifacts.ensure_artifact_installed("MPICH", MPICH_toml)

magemin_include_dir = normpath(MAGEMin_dir, "include") |> normpath
magemin_h = joinpath(magemin_include_dir, "MAGEMin.h")
@assert isfile(magemin_h)

nlopt_include_dir = normpath(NLopt_dir, "include") |> normpath
nlopt_h = joinpath(nlopt_include_dir, "nlopt.h")
@assert isfile(nlopt_h)

mpi_include_dir = joinpath(MPICH_dir, "include") |> normpath
mpi_h = joinpath(mpi_include_dir, "mpi.h")
@assert isfile(mpi_h)

# wrapper generator options
options = load_options(joinpath(@__DIR__, "generator.toml"))

# add compiler flags, e.g. "-DXXXXXXXXX"
args = get_default_args()
push!(args, "-I$magemin_include_dir")
push!(args, "-I$nlopt_include_dir")
push!(args, "-isystem$mpi_include_dir")

header_files = [magemin_h]


#headers = [joinpath(magemin_include_dir, header) for header in readdir(include_dir) if endswith(header, ".h")]
# there is also an experimental `detect_headers` function for auto-detecting top-level headers in the directory
# headers = detect_headers(clang_dir, args)

# create context
ctx = create_context(header_files, args, options)

# run generator
build!(ctx)