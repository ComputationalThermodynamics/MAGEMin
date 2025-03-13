using Clang.Generators
using Clang.LibClang.Clang_jll
using Pkg
using Pkg.Artifacts
#using MAGEMin_jll
using NLopt_jll
# using MPICH_jll

cd(@__DIR__)

# Exclude a few things

# headers 
#MAGEMin_toml = joinpath(dirname(pathof(MAGEMin_jll)), "..", "Artifacts.toml")
#MAGEMin_dir = Pkg.Artifacts.ensure_artifact_installed("MAGEMin", MAGEMin_toml)

# NOTE: for development, I use the local version of the header files (as they required some changes); 
#  Once those changes are merged into a MAGEMin_jll version, we can use the lines above  
MAGEMin_dir = joinpath(pwd(),"../src");  

NLopt_toml          = joinpath(dirname(pathof(NLopt_jll)), "..", "Artifacts.toml")
NLopt_dir           = Pkg.Artifacts.ensure_artifact_installed("NLopt", NLopt_toml)

# MPICH_toml          = joinpath(dirname(pathof(MPICH_jll)), "..", "Artifacts.toml")       # not for windows
# MPICH_dir           = Pkg.Artifacts.ensure_artifact_installed("MPICH", MPICH_toml)

magemin_include_dir = MAGEMin_dir
nlopt_include_dir   = normpath(NLopt_dir, "include") |> normpath
# mpi_include_dir     = joinpath(MPICH_dir, "include") |> normpath

# wrapper generator options
options = load_options(joinpath(@__DIR__, "generator.toml"))

# add compiler flags, e.g. "-DXXXXXXXXX"
args = get_default_args()
push!(args, "-I$magemin_include_dir")
push!(args, "-I$nlopt_include_dir")
# push!(args, "-isystem$mpi_include_dir")

# Process all header files in the magemin_include_dir directory:
header_files = [joinpath(magemin_include_dir, header) for header in readdir(magemin_include_dir) if endswith(header, ".h")]

# create context
ctx = create_context(header_files, args, options)

# run generator
build!(ctx)