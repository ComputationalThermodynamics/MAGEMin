# this tests the julia interface to MAGEMin
using Test

# MAGEMin_C loads a locally compiled version of the library if it detects
# a file with the appropriate name in the current working directory.
# To make use of this also for tests, we change the working directory to
# the parent directory of the package iff necessary.
cur_dir = pwd()

if endswith(cur_dir, "test")
    cd("../")           # change to main directory if we are in /test
end


@testset "serial" begin
    include(joinpath(@__DIR__, "tests.jl"))
end

@testset "threaded" begin
    # Do a dummy `@test true`:
    # If the process errors out the testset would error out as well
    @test true

    # Only run additional tests if we are running with a single thread right now
    if Threads.nthreads() == 1
        # We explicitly disable code coverage tracking with multiple threads since
        # this is expensive, see https://github.com/JuliaLang/julia/issues/36142
        run(`$(Base.julia_cmd()) --threads=2 --check-bounds=yes --code-coverage=none $(abspath(joinpath(@__DIR__, "tests.jl")))`)
    end
end


# Change back to the working directory we used when starting to run the tests
cd(cur_dir)
