# this tests the julia interface to MAGEMin
using Test

cur_dir = pwd();    

if  cur_dir[end-3:end]=="test"
    cd("../")           # change to main directory if we are in /test
end
using MAGEMin_C         # load MAGEMin (needs to be loaded from main directory to pick up correct library in case it is locally compiled)

# Initialize database 
gv, DB      = init_MAGEMin();

test        = 0;
sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
bulk_rock   = get_bulk_rock(gv, test)

# Call optimization routine for given P & T & bulk_rock
P           = 8.0
T           = 800.0
gv.verbose  = -1        # switch off any verbose
out         = point_wise_minimization(P,T, bulk_rock, gv, DB,sys_in);
@show out

@test out.G_system ≈ -797.7491824869334
@test out.ph == [ "opx", "cpx", "ol", "spn"]

@test all(abs.(out.ph_frac - [ 0.24226960158631541, 0.1416697366114075, 0.5880694152724345, 0.027991246529842587])  .< 1e-4)

# print more detailed info about this point:
print_info(out)

# Stores data of tests
mutable struct outP{ _T  } 
    P           ::  _T
    T           ::  _T 
    test        ::  Int64

    G           ::  _T
    ph          ::  Vector{String}
    ph_frac     ::  Vector{Float64}
end

print_error_msg(i,out) = println("ERROR for point $i with test=$(out.test); P=$(out.P); T=$(out.T); stable phases=$(out.ph), fractions=$(out.ph_frac)")
finalize_MAGEMin(gv,DB)

# Automatic testing of all points
function TestPoints(list, gv, DB)

    for i=1:size(list,1)
        bulk_rock   = get_bulk_rock(gv, list[i].test)
        out         = point_wise_minimization(list[i].P,list[i].T, bulk_rock, gv, DB)

        # We need to sort the phases (sometimes they are ordered differently)
        ind_sol = sortperm(list[i].ph)
        ind_out = sortperm(out.ph)
        
        result1 = @test out.G_system  ≈ list[i].G     rtol=1e-3
        result2 = @test out.ph[ind_out]        == list[i].ph[ind_sol]
        result3 = @test out.ph_frac[ind_out] ≈ list[i].ph_frac[ind_sol] atol=1.5e-2       # ok, this is really large (needs fixing for test6!)
        
        # print more info about the point if one of the tests above fails
        if isa(result1,Test.Fail) || isa(result2,Test.Fail) || isa(result3,Test.Fail)
            print_error_msg(i,list[i])
        end

    end
end

# load reference for built-in tests
println("Testing points from the reference diagrams:")
@testset verbose = true "Total tests" begin
    println("  Starting KLB-1 peridotite tests")
    gv, DB = init_MAGEMin();
    gv.verbose=-1;
    @testset "KLB-1 peridotite tests" begin
        include("test_diagram_test0.jl")
        TestPoints(list, gv, DB)
    end
    finalize_MAGEMin(gv,DB)

    println("  Starting RE-46 icelandic basalt tests")
    gv, DB = init_MAGEMin();
    gv.verbose=-1;
    @testset "RE-46 icelandic basalt tests" begin
        include("test_diagram_test1.jl")
        TestPoints(list, gv, DB)
    end
    finalize_MAGEMin(gv,DB)

    println("  Starting Wet MORB tests")
    gv, DB = init_MAGEMin();
    gv.verbose=-1;
    @testset "Wet MORB tests" begin
        include("test_diagram_test6.jl")
        TestPoints(list, gv, DB)
    end
    finalize_MAGEMin(gv,DB)
end



cd(cur_dir)
