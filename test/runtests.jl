# this tests the julia interface to MAGEMin
using Test

cur_dir = pwd();    

if  cur_dir[end-3:end]=="test"
    cd("../")           # change to main directory if we are in /test
end
using MAGEMin_C         # load MAGEMin (needs to be loaded from main directory to pick up correct library in case it is locally compiled)

# Initialize database 
gv, z_b, DB, splx_data      = init_MAGEMin();


sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
test        = 0         #KLB1
gv          = use_predefined_bulk_rock(gv, test)

# Call optimization routine for given P & T & bulk_rock
P           = 8.0
T           = 800.0
gv.verbose  = -1        # switch off any verbose
out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
@show out

@test out.G_system ≈ -797.7491947356159
@test out.ph == ["opx", "ol", "cpx", "spn"]
@test all(abs.(out.ph_frac - [ 0.24226960158631541, 0.5880694152724345, 0.1416697366114075,  0.027991246529842587])  .< 1e-4)

# print more detailed info about this point:
print_info(out)

@testset "pointwise tests  " begin

    for i=1:100
        # same but with int
        P           = 8
        T           = 800
        gv.verbose  = -1        # switch off any verbose
        out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)

        @test out.G_system ≈ -797.7491947356159
        @test out.ph == ["opx", "ol", "cpx", "spn"]
        @test all(abs.(out.ph_frac - [ 0.24226960158631541, 0.5880694152724345, 0.1416697366114075,  0.027991246529842587])  .< 1e-4)
    end
end
finalize_MAGEMin(gv,DB)


@testset "specify bulk rock" begin
    using MAGEMin_C
    gv, z_b, DB, splx_data      = init_MAGEMin();
    bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    bulk_in    = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in     = "wt"
    bulk_rock  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,sys_in);
    gv         = define_bulk_rock(gv, bulk_rock);
    P,T         = 10.0, 1100.0;
    gv.verbose  = -1;        # switch off any verbose
    out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
    finalize_MAGEMin(gv,DB)

    @test out.G_system ≈ -907.2788704076264
end

@testset "convert bulk rock" begin
    bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    bulk_in    = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    bulk_rock  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"wt");

    @test bulk_rock ≈ [45.322438151798686, 8.376385705825816, 11.59989542303507, 14.132959653999844, 7.5133873577293135, 0.3521517320409561, 1.696362856414661, 0.4786296255318463, 1.1547757294831036, 0.009999000099990002, 9.36301476404072]

end


@testset "test Seismic velocities & modulus" begin
    # Call optimization routine for given P & T & bulk_rock
    gv, z_b, DB, splx_data      = init_MAGEMin();
    test        = 0;
    sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
    gv          = use_predefined_bulk_rock(gv, test)
    
    P           = 8.0
    T           = 1200.0
    gv.verbose  = -1        # switch off any verbose
    out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
    
    
    tol = 1e-2;
    @test abs(out.bulkMod - 95.35222421341481           < tol)
    @test abs(out.shearMod - 29.907907390690557         < tol)
    @test abs(out.Vs - 3.056253320843246                < tol)
    @test abs(out.Vp - 6.498781717400121                < tol)
    @test abs(out.Vs_S - 4.30872049030154               < tol)
    @test abs(out.Vp_S - 7.392153167537697              < tol)
    @test abs(out.bulkModulus_M - 27.260603902167567    < tol)
    @test abs(out.bulkModulus_S - 95.74343528580735     < tol)
    @test abs(out.shearModulus_S - 59.4665150508297     < tol)

    finalize_MAGEMin(gv,DB)
end


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

# Automatic testing of all points
function TestPoints(list, gv, z_b, DB, splx_data)

    for i=1:size(list,1)
        gv          = use_predefined_bulk_rock(gv, list[i].test)
        out         = point_wise_minimization(list[i].P,list[i].T, gv, z_b, DB, splx_data, sys_in)

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
    gv, z_b, DB, splx_data      = init_MAGEMin();
    gv.verbose=-1;
    @testset "KLB-1 peridotite tests" begin
        include("test_diagram_test0.jl")
        TestPoints(list, gv, z_b, DB, splx_data)
    end
    finalize_MAGEMin(gv,DB)

    println("  Starting RE-46 icelandic basalt tests")
    gv, z_b, DB, splx_data      = init_MAGEMin();
    gv.verbose=-1;
    @testset "RE-46 icelandic basalt tests" begin
        include("test_diagram_test1.jl")
        TestPoints(list, gv, z_b, DB, splx_data)
    end
    finalize_MAGEMin(gv,DB)

    println("  Starting Wet MORB tests")
    gv, z_b, DB, splx_data      = init_MAGEMin();
    gv.verbose=-1;
    @testset "Wet MORB tests" begin
        include("test_diagram_test6.jl")
        TestPoints(list, gv, z_b, DB, splx_data)
    end
    finalize_MAGEMin(gv,DB)
end



cd(cur_dir)
