using Test

# Load MAGEMin (needs to be loaded from main directory to pick up correct
# library in case it is locally compiled). This is handled by the logic in
# runtests.jl
using MAGEMin_C

# Initialize database  - new way
data        =   Initialize_MAGEMin("ig", verbose=true);
test        =   0         #KLB1
data        =   use_predefined_bulk_rock(data, test);

# Call optimization routine for given P & T & bulk_rock
P           =   8.0
T           =   800.0
out         =   point_wise_minimization(P,T, data);

@show out

@test out.G_system ≈ -797.7491824683576
@test out.ph == ["opx", "ol", "cpx", "spn"]
@test all(abs.(out.ph_frac - [ 0.24226960158631541, 0.5880694152724345, 0.1416697366114075,  0.027991246529842587])  .< 1e-4)

# print more detailed info about this point:
print_info(out)
Finalize_MAGEMin(data)



# previous way we defined this (left here for backwards compatibility)
db          = "ig"
gv, z_b, DB, splx_data  = init_MAGEMin(db);
sys_in      =   "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
test        =   0         #KLB1
gv          =   use_predefined_bulk_rock(gv, test, db);
gv.verbose=-1
P           =   8.0
T           =   800.0
out         =   point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in);
@test out.G_system ≈ -797.7491824683576
@test out.ph == ["opx", "ol", "cpx", "spn"]
@test all(abs.(out.ph_frac - [ 0.24226960158631541, 0.5880694152724345, 0.1416697366114075,  0.027991246529842587])  .< 1e-4)
finalize_MAGEMin(gv,DB)


@testset "pointwise tests  " begin
    n       =   100;
    P       =   fill(8.0,n)
    T       =   fill(800.0,n)
    db      =   "ig"
    data    =   Initialize_MAGEMin(db, verbose=false);
    out     =   multi_point_minimization(P, T, data, test=0);
    @test out[end].G_system ≈ -797.7491824683576
    @test out[end].ph == ["opx", "ol", "cpx", "spn"]
    @test all(abs.(out[end].ph_frac - [ 0.24226960158631541, 0.5880694152724345, 0.1416697366114075,  0.027991246529842587])  .< 1e-4)

    Finalize_MAGEMin(data)
end


@testset "ultramafic database" begin
    P,T     = 10.0, 500.0
    data    = Initialize_MAGEMin("um", verbose=false);
    out     = single_point_minimization(P, T, data)
    @test abs(out.G_system + 556.243321)/abs(556.243321) < 2e-4
    @test all(abs.(out.ph_frac - [0.5753192548298169, 0.024582756340200455, 0.024184726745689656, 0.0554904375970385, 0.0016172391886508995, 0.3184197873654363, 0.0003857953985807551])  .< 1e-4)
    @test out.ph == ["atg", "ol", "spi", "chl", "po", "fluid", "pyr"]

end


@testset "specify bulk rock" begin

    data    = Initialize_MAGEMin("ig", verbose=false);

    # One bulk rock for all points
    P,T     = 10.0, 1100.0
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X       = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    sys_in  = "wt"
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

    @test abs(out.G_system + 917.008526)/abs(917.008526) < 2e-4


    # different bulk rock per point
    P       = [10.0, 10.0]
    T       = [1100.0, 1100.0]
    Xoxides = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    X1      = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X2      = [49.43; 14.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    X       = [X1,X2]
    sys_in  = "wt"
    out     = multi_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)

    @test out[1].G_system ≈ -917.0085258258146 rtol=2e-4
    @test out[2].G_system ≈ -912.7754865001111 rtol=2e-4

    Finalize_MAGEMin(data)
end

@testset "convert bulk rock" begin

    bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"];
    bulk_in    = [48.43; 15.19; 11.57; 10.13; 6.65; 1.64; 0.59; 1.87; 0.68; 0.0; 3.0];
    bulk_rock  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"wt","ig" );

    @test bulk_rock ≈ [46.12136551106488,8.52404156868422,11.80437413592006,14.382090296468109,6.470772496142926,0.3583593339444149,1.7262657202608955,0.487066733471891,0.5876026409473947,0.009999000099989998,9.528062562995213]

    bulk_in_ox = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "MnO"; "H2O"];
    bulk_in    = [69.64; 13.76; 1.77; 1.73; 4.32; 0.4; 2.61; 2.41; 0.80; 0.07; 0.0];
    bulk_rock  = convertBulk4MAGEMin(bulk_in,bulk_in_ox,"wt","mp" );

    @test bulk_rock ≈ [76.57038397179574,8.914984523583415,2.0849576977131403,2.835783318610597,4.30275071755529,1.8302970975627948,2.568605789798099,0.6615823604771729,0.16546809116073818,0.06518643174302832,0.0]
end


@testset "Seismic velocities & modulus" begin
    # Call optimization routine for given P & T & bulk_rock
    data         = Initialize_MAGEMin("ig", verbose=false);
    test        = 0;
    data         = use_predefined_bulk_rock(data, test)
    P           = 8.0
    T           = 1200.0
    out         = point_wise_minimization(P,T, data)

    tol = 5e-2;
    @test abs(out.bulkMod - 95.35222421341481           < tol)
    @test abs(out.shearMod - 29.907907390690557         < tol)
    @test abs(out.Vs - 3.056253320843246                < tol)
    @test abs(out.Vp - 6.498781717400121                < tol)
    @test abs(out.Vs_S - 4.30872049030154               < tol)
    @test abs(out.Vp_S - 7.392153167537697              < tol)
    # @test abs(out.bulkModulus_M - 27.260603902167567    < tol)
    @test abs(out.bulkModulus_S - 95.74343528580735     < tol)
    @test abs(out.shearModulus_S - 59.4665150508297     < tol)

    Finalize_MAGEMin(data)
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
function TestPoints(list, data::MAGEMin_Data)

    # Compute all points
    P = [ l.P for l in list]
    T = [ l.T for l in list]
    test = [ l.test for l in list]
    out_vec = multi_point_minimization(P, T, data, test = test[1]);

    # Check if the points this fit
    for (i,out) in enumerate(out_vec)
        VerifyPoint(out, list[i], i)
    end
    return nothing
end

# This checks whether a single point agrees with precomputed values & prints a message if not
function VerifyPoint(out, list, i)

     # We need to sort the phases (sometimes they are ordered differently)
     ind_sol = sortperm(list.ph)
     ind_out = sortperm(out.ph)

     result1 = @test out.G_system  ≈ list.G     rtol=1e-3
     result2 = @test out.ph[ind_out]        == list.ph[ind_sol]
     result3 = @test out.ph_frac[ind_out] ≈ list.ph_frac[ind_sol] atol=5e-2       # ok, this is really large (needs fixing for test6!)

     # print more info about the point if one of the tests above fails
     if isa(result1,Test.Fail) || isa(result2,Test.Fail) || isa(result3,Test.Fail)
         print_error_msg(i,list)
     end

     return nothing
end

# load reference for built-in tests
println("Testing points from the reference diagrams:")
@testset verbose = true "Total tests" begin
    println("  Starting KLB-1 peridotite tests")
    db  = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
    data = Initialize_MAGEMin(db, verbose=false);

    gv.verbose=-1;
    @testset "KLB-1 peridotite tests" begin
        include("test_diagram_test0.jl")
        TestPoints(list, data)
    end
    Finalize_MAGEMin(data)

    println("  Starting RE-46 icelandic basalt tests")
    db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
    data = Initialize_MAGEMin(db, verbose=false);
    gv.verbose=-1;
    @testset "RE-46 icelandic basalt tests" begin
        include("test_diagram_test1.jl")
        TestPoints(list, data)
    end


    println("  Starting Wet MORB tests")
    db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
    data = Initialize_MAGEMin(db, verbose=false);
    @testset "Wet MORB tests" begin
        include("test_diagram_test6.jl")
        TestPoints(list, data)
    end
    Finalize_MAGEMin(data)
end
