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
# Cases of the nullspace optimizer test suite (optimizer=1, solver=0), shared by gen_tests_ns.jl and test_nullspace.jl
# test >= 0 uses the built-in bulk-rock composition, test = -1 uses X/Xoxides (mol); oxides absent from Xoxides are exactly zero

struct ns_case
    name        :: String
    db          :: String
    test        :: Int64
    Xoxides     :: Vector{String}
    X           :: Vector{Float64}
    P           :: Vector{Float64}
    T           :: Vector{Float64}
end

ns_cases = [
    ns_case("ig_KLB1",      "ig",  0, String[], Float64[], [5.0, 10.0, 15.0, 20.0],      [900.0, 1050.0, 1200.0, 1350.0]),
    ns_case("mp_test0",     "mp",  0, String[], Float64[], [4.0, 6.0, 8.0, 10.0],        [550.0, 600.0, 650.0, 700.0]),
    ns_case("mb_test0",     "mb",  0, String[], Float64[], [5.0, 8.0, 11.0, 14.0],       [650.0, 750.0, 850.0, 950.0]),
    ns_case("um_test0",     "um",  0, String[], Float64[], [10.0, 15.0, 20.0, 25.0],     [500.0, 575.0, 650.0, 725.0]),
    ns_case("mtl_test0",    "mtl", 0, String[], Float64[], [150.0, 200.0, 250.0],        [1400.0, 1600.0, 1800.0]),
    ns_case("ig_wetMORB",   "ig",  6, String[], Float64[], [0.01, 5.01, 10.01, 15.01, 20.01, 25.01, 30.01], [800.0, 900.0, 1000.0, 1100.0, 1200.0, 1300.0, 1400.0]),
    ns_case("ig_KLB1_noNa", "ig", -1, ["SiO2", "Al2O3", "CaO", "MgO", "FeO", "K2O"], [38.494, 1.776, 2.824, 50.566, 5.886, 0.01], [5.0, 10.0, 15.0], [800.0, 1000.0, 1200.0]),
    ns_case("ig_MS",        "ig", -1, ["SiO2", "MgO"],                         [45.0, 35.0],             [5.0, 10.0, 20.0],     [900.0, 1200.0, 1500.0]),
    ns_case("ig_CMAS",      "ig", -1, ["SiO2", "Al2O3", "CaO", "MgO"],         [45.0, 10.0, 8.0, 35.0],  [5.0, 10.0, 20.0],     [900.0, 1200.0, 1500.0]),
    ns_case("mp_FASH",      "mp", -1, ["SiO2", "Al2O3", "FeO", "H2O"],         [60.0, 20.0, 10.0, 5.0],  [4.0, 6.0, 8.0],       [500.0, 600.0, 700.0]),
    ns_case("um_MSH",       "um", -1, ["SiO2", "MgO", "H2O"],                  [40.0, 50.0, 10.0],       [10.0, 20.0],          [400.0, 600.0, 800.0]),
    ns_case("ume_CMS",      "ume",-1, ["SiO2", "CaO", "MgO"],                  [45.0, 8.0, 35.0],        [10.0, 20.0],          [500.0, 700.0, 900.0]),
    ns_case("mtl_MF",       "mtl",-1, ["MgO", "FeO"],                          [35.0, 8.0],              [150.0, 250.0],        [1400.0, 1800.0]),
    ns_case("sb11_test0",   "sb11", 0, String[], Float64[], [50.0, 100.0, 150.0],        [1200.0, 1600.0, 2000.0]),
    ns_case("sb21_test0",   "sb21", 0, String[], Float64[], [50.0, 100.0, 150.0],        [1200.0, 1600.0, 2000.0]),
    ns_case("sb24_test0",   "sb24", 0, String[], Float64[], [50.0, 100.0, 150.0],        [1200.0, 1600.0, 2000.0]),
    ns_case("sb11_CMAS",    "sb11",-1, ["SiO2", "Al2O3", "CaO", "MgO"],        [45.0, 10.0, 8.0, 35.0],  [50.0, 100.0, 150.0],  [1200.0, 1600.0, 2000.0]),
    ns_case("sb21_MS",      "sb21",-1, ["SiO2", "MgO"],                        [50.0, 50.0],             [50.0, 100.0, 150.0],  [1200.0, 1600.0, 2000.0]),
    ns_case("sb24_noNaCr",  "sb24",-1, ["SiO2", "Al2O3", "CaO", "MgO", "Fe", "O"], [38.0, 2.8, 2.0, 50.0, 6.0, 0.1], [50.0, 100.0, 150.0], [1200.0, 1600.0, 2000.0]),
]

function ns_case_points(c :: ns_case)
    P = Float64[]; T = Float64[]
    for p in c.P, t in c.T
        push!(P, p); push!(T, t)
    end
    return P, T
end

function ns_case_run(c :: ns_case, P, T)
    data = Initialize_MAGEMin(c.db, verbose=-1, solver=0, optimizer=1)
    if c.test >= 0
        out = multi_point_minimization(P, T, data, test=c.test, progressbar=false)
    else
        out = multi_point_minimization(P, T, data, X=c.X, Xoxides=c.Xoxides, sys_in="mol", progressbar=false)
    end
    Finalize_MAGEMin(data)
    return out
end
