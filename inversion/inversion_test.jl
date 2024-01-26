#=
    NR, 14/07/2023
    Code prototype to use MAGEMin for calibration of database. 
    Note that the constraints are defined in the constraints.jl file
    To test the routine, you can generate a set of constraints using gen_test_constraints.jl (for production run, these would be experimental constraints)

    -> problem_definition.jl allow to define the set of variables to invert for

    --------------------------------------------    
                  Available database
    --------------------------------------------
    ig  igneous     (HP18 -> Green et al., 2023)
    igd igneous     (T21 -> Green et al., 2023) 
    alk alkaline    (Weller et al., 2023)
    mp  metapelite  (White et al 2014b) 
    mb  metabasite  (Green et al., 2016)
    um  ultramafic  (Evans & Frost, 2021)
    --------------------------------------------

    Missfit function:
        The missfit function being minimized is defined as:

            missfit = sum(mf[i]) (1)
            
                i is a constraint index (known phase equilibrium at given PTX conditions)

            mf[i] = sum( abs(cons_ph_fraction - est_ph_fraction)** )_PH_INTERSECTION + sum( abs(cons_ph_fraction OR est_ph_fraction)*10.0 )_PH_DIFFERENCE

                cons_ph_fraction is the phase fraction of reference (constraint)
                est_ph_fraction is the estimated (calculated) phase fraction
                _PH_INTERSECTION is the set of common phase between constraint and estimated stable phases
                _PH_DIFFERENCE is the set of phase that are different between constraint and estimated stable phases
             ** only if abs(cons_ph_fraction - est_ph_fraction) > 1e-4;

=#


#=
    Load libraries
=#
using MAGEMin_C
using NLopt
# using Plots



#=
    Load needed set of functions
=#
include("inversion_functions.jl")


#=
    Select options below
=#
database        = "ig";                 # select database here, ig, igd, alk, mp, mb, um
display_ss_info = 1;                    # shows the names of the solution phases of the database, and their sizes


#=
    Initialize MAGEMin with the selected database
=#
# global gv, z_b, DB, splx_data   = init_MAGEMin(database);


#=
    Here we retrieve the name  and size for the W's and g0's
    Note that  values have not been initialized yet, only the memory is allocated
=#
if display_ss_info == 1
    display_ss_infos( database );
end

#=
    Problem definition (defined in problem_definition.jl)
=#
include("problem_definition.jl");
# <- def_W

#=
    Include constraints (reference point to fit)
=#
include("constraints.jl")
# <- constraints

#=
    Get initial guess, lower bounds, upper bounds, id's and dimensionality of the problem
=#
x, lb, ub, id  = get_pb_dimensionality(def_W)



#=
    Call optimization routine (using Optim)
    -> note that gradient-free optimizer has to be choosen
=#
# opt = Opt(:LN_BOBYQA, length(x))
opt = Opt(:LN_PRAXIS, length(x))
# opt = Opt(:LN_COBYLA, length(x))
# opt = Opt(:LN_SBPLX, length(x))


opt.lower_bounds    = lb;               # lower bounds
opt.upper_bounds    = ub;               # upper bounds
opt.xtol_rel        = 1e-6;             # variable tolerance (stopping criteria)
g                   = [];               # here we don't use gradient based minimization, but a gradient empty variable still needs to be passed to NLopt

min_objective!(opt, (x,g) -> objective_function(    x,
                                                    id,
                                                    def_W,
                                                    constraints,

                                                    database,
                                                    0           ) )

(minf,minx,ret) = optimize(opt, x)

mfsol = objective_function(  minx,
                             id,
                             def_W,
                             constraints,

                             database,
                             1       )




#=
    save some tests results using default constraints
    
    -----------------------------------
            Missif       |  Algorithm
    -----------------------------------
    0.08410554083169643     LN_BOBYQA
    0.014417795440396328    LN_PRAXIS

=#

