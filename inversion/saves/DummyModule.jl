# module DummyModule
    # include("../gen/magemin_library.jl")
    # include("../julia/MAGEMin_wrappers.jl")

    function pwm(gv, z_b, DB, splx_data, P,T,db)
        
        # Initialize database 
        
        # gv, z_b, DB, splx_data      = init_MAGEMin(db);

        sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
        test        = 0         #KLB1
        gv          = use_predefined_bulk_rock(gv, test, db);

        # Call optimization routine for given P & T & bulk_rock
        P;
        T;
        gv.verbose  = -1        # switch off any verbose
        out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in);
        # @show out
        # finalize_MAGEMin(gv,DB)
        return out.G_system;
    end

# end