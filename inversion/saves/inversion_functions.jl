export  get_str_id, misfit_calc, display_ss_infos, objective_function


function objective_function(    x,
                                id,
                                def_W,
                                constraints,

                                database,
                                disp        )
    gv, z_b, DB, splx_data  = init_MAGEMin(database);

    gv.verbose  = -1        # switch off any verbose

    mf = 0.0;
    for c=1:length(constraints)    # loop through calibration points to compute total missift

        gv  = define_bulk_rock(gv, constraints[c].X, constraints[c].ox, constraints[c].sys, database);

        P          = constraints[c].P;
        T          = constraints[c].T;

        # initialize MAGEMin up to G0 and W's point
        gv, z_b, DB, splx_data = pwm_init(P, T, gv, z_b, DB, splx_data);

        # get the solution phase structure (size gv.len_ss)
        ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref}, DB.SS_ref_db,gv.len_ss);

        #here update W'S
        for i=1:3:length(x)
 
            ss      = id[i,1];      # get solution phase id
            w_n     = id[i,2];      # get W's id

            Wnew    = x[i] + x[i+1]*P + x[i+2]*T; 

            W       = unsafe_wrap(Vector{Cdouble},ss_struct[ss].W, ss_struct[ss].n_w)
            W[w_n]  = Wnew;
        end

        out       = pwm_run(gv, z_b, DB, splx_data);

        mf       += misfit_calc(constraints[c],out,disp);
    end
    print("mf: ",mf,"\n")
    return mf;
    finalize_MAGEMin(gv,DB)
end

function get_pb_dimensionality(def_W)
    n_x  = length(def_W)*3;
    x    = zeros(n_x)                           # define initial guess
    lb   = zeros(n_x)                           # define lower bound
    ub   = zeros(n_x)                           # define upper bound
    id   = zeros(Int8,n_x,2)                    # define id matrix, ss_n, W_n, variable (val, Pdep, Tdep)
    ix   = 1;
    for i=1:length(def_W)
        x[ix]   = def_W[i].val;             id[ix,1]   = def_W[i].ss_num; id[ix,2]   = def_W[i].W_num;
        x[ix+1] = def_W[i].P_dep;           id[ix+1,1] = def_W[i].ss_num; id[ix+1,2] = def_W[i].W_num;
        x[ix+2] = def_W[i].T_dep;           id[ix+2,1] = def_W[i].ss_num; id[ix+2,2] = def_W[i].W_num;

        lb[ix]   = def_W[i].val_range[1];
        lb[ix+1] = def_W[i].P_dep_range[1];
        lb[ix+2] = def_W[i].T_dep_range[1];

        ub[ix]   = def_W[i].val_range[2];
        ub[ix+1] = def_W[i].P_dep_range[2];
        ub[ix+2] = def_W[i].T_dep_range[2];
        ix += 3;
    end

    return x, lb, ub, id;
end


function misfit_calc(   constraints,
                        out,
                        disp             )
    # reference data for comparison
    ph_ref      = constraints.ph;           #stable phase names
    frac_ref    = constraints.ph_frac;      #stable phase fractions

    # newly calculated point
    ph_in       = out.ph;
    frac_in     = out.ph_frac;

    # sort points
    ind_ref     = sortperm(ph_ref);
    ind_in      = sortperm(ph_in);

    a           = ph_ref[ind_ref];
    b           = ph_in[ind_in];
    af          = frac_ref[ind_ref];
    bf          = frac_in[ind_in];

    # find phases that differ from one set to another
    ph_mf = 0.0;

    k = intersect(a,b)
    l = symdiff(a,b)

    if ~isempty(l)
        for i=1:length(l)
            if issubset([l[i]],a)
                id      = get_str_id(l[i],a)
                ph_mf  += frac_ref[ind_ref[id]]*10.0;
            else
                id      = get_str_id(l[i],b)
                ph_mf  += frac_in[ind_in[id]]*10.0;
            end
        end
    end

    if ~isempty(k)

        for i=1:length(k)
            ida      = get_str_id(k[i],a);
            idb      = get_str_id(k[i],b);

            d        = abs(bf[idb] - af[ida]);

            if (d > 1e-4)
                ph_mf   += d;
            end
        end
    end

    if disp == 1
        print(a," ",af,"\n")
        print(b," ",bf,"\n\n")
    end
    return ph_mf;
end


function display_ss_infos( database )

    gv, z_b, DB, splx_data   = init_MAGEMin(database);

    # here the database is initialized with default test value and random PT conditions, this to access W's and G's
    gv          = use_predefined_bulk_rock(gv, 0, database);
    gv, z_b, DB, splx_data = pwm_init(8.0,800.0, gv, z_b, DB, splx_data);

    ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss));
    ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);
    
    print("\n---------------------------------\n");
    print("    Database information (",database,")\n");
    print("----------------------------------\n");
    for i=1:gv.len_ss
        print("   ",ss_names[i],": ",i,"; n_W's, ",ss_struct[i].n_w,"; n_em's, ",ss_struct[i].n_em,"\n")
    end
    print("----------------------------------\n\n");
    finalize_MAGEMin(gv,DB)
end


function get_str_id(str,str_array)
    id = -1;
    for i=1:length(str_array)
        if str == str_array[i];
            id = i;
            break;
        end
    end

    return id
end
