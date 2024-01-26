
function julia_initialize_MAGEMin()

    # Initialize database 
    db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
    gv, z_b, DB, splx_data = init_MAGEMin(db);

    sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
    test        = 0         #KLB1
    gv          = use_predefined_bulk_rock(gv, test, db);

    return gv, z_b, DB, splx_data, sys_in
end

function minimize_point(gv, z_b, DB, splx_data, sys_in, P, T)

    gv.verbose  = -1        # switch off any verbose
    out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in);

    return out.G_system;
end

function julia_finalize_MAGEMin(gv,DB)
    finalize_MAGEMin(gv,DB)
end


function fake_min_function(Pr,Tr,   gv, z_b, DB, splx_data, sys_in)

    # gv, z_b, DB, splx_data, sys_in  = julia_initialize_MAGEMin()
    G_system                        = minimize_point(gv, z_b, DB, splx_data, sys_in,Pr,Tr)

    # julia_finalize_MAGEMin(gv,DB)
    
    return G_system;
end

# the goal is to run this function on #0 and send point on other # then reduce to #0
function fake_obj(Pr,Tr,size,ppp,   gv, z_b, DB, splx_data, sys_in, x)

    # msg to continue parallel computation
    for i=2:size
        dst                 = i-1;                                          # target processor
        send_cont           = Array{Int64}(undef, 1)                        # file right point list
        send_cont[1]        = 1;
        creq                = MPI.Isend(send_cont, dst, rank+32, comm)      # send list
    end 
    MPI.Barrier(comm)

    # send initial guess to the sidekicks
    for i=2:size
        dst                 = i-1;                                          # target processor
        send_x              = Array{Float64}(undef, length(x))              # initialize right size
        send_x             .= x;                                            # file right point list
        sx                  = MPI.Isend(send_x, dst, rank+32, comm)         # send list
    end
    MPI.Barrier(comm)

    # send point list to sidekick processors
    for i=2:size
        dst                 = i-1;                                          # target processor
        send_mesg           = Array{Int64}(undef, length(ppp[i]))       # initialize right size
        send_mesg          .= ppp[i];   
        # send_mesg[1]        = 1;                                            # file right point list
        sreq                = MPI.Isend(send_mesg, dst, rank+32, comm)      # send list
    end

    # run other point on main processor
    missfit = 0.0;    
    for i=1:length(ppp[1])
        G = fake_min_function(Pr[ppp[1][i]],Tr[ppp[1][i]],   gv, z_b, DB, splx_data, sys_in)
        # print("Gsys = $G; on rank #$rank\n")
        missfit += G
    end

    MPI.Barrier(comm)

    return missfit
end

# the goal is to run this function on #0 and send point on other # then reduce to #0
function fake_NLopt(Pr,Tr,size,ppp,   gv, z_b, DB, splx_data, sys_in, x)

    n_ite = 16;
    for i=1:n_ite
        print("iteration #$i\n")
        missfit     = fake_obj(Pr,Tr,size,ppp,   gv, z_b, DB, splx_data, sys_in, x)
        fake_mf     = MPI.Reduce(missfit, +, 0, comm)
        print("   fake missfit (sum G) $fake_mf\n\n")
    end

    # msg to kill parallel computation
    for i=2:size
        dst                 = i-1;                                          # target processor
        send_cont           = Array{Int64}(undef, 1)                        # file right point list
        send_cont[1]        = 0;
        creq                = MPI.Isend(send_cont, dst, rank+32, comm)      # send list
    end 
       
    MPI.Barrier(comm)
end


function fake_side(Pr,Tr,ppp,   gv, z_b, DB, splx_data, sys_in, len_x)

    comm        = MPI.COMM_WORLD
    rank        = MPI.Comm_rank(comm)
    recv_mesg   = Array{Int64}(undef, length(ppp[rank+1]))  
    recv_x      = Array{Float64}(undef, len_x)  
    recv_cont   = Array{Int64}(undef, 1)  
    compute     = 1

    while compute == 1

        creq        = MPI.Irecv!(recv_cont, 0,  32, comm)
        MPI.Barrier(comm) 
        print(" received continue on rank $rank is $creq[1]\n")
        compute = recv_cont[1];

        if compute == 1
            xreq        = MPI.Irecv!(recv_x, 0,  32, comm)
            MPI.Barrier(comm) 
            # print(" received guess on rank $rank is $recv_x\n")

            rreq        = MPI.Irecv!(recv_mesg, 0,  32, comm)
            points      = recv_mesg;  
            stats       = MPI.Wait(rreq)

            missfit     = 0.0;
            print(" compute = $compute\n")

                for i=1:length(ppp[rank+1])
                    G = fake_min_function(Pr[ppp[rank+1][i]],Tr[ppp[rank+1][i]],   gv, z_b, DB, splx_data, sys_in)
                    # print("Gsys = $G; on rank #$rank\n")
                    missfit += G
                end

            MPI.Barrier(comm)
            fake_mf = MPI.Reduce(missfit, +, 0, comm)

        end
    end
end


