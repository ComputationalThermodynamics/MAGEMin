using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
size = MPI.Comm_size(comm)

# load MAGEMin libraries
using NLopt
include("../gen/magemin_library_mpi.jl")
include("../julia/MAGEMin_wrappers.jl")

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


function fake_min_function(id,Pr,Tr)

    gv, z_b, DB, splx_data, sys_in  = julia_initialize_MAGEMin()
    G_system                        = minimize_point(gv, z_b, DB, splx_data, sys_in,Pr[id],Tr[id])

    julia_finalize_MAGEMin(gv,DB)
    
    return G_system;
end


# the goal is to run this function on #0 and send point on other # then reduce to #0
function fake_NLop(mf_array,Pr,Tr,npoints)
    p = 1;
    
    # nloop = 4;
    # for loop=1:nloop    

    for point=1:npoints
        for i=1:size-1

            send_mesg = Array{Int64}(undef, 2)
            send_mesg[1] = i; 
        
            if loop==nloop
                send_mesg[2] = 0;
            else
                send_mesg[2] = 1;
            end

            dst=i;
            sreq = MPI.Isend(send_mesg, dst, rank+32, comm)
        end

        G = fake_min_function(rank  + 1,Pr,Tr)
        mf_array[p] = G
        p += 1;
        print("Gsys = $G; on rank #$rank\n")
    end

end


function fake_side(mf_array,Pr,Tr)
    compute = 1
    p       = 1;
    while compute == 1 

        recv_mesg   = Array{Int64}(undef, 2)
        rreq        = MPI.Irecv!(recv_mesg, 0,  32, comm)
        stats       = MPI.Wait(rreq)

        G = fake_min_function(recv_mesg[1]  + 1,Pr,Tr)
        mf_array[p] = G
        p          += 1;
        print("Gsys = $G; on rank #$rank\n")

        if recv_mesg[2] == 0
            compute = 0
        end

    end
end



# distributed point informations
Pr      =  [2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24];
Tr      =  [700.,800.,900.,1000.,700.,800.,900.,1000.,700.,800.,900.,1000.];


npoints     = length(Pr);
# npoints     = 16;
mf_array    = zeros(npoints)

if rank == 0
    fake_NLop(mf_array,Pr,Tr,npoints)
else
    fake_side(mf_array,Pr,Tr)
end

print(mf_array,"\n")


MPI.Barrier(comm)



