


# # @everywhere function init_db()

# #     # Initialize database 
# #     db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
# #     gv, z_b, DB, splx_data      = init_MAGEMin(db);

# #     # sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
# #     test        = 0         #KLB1
# #     gv          = use_predefined_bulk_rock(gv, test, db);
# #     return (gv, z_b, DB, splx_data,db)
# #     # return gv
# # end

# # for i=1:nworkers()
# #         @spawnat i gv, z_b, DB, splx_data,db = init_db()
# # end


# # gv, z_b, DB, splx_data,db = init_db()
# # sendto(workers(), gv=gv)
# db          = "ig"
# gv, z_b, DB, splx_data      = init_MAGEMin(db);
# test        = 0         #KLB1
# gv          = use_predefined_bulk_rock(gv, test, db);
# sendto(workers(), gv=gv)
# sendto(workers(), z_b=z_b)
# sendto(workers(), DB=DB)
# sendto(workers(), splx_data=splx_data)


# @distributed (+) for i = 1:4
#     # gv, z_b, DB, splx_data,db = init_db()
#     sys_in      = "mol"
#     P           = 8.0
#     T           = 800.0
#     out  = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in);
#     G    = out.G_system;
#     # GC.gc()
# end


# out = @spawnat 2 point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)






# # @everywhere gv, z_b, DB, splx_data,db = init_db()

# # gv, z_b, DB, splx_data,db = fetch(1)
# @distributed (+) for i = 1:8
#     # gv, z_b, DB, splx_data,db = init_db()
#     sys_in      = "mol"
#     P           = 8.0
#     T           = 800.0
#     out  = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)
#     G    = out.G_system;
#     # GC.gc()
# end

# @everywhere GC.gc()

# @spawnat 1 begin
#     GC.gc()
# end










    # # Call optimization routine for given P & T & bulk_rock
    # P           = 8.0
    # T           = 800.0
    # # n = 1e2
    # # for i=1:n
    # # gv.verbose  = -1        # switch off any verbose
    # out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)

    # print("point $i G -> $(out.G_system)\n")
    #     if i % 100 == 0
    #         GC.gc()
    #     end
    # end












# using Distributed
# addprocs(6)
# @everywhere include("gen/magemin_library.jl")
# @everywhere include("julia/MAGEMin_wrappers.jl")

# np = nworkers();

# @everywhere include("inversion/DummyModule.jl")

# @everywhere db = "ig"
# @everywhere gv, z_b, DB, splx_data      = init_MAGEMin(db);

# @everywhere pwm(gv, z_b, DB, splx_data,8.0,800.0,db)

# wp = WorkerPool([2, 3]);
# A = rand(3000);
# f = remotecall(maximum, wp, A)

# Pkg.activate(joinpath(@__DIR__, ".", "examples"))
# # this might be a way to initialize and reuse workers...
# for i=1:np

#         # gv = @fetchfrom i gv
#         # z_b = @fetchfrom i z_b
#         # DB = @fetchfrom i DB
#         # splx_data = @fetchfrom i splx_data
#         # db = @fetchfrom i db
#         # fetch(@spawnat i pwm(gv, z_b, DB, splx_data,8.0,800.0,db))

#     fetch(@spawnat i pwm(@fetchfrom i gv, @fetchfrom i z_b, @fetchfrom i DB, @fetchfrom i splx_data,8.0,800.0,@fetchfrom i db))
# end





# db = "ig"
# gv, z_b, DB, splx_data      = init_MAGEMin(db);

# P = [1,2,3,4,5,6,7,8,9,10,11,12]
# T = [800,800,800,800,800,800,800,800,800,800,800,800]

# np = 12

# @time G  = @distributed (+) for i = 1:np
#     pwm(gv, z_b, DB, splx_data, 8.0,800.0,db)
# end



# global G = 0.0

# @time for i = 1:np
#     global G +=fetch(@spawnat :any pwm(P[i],T[i],db))
# end


# using MPI


# # using Plots

# #=
#     Initialize MPI (run single point minimization in parallel to speed-up calculation)
# =#
# MPI.Init()
# comm = MPI.COMM_WORLD
# print("Hello world, I am rank $(MPI.Comm_rank(comm)) of $(MPI.Comm_size(comm))\n")


# rank = MPI.Comm_rank(comm)
# size = MPI.Comm_size(comm)

# using MAGEMin_C
# # using NLopt

# Pr = [2.,4.,6.,8.];
# Tr = [700.,800.,900.,1000.];

# # Initialize database 
# db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
# gv, z_b, DB, splx_data      = init_MAGEMin(db);

# sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
# test        = 0         #KLB1
# gv          = use_predefined_bulk_rock(gv, test, db);

# # Call optimization routine for given P & T & bulk_rock
# P           = Pr[rank+1];
# T           = Tr[rank+1];
# gv.verbose  = -1        # switch off any verbose
# out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in);
# @show out

# finalize_MAGEMin(gv,DB)

# MPI.Barrier(comm)




# nheads = 0;

# # @everywhere include_string(Main, $(read("inversion/count_heads.jl", String)), "count_heads.jl")

# @everywhere function count_heads(n)
#     c::Int = 0
#     for i = 1:n
#         c += rand(Bool)
#     end
#     c
# end	

# nheads = @distributed (+) for i = 1:10
#     count_heads(1e9)
# end
# for ite=1:10
# a = @spawnat :any count_heads(1e9)
# b = @spawnat :any count_heads(1e9)
# c = @spawnat :any count_heads(1e9)
# d = @spawnat :any count_heads(1e9)
# e = @spawnat :any count_heads(1e9)
# f = @spawnat :any count_heads(1e9)

# nheads += fetch(a)+fetch(b)+fetch(c)+fetch(d)+fetch(e)+fetch(f)
# end
# print(nheads)

# julia -p 12


# @everywhere using MAGEMin_C
# @everywhere function pwm()
    
#     # Initialize database 
#     db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
#     gv, z_b, DB, splx_data      = init_MAGEMin(db);

#     sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
#     test        = 0         #KLB1
#     gv          = use_predefined_bulk_rock(gv, test, db);

#     # Call optimization routine for given P & T & bulk_rock
#     P           = 5.0;
#     T           = 800.0;
#     gv.verbose  = -1        # switch off any verbose
#     out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in);
#     @show out

#     return out.G_system;

#     finalize_MAGEMin(gv,DB)
# end

# np = 10
# @distributed (+) for i = 1:np
#     DummyModule.pwm()
# end

# print(G/np)


# # where i is the process!
# put!(RemoteChannel(2), pwm(5.0,800.0))

# np = 50
# @time for i = 1:np
#     put!(RemoteChannel(), pwm(5.0,800.0))
#     # pwm(5.0,800.0)
# end
# # @spawnat :any pwm(5.0,800.0)

# fetch(@spawnat 2 pwm(5.0,800.0))



# using MAGEMin_C         # load MAGEMin (needs to be loaded from main directory to pick up correct library in case it is locally compiled)

# # Initialize database 
# db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
# gv, z_b, DB, splx_data      = init_MAGEMin(db);

# sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
# test        = 0         #KLB1
# gv          = use_predefined_bulk_rock(gv, test, db);


include("gen/magemin_library.jl") 
include("julia/MAGEMin_wrappers.jl")
# Initialize database 
db          = "ig"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
gv, z_b, DB, splx_data      = init_MAGEMin(db);

sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
test        = 0         #KLB1
gv          = use_predefined_bulk_rock(gv, test, db);


# Call optimization routine for given P & T & bulk_rock
P           = 8.0
T           = 800.0
n = 1e4
for i=1:n
gv.verbose  = -1        # switch off any verbose
out         = point_wise_minimization(P,T, gv, z_b, DB, splx_data, sys_in)

print("point $i G -> $(out.G_system)\n")
#     if i % 100 == 0
#         GC.gc()
#     end
end
