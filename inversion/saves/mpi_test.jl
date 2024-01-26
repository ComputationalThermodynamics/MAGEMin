using MPI
MPI.Init()


# infos:

# add MPI
# julia --project -e 'using Pkg; Pkg.add("MPIPreferences")'
# julia --project -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
# using MPI
# MPI.install_mpiexecjl()
# alias mpiexecjl="/home/nico/.julia/bin/mpiexecjl"

## julia --project -e 'ENV["JULIA_MPI_BINARY"]="system"; using Pkg; Pkg.build("MPI"; verbose=true)'
## using Pkg; Pkg.build("MPI"; verbose=true)

# locate mpiexec of system to export it:
# export PATH="$PATH:/home/nico/julia-1.9.2-linux-x86_64/julia-1.9.2/bin"
# export JULIA_MPI_PATH=/usr/bin/



# load MAGEMin libraries
using NLopt
include("../gen/magemin_library_mpi.jl")
include("../julia/MAGEMin_wrappers.jl")
include("./mpi_test_functions.jl")

# distributed point informations
Pr          = [2.,4.,6.,8.,10.,12.,14.,16.,18.,20.,22.,24.,28.,30.,32.,34.];
Tr          = [700.,800.,900.,1000.,700.,800.,900.,1000.,700.,800.,900.,1000.,800.,900.,1000.,700.];
npoints     = length(Pr);


comm        = MPI.COMM_WORLD
rank        = MPI.Comm_rank(comm)
size        = MPI.Comm_size(comm)

# here we compute the list of point of each rank, on rank 0 (ppp = points per processor)
ppp     = [];
for p=1:size
    pp  = [];
    for i=1:length(Pr)
        if mod(i,size) == (p-1)
            push!(pp,i)
        end
    end
    push!(ppp,pp)
end

gv, z_b, DB, splx_data, sys_in  = julia_initialize_MAGEMin()

x               = zeros(12) .+ rand()                       # fake example of variables to be send to all workers

if rank == 0
    fake_NLopt(Pr,Tr,size,ppp,   gv, z_b, DB, splx_data, sys_in, x)
else
    fake_side(Pr,Tr,ppp,   gv, z_b, DB, splx_data, sys_in, length(x))
end
# print(" Fake missfit is $missfit\n")

MPI.Barrier(comm)

julia_finalize_MAGEMin(gv,DB)
    

