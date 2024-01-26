using MPI
MPI.Init()

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
println("Hello world, I am the MPI process no $rank.")
MPI.Barrier(comm)

MPI.Finalize()