# This script helps to generate a list of points for testing MAGEMin using reference built-in bulk-rock compositions

cur_dir = pwd();    
if  cur_dir[end-3:end]=="test"
    cd("../")           # change to main directory if we are in /test
end
using MAGEMin_C

# Initialize database 
db          = "igad"  # database: ig, igneous (Holland et al., 2018); mp, metapelite (White et al 2014b)
gv, z_b, DB, splx_data      = init_MAGEMin(db);


sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
test        = 0
gv          = use_predefined_bulk_rock(gv, test, db);

mutable struct outP{ _T  } 
    P           ::  _T
    T           ::  _T 
    test        ::  Int64

    G           ::  _T
    ph          ::  Vector{String}
    ph_frac     ::  Vector{Float64}
end

outList       = Vector{outP}(undef, 0)

# test 0
gv.verbose  = -1    # switch off any verbose

Tmin        = 600.0
Tmax        = 1200.0
Tstep       = 100.0
Pmin        = 0.01
Pmax        = 10.01
Pstep       = 2.0

for i=Pmin:Pstep:Pmax
    for j=Tmin:Tstep:Tmax

        out         = point_wise_minimization(i,j, gv, z_b, DB, splx_data, sys_in)
        push!(outList,outP(i,j,test,out.G_system,out.ph,out.ph_frac[:]))
        print(out)
    end
end

for i=1:length(outList)
    print(outList[i],",\n")
end

