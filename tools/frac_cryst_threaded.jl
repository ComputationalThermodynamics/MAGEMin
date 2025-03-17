#=
Test file to run threaded fractional crystallization model 
NR - 5/3/2025
=#
using ProgressMeter
using MAGEMin_C
using Base.Threads: @threads

function get_data_thread( MAGEMin_db :: MAGEMin_Data )

    id          = Threads.threadid()
    gv          = MAGEMin_db.gv[id]
    z_b         = MAGEMin_db.z_b[id]
    DB          = MAGEMin_db.DB[id]
    splx_data   = MAGEMin_db.splx_data[id]
    
   return (gv, z_b, DB, splx_data)
end

function example_of_threaded_MAGEMin_calc(  data_thread :: Tuple{Any, Any, Any, Any}, dtb :: String,

                                            starting_P :: Float64,
                                            starting_T :: Float64,
                                            ending_T   :: Float64,
                                            n_steps    :: Int64,

                                            sys_in     :: String,
                                            bulk       :: Vector{Float64},
                                            Xoxides    :: Vector{String}           )

    gv, z_b, DB, splx_data = data_thread        # Unpack the MAGEMin data

    Out_PT = Vector{MAGEMin_C.gmin_struct{Float64, Int64}}(undef, n_steps)
    gv      = define_bulk_rock(gv, bulk, Xoxides, sys_in, dtb);

    for i = 1:n_steps


        P       = Float64(starting_P)
        T       = Float64(starting_T - (starting_T - ending_T) * (i-1)/(n_steps-1))

        out     = point_wise_minimization(  P, T, gv, z_b, DB, splx_data;
                                                name_solvus=true)                 #Gi=Gi, scp=scp, rm_list=rm_list

        Out_PT[i] = deepcopy(out)

        if "liq" in out.ph 
            bulk    = out.bulk_M
            oxides  = out.oxides

            gv      = define_bulk_rock(gv, bulk, oxides, "mol", dtb);
        end

    end

    return Out_PT

end

function perform_threaded_calc( Out_all     :: Vector{Vector{MAGEMin_C.gmin_struct{Float64, Int64}}},
                                data        :: MAGEMin_Data,
                                dtb         :: String,
                                n_starting_points :: Int64,
                                starting_P  :: Vector{Float64},
                                starting_T  :: Vector{Float64},
                                ending_T    :: Vector{Float64},
                                n_steps     :: Int64,
                                sys_in      :: String,
                                bulk        :: Matrix{Float64},
                                Xoxides     :: Vector{String} )

    progr = Progress(n_starting_points, desc="Computing $n_starting_points examples of threaded MAGEMin_calc...") # progress meter
    @threads :static for i=1:n_starting_points

        data_thread = get_data_thread(data)
        starting_P_  = starting_P[i]
        starting_T_  = starting_T[i]
        ending_T_    = ending_T[i]
        n_steps_     = n_steps
        bulk_        = bulk[i,:]

        Out_PT = example_of_threaded_MAGEMin_calc(  data_thread, dtb,

                                                    starting_P_,
                                                    starting_T_,
                                                    ending_T_,
                                                    n_steps_,
                                                    sys_in,
                                                    bulk_,
                                                    Xoxides )

        Out_all[i] = Out_PT
        next!(progr)

    end
    finish!(progr)

    return Out_all
end



# first initialize MAGEMin
dtb     = "mp"
data    = Initialize_MAGEMin(dtb, verbose=-1; solver=2);

n_starting_points  = 64

# Allocate memory for the output (Nested_structure where each element is a vector of gmin_struct)
Out_all  =  Vector{Vector{MAGEMin_C.gmin_struct{Float64, Int64}}}(undef, n_starting_points);


starting_P  = [range(1.0,10.0,n_starting_points);]      # 10 starting points
starting_T  = ones(n_starting_points) .* 1300.0      
ending_T    = ones(n_starting_points) .* 600.0  
n_steps     = 128

sys_in      = "wt"
bulk        = repeat([58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.0]', n_starting_points)
Xoxides     = ["SiO2", "TiO2", "Al2O3", "FeO",  "MnO",  "MgO",  "CaO",  "Na2O", "K2O","H2O","O"]

Out_all     = perform_threaded_calc(Out_all, data, dtb, n_starting_points, starting_P, starting_T, ending_T, n_steps, sys_in, bulk, Xoxides);

Finalize_MAGEMin(data)