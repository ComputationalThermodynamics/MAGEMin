using MAGEMin_C

dtb                     = "mp"
gv, z_b, DB, splx_data  = init_MAGEMin(dtb);

# here the database is initialized with default test value and random PT conditions, this to access W's and G's
gv                      = use_predefined_bulk_rock(gv, 0, dtb);
gv, z_b, DB, splx_data  = pwm_init(8.0, 800.0, gv, z_b, DB, splx_data; G0 = false);

# Here you need to provide the right tc-ds number for the database  you want to access
arr_ptr                 = LibMAGEMin.get_arr_em_db_tc(gv.EM_dataset)   # Returns Ptr{EM_db}  # for tc-ds62

# The following loads all structs into a Julia array
arr                     = [unsafe_load(arr_ptr, i) for i in 1:gv.n_em_db] 

# The following extracts the names as strings
em_names = [let nb=collect(arr[i].Name); np=findfirst(==(0x00), reinterpret(UInt8, nb)); String(reinterpret(UInt8, nb)[1:(np === nothing ? length(nb) : np-1)]) end for i in 1:gv.n_em_db]

#= val content:
    enthalpy 		    = input_1[0];
    entropy  		    = input_1[1];
    volume   		    = input_1[2];

    cpa      		    = input_2[0];
    cpb      		    = input_2[1];
    cpc      		    = input_2[2];
    cpd      		    = input_2[3];

    alpha0   		    = input_3[0];
    kappa0   		    = input_3[1];
    kappa0p  		    = input_3[2];
    kappa0pp 		    = input_3[3];	
=#

id_fa                   = findfirst(isequal("fa"), em_names) -1  # find index of olivine; -1 for the C memory starting at 0

val                     = unsafe_load(arr_ptr, id_fa)  # Get the struct at index 1
HSV                     = val.input_1 # store H, S and V just for clarity

fac                     = 0.99 # modify H, S and V, e.g. reduce by 1%
setfield!(val, :input_1, tuple((HSV[i]*fac for i in 1:3)...))  # Modify the input_1 field with a new tuple

# Now store the modified struct back to the array. Note that id_fa is used to store it back in the same position
# Also note that saving the original values of the variables is important, as unsafe_store! overwrites the entire struct
unsafe_store!(arr_ptr, val, id_fa)

gv      = LibMAGEMin.ComputeG0_point(gv.EM_database, z_b, gv, DB.PP_ref_db,DB.SS_ref_db);


# then: pwm_run
