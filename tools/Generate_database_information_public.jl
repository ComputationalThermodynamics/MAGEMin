using MAGEMin_C


function get_database_infos()
    
    db_details      = [ "Metapelite (White et al., 2014)",
                        "Metabasite (Green et al., 2016)",
                        "Igneous (Green et al., 2025, corrected after Holland et al., 2018)",
                        "Igneous alkaline dry (Weller et al., 2024)",
                        "Ultramafic (Evans & Frost., 2021)",
                        "Ultramafic extended (Evans & Frost., 2021) with pl, amp and aug from Green et al., 2016",
                        "Mantle (Holland et al., 2013)",
                        "Metapelite extended (White et al., 2014 with po from Evans & Frost., 2021, amp dio and aug from Green et al., 2016)",
                        "Stixrude & Lithgow-Bertelloni (2011)" ]

    database_list   = ["mp","mb","ig","igad","um","ume","mtl","mpe","sb11"]

    db_inf          = Array{db_infos, 1}(undef, length(database_list))

    for k in eachindex(database_list)
        datab         = database_list
        gv, z_b, DB, splx_data  = init_MAGEMin(datab; mbCpx = 1);
        gv          =   use_predefined_bulk_rock(gv, 0, datab);
        gv.verbose  =  -1
        P, T        =   8.0,800.0
        gv, z_b, DB, splx_data = pwm_init(P, T, gv, z_b, DB, splx_data);
    
        ss_struct  = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);
        ss_names   = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss));
        pp_names   = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.PP_list, gv.len_pp));
    
        ss = Array{ss_infos, 1}(undef, gv.len_ss)
    
        for i=1:gv.len_ss
            n_em 	= ss_struct[i].n_em
            n_xeos 	= ss_struct[i].n_xeos
            fName   = unsafe_string.(ss_struct[i].fName)
            
            em_names   = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, ss_struct[i].EM_list, n_em))
    
            em_names2  = Vector{String}(undef,n_em+1)
            em_names2[1] = "none"
            em_names2[2:end] = em_names
    
            xeos_names = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, ss_struct[i].CV_list, n_xeos))
            xeos_names2  = Vector{String}(undef,n_xeos+1)
            xeos_names2[1] = "none"
            xeos_names2[2:end] = xeos_names
    
            ss[i]   = ss_infos(fName, ss_names[i], n_em, n_xeos, em_names2, xeos_names2)
        end
    
        db_inf[k] = db_infos(database_list[k],db_details[k],ss,ss_names,pp_names)
        # finalize_MAGEMin(gv,DB, z_b)
    end
    
    return db_inf
end


mutable struct ss_infos
    ss_fName:: String
    ss_name :: String
    n_em    :: Int64
    n_xeos  :: Int64
    ss_em   :: Vector{String}
    ss_xeos :: Vector{String}
end


mutable struct db_infos
    db_name :: String
    db_info :: String
    data_ss :: Array{ss_infos}
    ss_name :: Array{String}
    data_pp :: Array{String}
end


db_inf = get_database_infos()

print(db_inf)

