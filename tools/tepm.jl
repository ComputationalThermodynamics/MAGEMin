# if tepm == "true"
#     Out_TE_XY       = Vector{MAGEMin_C.tepm_struct{Float64}}(undef,length(data.x))
#     Out_TE_XY_new   = Vector{MAGEMin_C.tepm_struct{Float64}}(undef,n_new_points)
# end
    if tepm == 1
        Out_PT_TE = Vector{tepm_struct{Float64}}(undef, length(P))
    end
    # here we compute trace element partitioning and zircon saturation
    if (tepm == 1)
        if (out.frac_M > 0.0 && out.frac_S > 0.0)
            Cliq, Cmin, ph_TE, ph_wt_norm, liq_wt_norm, Cliq_Zr, te_names  = compute_TE_partitioning(   te_X,
                                                                                                        out,
                                                                                                        dtb;
                                                                                                        TE_db = te_db)

            # Then we compute zirconium saturation
            Sat_zr_liq  = zirconium_saturation( out; 
                                                model = zr_sat)   

            if Cliq_Zr > Sat_zr_liq
                zircon_wt, SiO2_wt, O_wt  = adjust_bulk_4_zircon(Cliq_Zr, Sat_zr_liq)
                SiO2_id     = findall(out.oxides .== "SiO2")[1]

                bulk_act    = copy(out.bulk_wt)
                bulk_act[SiO2_id]    = out.bulk_wt[SiO2_id] - SiO2_wt 
                bulk_act  ./= sum(bulk_act)
                gv          = define_bulk_rock(gv, bulk_act, out.oxides, "wt", dtb);
                mSS_vec     = deepcopy(out.mSS_vec)
                out_cor     = point_wise_minimization_with_guess(mSS_vec, P, T, gv, z_b, DB, splx_data)

                Cliq, Cmin, ph_TE, ph_wt_norm, liq_wt_norm, Cliq_Zr, te_names = compute_TE_partitioning(    te_X,
                                                                                                            out_cor,
                                                                                                            dtb;
                                                                                                            TE_db = te_db)

                # Then we compute zirconium saturation
                Sat_zr_liq  = zirconium_saturation( out; 
                                                    model = zr_sat)     

                zircon_wt, SiO2_wt, O_wt  = adjust_bulk_4_zircon(Cliq_Zr, Sat_zr_liq)
            else
                zircon_wt = 0.0;
            end

        elseif out.frac_M == 1.0
            TE_dtb      =  get_TE_database("TE_OL_felsic")

            Sat_zr_liq  = zirconium_saturation( out; 
                                                model = zr_sat)     

            zircon_wt, SiO2_wt, O_wt  = adjust_bulk_4_zircon(Cliq_Zr, Sat_zr_liq)

            te_names    = TE_dtb.element_name
            ph_TE       = nothing
            Cliq        = te_X 
            Cmin        = nothing
            te_db       = te_db
            zr_sat      = zr_sat
        else 
            TE_dtb      =  get_TE_database("TE_OL_felsic")

            te_names    = TE_dtb.element_name
            ph_TE       = nothing
            Cliq        = nothing 
            Cmin        = nothing
            te_db       = te_db
            zr_sat      = zr_sat
            Sat_zr_liq  = nothing
            zircon_wt   = 0.0
        end

        out_te = tepm_struct{Float64}(  te_names, ph_TE, te_X, Cliq, Cmin,
                                        te_db, zr_sat, 
                                        Sat_zr_liq, zircon_wt)

        return out, out_te
    else
        return out
    end