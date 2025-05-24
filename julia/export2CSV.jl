#=~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#   Project      : MAGEMin_C
#   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
#   Developers   : Nicolas Riel, Boris Kaus
#   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
#   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
#   Contact      : nriel[at]uni-mainz.de
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ =#

"""
MAGEMin_dataTE2dataframe( out:: Union{Vector{out_tepm}, out_tepm},dtb,fileout)

    Transform MAGEMin trace-element output into a dataframe for quick(ish) save

"""
function MAGEMin_dataTE2dataframe(  out     :: Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}},
                                    out_te  :: Union{Vector{out_tepm}, out_tepm},
                                    dtb,
                                    fileout)

    # here we fill the dataframe with the all minimized point entries
    if typeof(out) == MAGEMin_C.gmin_struct{Float64, Int64}
        out = [out]
    end
    if typeof(out_te) == MAGEMin_C.out_tepm
        out_te = [out_te]
    end
    np        = length(out)

    db_in     = retrieve_solution_phase_information(dtb)
    datetoday = string(Dates.today())
    rightnow  = string(Dates.Time(Dates.now()))

    metadata   = "# MAGEMin_TE " * " $(out[1].MAGEMin_ver);" * datetoday * ", " * rightnow * "; using database " * db_in.db_info * "\n"

    # Here we create the dataframe's header:
    MAGEMin_db = DataFrame(         Symbol("point[#]")          => Int64[],
                                    Symbol("X[0.0-1.0]")        => Float64[],
                                    Symbol("P[kbar]")           => Float64[],
                                    Symbol("T[°C]")             => Float64[],
                                    Symbol("phase")             => String[],
                                    Symbol("mode[wt%]")         => Float64[],
                                    Symbol("Zr_sat[μg/g]")      => Float64[],
    )
    for i in out[1].oxides
        col = i*"_cor[wt%]"
        MAGEMin_db[!, col] = Float64[] 
    end

    for i in out_te[1].elements
        col = i*"_[μg/g]"
        MAGEMin_db[!, col] = Float64[] 
    end

    print("\noutput path: $(pwd())\n")
    @showprogress "Saving data to csv..." for k=1:np

        part_1 = Dict(  "point[#]"      => k,
                        "X[0.0-1.0]"    => out[k].X[1],
                        "P[kbar]"       => out[k].P_kbar,
                        "T[°C]"         => out[k].T_C,
                        "phase"         => "system",
                        "mode[wt%]"     => 100.0,
                        "Zr_sat[μg/g]"  => "-")

        if ~isnothing(out_te[k].bulk_cor_wt)
            part_2 = Dict(  (out[1].oxides[j]*"_cor[wt%]" => out_te[k].bulk_cor_wt[j]*100.0)
                            for j in eachindex(out[1].oxides))
        else
            part_2 = Dict(  (out[1].oxides[j]*"_cor[wt%]" => "-")
                            for j in eachindex(out[1].oxides))
        end

        part_3 = Dict(  ( out_te[1].elements[j]*"_[μg/g]" => out_te[k].C0[j])
                        for j in eachindex(out_te[1].elements))

        row    = merge(part_1,part_2,part_3)   

        push!(MAGEMin_db, row, cols=:union)
        
        # liquid
        if isnothing(out_te[k].Sat_zr_liq)
            Sat_zr_liq = "-"
        else 
            Sat_zr_liq = out_te[k].Sat_zr_liq
        end

        if ~isnothing(out_te[k].liq_wt_norm)
            part_1 = Dict(  "point[#]"      => k,
                            "X[0.0-1.0]"    => out[k].X[1],
                            "P[kbar]"       => out[k].P_kbar,
                            "T[°C]"         => out[k].T_C,
                            "phase"         => "liq",
                            "mode[wt%]"     => out_te[k].liq_wt_norm .* 100.0,
                            "Zr_sat[μg/g]"  => Sat_zr_liq)

            part_2 = Dict(  (out[1].oxides[j]*"_cor[wt%]" => "-")
                            for j in eachindex(out[1].oxides))

            part_3 = Dict(  ( out_te[1].elements[j]*"_[μg/g]" => out_te[k].Cliq[j])
                            for j in eachindex(out_te[1].elements))

            row    = merge(part_1,part_2,part_3)   

            push!(MAGEMin_db, row, cols=:union)
        end

        # solid
        if ~isnothing(out_te[k].Csol)
            if ~isnothing(out_te[k].liq_wt_norm)
                sol_mode = (1.0-out_te[k].liq_wt_norm)*100.0
            else
                sol_mode = 100.0
            end
            part_1 = Dict(  "point[#]"      => k,
                            "X[0.0-1.0]"    => out[k].X[1],
                            "P[kbar]"       => out[k].P_kbar,
                            "T[°C]"         => out[k].T_C,
                            "phase"         => "sol",
                            "mode[wt%]"     => sol_mode,
                            "Zr_sat[μg/g]"  => "-")

            part_2 = Dict(  (out[1].oxides[j]*"_cor[wt%]" => "-")
                            for j in eachindex(out[1].oxides))

            part_3 = Dict(  ( out_te[1].elements[j]*"_[μg/g]" => out_te[k].Csol[j])
                            for j in eachindex(out_te[1].elements))

            row    = merge(part_1,part_2,part_3)   

            push!(MAGEMin_db, row, cols=:union)
        end
        
        if ~isnothing(out_te[k].ph_TE)
            nph  = length(out_te[k].ph_TE)
            for i=1:nph
                part_1 = Dict(  "point[#]"      => k,
                                "X[0.0-1.0]"    => out[k].X[1],
                                "P[kbar]"       => out[k].P_kbar,
                                "T[°C]"         => out[k].T_C,
                                "phase"         => out_te[k].ph_TE[i],
                                "mode[wt%]"     => out_te[k].ph_wt_norm[i].*100.0,
                                "Zr_sat[μg/g]"  => "-")

                part_2 = Dict(  (out[1].oxides[j]*"_cor[wt%]" => "-")
                                for j in eachindex(out[1].oxides))

                part_3 = Dict(  ( out_te[1].elements[j]*"_[μg/g]" => out_te[k].Cmin[i,j])
                                for j in eachindex(out_te[1].elements))

                row    = merge(part_1,part_2,part_3)   

                push!(MAGEMin_db, row, cols=:union)

            end
        end

    end


    meta = fileout*"_metadata.txt"
    filename = fileout*".csv"
    CSV.write(filename, MAGEMin_db)

    open(meta, "w") do file
        write(file, metadata)
    end

    return nothing
end

"""
MAGEMin_data2dataframe( out:: Union{Vector{MAGEMin_C.gmin_struct{Float64, Int64}}, MAGEMin_C.gmin_struct{Float64, Int64}})

    Transform MAGEMin output into a dataframe for quick(ish) save

"""
function MAGEMin_data2dataframe( out:: Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}},dtb,fileout)

    # here we fill the dataframe with the all minimized point entries
    if typeof(out) == MAGEMin_C.gmin_struct{Float64, Int64}
        out = [out]
    end
    np        = length(out)

    db_in     = retrieve_solution_phase_information(dtb)
    datetoday = string(Dates.today())
    rightnow  = string(Dates.Time(Dates.now()))

    metadata   = "# MAGEMin " * " $(out[1].MAGEMin_ver);" * datetoday * ", " * rightnow * "; using database " * db_in.db_info * "\n"

    # Here we create the dataframe's header:
    MAGEMin_db = DataFrame(         Symbol("point[#]")      => Int64[],
                                    Symbol("X[0.0-1.0]")    => Float64[],
                                    Symbol("P[kbar]")       => Float64[],
                                    Symbol("T[°C]")         => Float64[],
                                    Symbol("phase")         => String[],
                                    Symbol("mode[mol%]")    => Float64[],
                                    Symbol("mode[wt%]")     => Float64[],
                                    Symbol("log10(fO2)")    => Float64[],
                                    Symbol("log10(dQFM)")   => Float64[],
                                    Symbol("aH2O")          => Float64[],
                                    Symbol("aSiO2")         => Float64[],
                                    Symbol("aTiO2")         => Float64[],
                                    Symbol("aAl2O3")        => Float64[],
                                    Symbol("aMgO")          => Float64[],
                                    Symbol("aFeO")          => Float64[],
                                    Symbol("density[kg/m3]")    => Float64[],
                                    Symbol("volume[cm3/mol]")   => Float64[],
                                    Symbol("heatCapacity[J/kg/K]")=> Float64[],
                                    Symbol("alpha[1/K]")    => Float64[],
                                    Symbol("Entropy[J/K]")  => Float64[],
                                    Symbol("Enthalpy[J]")   => Float64[],
                                    Symbol("Vp[km/s]")      => Float64[],
                                    Symbol("Vs[km/s]")      => Float64[],
                                    Symbol("Vp_S[km/s]")      => Float64[],
                                    Symbol("Vs_S[km/s]")      => Float64[],                                    
                                    Symbol("BulkMod[GPa]")  => Float64[],
                                    Symbol("ShearMod[GPa]") => Float64[],
    )

    for i in out[1].oxides
        col = i*"[mol%]"
        MAGEMin_db[!, col] = Float64[] 
    end

    for i in out[1].oxides
        col = i*"[wt%]"
        MAGEMin_db[!, col] = Float64[] 
    end


    print("\noutput path: $(pwd())\n")
    @showprogress "Saving data to csv..." for k=1:np
        np  = length(out[k].ph)
        nss = out[k].n_SS
        npp = out[k].n_PP

        part_1 = Dict(  "point[#]"      => k,
                        "X[0.0-1.0]"    => out[k].X[1],
                        "P[kbar]"       => out[k].P_kbar,
                        "T[°C]"         => out[k].T_C,
                        "phase"         => "system",
                        "mode[mol%]"    => 100.0,
                        "mode[wt%]"     => 100.0,
                        "log10(fO2)"    => out[k].fO2[1],
                        "log10(dQFM)"   => out[k].dQFM[1],
                        "aH2O"          => out[k].aH2O,
                        "aSiO2"         => out[k].aSiO2,
                        "aTiO2"         => out[k].aTiO2,
                        "aAl2O3"        => out[k].aAl2O3,
                        "aMgO"          => out[k].aMgO,
                        "aFeO"          => out[k].aFeO,
                        "density[kg/m3]"    => out[k].rho,
                        "volume[cm3/mol]"   => out[k].V,
                        "heatCapacity[J/kg/K]"=> out[k].s_cp,
                        "alpha[1/K]"    => out[k].alpha,
                        "Entropy[J/K]"  => out[k].entropy,
                        "Enthalpy[J]"   => out[k].enthalpy,
                        "Vp[km/s]"      => out[k].Vp,
                        "Vs[km/s]"      => out[k].Vs,
                        "Vp_S[km/s]"    => out[k].Vp_S,
                        "Vs_S[km/s]"    => out[k].Vs_S,
                        "BulkMod[GPa]"  => out[k].bulkMod,
                        "ShearMod[GPa]" => out[k].shearMod )          

        part_2 = Dict(  (out[1].oxides[j]*"[mol%]" => out[k].bulk[j]*100.0)
                        for j in eachindex(out[1].oxides))

        part_3 = Dict(  (out[1].oxides[j]*"[wt%]" => out[k].bulk_wt[j]*100.0)
                        for j in eachindex(out[1].oxides))
   
        row    = merge(part_1,part_2,part_3)   

        push!(MAGEMin_db, row, cols=:union)

        for i=1:nss
            part_1 = Dict(  "point[#]"      => k,
                            "X[0.0-1.0]"    => out[k].X[1],
                            "P[kbar]"       => out[k].P_kbar,
                            "T[°C]"         => out[k].T_C,
                            "phase"         => out[k].ph[i],
                            "mode[mol%]"    => out[k].ph_frac[i].*100.0,
                            "mode[wt%]"     => out[k].ph_frac_wt[i].*100.0,
                            "log10(fO2)"    => "-",
                            "log10(dQFM)"   => "-",
                            "aH2O"          => "-",
                            "aSiO2"         => "-",
                            "aTiO2"         => "-",
                            "aAl2O3"        => "-",
                            "aMgO"          => "-",
                            "aFeO"          => "-",
                            "density[kg/m3]"    => out[k].SS_vec[i].rho,
                            "volume[cm3/mol]"   => out[k].SS_vec[i].V,
                            "heatCapacity[J/kg/K]"=> out[k].SS_vec[i].cp,
                            "alpha[1/K]"    => out[k].SS_vec[i].alpha,
                            "Entropy[J/K]"  => out[k].SS_vec[i].entropy,
                            "Enthalpy[J]"   => out[k].SS_vec[i].enthalpy,
                            "Vp[km/s]"      => out[k].SS_vec[i].Vp,
                            "Vs[km/s]"      => out[k].SS_vec[i].Vs,
                            "Vp_S[km/s]"      => "-",
                            "Vs_S[km/s]"      => "-",                            
                            "BulkMod[GPa]"  => out[k].SS_vec[i].bulkMod,
                            "ShearMod[GPa]" => out[k].SS_vec[i].shearMod )  

            part_2 = Dict(  (out[1].oxides[j]*"[mol%]" => out[k].SS_vec[i].Comp[j]*100.0)
                            for j in eachindex(out[1].oxides))

            part_3 = Dict(  (out[1].oxides[j]*"[wt%]" => out[k].SS_vec[i].Comp_wt[j]*100.0)
                            for j in eachindex(out[1].oxides))

            row    = merge(part_1,part_2,part_3)   

            push!(MAGEMin_db, row, cols=:union)
            
        end

        if npp > 0
            for i=1:npp
                pos = i + nss

                part_1 = Dict(  "point[#]"      => k,
                                "X[0.0-1.0]"    => out[k].X[1],
                                "P[kbar]"       => out[k].P_kbar,
                                "T[°C]"         => out[k].T_C,
                                "phase"         => out[k].ph[pos],
                                "mode[mol%]"    => out[k].ph_frac[pos].*100.0,
                                "mode[wt%]"     => out[k].ph_frac_wt[pos].*100.0,
                                "log10(fO2)"    => "-",
                                "log10(dQFM)"   => "-",
                                "aH2O"          => "-",
                                "aSiO2"         => "-",
                                "aTiO2"         => "-",
                                "aAl2O3"        => "-",
                                "aMgO"          => "-",
                                "aFeO"          => "-",
                                "density[kg/m3]"    => out[k].PP_vec[i].rho,
                                "volume[cm3/mol]"   => out[k].PP_vec[i].V,
                                "heatCapacity[J/kg/K]"=> out[k].PP_vec[i].cp,
                                "alpha[1/K]"    => out[k].PP_vec[i].alpha,
                                "Entropy[J/K]"  => out[k].PP_vec[i].entropy,
                                "Enthalpy[J]"   => out[k].PP_vec[i].enthalpy,
                                "Vp[km/s]"      => out[k].PP_vec[i].Vp,
                                "Vs[km/s]"      => out[k].PP_vec[i].Vs,
                                "Vp_S[km/s]"      => "-",
                                "Vs_S[km/s]"      => "-",           
                                "BulkMod[GPa]"  => out[k].PP_vec[i].bulkMod,
                                "ShearMod[GPa]" => out[k].PP_vec[i].shearMod )  

                part_2 = Dict(  (out[1].oxides[j]*"[mol%]" => out[k].PP_vec[i].Comp[j]*100.0)
                                for j in eachindex(out[1].oxides))

                part_3 = Dict(  (out[1].oxides[j]*"[wt%]" => out[k].PP_vec[i].Comp_wt[j]*100.0)
                                for j in eachindex(out[1].oxides))

                row    = merge(part_1,part_2,part_3)   

                push!(MAGEMin_db, row, cols=:union)

            end
        end

    end

    meta = fileout*"_metadata.txt"
    filename = fileout*".csv"
    CSV.write(filename, MAGEMin_db)

    open(meta, "w") do file
        write(file, metadata)
    end

    return nothing
end


"""
    get_all_stable_phases(out:: Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}})
    return ph_name

    The function receives as an input a single/Vector of MAGEMin_C output structure and returns the list (Vector{String}) of unique stable phases
        - Note that first the sorted solution phase names are provided, followed by the sorted pure phase names
          e.g., ["amp", "bi", "chl", "cpx", "ep", "fl", "fsp", "liq", "opx", "sph"]
"""
function get_all_stable_phases(out:: Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}})
     np     = length(out)
    if typeof(out) == MAGEMin_C.gmin_struct{Float64, Int64}
        out = [out]
    end
    ss_ph_names = Vector{String}()
    pp_ph_names = Vector{String}()
    for k = 1:np
        for l=1:length(out[k].ph)

            if l > out[k].n_SS
                if ~(out[k].ph[l] in pp_ph_names)
                    push!(pp_ph_names,out[k].ph[l])
                end
            else
                if ~(out[k].ph[l] in ss_ph_names)
                    push!(ss_ph_names,out[k].ph[l])
                end
            end
        end
    end

    ph_names = vcat(sort(ss_ph_names),sort(pp_ph_names))


    return ph_names, length(ss_ph_names), length(pp_ph_names)
end


"""
MAGEMin_data2dataframe( out:: Union{Vector{MAGEMin_C.gmin_struct{Float64, Int64}}, MAGEMin_C.gmin_struct{Float64, Int64}})

    Transform MAGEMin output into a dataframe for quick(ish) save

"""
function MAGEMin_data2dataframe_inlined( out:: Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}},dtb,fileout)
    print("\noutput path: $(pwd())\n")
  
    # here we fill the dataframe with the all minimized point entries
    if typeof(out) == MAGEMin_C.gmin_struct{Float64, Int64}
        out = [out]
    end
    np        = length(out)
    db_in     = retrieve_solution_phase_information(dtb)
    datetoday = string(Dates.today())
    rightnow  = string(Dates.Time(Dates.now()))

    metadata  = "# MAGEMin " * " $(out[1].MAGEMin_ver);" * datetoday * ", " * rightnow * "; using database " * db_in.db_info * "\n"

    ph_names, n_ss, n_pp  = get_all_stable_phases(out)

    # Here we create the dataframe's header:
    sys_db = DataFrame(             Symbol("point[#]")      => Int64[],
                                    Symbol("X[0.0-1.0]")    => Float64[],
                                    Symbol("P[kbar]")       => Float64[],
                                    Symbol("T[°C]")         => Float64[],
                                    Symbol("sys_log10(fO2)")    => Float64[],
                                    Symbol("sys_log10(dQFM)")   => Float64[],
                                    Symbol("sys_aH2O")          => Float64[],
                                    Symbol("sys_aSiO2")         => Float64[],
                                    Symbol("sys_aTiO2")         => Float64[],
                                    Symbol("sys_aAl2O3")        => Float64[],
                                    Symbol("sys_aMgO")          => Float64[],
                                    Symbol("sys_aFeO")          => Float64[],
                                    Symbol("sys_density[kg/m3]")    => Float64[],
                                    Symbol("sys_volume[cm3/mol]")   => Float64[],
                                    Symbol("sys_heatCapacity[J/kg/K]")=> Float64[],
                                    Symbol("sys_alpha[1/K]")    => Float64[],
                                    Symbol("sys_Entropy[J/K]")  => Float64[],
                                    Symbol("sys_Enthalpy[J]")   => Float64[],
                                    Symbol("sys_Vp[km/s]")      => Float64[],
                                    Symbol("sys_Vs[km/s]")      => Float64[],
                                    Symbol("Vp_S[km/s]")      => Float64[],
                                    Symbol("Vs_S[km/s]")      => Float64[],                                    
                                    Symbol("sys_BulkMod[GPa]")  => Float64[],
                                    Symbol("sys_ShearMod[GPa]") => Float64[])

    for i in out[1].oxides
        col = "sys_"*i*"[mol%]"
        sys_db[!, col] = Float64[] 
    end

    for i in out[1].oxides
        col = "sys_"*i*"[wt%]"
        sys_db[!, col] = Float64[] 
    end

    @showprogress "Saving data to csv..." for k=1:np
        part_1 = Dict(  "point[#]"      => k,
                        "X[0.0-1.0]"    => out[k].X[1],
                        "P[kbar]"       => out[k].P_kbar,
                        "T[°C]"         => out[k].T_C,
                        "sys_log10(fO2)"    => out[k].fO2[1],
                        "sys_log10(dQFM)"   => out[k].dQFM[1],
                        "sys_aH2O"          => out[k].aH2O,
                        "sys_aSiO2"         => out[k].aSiO2,
                        "sys_aTiO2"         => out[k].aTiO2,
                        "sys_aAl2O3"        => out[k].aAl2O3,
                        "sys_aMgO"          => out[k].aMgO,
                        "sys_aFeO"          => out[k].aFeO,
                        "sys_density[kg/m3]"    => out[k].rho,
                        "sys_volume[cm3/mol]"   => out[k].V,
                        "sys_heatCapacity[J/kg/K]"=> out[k].s_cp,
                        "sys_alpha[1/K]"    => out[k].alpha,
                        "sys_Entropy[J/K]"  => out[k].entropy,
                        "sys_Enthalpy[J]"   => out[k].enthalpy,
                        "sys_Vp[km/s]"      => out[k].Vp,
                        "sys_Vs[km/s]"      => out[k].Vs,
                        "Vp_S[km/s]"    => out[k].Vp_S,
                        "Vs_S[km/s]"    => out[k].Vs_S,
                        "sys_BulkMod[GPa]"  => out[k].bulkMod,
                        "sys_ShearMod[GPa]" => out[k].shearMod )          

        part_2 = Dict(  ("sys_"*out[1].oxides[j]*"[mol%]" => out[k].bulk[j]*100.0)
                        for j in eachindex(out[1].oxides))

        part_3 = Dict(  ("sys_"*out[1].oxides[j]*"[wt%]" => out[k].bulk_wt[j]*100.0)
                        for j in eachindex(out[1].oxides))
   
        row    = merge(part_1,part_2,part_3)   

        push!(sys_db, row, cols=:union)
    end


    ss_db   = DataFrame()
    for i=1:n_ss
        ph = ph_names[i]
        ss_db_ =     DataFrame(     Symbol("$(ph)_mode[mol%]")    => Float64[],
                                    Symbol("$(ph)_mode[wt%]")     => Float64[],
                                    Symbol("$(ph)_mode[vol%]")     => Float64[],
                                    Symbol("$(ph)_density[kg/m3]")    => Float64[],
                                    Symbol("$(ph)_volume[cm3/mol]")   => Float64[],
                                    Symbol("$(ph)_heatCapacity[J/kg/K]")=> Float64[],
                                    Symbol("$(ph)_alpha[1/K]")    => Float64[],
                                    Symbol("$(ph)_Entropy[J/K]")  => Float64[],
                                    Symbol("$(ph)_Enthalpy[J]")   => Float64[],
                                    Symbol("$(ph)_Vp[km/s]")      => Float64[],
                                    Symbol("$(ph)_Vs[km/s]")      => Float64[],                                  
                                    Symbol("$(ph)_BulkMod[GPa]")  => Float64[],
                                    Symbol("$(ph)_ShearMod[GPa]") => Float64[])

        for i in out[1].oxides
            col = ph*"_"*i*"[mol%]"
            ss_db_[!, col] = Float64[] 
        end

        for i in out[1].oxides
            col = ph*"_"*i*"[wt%]"
            ss_db_[!, col] = Float64[] 
        end
        ss_db = hcat(ss_db,ss_db_)
    end
  

    @showprogress "Saving data to csv..." for k=1:np
        row     = Dict{String, Float64}()
        for j=1:n_ss
            
            ph      = ph_names[j]
            if ph in out[k].ph
                i = findfirst(out[k].ph .== ph)
                ss_part_1 = Dict(   "$(ph)_mode[mol%]"    => out[k].ph_frac[i].*100.0,
                        "$(ph)_mode[wt%]"     => out[k].ph_frac_wt[i].*100.0,
                        "$(ph)_mode[vol%]"     => out[k].ph_frac_vol[i].*100.0,
                        "$(ph)_density[kg/m3]"    => out[k].SS_vec[i].rho,
                        "$(ph)_volume[cm3/mol]"   => out[k].SS_vec[i].V,
                        "$(ph)_heatCapacity[J/kg/K]"=> out[k].SS_vec[i].cp,
                        "$(ph)_alpha[1/K]"    => out[k].SS_vec[i].alpha,
                        "$(ph)_Entropy[J/K]"  => out[k].SS_vec[i].entropy,
                        "$(ph)_Enthalpy[J]"   => out[k].SS_vec[i].enthalpy,
                        "$(ph)_Vp[km/s]"      => out[k].SS_vec[i].Vp,
                        "$(ph)_Vs[km/s]"      => out[k].SS_vec[i].Vs,                         
                        "$(ph)_BulkMod[GPa]"  => out[k].SS_vec[i].bulkMod,
                        "$(ph)_ShearMod[GPa]" => out[k].SS_vec[i].shearMod )  

                ss_part_2 = Dict(  (ph*"_"*out[1].oxides[j]*"[mol%]" => out[k].SS_vec[i].Comp[j]*100.0)
                                for j in eachindex(out[1].oxides))

                ss_part_3 = Dict(  (ph*"_"*out[1].oxides[j]*"[wt%]" => out[k].SS_vec[i].Comp_wt[j]*100.0)
                                for j in eachindex(out[1].oxides))
            else
                ss_part_1 = Dict(   "$(ph)_mode[mol%]"    => NaN,
                                    "$(ph)_mode[wt%]"     => NaN,
                                    "$(ph)_mode[vol%]"     => NaN,
                                    "$(ph)_density[kg/m3]"    => NaN,
                                    "$(ph)_volume[cm3/mol]"   => NaN,
                                    "$(ph)_heatCapacity[J/kg/K]"=> NaN,
                                    "$(ph)_alpha[1/K]"    => NaN,
                                    "$(ph)_Entropy[J/K]"  => NaN,
                                    "$(ph)_Enthalpy[J]"   => NaN,
                                    "$(ph)_Vp[km/s]"      => NaN,
                                    "$(ph)_Vs[km/s]"      => NaN,                         
                                    "$(ph)_BulkMod[GPa]"  => NaN,
                                    "$(ph)_ShearMod[GPa]" => NaN )  

                ss_part_2 = Dict(  (ph*"_"*out[1].oxides[j]*"[mol%]" => NaN)
                                for j in eachindex(out[1].oxides))

                ss_part_3 = Dict(  (ph*"_"*out[1].oxides[j]*"[wt%]" => NaN)
                                for j in eachindex(out[1].oxides))

            end
            row_    = merge(ss_part_1,ss_part_2,ss_part_3)  
            row = merge(row,row_) 
        end
        push!(ss_db, row, cols=:union)
    end

    if n_pp > 0
        pp_db   = DataFrame()
        for i=1:n_pp
            ph = ph_names[i+n_ss]
            pp_db_ =     DataFrame(     Symbol("$(ph)_mode[mol%]")    => Float64[],
                                        Symbol("$(ph)_mode[wt%]")     => Float64[],
                                        Symbol("$(ph)_mode[vol%]")     => Float64[],
                                        Symbol("$(ph)_density[kg/m3]")    => Float64[],
                                        Symbol("$(ph)_volume[cm3/mol]")   => Float64[],
                                        Symbol("$(ph)_heatCapacity[J/kg/K]")=> Float64[],
                                        Symbol("$(ph)_alpha[1/K]")    => Float64[],
                                        Symbol("$(ph)_Entropy[J/K]")  => Float64[],
                                        Symbol("$(ph)_Enthalpy[J]")   => Float64[],
                                        Symbol("$(ph)_Vp[km/s]")      => Float64[],
                                        Symbol("$(ph)_Vs[km/s]")      => Float64[],                                  
                                        Symbol("$(ph)_BulkMod[GPa]")  => Float64[],
                                        Symbol("$(ph)_ShearMod[GPa]") => Float64[])

            for i in out[1].oxides
                col = ph*"_"*i*"[mol%]"
                pp_db_[!, col] = Float64[] 
            end

            for i in out[1].oxides
                col = ph*"_"*i*"[wt%]"
                pp_db_[!, col] = Float64[] 
            end
            pp_db = hcat(pp_db,pp_db_)
        end
    

        @showprogress "Saving data to csv..." for k=1:np
            row     = Dict{String, Float64}()
            for j=1:n_pp
                ph      = ph_names[j+n_ss]
                if ph in out[k].ph
                    i = findfirst(out[k].ph .== ph)
                    s = out[k].n_SS
                    pp_part_1 = Dict(   "$(ph)_mode[mol%]"    => out[k].ph_frac[i].*100.0,
                                        "$(ph)_mode[wt%]"     => out[k].ph_frac_wt[i].*100.0,
                                        "$(ph)_mode[vol%]"     => out[k].ph_frac_vol[i].*100.0,
                                        "$(ph)_density[kg/m3]"    => out[k].PP_vec[i-s].rho,
                                        "$(ph)_volume[cm3/mol]"   => out[k].PP_vec[i-s].V,
                                        "$(ph)_heatCapacity[J/kg/K]"=> out[k].PP_vec[i-s].cp,
                                        "$(ph)_alpha[1/K]"    => out[k].PP_vec[i-s].alpha,
                                        "$(ph)_Entropy[J/K]"  => out[k].PP_vec[i-s].entropy,
                                        "$(ph)_Enthalpy[J]"   => out[k].PP_vec[i-s].enthalpy,
                                        "$(ph)_Vp[km/s]"      => out[k].PP_vec[i-s].Vp,
                                        "$(ph)_Vs[km/s]"      => out[k].PP_vec[i-s].Vs,                         
                                        "$(ph)_BulkMod[GPa]"  => out[k].PP_vec[i-s].bulkMod,
                                        "$(ph)_ShearMod[GPa]" => out[k].PP_vec[i-s].shearMod )  

                    pp_part_2 = Dict(  (ph*"_"*out[1].oxides[j]*"[mol%]" => out[k].PP_vec[i-s].Comp[j]*100.0)
                                    for j in eachindex(out[1].oxides))

                    pp_part_3 = Dict(  (ph*"_"*out[1].oxides[j]*"[wt%]" => out[k].PP_vec[i-s].Comp_wt[j]*100.0)
                                    for j in eachindex(out[1].oxides))
                else
                    pp_part_1 = Dict(   "$(ph)_mode[mol%]"    => NaN,
                                        "$(ph)_mode[wt%]"     => NaN,
                                        "$(ph)_mode[vol%]"     => NaN,
                                        "$(ph)_density[kg/m3]"    => NaN,
                                        "$(ph)_volume[cm3/mol]"   => NaN,
                                        "$(ph)_heatCapacity[J/kg/K]"=> NaN,
                                        "$(ph)_alpha[1/K]"    => NaN,
                                        "$(ph)_Entropy[J/K]"  => NaN,
                                        "$(ph)_Enthalpy[J]"   => NaN,
                                        "$(ph)_Vp[km/s]"      => NaN,
                                        "$(ph)_Vs[km/s]"      => NaN,                         
                                        "$(ph)_BulkMod[GPa]"  => NaN,
                                        "$(ph)_ShearMod[GPa]" => NaN )  

                    pp_part_2 = Dict(  (ph*"_"*out[1].oxides[j]*"[mol%]" => NaN)
                                    for j in eachindex(out[1].oxides))

                    pp_part_3 = Dict(  (ph*"_"*out[1].oxides[j]*"[wt%]" => NaN)
                                    for j in eachindex(out[1].oxides))

                end
                row_    = merge(pp_part_1,pp_part_2,pp_part_3)  
                row = merge(row,row_) 
            end
            push!(pp_db, row, cols=:union)
        end
    end

    if n_pp > 0
        MAGEMin_db = hcat(sys_db,ss_db,pp_db)
    else
        MAGEMin_db = hcat(sys_db,ss_db)
    end

    meta = fileout*"_metadata.txt"
    filename = fileout*"_inlined.csv"
    CSV.write(filename, MAGEMin_db)

    open(meta, "w") do file
        write(file, metadata)
    end

    return nothing
end