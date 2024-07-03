
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
        if ~isnothing(out_te[k].liq_wt_norm)
            part_1 = Dict(  "point[#]"      => k,
                            "X[0.0-1.0]"    => out[k].X[1],
                            "P[kbar]"       => out[k].P_kbar,
                            "T[°C]"         => out[k].T_C,
                            "phase"         => "liq",
                            "mode[wt%]"     => out_te[k].liq_wt_norm,
                            "Zr_sat[μg/g]"  => out_te[k].Sat_zr_liq)

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
                                "phase"         => out_te[k].ph_TE,
                                "mode[wt%]"     => out_te[k].ph_wt_norm,
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
                                    Symbol("heatCapacity[kJ/K]")=> Float64[],
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
                        "heatCapacity[kJ/K]"=> out[k].cp,
                        "alpha[1/K]"    => out[k].alpha,
                        "Entropy[J/K]"  => out[k].entropy,
                        "Enthalpy[J]"   => out[k].enthalpy,
                        "Vp[km/s]"      => out[k].Vp,
                        "Vs[km/s]"      => out[k].Vs,
                        "Vp_S[km/s]"      => out[k].Vp_S,
                        "Vs_S[km/s]"      => out[k].Vs_S,
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
                            "heatCapacity[kJ/K]"=> out[k].SS_vec[i].cp,
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
                                "heatCapacity[kJ/K]"=> out[k].PP_vec[i].cp,
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