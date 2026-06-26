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
    MAGEMin_dataTE2dataframe(out, out_te, dtb, fileout)

    Export MAGEMin minimization results together with trace element partitioning
    output to a CSV file, writing one row per phase per point.

    The CSV contains columns for P, T, X, phase name, mode [wt%], Zr saturation
    [μg/g], corrected bulk oxide composition [wt%], and per-element concentrations
    [μg/g]. A companion metadata text file recording the MAGEMin version, database,
    date, and time is written alongside the CSV.

    Rows are written for the bulk system (`"system"`), the melt (`"liq"`), the
    bulk solid aggregate (`"sol"`), and each individual mineral phase.

    Parameters
    ----------
    out : Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}}
        MAGEMin minimization output for one or more points.
    out_te : Union{Vector{out_tepm}, out_tepm}
        Trace element partitioning output corresponding to each point in `out`.
    dtb : String
        Database identifier used to retrieve database metadata (e.g., `"ig"`, `"mp"`).
    fileout : String
        Base path for output files. Two files are created: `fileout.csv` and
        `fileout_metadata.txt`.
    use_Warr2021 : Bool, optional
        If `true`, mineral phase names are converted to IMA-CNMNC approved symbols
        following Warr (2021). Names with no official symbol are returned with a
        trailing `"*"`. Default is `false`.
    use_GPA : Bool, optional
        If `true`, pressure is reported in GPa (column `"P[GPa]"`) instead of kbar
        (column `"P[kbar]"`). Default is `false`.

    Returns
    -------
    nothing
"""
function MAGEMin_dataTE2dataframe(  out     :: Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}},
                                    out_te  :: Union{Vector{out_tepm}, out_tepm},
                                    dtb,
                                    fileout;
                                    use_Warr2021 :: Bool = false,
                                    use_GPA :: Bool = false)

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

    P_colname  = use_GPA ? "P[GPa]" : "P[kbar]"

    # Here we create the dataframe's header:
    MAGEMin_db = DataFrame(         Symbol("point[#]")          => Int64[],
                                    Symbol("X[0.0-1.0]")        => Float64[],
                                    Symbol(P_colname)            => Float64[],
                                    Symbol("T[°C]")             => Float64[],
                                    Symbol("phase")             => String[],
                                    Symbol("mode[wt%]")         => Float64[],
                                    Symbol("Zr_sat[μg/g]")      => Float64[],
                                    Symbol("S_sat[μg/g]")       => Float64[],
                                    Symbol("P2O5_sat[wt%]")     => Float64[],
                                    Symbol("CO2_sat[wt%]")      => Float64[],
                                    Symbol("zrc_wt[wt%]")       => Float64[],
                                    Symbol("sulf_wt[wt%]")      => Float64[],
                                    Symbol("fapt_wt[wt%]")      => Float64[],
                                    Symbol("fl_CO2_wt[wt%]")    => Float64[],
                                    Symbol("bulk_D[-]")         => Float64[],
    )
    for i in out[1].oxides
        col = i*"_cor[wt%]"
        MAGEMin_db[!, col] = Float64[]
    end

    for i in out[1].oxides
        col = i*"_cor[mol%]"
        MAGEMin_db[!, col] = Float64[]
    end

    for i in out_te[1].elements
        col = i*"_[μg/g]"
        MAGEMin_db[!, col] = Float64[]
    end

    print("\noutput path: $(pwd())\n")
    @showprogress "Saving data to csv..." for k=1:np
        P_val = use_GPA ? out[k].P_kbar * 0.1 : out[k].P_kbar

        # system
        part_1 = Dict(  "point[#]"      => k,
                        "X[0.0-1.0]"    => out[k].X[1],
                        P_colname       => P_val,
                        "T[°C]"         => out[k].T_C,
                        "phase"         => "system",
                        "mode[wt%]"     => 100.0,
                        "Zr_sat[μg/g]"  => "-",
                        "S_sat[μg/g]"   => "-",
                        "P2O5_sat[wt%]" => "-",
                        "CO2_sat[wt%]"  => "-",
                        "zrc_wt[wt%]"   => "-",
                        "sulf_wt[wt%]"  => "-",
                        "fapt_wt[wt%]"  => "-",
                        "fl_CO2_wt[wt%]"=> "-",
                        "bulk_D[-]"     => out_te[k].bulk_D)

        if ~isnothing(out_te[k].bulk_cor_wt) && !any(isnan, out_te[k].bulk_cor_wt)
            part_2 = Dict(  (out[1].oxides[j]*"_cor[wt%]" => out_te[k].bulk_cor_wt[j]*100.0)
                            for j in eachindex(out[1].oxides))
        else
            part_2 = Dict(  (out[1].oxides[j]*"_cor[wt%]" => "-")
                            for j in eachindex(out[1].oxides))
        end

        if ~isnothing(out_te[k].bulk_cor_mol) && !any(isnan, out_te[k].bulk_cor_mol)
            part_2b = Dict( (out[1].oxides[j]*"_cor[mol%]" => out_te[k].bulk_cor_mol[j]*100.0)
                            for j in eachindex(out[1].oxides))
        else
            part_2b = Dict( (out[1].oxides[j]*"_cor[mol%]" => "-")
                            for j in eachindex(out[1].oxides))
        end

        part_3 = Dict(  ( out_te[1].elements[j]*"_[μg/g]" => out_te[k].C0[j])
                        for j in eachindex(out_te[1].elements))

        row    = merge(part_1,part_2,part_2b,part_3)

        push!(MAGEMin_db, row, cols=:union)
        
        # liquid
        if isnothing(out_te[k].Sat_Zr_liq)
            Sat_Zr_liq = "-"
        else
            Sat_Zr_liq = out_te[k].Sat_Zr_liq
        end

        if isnothing(out_te[k].Sat_S_liq)
            Sat_S_liq = "-"
        else
            Sat_S_liq = out_te[k].Sat_S_liq
        end

        if isnothing(out_te[k].Sat_P2O5_liq)
            Sat_P2O5_liq = "-"
        else
            Sat_P2O5_liq = out_te[k].Sat_P2O5_liq
        end

        if isnothing(out_te[k].Sat_CO2_liq)
            Sat_CO2_liq = "-"
        else
            Sat_CO2_liq = out_te[k].Sat_CO2_liq
        end

        if isnothing(out_te[k].zrc_wt)
            zrc_wt = "-"
        else
            zrc_wt = out_te[k].zrc_wt
        end

        if isnothing(out_te[k].sulf_wt)
            sulf_wt = "-"
        else
            sulf_wt = out_te[k].sulf_wt
        end

        if isnothing(out_te[k].fapt_wt)
            fapt_wt = "-"
        else
            fapt_wt = out_te[k].fapt_wt
        end

        if isnothing(out_te[k].fl_CO2_wt)
            fl_CO2_wt = "-"
        else
            fl_CO2_wt = out_te[k].fl_CO2_wt
        end

        if ~isnothing(out_te[k].liq_wt_norm) && !any(isnan, out_te[k].liq_wt_norm)
            part_1 = Dict(  "point[#]"      => k,
                            "X[0.0-1.0]"    => out[k].X[1],
                            P_colname       => P_val,
                            "T[°C]"         => out[k].T_C,
                            "phase"         => "liq",
                            "mode[wt%]"     => out_te[k].liq_wt_norm .* 100.0,
                            "Zr_sat[μg/g]"  => Sat_Zr_liq,
                            "S_sat[μg/g]"   => Sat_S_liq,
                            "P2O5_sat[wt%]" => Sat_P2O5_liq,
                            "CO2_sat[wt%]"  => Sat_CO2_liq,
                            "zrc_wt[wt%]"   => zrc_wt,
                            "sulf_wt[wt%]"  => sulf_wt,
                            "fapt_wt[wt%]"  => fapt_wt,
                            "fl_CO2_wt[wt%]"=> fl_CO2_wt,
                            "bulk_D[-]"     => "-")

            part_2 = Dict(  (out[1].oxides[j]*"_cor[wt%]" => "-")
                            for j in eachindex(out[1].oxides))

            part_2b = Dict( (out[1].oxides[j]*"_cor[mol%]" => "-")
                            for j in eachindex(out[1].oxides))

            part_3 = Dict(  ( out_te[1].elements[j]*"_[μg/g]" => out_te[k].Cliq[j])
                            for j in eachindex(out_te[1].elements))

            row    = merge(part_1,part_2,part_2b,part_3)

            push!(MAGEMin_db, row, cols=:union)
        end

        # solid
        if ~isnothing(out_te[k].Csol) && !any(isnan, out_te[k].Csol)
            if ~isnothing(out_te[k].liq_wt_norm) && !isnan(out_te[k].liq_wt_norm)
                sol_mode = (1.0-out_te[k].liq_wt_norm)*100.0
            else
                sol_mode = 100.0
            end
            part_1 = Dict(  "point[#]"      => k,
                            "X[0.0-1.0]"    => out[k].X[1],
                            P_colname       => P_val,
                            "T[°C]"         => out[k].T_C,
                            "phase"         => "sol",
                            "mode[wt%]"     => sol_mode,
                            "Zr_sat[μg/g]"  => "-",
                            "S_sat[μg/g]"   => "-",
                            "P2O5_sat[wt%]" => "-",
                            "CO2_sat[wt%]"  => "-",
                            "zrc_wt[wt%]"   => "-",
                            "sulf_wt[wt%]"  => "-",
                            "fapt_wt[wt%]"  => "-",
                            "fl_CO2_wt[wt%]"=> "-",
                            "bulk_D[-]"     => "-")

            part_2 = Dict(  (out[1].oxides[j]*"_cor[wt%]" => "-")
                            for j in eachindex(out[1].oxides))

            part_2b = Dict( (out[1].oxides[j]*"_cor[mol%]" => "-")
                            for j in eachindex(out[1].oxides))

            part_3 = Dict(  ( out_te[1].elements[j]*"_[μg/g]" => out_te[k].Csol[j])
                            for j in eachindex(out_te[1].elements))

            row    = merge(part_1,part_2,part_2b,part_3)

            push!(MAGEMin_db, row, cols=:union)
        end

        if ~isnothing(out_te[k].ph_TE)
            nph  = length(out_te[k].ph_TE)
            for i=1:nph
                ph_label = use_Warr2021 ? get_Warr_name(out_te[k].ph_TE[i]) : out_te[k].ph_TE[i]
                part_1 = Dict(  "point[#]"      => k,
                                "X[0.0-1.0]"    => out[k].X[1],
                                P_colname       => P_val,
                                "T[°C]"         => out[k].T_C,
                                "phase"         => ph_label,
                                "mode[wt%]"     => out_te[k].ph_wt_norm[i].*100.0,
                                "Zr_sat[μg/g]"  => "-",
                                "S_sat[μg/g]"   => "-",
                                "P2O5_sat[wt%]" => "-",
                                "CO2_sat[wt%]"  => "-",
                                "zrc_wt[wt%]"   => "-",
                                "sulf_wt[wt%]"  => "-",
                                "fapt_wt[wt%]"  => "-",
                                "fl_CO2_wt[wt%]"=> "-",
                                "bulk_D[-]"     => "-")

                part_2 = Dict(  (out[1].oxides[j]*"_cor[wt%]" => "-")
                                for j in eachindex(out[1].oxides))

                part_2b = Dict( (out[1].oxides[j]*"_cor[mol%]" => "-")
                                for j in eachindex(out[1].oxides))

                part_3 = Dict(  ( out_te[1].elements[j]*"_[μg/g]" => out_te[k].Cmin[i,j])
                                for j in eachindex(out_te[1].elements))

                row    = merge(part_1,part_2,part_2b,part_3)

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
    MAGEMin_data2dataframe(out, dtb, fileout)

    Export MAGEMin minimization results to a CSV file in long (stacked) format,
    writing one row per phase per point.

    Each row contains point index, X fraction, P [kbar], T [°C], phase name,
    modal abundances (mol%, wt%, vol%), system thermodynamic properties
    (ρ, Vp, Vs, Cp, α, S, H, K, G), activity/fugacity variables, and oxide
    compositions (mol% and wt%) plus apfu element compositions for each phase.
    System-level rows (`"system"`) additionally carry fO₂, ΔQFM, melt
    viscosity, and the seismically-corrected solid-aggregate velocities
    (`Vp_cor`, `Vs_cor`, NaN if seismic corrections were not requested).
    A companion metadata file is written at `fileout_metadata.txt`.

    Parameters
    ----------
    out : Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}}
        MAGEMin minimization output for one or more points.
    dtb : String
        Database identifier used to retrieve database metadata (e.g., `"ig"`, `"mp"`).
    fileout : String
        Base path for output files. Two files are created: `fileout.csv` and
        `fileout_metadata.txt`.
    use_Warr2021 : Bool, optional
        If `true`, mineral phase names are converted to IMA-CNMNC approved symbols
        following Warr (2021). Names with no official symbol are returned with a
        trailing `"*"`. Default is `false`.
    use_GPA : Bool, optional
        If `true`, pressure is reported in GPa (column `"P[GPa]"`) instead of kbar
        (column `"P[kbar]"`). Default is `false`.

    Returns
    -------
    nothing
"""
function MAGEMin_data2dataframe( out:: Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}},dtb,fileout; use_Warr2021 :: Bool = false, use_GPA :: Bool = false)

    # here we fill the dataframe with the all minimized point entries
    if typeof(out) == MAGEMin_C.gmin_struct{Float64, Int64}
        out = [out]
    end
    np        = length(out)

    db_in     = retrieve_solution_phase_information(dtb)
    datetoday = string(Dates.today())
    rightnow  = string(Dates.Time(Dates.now()))

    metadata   = "# MAGEMin " * " $(out[1].MAGEMin_ver);" * datetoday * ", " * rightnow * "; using database " * db_in.db_info * "\n"

    P_colname  = use_GPA ? "P[GPa]" : "P[kbar]"

    # Here we create the dataframe's header:
    MAGEMin_db = DataFrame(         Symbol("point[#]")      => Int64[],
                                    Symbol("X[0.0-1.0]")    => Float64[],
                                    Symbol(P_colname)       => Float64[],
                                    Symbol("T[°C]")         => Float64[],
                                    Symbol("phase")         => String[],
                                    Symbol("mode[mol%]")    => Float64[],
                                    Symbol("mode[wt%]")     => Float64[],
                                    Symbol("mode[vol%]")    => Float64[],
                                    Symbol("log10(fO2)")    => Float64[],
                                    Symbol("log10(dQFM)")   => Float64[],
                                    Symbol("aH2O")          => Float64[],
                                    Symbol("aSiO2")         => Float64[],
                                    Symbol("aTiO2")         => Float64[],
                                    Symbol("aAl2O3")        => Float64[],
                                    Symbol("aMgO")          => Float64[],
                                    Symbol("aFeO")          => Float64[],
                                    Symbol("meltViscosity[Pa.s]")    => Float64[],
                                    Symbol("density[kg/m3]")    => Float64[],
                                    Symbol("volume[cm3/mol]")   => Float64[],
                                    Symbol("heatCapacity[J/kg/K]")=> Float64[],
                                    Symbol("alpha[1/K]")    => Float64[],
                                    Symbol("Entropy[kJ/K]")  => Float64[],
                                    Symbol("Enthalpy[kJ/mol]")   => Float64[],
                                    Symbol("Vp[km/s]")      => Float64[],
                                    Symbol("Vs[km/s]")      => Float64[],
                                    Symbol("Vp_S[km/s]")      => Float64[],
                                    Symbol("Vs_S[km/s]")      => Float64[],
                                    Symbol("Vp_cor[km/s]")      => Float64[],
                                    Symbol("Vs_cor[km/s]")      => Float64[],
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
        P_val = use_GPA ? out[k].P_kbar * 0.1 : out[k].P_kbar

        part_1 = Dict(  "point[#]"      => k,
                        "X[0.0-1.0]"    => out[k].X[1],
                        P_colname       => P_val,
                        "T[°C]"         => out[k].T_C,
                        "phase"         => "system",
                        "mode[mol%]"    => 100.0,
                        "mode[wt%]"     => 100.0,
                        "mode[vol%]"    => 100.0,
                        "log10(fO2)"    => out[k].fO2[1],
                        "log10(dQFM)"   => out[k].dQFM[1],
                        "aH2O"          => out[k].aH2O,
                        "aSiO2"         => out[k].aSiO2,
                        "aTiO2"         => out[k].aTiO2,
                        "aAl2O3"        => out[k].aAl2O3,
                        "aMgO"          => out[k].aMgO,
                        "aFeO"          => out[k].aFeO,
                        "meltViscosity[Pa.s]"    => out[k].eta_M,
                        "density[kg/m3]"    => out[k].rho,
                        "volume[cm3/mol]"   => out[k].V,
                        "heatCapacity[J/kg/K]"=> out[k].s_cp,
                        "alpha[1/K]"    => out[k].alpha,
                        "Entropy[kJ/K]"  => out[k].entropy,
                        "Enthalpy[kJ/mol]"   => out[k].enthalpy,
                        "Vp[km/s]"      => out[k].Vp,
                        "Vs[km/s]"      => out[k].Vs,
                        "Vp_S[km/s]"    => out[k].Vp_S,
                        "Vs_S[km/s]"    => out[k].Vs_S,
                        "Vp_cor[km/s]"  => out[k].Vp_cor,
                        "Vs_cor[km/s]"  => out[k].Vs_cor,
                        "BulkMod[GPa]"  => out[k].bulkMod,
                        "ShearMod[GPa]" => out[k].shearMod )

        part_2 = Dict(  (out[1].oxides[j]*"[mol%]" => out[k].bulk[j]*100.0)
                        for j in eachindex(out[1].oxides))

        part_3 = Dict(  (out[1].oxides[j]*"[wt%]" => out[k].bulk_wt[j]*100.0)
                        for j in eachindex(out[1].oxides))
   
        row    = merge(part_1,part_2,part_3)   

        push!(MAGEMin_db, row, cols=:union)

        for i=1:nss
            ph_label = use_Warr2021 ? get_Warr_name(out[k].ph[i]) : out[k].ph[i]
            part_1 = Dict(  "point[#]"      => k,
                            "X[0.0-1.0]"    => out[k].X[1],
                            P_colname       => P_val,
                            "T[°C]"         => out[k].T_C,
                            "phase"         => ph_label,
                            "mode[mol%]"    => out[k].ph_frac[i].*100.0,
                            "mode[wt%]"     => out[k].ph_frac_wt[i].*100.0,
                            "mode[vol%]"    => out[k].ph_frac_vol[i].*100.0,
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
                            "Entropy[kJ/K]"  => out[k].SS_vec[i].entropy,
                            "Enthalpy[kJ/mol]"   => out[k].SS_vec[i].enthalpy,
                            "Vp[km/s]"      => out[k].SS_vec[i].Vp,
                            "Vs[km/s]"      => out[k].SS_vec[i].Vs,
                            "Vp_S[km/s]"      => "-",
                            "Vs_S[km/s]"      => "-",
                            "Vp_cor[km/s]"    => "-",
                            "Vs_cor[km/s]"    => "-",
                            "BulkMod[GPa]"  => out[k].SS_vec[i].bulkMod,
                            "ShearMod[GPa]" => out[k].SS_vec[i].shearMod )  

            part_2 = Dict(  (out[1].oxides[j]*"[mol%]" => out[k].SS_vec[i].Comp[j]*100.0)
                            for j in eachindex(out[1].oxides))

            part_3 = Dict(  (out[1].oxides[j]*"[wt%]" => out[k].SS_vec[i].Comp_wt[j]*100.0)
                            for j in eachindex(out[1].oxides))

            part_4 = Dict(  (out[1].elements[j]*"[apfu]" => out[k].SS_vec[i].Comp_apfu[j])
                            for j in eachindex(out[1].elements))

            row    = merge(part_1,part_2,part_3,part_4)   

            push!(MAGEMin_db, row, cols=:union)
            
        end

        if npp > 0
            for i=1:npp
                pos = i + nss
                ph_label = use_Warr2021 ? get_Warr_name(out[k].ph[pos]) : out[k].ph[pos]
                part_1 = Dict(  "point[#]"      => k,
                                "X[0.0-1.0]"    => out[k].X[1],
                                P_colname       => P_val,
                                "T[°C]"         => out[k].T_C,
                                "phase"         => ph_label,
                                "mode[mol%]"    => out[k].ph_frac[pos].*100.0,
                                "mode[wt%]"     => out[k].ph_frac_wt[pos].*100.0,
                                "mode[vol%]"    => out[k].ph_frac_vol[pos].*100.0,
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
                                "Entropy[kJ/K]"  => out[k].PP_vec[i].entropy,
                                "Enthalpy[kJ/mol]"   => out[k].PP_vec[i].enthalpy,
                                "Vp[km/s]"      => out[k].PP_vec[i].Vp,
                                "Vs[km/s]"      => out[k].PP_vec[i].Vs,
                                "Vp_S[km/s]"      => "-",
                                "Vs_S[km/s]"      => "-",
                                "Vp_cor[km/s]"    => "-",
                                "Vs_cor[km/s]"    => "-",
                                "BulkMod[GPa]"  => out[k].PP_vec[i].bulkMod,
                                "ShearMod[GPa]" => out[k].PP_vec[i].shearMod )  

                part_2 = Dict(  (out[1].oxides[j]*"[mol%]" => out[k].PP_vec[i].Comp[j]*100.0)
                                for j in eachindex(out[1].oxides))

                part_3 = Dict(  (out[1].oxides[j]*"[wt%]" => out[k].PP_vec[i].Comp_wt[j]*100.0)
                                for j in eachindex(out[1].oxides))

                part_4 = Dict(  (out[1].elements[j]*"[apfu]" => out[k].PP_vec[i].Comp_apfu[j])
                                for j in eachindex(out[1].elements))

                row    = merge(part_1,part_2,part_3,part_4)   

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
    get_all_stable_phases(out)

    Collect the unique set of stable phase names across all minimization points.

    Solution phases are sorted alphabetically and listed first; pure phases follow
    in their own sorted block. For example:
    `["amp", "bi", "chl", "cpx", "ep", "fl", "fsp", "liq", "opx", "sph", "q", "ru"]`

    Parameters
    ----------
    out : Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}}
        MAGEMin minimization output for one or more points.

    Returns
    -------
    ph_names : Vector{String}
        Sorted unique phase names (solution phases first, then pure phases).
    n_ss : Int64
        Number of unique solution phases.
    n_pp : Int64
        Number of unique pure phases.
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
    MAGEMin_data2dataframe_inlined(out, dtb, fileout)

    Export MAGEMin minimization results to a CSV file in wide (inlined) format,
    with one row per point and all phase data spread across columns.

    The output DataFrame is built in three blocks that are horizontally
    concatenated: system-level properties (`sys_*` prefix), solution-phase
    columns (`<ph>_*` prefix, `NaN`-filled when a phase is absent at a given
    point), and pure-phase columns (same convention). System-level properties
    additionally include the solid-aggregate seismic velocities (`Vp_S`,
    `Vs_S`) and their seismically-corrected counterparts (`Vp_cor`, `Vs_cor`,
    NaN if seismic corrections were not requested). Column groups per phase
    include modal fractions (mol%, wt%, vol%), thermodynamic properties
    (ρ, Vp, Vs, Cp, α, S, H, K, G), and oxide compositions (mol%, wt%) plus
    apfu element compositions. A companion metadata file is written at
    `fileout_metadata.txt`; the CSV is written to `fileout_inlined.csv`.

    Parameters
    ----------
    out : Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}}
        MAGEMin minimization output for one or more points.
    dtb : String
        Database identifier used to retrieve database metadata (e.g., `"ig"`, `"mp"`).
    fileout : String
        Base path for output files. Two files are created: `fileout_inlined.csv`
        and `fileout_metadata.txt`.
    use_Warr2021 : Bool, optional
        If `true`, mineral phase names used as column prefixes are converted to
        IMA-CNMNC approved symbols following Warr (2021). Names with no official
        symbol are returned with a trailing `"*"`. Default is `false`.
    use_GPA : Bool, optional
        If `true`, pressure is reported in GPa (column `"P[GPa]"`) instead of kbar
        (column `"P[kbar]"`). Default is `false`.

    Returns
    -------
    nothing
"""
function MAGEMin_data2dataframe_inlined( out:: Union{Vector{gmin_struct{Float64, Int64}}, gmin_struct{Float64, Int64}},dtb,fileout; use_Warr2021 :: Bool = false, use_GPA :: Bool = false)
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

    P_colname  = use_GPA ? "P[GPa]" : "P[kbar]"

    # Here we create the dataframe's header:
    sys_db = DataFrame(             Symbol("point[#]")      => Int64[],
                                    Symbol("X[0.0-1.0]")    => Float64[],
                                    Symbol(P_colname)       => Float64[],
                                    Symbol("T[°C]")         => Float64[],
                                    Symbol("sys_log10(fO2)")    => Float64[],
                                    Symbol("sys_log10(dQFM)")   => Float64[],
                                    Symbol("sys_aH2O")          => Float64[],
                                    Symbol("sys_aSiO2")         => Float64[],
                                    Symbol("sys_aTiO2")         => Float64[],
                                    Symbol("sys_aAl2O3")        => Float64[],
                                    Symbol("sys_aMgO")          => Float64[],
                                    Symbol("sys_aFeO")          => Float64[],
                                    Symbol("meltViscosity[Pa.s]")    => Float64[],
                                    Symbol("sys_density[kg/m3]")    => Float64[],
                                    Symbol("sys_volume[cm3/mol]")   => Float64[],
                                    Symbol("sys_heatCapacity[J/kg/K]")=> Float64[],
                                    Symbol("sys_alpha[1/K]")    => Float64[],
                                    Symbol("sys_Entropy[kJ/K]")  => Float64[],
                                    Symbol("sys_Enthalpy[kJ/mol]")   => Float64[],
                                    Symbol("sys_Vp[km/s]")      => Float64[],
                                    Symbol("sys_Vs[km/s]")      => Float64[],
                                    Symbol("Vp_S[km/s]")      => Float64[],
                                    Symbol("Vs_S[km/s]")      => Float64[],
                                    Symbol("Vp_cor[km/s]")      => Float64[],
                                    Symbol("Vs_cor[km/s]")      => Float64[],
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
        P_val = use_GPA ? out[k].P_kbar * 0.1 : out[k].P_kbar

        part_1 = Dict(  "point[#]"      => k,
                        "X[0.0-1.0]"    => out[k].X[1],
                        P_colname       => P_val,
                        "T[°C]"         => out[k].T_C,
                        "sys_log10(fO2)"    => out[k].fO2[1],
                        "sys_log10(dQFM)"   => out[k].dQFM[1],
                        "sys_aH2O"          => out[k].aH2O,
                        "sys_aSiO2"         => out[k].aSiO2,
                        "sys_aTiO2"         => out[k].aTiO2,
                        "sys_aAl2O3"        => out[k].aAl2O3,
                        "sys_aMgO"          => out[k].aMgO,
                        "sys_aFeO"          => out[k].aFeO,
                        "meltViscosity[Pa.s]"    => out[k].eta_M,
                        "sys_density[kg/m3]"    => out[k].rho,
                        "sys_volume[cm3/mol]"   => out[k].V,
                        "sys_heatCapacity[J/kg/K]"=> out[k].s_cp,
                        "sys_alpha[1/K]"    => out[k].alpha,
                        "sys_Entropy[kJ/K]"  => out[k].entropy,
                        "sys_Enthalpy[kJ/mol]"   => out[k].enthalpy,
                        "sys_Vp[km/s]"      => out[k].Vp,
                        "sys_Vs[km/s]"      => out[k].Vs,
                        "Vp_S[km/s]"        => out[k].Vp_S,
                        "Vs_S[km/s]"        => out[k].Vs_S,
                        "Vp_cor[km/s]"      => out[k].Vp_cor,
                        "Vs_cor[km/s]"      => out[k].Vs_cor,
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
        ph     = ph_names[i]
        ph_col = use_Warr2021 ? get_Warr_name(ph) : ph
        ss_db_ =     DataFrame(     Symbol("$(ph_col)_mode[mol%]")    => Float64[],
                                    Symbol("$(ph_col)_mode[wt%]")     => Float64[],
                                    Symbol("$(ph_col)_mode[vol%]")     => Float64[],
                                    Symbol("$(ph_col)_density[kg/m3]")    => Float64[],
                                    Symbol("$(ph_col)_volume[cm3/mol]")   => Float64[],
                                    Symbol("$(ph_col)_heatCapacity[J/kg/K]")=> Float64[],
                                    Symbol("$(ph_col)_alpha[1/K]")    => Float64[],
                                    Symbol("$(ph_col)_Entropy[kJ/K]")  => Float64[],
                                    Symbol("$(ph_col)_Enthalpy[kJ/mol]")   => Float64[],
                                    Symbol("$(ph_col)_Vp[km/s]")      => Float64[],
                                    Symbol("$(ph_col)_Vs[km/s]")      => Float64[],
                                    Symbol("$(ph_col)_BulkMod[GPa]")  => Float64[],
                                    Symbol("$(ph_col)_ShearMod[GPa]") => Float64[])

        for i in out[1].oxides
            col = ph_col*"_"*i*"[mol%]"
            ss_db_[!, col] = Float64[]
        end

        for i in out[1].oxides
            col = ph_col*"_"*i*"[wt%]"
            ss_db_[!, col] = Float64[]
        end

        for i in out[1].elements
            col = ph_col*"_"*i*"[apfu]"
            ss_db_[!, col] = Float64[]
        end

        ss_db = hcat(ss_db,ss_db_)
    end
  

    @showprogress "Saving data to csv..." for k=1:np
        row     = Dict{String, Float64}()
        for j=1:n_ss
            ph      = ph_names[j]
            ph_col  = use_Warr2021 ? get_Warr_name(ph) : ph
            if ph in out[k].ph
                i = findfirst(out[k].ph .== ph)
                ss_part_1 = Dict(   "$(ph_col)_mode[mol%]"    => out[k].ph_frac[i].*100.0,
                        "$(ph_col)_mode[wt%]"     => out[k].ph_frac_wt[i].*100.0,
                        "$(ph_col)_mode[vol%]"     => out[k].ph_frac_vol[i].*100.0,
                        "$(ph_col)_density[kg/m3]"    => out[k].SS_vec[i].rho,
                        "$(ph_col)_volume[cm3/mol]"   => out[k].SS_vec[i].V,
                        "$(ph_col)_heatCapacity[J/kg/K]"=> out[k].SS_vec[i].cp,
                        "$(ph_col)_alpha[1/K]"    => out[k].SS_vec[i].alpha,
                        "$(ph_col)_Entropy[kJ/K]"  => out[k].SS_vec[i].entropy,
                        "$(ph_col)_Enthalpy[kJ/mol]"   => out[k].SS_vec[i].enthalpy,
                        "$(ph_col)_Vp[km/s]"      => out[k].SS_vec[i].Vp,
                        "$(ph_col)_Vs[km/s]"      => out[k].SS_vec[i].Vs,
                        "$(ph_col)_BulkMod[GPa]"  => out[k].SS_vec[i].bulkMod,
                        "$(ph_col)_ShearMod[GPa]" => out[k].SS_vec[i].shearMod )

                ss_part_2 = Dict(  (ph_col*"_"*out[1].oxides[j]*"[mol%]" => out[k].SS_vec[i].Comp[j]*100.0)
                                for j in eachindex(out[1].oxides))

                ss_part_3 = Dict(  (ph_col*"_"*out[1].oxides[j]*"[wt%]" => out[k].SS_vec[i].Comp_wt[j]*100.0)
                                for j in eachindex(out[1].oxides))

                ss_part_4 = Dict(  (ph_col*"_"*out[1].elements[j]*"[apfu]" => out[k].SS_vec[i].Comp_apfu[j])
                                for j in eachindex(out[1].elements))
            else
                ss_part_1 = Dict(   "$(ph_col)_mode[mol%]"    => NaN,
                                    "$(ph_col)_mode[wt%]"     => NaN,
                                    "$(ph_col)_mode[vol%]"     => NaN,
                                    "$(ph_col)_density[kg/m3]"    => NaN,
                                    "$(ph_col)_volume[cm3/mol]"   => NaN,
                                    "$(ph_col)_heatCapacity[J/kg/K]"=> NaN,
                                    "$(ph_col)_alpha[1/K]"    => NaN,
                                    "$(ph_col)_Entropy[kJ/K]"  => NaN,
                                    "$(ph_col)_Enthalpy[kJ/mol]"   => NaN,
                                    "$(ph_col)_Vp[km/s]"      => NaN,
                                    "$(ph_col)_Vs[km/s]"      => NaN,
                                    "$(ph_col)_BulkMod[GPa]"  => NaN,
                                    "$(ph_col)_ShearMod[GPa]" => NaN )

                ss_part_2 = Dict(  (ph_col*"_"*out[1].oxides[j]*"[mol%]" => NaN)
                                for j in eachindex(out[1].oxides))

                ss_part_3 = Dict(  (ph_col*"_"*out[1].oxides[j]*"[wt%]" => NaN)
                                for j in eachindex(out[1].oxides))

                ss_part_4 = Dict(  (ph_col*"_"*out[1].elements[j]*"[apfu]" => NaN)
                                for j in eachindex(out[1].elements))
            end
            row_    = merge(ss_part_1,ss_part_2,ss_part_3,ss_part_4)
            row = merge(row,row_)
        end
        push!(ss_db, row, cols=:union)
    end

    if n_pp > 0
        pp_db   = DataFrame()
        for i=1:n_pp
            ph     = ph_names[i+n_ss]
            ph_col = use_Warr2021 ? get_Warr_name(ph) : ph
            pp_db_ =     DataFrame(     Symbol("$(ph_col)_mode[mol%]")    => Float64[],
                                        Symbol("$(ph_col)_mode[wt%]")     => Float64[],
                                        Symbol("$(ph_col)_mode[vol%]")     => Float64[],
                                        Symbol("$(ph_col)_density[kg/m3]")    => Float64[],
                                        Symbol("$(ph_col)_volume[cm3/mol]")   => Float64[],
                                        Symbol("$(ph_col)_heatCapacity[J/kg/K]")=> Float64[],
                                        Symbol("$(ph_col)_alpha[1/K]")    => Float64[],
                                        Symbol("$(ph_col)_Entropy[kJ/K]")  => Float64[],
                                        Symbol("$(ph_col)_Enthalpy[kJ/mol]")   => Float64[],
                                        Symbol("$(ph_col)_Vp[km/s]")      => Float64[],
                                        Symbol("$(ph_col)_Vs[km/s]")      => Float64[],
                                        Symbol("$(ph_col)_BulkMod[GPa]")  => Float64[],
                                        Symbol("$(ph_col)_ShearMod[GPa]") => Float64[])

            for i in out[1].oxides
                col = ph_col*"_"*i*"[mol%]"
                pp_db_[!, col] = Float64[]
            end

            for i in out[1].oxides
                col = ph_col*"_"*i*"[wt%]"
                pp_db_[!, col] = Float64[]
            end

            for i in out[1].elements
                col = ph_col*"_"*i*"[apfu]"
                pp_db_[!, col] = Float64[]
            end
            pp_db = hcat(pp_db,pp_db_)
        end


        @showprogress "Saving data to csv..." for k=1:np
            row     = Dict{String, Float64}()
            for j=1:n_pp
                ph      = ph_names[j+n_ss]
                ph_col  = use_Warr2021 ? get_Warr_name(ph) : ph
                if ph in out[k].ph
                    i = findfirst(out[k].ph .== ph)
                    s = out[k].n_SS
                    pp_part_1 = Dict(   "$(ph_col)_mode[mol%]"    => out[k].ph_frac[i].*100.0,
                                        "$(ph_col)_mode[wt%]"     => out[k].ph_frac_wt[i].*100.0,
                                        "$(ph_col)_mode[vol%]"     => out[k].ph_frac_vol[i].*100.0,
                                        "$(ph_col)_density[kg/m3]"    => out[k].PP_vec[i-s].rho,
                                        "$(ph_col)_volume[cm3/mol]"   => out[k].PP_vec[i-s].V,
                                        "$(ph_col)_heatCapacity[J/kg/K]"=> out[k].PP_vec[i-s].cp,
                                        "$(ph_col)_alpha[1/K]"    => out[k].PP_vec[i-s].alpha,
                                        "$(ph_col)_Entropy[kJ/K]"  => out[k].PP_vec[i-s].entropy,
                                        "$(ph_col)_Enthalpy[kJ/mol]"   => out[k].PP_vec[i-s].enthalpy,
                                        "$(ph_col)_Vp[km/s]"      => out[k].PP_vec[i-s].Vp,
                                        "$(ph_col)_Vs[km/s]"      => out[k].PP_vec[i-s].Vs,
                                        "$(ph_col)_BulkMod[GPa]"  => out[k].PP_vec[i-s].bulkMod,
                                        "$(ph_col)_ShearMod[GPa]" => out[k].PP_vec[i-s].shearMod )

                    pp_part_2 = Dict(  (ph_col*"_"*out[1].oxides[j]*"[mol%]" => out[k].PP_vec[i-s].Comp[j]*100.0)
                                    for j in eachindex(out[1].oxides))

                    pp_part_3 = Dict(  (ph_col*"_"*out[1].oxides[j]*"[wt%]" => out[k].PP_vec[i-s].Comp_wt[j]*100.0)
                                    for j in eachindex(out[1].oxides))

                    pp_part_4 = Dict(  (ph_col*"_"*out[1].elements[j]*"[apfu]" => out[k].PP_vec[i-s].Comp_apfu[j])
                                    for j in eachindex(out[1].elements))
                else
                    pp_part_1 = Dict(   "$(ph_col)_mode[mol%]"    => NaN,
                                        "$(ph_col)_mode[wt%]"     => NaN,
                                        "$(ph_col)_mode[vol%]"     => NaN,
                                        "$(ph_col)_density[kg/m3]"    => NaN,
                                        "$(ph_col)_volume[cm3/mol]"   => NaN,
                                        "$(ph_col)_heatCapacity[J/kg/K]"=> NaN,
                                        "$(ph_col)_alpha[1/K]"    => NaN,
                                        "$(ph_col)_Entropy[kJ/K]"  => NaN,
                                        "$(ph_col)_Enthalpy[kJ/mol]"   => NaN,
                                        "$(ph_col)_Vp[km/s]"      => NaN,
                                        "$(ph_col)_Vs[km/s]"      => NaN,
                                        "$(ph_col)_BulkMod[GPa]"  => NaN,
                                        "$(ph_col)_ShearMod[GPa]" => NaN )

                    pp_part_2 = Dict(  (ph_col*"_"*out[1].oxides[j]*"[mol%]" => NaN)
                                    for j in eachindex(out[1].oxides))

                    pp_part_3 = Dict(  (ph_col*"_"*out[1].oxides[j]*"[wt%]" => NaN)
                                    for j in eachindex(out[1].oxides))

                    pp_part_4 = Dict(  (ph_col*"_"*out[1].elements[j]*"[apfu]" => NaN)
                                    for j in eachindex(out[1].elements))
                end
                row_    = merge(pp_part_1,pp_part_2,pp_part_3,pp_part_4)
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