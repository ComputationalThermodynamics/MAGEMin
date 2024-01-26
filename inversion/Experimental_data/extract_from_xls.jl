#=
    NR 21-07-2023
    Script to extract experimental data from xls file (copy/paste from pdf articles)

    The goal here is to copy and paste table from published experimental work into an excel spreadsheet and extract/sort all needed data
    Although excel is evil, it works decently to copy and paste from pdf

    *xls file must be formatted correctly! (check xls files available)
    *one spreadsheet for mode
    *one spreadsheet for mineral composition

    NB: I did not bother with type stability at all, as reading a few dozens of files is not a bottleneck in term of performance 
=#

using XLSX
include("extract_function.jl")

#=
    Below is where articles can be added
=#
include("publication_list.jl")


#=
    Effectively builds the dabase
=#
exp_data = [];
for p in eachindex(list_articles)

    xf      = XLSX.readxlsx(list_articles[p].file);
    sh      = XLSX.sheetnames(xf);
    wox     = list_articles[p].wox;
    X_unit  = list_articles[p].X_unit;
    P_unit  = list_articles[p].P_unit;
    T_unit  = list_articles[p].T_unit;

    # retrieve max stable phase list
    ph_list   = [];
    all_ph    = xf["modes"][:][:,1];
    for i in eachindex(all_ph)
        if all_ph[i] != "Run" && all_ph[i] != "P" && all_ph[i] != "T" && all_ph[i] != "Bulk" && all_ph[i] != "dFMQ"
        push!(ph_list,all_ph[i])
        end
    end

    nExp      = length(xf["modes"][:][1,:])-1;
    run_names = string.(xf["modes"][:][1,2:end]);
    modInfo   = string.(xf["modes"][:][:,1]);

    Pid       = get_id_string(modInfo,"P");
    Tid       = get_id_string(modInfo,"T");
    Xid       = get_id_string(modInfo,"Bulk");
    dFMQid    = get_id_string(modInfo,"dFMQ");

    for i=1:nExp
        print(i)
        authorID    = xf["info"][:][1,2];
        publication = xf["info"][:][2,2];
        P           = xf["modes"][:][Pid,i+1];
        T           = xf["modes"][:][Tid,i+1];
        dFMQ        = parse(Float64,string.(xf["modes"][:][dFMQid,i+1]));
        X           = get_bulk(xf["modes"][:][Xid,i+1],xf["info"][:],wox);
        X_unit      = X_unit;                       # just for clarity
        P_unit      = P_unit;                       # just for clarity
        T_unit      = T_unit;                       # just for clarity
        run         = run_names[i];
        wox         = wox;

        din         = parse.(Float64,string.(xf["modes"][:][6:end,i+1]));
        ph, ph_frac = get_ph_ph_frac(din,ph_list);
        ph_comp     = get_phase_comp(run,ph,xf,wox)

        push!(exp_data,data{Float64}(authorID,publication,X_unit,P_unit,T_unit,run,wox,X,P,T,dFMQ,ph,ph_frac,ph_comp));

    end

end

#=
    print database
=#
for i in eachindex(exp_data)
    print(exp_data[i],"\n")
end