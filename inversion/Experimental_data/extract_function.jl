mutable struct data{ _T  } 
    authorID    ::  String
    publication ::  String
    X_unit      ::  String
    P_unit      ::  String
    T_unit      ::  String

    run         ::  String
    wox         ::  Vector{String}
    X           ::  Vector{Float64}
    P           ::  _T
    T           ::  _T 
    dFMQ        ::  _T

    ph          ::  Vector{String}
    ph_frac     ::  Vector{Float64}
    ph_comp     ::  Matrix{Float64}
end

function get_id_string(dataArray,data)
    index = -1;
    for i in eachindex(dataArray)
        if dataArray[i] == data
            index = i;
            break;
        end
    end

    return index;
end


function get_ph_ph_frac(all_ph_frac,ph_list)
    ph = [];
    ph_frac = [];

    for i in eachindex(all_ph_frac)
        if (all_ph_frac[i] != 0.0)
            push!(ph,ph_list[i]);
            push!(ph_frac,all_ph_frac[i])
        end
    end

    return ph, ph_frac;
end


function get_phase_comp(run,ph,xf,wox)
    ix      = zeros(Int64,length(ph)) .-1;
    locWox  = zeros(Int64,length(wox));
    sheets  = XLSX.sheetnames(xf);
    # get position of compositions

    act     = 0;
    for i in eachindex(ph)

        phOn     = get_id_string(sheets,ph[i])
        if (phOn != -1)

            ph_run   = string.(xf[ph[i]][:][1,1:end]);
            id       = get_id_string(ph_run,run);
            if (id != -1)
                ix[i] = id;
                act = 1;
            end

        end
    end

    # map oxides indices, note that if id = -1, content will be set to 0.0

    if act == 1
        ph_comp = zeros(Float64,length(wox),length(ph));
        for i in eachindex(ph)
            # print(i)
            if ix[i] == -1
                vals = zeros(Float64,length(wox)) .+ NaN;
                ph_comp[:,i] = vals;
            else

                oxlist = xf[ph[i]][:][:,1];
                for i in eachindex(wox)
                    locWox[i] = get_id_string(oxlist,wox[i])
                end            

                vals = zeros(Float64,length(wox));
                tmp = xf[ph[i]][:][:,ix[i]];
                for j in eachindex(wox)
                    if locWox[j] == -1
                        vals[j] = 0.0;
                    else
                        vals[j] = tmp[locWox[j]];
                    end
                end
            end

            ph_comp[:,i] = vals;
        end
    else
        ph_comp = [;;];
    end

    return ph_comp;
end


function get_bulk(id,infos,wox)
    X = zeros(length(wox));

    ox  = infos[:,1];
    col = id + 1;
    for i in eachindex(wox)
        ix = get_id_string(ox,wox[i])
        if ix == -1
            X[i] = 0.0;
        else
            X[i] = infos[ix,col];
        end
    end

    return X    
end