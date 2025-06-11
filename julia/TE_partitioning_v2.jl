

#= generalization of the KDs database =#


"""
    Holds the partitioning coefficient database
"""
struct custom_KDs_database
    infos           #:: String
    element_name    #:: Vector{String}
    phase_name      #:: Vector{String}, String
    KDs_expr        #:: Matrix{Expr}, Vector{Expr}
end


function retrieve_eval_rules_TE()

    in_eval_TE     = ["T_C","P_kbar"]
    out_eval_TE    = ["out.T_C","out.P_kbar"]

    return in_eval_TE, out_eval_TE
end



function convert_SS_eval_TE(str)
    pattern = r"\[:([A-Za-z_][A-Za-z0-9_]*)\]"
    matches = collect(eachmatch(pattern, str))

    if !isempty(matches)
        for i in eachindex(matches)
            name = matches[i].captures[1]
            println(name)

            str = replace(str, matches[i].match => "out.SS_vec[out.SS_syms[:$name]]")
        end
    end

    return str
end


"""
    KDs_database = custom_KDs_database(infos::String, 
                        element_name::Vector{String}, 
                        phase_name::Vector{String}, 
                        KDs_expr::Matrix{Expr})

Create a custom KDs database from the given information.

returns a custom_KDs_database object that can be used in the TE_partitioning.jl module.

"""
function create_custom_KDs_database(el_name         :: Vector{String},
                                    phase_name      :: Vector{String},
                                    KDs_expr_str    :: Vector{String} ;
                                    info            :: String = "Custom KDs database")

    n_el    = length(el_name)
    n_ph    = length(phase_name)

    KDs_expr = Matrix{Expr}(undef, n_ph, n_el)
    in_eval_TE, out_eval_TE = retrieve_eval_rules_TE()

    for i = 1:n_ph
        for j = 1:n_el

            expr_str      = KDs_expr_str[i]
            for i=1:length(in_eval_TE)
                expr_str = replace(expr_str, in_eval_TE[i] => out_eval_TE[i])
            end

            expr_str      = convert_SS_eval_TE(expr_str)

            KDs_expr[i,j] = Meta.parse("$(expr_str) +0")
        end
    end

    return custom_KDs_database(info, el_name, phase_name, KDs_expr)
end





if 0==0         #run test
    using MAGEMin_C
    data    = Initialize_MAGEMin("mp", verbose=1, solver=0);
    P,T     = 6.0, 699.0
    Xoxides = ["SiO2";  "TiO2";  "Al2O3";  "FeO";   "MnO";   "MgO";   "CaO";   "Na2O";  "K2O"; "H2O"; "O"];
    X       = [58.509,  1.022,   14.858, 4.371, 0.141, 4.561, 5.912, 3.296, 2.399, 10.0, 0.2];
    sys_in  = "wt"
    out     = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, name_solvus=true)
    Finalize_MAGEMin(data)



    # create database on the fly
    el      = ["Li"]
    ph      = ["q","afs","pl","bi","opx","cd","mu","amp","fl","cpx","g"]
    KDs     = ["0.17";"0.14 * T_C/1000.0 + [:bi].compVariables[1]";"0.33 + 0.01*P_kbar";"1.67 * P_kbar / 10.0 + T_C/1000.0";"0.2";"125";"0.82";"0.2";"0.65";"0.26";"0.01"] 

    KDs_database = create_custom_KDs_database(el, ph, KDs)


    # input data
    liq_wt      = out.frac_M_wt
    C0          = 100.0

    KDs         = KDs_database.KDs_expr
    phase_name  = KDs_database.phase_name

    ph, ph_wt   =  MAGEMin_C.mineral_classification(out, "mp");
    TE_ph       =  intersect(ph, phase_name);

    # get indexes of the phase with respect to MAGEMin output and TE_database
    MM_ph_idx   = [findfirst(isequal(x), ph) for x in TE_ph]
    TE_ph_idx   = [findfirst(isequal(x), phase_name) for x in TE_ph]

    # normalize phase fractions
    sum_ph_frac = sum(ph_wt[MM_ph_idx]);
    liq_wt_norm = liq_wt/(sum_ph_frac+liq_wt);
    ph_wt_norm  = ph_wt[MM_ph_idx]./sum_ph_frac;
    ph_TE       = ph[MM_ph_idx];

    P_kbar = out.P_kbar
    nel         = size(KDs_database.KDs_expr)[2]
    np          = length(ph_TE)
    Cmin        = zeros(size(KDs[TE_ph_idx,:])) 
    D           = zeros(Float64, nel)


    for j=1:nel
        for i=1:np
            expr          = KDs_database.KDs_expr[TE_ph_idx[i],j]
            val_expr      = eval(expr)
            D           .+= val_expr* ph_wt_norm[i]
            Cmin[i,j]     = val_expr 
        end
    end

    Cliq        = C0 ./ (D .+ liq_wt_norm.*(1.0 .- D));
    Csol        = (C0 .- Cliq .*  liq_wt_norm) ./ (1.0 .- liq_wt_norm)

    for i = 1:np
        for j=1:nel
            Cmin[i,j] =  Cmin[i,j] * Cliq[j];
        end
    end

end

