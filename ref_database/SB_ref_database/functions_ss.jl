# how to compute chear modulus:
# mu = mu0 + (P-Pr)*dmu/dP + (T-Tr)*dmu/dT
# where r is the reference pressure and temperature usually, 1bar and 298.15K
using LinearAlgebra

global tab                  = "    "
global tab2                  = "        "
global tab3                  = "            "

struct Phase
    id::String                  # Name id
    abbrev::String              # Abbreviation
    fml::String                 # Chemical formula
    oxides::Dict{String, Float64}  # Oxides
    F0::Float64                 # Helmoltz energy (F0, J/mol)
    n::Float64                  # negative of the number of atoms per formula unit (-n)
    V0::Float64                 # negative of the volume (-V0)
    K0::Float64                 # c1: isothermal bulk modulus (K0, bar)
    Kp::Float64                 # c2: pressure derivative of the isothermal bulk modulus (K')
    Θ0::Float64                 # c3: Debye Temperature (Θ0, K)
    γ0::Float64                 # c4: Gruneisen thermal parameter (γ0)
    q0::Float64                 # c5: Mie-Gruneisen exponent (q0)
    ηS0::Float64                # c6: Shear strain derivative of the tensorial Gruneisen parameter (ηS0)
    cme::Float64                # c7: Configurational (and magnetic) entropy (J/mol/K)
    aSm::Float64                # Ambient shear modulus (GPa)
    pd::Float64                 # Pressure derivative
    td::Float64                 # Temperature derivative
end

function read_data(fname::String)
    return JSON3.read(fname, Vector{Phase}) |> DataFrame
end

data2 = read_data("stx11_data.json")
@save "STIX11.jld2" data2

# Function to generate a simplex in n dimensions with given step size
function generate_simplex(n, step)
    # Generate all possible points with the given step size
    points = collect(0:step:1)
    
    # Generate all combinations of points in n dimensions
    combinations = Iterators.product(ntuple(_ -> points, n)...)
    
    # Filter combinations to keep only those where the sum of coordinates is 1
    simplex = [collect(comb) for comb in combinations if sum(comb) ≈ 1.0]
    
    return simplex
end


function get_w_ids(ss,ii,data2)
    ks          = keys(ss[ii].margules)
    # ks_strings  = [String(s) for s in ks]
    n_em        = length(ss[ii].endmembers)
    em          = keys(ss[ii].endmembers)
    em_strings  = [String(s) for s in em]
    marg_out    = []

    for i=1:n_em-1
        for j=i+1:n_em
            push!(marg_out, sort([em_strings[i],em_strings[j]]))
        end
    end

    val_in = []
    for (key, value) in ss[ii].margules
        push!(val_in, value)
    end

    em_out      = keys(ss[ii].endmembers)
    em_names    = data2[:,:id]
    
    em_ids      = []
    for i in em_out
        id = findfirst(x->x==i,em_names)
        push!(em_ids,id)
    end
    sorted_indexes = sortperm(em_ids)
    em_in       = em_names[em_ids[sorted_indexes]]

    ids_v =  []
    for i in em_out
        id = findfirst(x->x==i,em_in)
        push!(ids_v,id)
    end
    # println(em_in)

    # em_in is the same order as v


    marg_in  = []
    em_list  = []
    j = 1
    for i in ks
        w_em         = split(i,",")
        w_em_strings = [String(s) for s in w_em]
        push!(marg_in, sort(w_em_strings))
        if w_em_strings[1] in em_list
        else
            push!(em_list, w_em_strings[1])
        end
    end

    ids_W = []
    for i in marg_out
        if i in marg_in
            idx = findfirst(x->x==i, marg_in)
            push!(ids_W, idx)
        else
            push!(ids_W, 0)
        end
    end
    return ids_W,ids_v
end


"""
    get_mu_Gex(W, v, n_em, sym)
"""
function get_mu_Gex(W, v, n_em, sym)
    @variables p[1:n_em]
    mu_Gex      = Vector{Num}(undef,n_em)
    phi_p       = Float64.(I(n_em));

    if sym == 0
        mat_phi      = Vector{Num}(undef,n_em)

        sum_v = 0.0
        for i=1:n_em
            sum_v  += p[i]*v[i];
        end
        for i=1:n_em
            mat_phi[i]   =p[i]*v[i]/sum_v;
        end

        mu = 0.
        for i = 1:n_em
            mu = 0;
            it = 0;
            for j = 1:(n_em-1)
                for k = j+1:n_em
                    it = it + 1;
                    mu -= (phi_p[i,j]-mat_phi[j])*(phi_p[i,k]-mat_phi[k])*(W[it]*2*v[i]/(v[j]+v[k]));
                end
            end
            mu_Gex[i]= mu 
        end
        
    elseif sym == 1

        for i = 1:n_em
            mu = 0;
            it = 0;
            for j = 1:n_em-1
                for k = j+1:n_em
                    it  = it + 1;
                    mu -= (phi_p[i,j] - p[j])*(phi_p[i,k]-p[k])*W[it];
              
                end
            end
           
            mu_Gex[i]       = mu; 
        end
    end

    return mu_Gex
end

# Function to adjust the simplex points
function adjust_simplex(simplex, eps)
    adjusted_simplex = []
    for point in simplex
        adjusted_point = copy(point)
        for i in 1:length(point)
            if point[i] == 0.0
                adjusted_point[i] = eps
            elseif point[i] == 1.0
                adjusted_point[i] = 1.0 - eps
            end
        end
        # Adjust the remaining coordinates to ensure the sum is still 1
        sum_adjusted = sum(adjusted_point)
        if sum_adjusted ≈ 1.0
            push!(adjusted_simplex, adjusted_point)
        else
            # Find the index of the largest coordinate
            max_index = argmax(adjusted_point)
            adjusted_point[max_index] += 1.0 - sum_adjusted
            push!(adjusted_simplex, adjusted_point)
        end
    end
    return adjusted_simplex
end

struct ModelJSON
    name            :: String
    abbrev          :: String
    endmembers      :: Dict{String, Vector{String}}
    margules        :: Dict{String, Float64}
    van_laar        :: Vector{Float64}
end

function retrieve_site_cmp(ss, i)

    em          = String.(keys(ss[i].endmembers))
    n_em        = length(em)

    fml         = ss[i].endmembers[em[1]][2]

    matches     = eachmatch(r"\[([^\[\]]+)\]", fml)
    contents    = [m.captures[1] for m in matches]

    n_sf        = length(contents)

    mul         = zeros(Float64, n_sf)
    site_cmp    = zeros(Float64, n_sf,n_el, n_em)

    for j=1:n_em
        fml         = ss[i].endmembers[em[j]][2]
        matches     = eachmatch(r"\[([^\[\]]+)\]", fml)
        contents    = [m.captures[1] for m in matches]

        for k=1:length(contents)
            matches = eachmatch(r"_(\d+)", contents[k])
            n_atoms = [parse(Int, m.captures[1]) for m in matches]
            matches = eachmatch(r"[A-Za-z]+", contents[k])
            elements= [m.match for m in matches]
            mul[k]  = sum(n_atoms)

            if length(n_atoms) == 1
                id              = findfirst(elems .== elements[1])
                site_cmp[k,id, j] = Float64(n_atoms[1])
            else
                for l=1:length(n_atoms)
                    id              = findfirst(elems .== elements[l])
                    site_cmp[k,id, j] = Float64(n_atoms[l])
                end
            end
        end
    end


    return mul, site_cmp
end



function get_sb_gss_init_function(sb_ver,ss)
    sb_gss_init_function    = ""

    n_ss = length(ss)
    for i = 1:n_ss

        mul, site_cmp = retrieve_site_cmp(ss, i)

        W   = [ ss[i].margules[j] for j in  keys(ss[i].margules)]

        v   = [ ss[i].van_laar[j] for j in  keys(ss[i].van_laar)]

        em  = [ ss[i].endmembers[j][1] for j in  keys(ss[i].endmembers)]

        n_sf    = size(site_cmp)[1]
        n_ox    = size(site_cmp)[2]
        n_em    = size(site_cmp)[3]
        M   = Float64[]
        C   = Vector{Float64}[]
        for k=1:n_sf
            for l=1:n_ox
                if ~all(site_cmp[k,l,:] .== 0.0)
                    push!(C,site_cmp[k,l,:]./mul[k])
                    push!(M,mul[k])
                end
            end
        end
        C   = hcat(C...)'

        # println("C: ", C)
        # println(mul, site_cmp, W, v, em)

        sym = 1
        if ~isempty(v)
            sym = 0
        end
        
        sb_gss_init_function *= "/**\n"
        sb_gss_init_function *= "    allocate memory for $(sb_ver)_$(ss[i].abbrev)\n"
        sb_gss_init_function *= "*/\n"
        sb_gss_init_function *= "SS_ref G_SS_$(sb_ver)_$(ss[i].abbrev)_init_function(SS_ref SS_ref_db,  global_variable gv){\n\n"
        sb_gss_init_function *= "    SS_ref_db.is_liq    = 0;\n"
        sb_gss_init_function *= "    SS_ref_db.symmetry  = $sym;\n"
        sb_gss_init_function *= "    SS_ref_db.n_cat     = $(size(C,1));\n"
        sb_gss_init_function *= "    SS_ref_db.n_xeos    = $(length(em));\n"
        sb_gss_init_function *= "    SS_ref_db.n_em      = $(length(em));\n"
        sb_gss_init_function *= "    SS_ref_db.n_sf      = $(length(mul));\n"
        sb_gss_init_function *= "    SS_ref_db.n_w       = $(length(W));\n"
        if sym == 0
            sb_gss_init_function *= "    SS_ref_db.n_v       = $(length(v));\n"
        end
        sb_gss_init_function *= "\n"
        sb_gss_init_function *= "     return SS_ref_db;\n"
        sb_gss_init_function *= "}\n\n"
    end

    sb_gss_init_function *= "void SS_init_sb11(	    SS_init_type 		*SS_init,\n"
    sb_gss_init_function *= "                            global_variable 	 gv				){\n\n"
    sb_gss_init_function *= "$(tab)for (int iss = 0; iss < gv.len_ss; iss++){\n"
    for i = 1:n_ss
        if i == 1
            sb_gss_init_function *= "$(tab)$(tab)if      (strcmp( gv.SS_list[iss], \"$(ss[i].abbrev)\")  == 0 ){\n"
            sb_gss_init_function *= "$(tab)$(tab)$(tab)SS_init[iss]  = G_SS_$(sb_ver)_$(ss[i].abbrev)_init_function; 		}\n"
        else
            sb_gss_init_function *= "$(tab)$(tab)else if (strcmp( gv.SS_list[iss], \"$(ss[i].abbrev)\")  == 0 ){\n"
            sb_gss_init_function *= "$(tab)$(tab)$(tab)SS_init[iss]  = G_SS_$(sb_ver)_$(ss[i].abbrev)_init_function; 		}\n"
        end
    end
    sb_gss_init_function *= "$(tab)$(tab)else{\n"
    sb_gss_init_function *= "$(tab)$(tab)$(tab)printf(\"\\nsolid solution '%s' is not in the database, cannot be initiated\\n\", gv.SS_list[iss]);\n"
    sb_gss_init_function *= "$(tab)$(tab)}\n"
    sb_gss_init_function *= "$(tab)}\n"
    sb_gss_init_function *= "}\n"
    return sb_gss_init_function
end
    

function get_sb_objective_functions(sb_ver,ss)
    sb_objective_functions    = ""

    n_ss = length(ss)
    for ii = 1:n_ss

        mul, site_cmp = retrieve_site_cmp(ss, ii)

        W   = [ ss[ii].margules[j] for j in  keys(ss[ii].margules)]

        v   = [ ss[ii].van_laar[j] for j in  keys(ss[ii].van_laar)]

        em  = [ ss[ii].endmembers[j][1] for j in  keys(ss[ii].endmembers)]

        sym = 1
        if ~isempty(v)
            sym = 0
        end

        n_sf    = size(site_cmp)[1]
        n_ox    = size(site_cmp)[2]
        n_em    = size(site_cmp)[3]
        M   = Float64[]
        C   = Vector{Float64}[]
        for k=1:n_sf
            for l=1:n_ox
                if ~all(site_cmp[k,l,:] .== 0.0)
                    push!(C,site_cmp[k,l,:]./mul[k])
                    push!(M,mul[k])
                end
            end
        end
        C       = hcat(C...)'
        n_cat   = size(C,1)

        @variables p[0:n_em-1] R T gb[0:n_em-1] ape[0:n_em-1] fbc
        X               = Symbolics.scalarize(p)
        Gref            = Symbolics.scalarize(gb)
        A               = Symbolics.scalarize(ape)
        Xo              = C*X
        fac             = fbc/sum(A.*X)

        config          = R * T * (M' * Diagonal(Xo) * log.(Xo))
        grad_config     = Symbolics.gradient(config, X)
        mu_Gex          = get_mu_Gex(W, v, n_em, sym)

  

        G = fac*(Gref'*X + mu_Gex'*X + config);
        grad_G     = Symbolics.gradient(G, X)
        
        grad_fac = Symbolics.gradient(fac, X)
        # for i=1:n_em
        #     println("grad_fac $(grad_fac[i])")
        # end

        # -ape[0]*(fbc / ((ape[0]*p[0] + ape[1]*p[1])^2))
        # grad_muGex = Symbolics.gradient(mu_Gex', X)
        # println("mu_Gex $(mu_Gex)")
        # for i = 1:n_em
        #     println("grad_config $(grad_config[i])")
        #     println("grad_G $(grad_G[i])")
        # end
        # println("")


        sb_objective_functions *= "/**\n"
        sb_objective_functions *= "    Objective function for $(sb_ver)_$(ss[ii].abbrev)\n"
        sb_objective_functions *= "*/\n"
        sb_objective_functions *= "double obj_$(sb_ver)_$(ss[ii].abbrev)(unsigned n, const double *x, double *grad, void *SS_ref_db){\n"
        sb_objective_functions *= "$(tab)SS_ref *d         = (SS_ref *) SS_ref_db;\n\n"
        sb_objective_functions *= "$(tab)int n_em          = d->n_em;\n"
        sb_objective_functions *= "$(tab)double P          = d->P;\n"
        sb_objective_functions *= "$(tab)double T          = d->T;\n"
        sb_objective_functions *= "$(tab)double R          = d->R;\n"

        sb_objective_functions *= "$(tab)double *dfx       = d->dfx;\n"
        sb_objective_functions *= "$(tab)double **N        = d->N;\n"
        sb_objective_functions *= "$(tab)double *Vec1      = d->Vec1;\n"
        sb_objective_functions *= "$(tab)double *Vec2      = d->Vec2;\n"
        sb_objective_functions *= "$(tab)double *p         = d->p;\n"
        sb_objective_functions *= "$(tab)double *gb        = d->gb_lvl;\n"
        sb_objective_functions *= "$(tab)double *mu_Gex    = d->mu_Gex;\n"

        if sym == 0
            sb_objective_functions *= "$(tab)double *mat_phi   = d->mat_phi;\n"
        end
        sb_objective_functions *= "\n"
        sb_objective_functions *= "$(tab)for (int i = 0; i < n_em; i++){\n"
        sb_objective_functions *= "$(tab)$(tab)p[i]   = x[i];\n"
        sb_objective_functions *= "$(tab)}\n\n"
        if sym == 0
            sb_objective_functions *= "$(tab)d->sum_v = 0.0;\n"
            sb_objective_functions *= "$(tab)for (int i = 0; i < n_em; i++){\n"
            sb_objective_functions *= "$(tab)$(tab)d->sum_v += p[i]*d->v[i];\n"
            sb_objective_functions *= "$(tab)}\n"
            sb_objective_functions *= "$(tab)for (int i = 0; i < n_em; i++){\n"
            sb_objective_functions *= "$(tab)$(tab)d->mat_phi[i] = (p[i]*d->v[i])/d->sum_v;\n"
            sb_objective_functions *= "$(tab)}\n\n"
            
            sb_objective_functions *= "$(tab)double tmp = 0.0;\n"
            sb_objective_functions *= "$(tab)double Gex = 0.0;\n"
            sb_objective_functions *= "$(tab)for (int i = 0; i < d->n_em; i++){\n"
            sb_objective_functions *= "$(tab)$(tab)Gex = 0.0;\n"
            sb_objective_functions *= "$(tab)$(tab)int it = 0;\n"
            sb_objective_functions *= "$(tab)$(tab)for (int j = 0; j < d->n_xeos; j++){\n"
            sb_objective_functions *= "$(tab)$(tab)$(tab)tmp = (d->eye[i][j] - d->mat_phi[j]);\n"
            sb_objective_functions *= "$(tab)$(tab)$(tab)for (int k = j+1; k < d->n_em; k++){\n"
            sb_objective_functions *= "$(tab)$(tab)$(tab)$(tab)Gex -= tmp*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));\n"
            sb_objective_functions *= "$(tab)$(tab)$(tab)$(tab)it += 1;\n"
            sb_objective_functions *= "$(tab)$(tab)$(tab)}\n"
            sb_objective_functions *= "$(tab)$(tab)}\n"
            sb_objective_functions *= "$(tab)$(tab)mu_Gex[i] = Gex/1000.0;\n"
            sb_objective_functions *= "$(tab)}\n"
        else
            sb_objective_functions *= "$(tab)double tmp = 0.0;\n"
            sb_objective_functions *= "$(tab)double Gex = 0.0;\n"
            sb_objective_functions *= "$(tab)for (int i = 0; i < n_em; i++){\n"
            sb_objective_functions *= "$(tab)$(tab)Gex = 0.0;\n"
            sb_objective_functions *= "$(tab)$(tab)int it    = 0;\n"
            sb_objective_functions *= "$(tab)$(tab)for (int j = 0; j < d->n_xeos; j++){\n"
            sb_objective_functions *= "$(tab)$(tab)$(tab)tmp = (d->eye[i][j] - p[j]);\n"
            sb_objective_functions *= "$(tab)$(tab)$(tab)for (int k = j+1; k < n_em; k++){\n"
            sb_objective_functions *= "$(tab)$(tab)$(tab)$(tab)Gex -= tmp*(d->eye[i][k] - p[k])*(d->W[it]);\n"
            sb_objective_functions *= "$(tab)$(tab)$(tab)$(tab)it += 1;\n"
            sb_objective_functions *= "$(tab)$(tab)$(tab)}\n"
            sb_objective_functions *= "$(tab)$(tab)}\n"
            sb_objective_functions *= "$(tab)$(tab)mu_Gex[i] = Gex/1000.0;\n"
            sb_objective_functions *= "$(tab)}\n\n"
        end

        sb_objective_functions *= "$(tab)d->sum_apep = 0.0;\n"
        sb_objective_functions *= "$(tab)for (int i = 0; i < n_em; i++){\n"
        sb_objective_functions *= "$(tab)$(tab)d->sum_apep += d->ape[i]*p[i];\n"
        sb_objective_functions *= "$(tab)}\n"
        sb_objective_functions *= "$(tab)d->factor = d->fbc/d->sum_apep;\n\n"

        Sconfig = replace(string(config), r"(\d)(?=[\[\(a-zA-Z])" => s"\1*")
        sb_objective_functions *= "$(tab)double Sconfig    = $Sconfig;\n\n"

        # sb_objective_functions *= "$(tab)d->df_raw = Sconfig;\n"
        # sb_objective_functions *= "$(tab)for (int i = 0; i < n_em; i++){\n"
        # sb_objective_functions *= "$(tab)$(tab)d->df_raw += (mu_Gex[i] + gb[i])*p[i];\n"
        # sb_objective_functions *= "$(tab)}\n"
        # sb_objective_functions *= "$(tab)d->df = d->df_raw * d->factor;\n\n"
        sb_objective_functions *= "$(tab)d->df_raw = 0.0;\n"
        sb_objective_functions *= "$(tab)for (int i = 0; i < n_em; i++){\n"
        sb_objective_functions *= "$(tab)$(tab)d->df_raw += (mu_Gex[i] + gb[i] + Sconfig)*p[i];\n"
        sb_objective_functions *= "$(tab)}\n"
        sb_objective_functions *= "$(tab)d->df = d->df_raw * d->factor;\n\n"

        sb_objective_functions *= "$(tab)if (grad){\n"
        for i=1:n_em
            dS = replace(string(grad_config[i]), r"(\d)(?=[\[\(a-zA-Z])" => s"\1*")
            sb_objective_functions *= "$(tab2)grad[$(i-1)] = ($(dS) + mu_Gex[$(i-1)] + gb[$(i-1)])* d->factor - (d->df_raw * d->factor * (d->ape[$(i-1)]/d->sum_apep));\n"
            # PREV
            # dS = replace(string(grad_G[i]), r"(\d)(?=[\[\(a-zA-Z])" => s"\1*")
            # sb_objective_functions *= "$(tab2)grad[$(i-1)] = ($dS)* d->factor;\n"
        end
        sb_objective_functions *= "$(tab)}\n"


        sb_objective_functions *= "$(tab)return d->df;\n"
        sb_objective_functions *= "}\n\n"
    end

    sb_objective_functions *= "/**\n"
    sb_objective_functions *= "$(tab)associate the array of pointer with the right solution phase\n"
    sb_objective_functions *= "*/\n"
    sb_objective_functions *= "void SB_$(sb_ver)_objective_init_function(	obj_type 			*SS_objective,\n"
    sb_objective_functions *= "$(tab3)$(tab3)$(tab3)$(tab)global_variable 	 gv				){	\n"
                            
    sb_objective_functions *= "$(tab)for (int iss = 0; iss < gv.len_ss; iss++){\n"
    for i=1:n_ss
        if i == 1
            sb_objective_functions *= "$(tab2)if      (strcmp( gv.SS_list[iss], \"$(ss[i].abbrev)\")  == 0 ){\n"
            sb_objective_functions *= "$(tab3)SS_objective[iss]  = obj_$(sb_ver)_$(ss[i].abbrev); 		}\n"
        else
            sb_objective_functions *= "$(tab2)else if (strcmp( gv.SS_list[iss], \"$(ss[i].abbrev)\")  == 0 ){\n"
            sb_objective_functions *= "$(tab3)SS_objective[iss]  = obj_$(sb_ver)_$(ss[i].abbrev); 		}\n"
        end
    end
     
    sb_objective_functions *= "$(tab2)else{\n"
    sb_objective_functions *= "$(tab3)printf(\"\\nsolid solution '%s' is not in the database, cannot be initiated\\n\", gv.SS_list[iss]);\n"
    sb_objective_functions *= "$(tab2)}	\n"
    sb_objective_functions *= "$(tab)};\n"	
    sb_objective_functions *= "}\n"



    sb_objective_functions *= "/**\n"
    sb_objective_functions *= "$(tab)associate the array of pointer with the right solution phase\n"
    sb_objective_functions *= "*/\n"
    sb_objective_functions *= "void SB_$(sb_ver)_PC_init(	PC_type 			*PC_read,\n"
    sb_objective_functions *= "$(tab3)$(tab3)global_variable 	 gv				){	\n"
                            
    sb_objective_functions *= "$(tab)for (int iss = 0; iss < gv.len_ss; iss++){\n"
    for i=1:n_ss
        if i == 1
            sb_objective_functions *= "$(tab2)if      (strcmp( gv.SS_list[iss], \"$(ss[i].abbrev)\")  == 0 ){\n"
            sb_objective_functions *= "$(tab3)PC_read[iss]   = obj_$(sb_ver)_$(ss[i].abbrev); 		}\n"
        else
            sb_objective_functions *= "$(tab2)else if (strcmp( gv.SS_list[iss], \"$(ss[i].abbrev)\")  == 0 ){\n"
            sb_objective_functions *= "$(tab3)PC_read[iss]   = obj_$(sb_ver)_$(ss[i].abbrev); 		}\n"
        end
    end
     
    sb_objective_functions *= "$(tab2)else{\n"
    sb_objective_functions *= "$(tab3)printf(\"\\nsolid solution '%s' is not in the database, cannot be initiated\\n\", gv.SS_list[iss]);\n"
    sb_objective_functions *= "$(tab2)}	\n"
    sb_objective_functions *= "$(tab)};\n"	
    sb_objective_functions *= "}\n\n\n"

    sb_objective_functions *= "//headers for objectives functions\n"
    for ii = 1:n_ss
        sb_objective_functions *= "double obj_$(sb_ver)_$(ss[ii].abbrev)(unsigned n, const double *x, double *grad, void *SS_ref_db);\n"
    end


    return sb_objective_functions
end
    


function get_sb_gss_function(sb_ver,ss,data2)
    sb_gss_function    = ""

    n_ss = length(ss)
    for ii = 1:n_ss
        println(ii)
        mul, site_cmp = retrieve_site_cmp(ss, ii)

        W       = [ ss[ii].margules[j] for j in  keys(ss[ii].margules)]

        v       = [ ss[ii].van_laar[j] for j in  keys(ss[ii].van_laar)]

        em      = [ ss[ii].endmembers[j][1] for j in  keys(ss[ii].endmembers)]
        ids_w,ids_v   = get_w_ids(ss,ii,data2)


        # println(mul, site_cmp, W, v, em)

        sym = 1
        if ~isempty(v)
            sym = 0
        end

        n_sf    = size(site_cmp)[1]
        n_ox    = size(site_cmp)[2]
        n_em    = size(site_cmp)[3]
        M   = Float64[]
        C   = Vector{Float64}[]
        for k=1:n_sf
            for l=1:n_ox
                if ~all(site_cmp[k,l,:] .== 0.0)
                    push!(C,site_cmp[k,l,:]./mul[k])
                    push!(M,mul[k])
                end
            end
        end
        C       = hcat(C...)'
        n_cat   = size(C,1)

        sb_gss_function *= "/**\n"
        sb_gss_function *= "    Solution phase data for $(sb_ver)_$(ss[ii].abbrev)\n"
        sb_gss_function *= "*/\n"
        sb_gss_function *= "SS_ref G_SS_$(sb_ver)_$(ss[ii].abbrev)_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){\n"
        sb_gss_function *= "\n"
        sb_gss_function *= "$(tab)int i, j;\n"
        sb_gss_function *= "$(tab)int n_em = SS_ref_db.n_em;\n"
        sb_gss_function *= "\n"
        sb_gss_function *= "$(tab)char   *EM_tmp[] 		= {"
        for i in em
            sb_gss_function *= "\"$(i)\","
        end
        sb_gss_function *= "};\n"
        sb_gss_function *= "$(tab)for (int i = 0; i < SS_ref_db.n_em; i++){\n"
        sb_gss_function *= "$(tab)$(tab)strcpy(SS_ref_db.EM_list[i],EM_tmp[i]);\n"
        sb_gss_function *= "$(tab)};\n\n"

        sb_gss_function *= "$(tab)// Site mixing composition;\n"
        for i = 1:n_cat
            for j = 1:n_em
                sb_gss_function *= "$(tab)SS_ref_db.C[$(i-1)][$(j-1)] = $(C[i,j]);"
            end
            sb_gss_function *= "\n"
        end
        sb_gss_function *= "\n"
        N = nullspace(ones(n_em)')
        sb_gss_function *= "$(tab)// pre-computed Nullspace;\n"
        for i = 1:n_em
            for j = 1:n_em-1
                sb_gss_function *= "$(tab)SS_ref_db.N[$(i-1)][$(j-1)] = $(N[i,j]);"
            end
            sb_gss_function *= "\n"
        end

        sb_gss_function *= "\n"
        for i=1:length(W)
            sb_gss_function *= "$(tab)SS_ref_db.W[$(i-1)] = $(W[ids_w[i]]);\n"
        end
        sb_gss_function *= "\n"
        if ~isempty(v)
            for i=1:length(v)
                sb_gss_function *= "$(tab)SS_ref_db.v[$(i-1)] = $(v[ids_v[i]]);\n"
            end
        end
        sb_gss_function *= "\n"
        for i in em
            sb_gss_function *= "$(tab)em_data $i$(tab)$(tab)$(tab)= get_em_data($(tab)research_group, EM_dataset,\n"
            sb_gss_function *= "$(tab)"^11*"len_ox,\n"
            sb_gss_function *= "$(tab)"^11*"z_b,\n"
            sb_gss_function *= "$(tab)"^11*"SS_ref_db.P,\n"
            sb_gss_function *= "$(tab)"^11*"SS_ref_db.T,\n"
            sb_gss_function *= "$(tab)"^11*"\"$i\""*",\n"
            sb_gss_function *= "$(tab)"^11*"\"equilibrium\");\n\n"
        end

        for i = 1:length(em) 
            sb_gss_function *= "$(tab)SS_ref_db.gbase[$(i-1)]"*"$(tab)"^2*"= $(em[i]).gb;\n"
        end
        sb_gss_function *= "\n"
 
        for i = 1:length(em) 
            sb_gss_function *= "$(tab)SS_ref_db.ElShearMod[$(i-1)]"*"$(tab)"^2*"= $(em[i]).ElShearMod;\n"
        end
        sb_gss_function *= "\n"    

        for i = 1:length(em) 
            sb_gss_function *= "$(tab)SS_ref_db.ElBulkMod[$(i-1)]"*"$(tab)"^2*"= $(em[i]).ElBulkMod;\n"
        end
        sb_gss_function *= "\n"       

        sb_gss_function *= "$(tab)for (i = 0; i < len_ox; i++){\n"
        for i = 1:length(em) 
            sb_gss_function *= "$(tab)$(tab)SS_ref_db.Comp[$(i-1)][i] 	= $(em[i]).C[i];\n"
        end
        sb_gss_function *= "$(tab)}\n\n"   

        sb_gss_function *= "$(tab)for (i = 0; i < n_em; i++){\n"       
        sb_gss_function *= "$(tab)$(tab)SS_ref_db.z_em[i] = 1.0;\n"       
        sb_gss_function *= "$(tab)};\n\n"       

        for i = 1:length(em) 
            sb_gss_function *= "$(tab)SS_ref_db.bounds_ref[$(i-1)][0] = 0.0+eps;  SS_ref_db.bounds_ref[$(i-1)][1] = 1.0-eps;\n"
        end
        sb_gss_function *= "\n\n$(tab)return SS_ref_db;\n"
        sb_gss_function *= "}\n\n"
    end

    sb_gss_function *= "$(tab)SS_ref G_SS_$(sb_ver)_EM_function($(tab)$(tab)global_variable$(tab)$(tab)gv,\n"
    sb_gss_function *= "$(tab)"^10*"SS_ref$(tab)$(tab)$(tab)SS_ref_db,\n"
    sb_gss_function *= "$(tab)"^10*"int$(tab)$(tab)$(tab)EM_dataset,\n"
    sb_gss_function *= "$(tab)"^10*"bulk_info$(tab)$(tab)$(tab)z_b,\n"
    sb_gss_function *= "$(tab)"^10*"char$(tab)$(tab)$(tab)*name){$(tab)$(tab)$(tab)\n"
    sb_gss_function *= "\n"    
    sb_gss_function *= "$(tab)double eps 		   	= gv.bnd_val;\n"    
    sb_gss_function *= "$(tab)double P 			= SS_ref_db.P;\n"    
    sb_gss_function *= "$(tab)double T 			= SS_ref_db.T;\n\n"     
    sb_gss_function *= "$(tab)SS_ref_db.ss_flags[0]  = 1;\n\n"    

    sb_gss_function *= "$(tab)/* Associate the right solid-solution data */\n"    
    sb_gss_function *= "$(tab)for (int FD = 0; FD < gv.n_Diff; FD++){	\n\n"  
    sb_gss_function *= "$(tab)$(tab)if (FD == 8 || FD == 9){\n"  
    sb_gss_function *= "$(tab)$(tab)$(tab)SS_ref_db.P = 1.+ gv.gb_P_eps*gv.pdev[0][FD];\n"  
    sb_gss_function *= "$(tab)$(tab)$(tab)SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];\n"  
    sb_gss_function *= "$(tab)$(tab)}\n"  
    sb_gss_function *= "$(tab)$(tab)else{\n"  
    sb_gss_function *= "$(tab)$(tab)$(tab)SS_ref_db.P = P + gv.gb_P_eps*gv.pdev[0][FD];\n"  
    sb_gss_function *= "$(tab)$(tab)$(tab)SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];\n"  
    sb_gss_function *= "$(tab)$(tab)}\n\n"  

    for ii = 1:n_ss
        if ii == 1
            sb_gss_function *= "$(tab2)if (strcmp( name, \"$(ss[ii].abbrev)\") == 0 ){\n"    
        else
            sb_gss_function *= "$(tab2)else if (strcmp( name, \"$(ss[ii].abbrev)\") == 0 ){\n"    
        end
        sb_gss_function *= "$(tab3)SS_ref_db  = G_SS_$(sb_ver)_$(ss[ii].abbrev)_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);	}\n"  
    end
    sb_gss_function *= "$(tab2)else{\n"
    sb_gss_function *= "$(tab3)printf(\"\\nsolid solution '%s' is not in the database\\n\",name);	}\n\n"

    sb_gss_function *= "$(tab2)for (int j = 0; j < SS_ref_db.n_em; j++){\n"  
    sb_gss_function *= "$(tab3)SS_ref_db.mu_array[FD][j] = SS_ref_db.gbase[j];\n"  
    sb_gss_function *= "$(tab2)}\n"  
    sb_gss_function *= "$(tab)}\n\n"  

    sb_gss_function *= "$(tab)for (int j = 0; j < SS_ref_db.n_em; j++){\n"  
    sb_gss_function *= "$(tab2)SS_ref_db.bounds[j][0] = SS_ref_db.bounds_ref[j][0];\n"  
    sb_gss_function *= "$(tab2)SS_ref_db.bounds[j][1] = SS_ref_db.bounds_ref[j][1];\n"  
    sb_gss_function *= "$(tab)}\n\n"  

    sb_gss_function *= "$(tab)/* Calculate the number of atoms in the bulk-rock composition */\n"  
    sb_gss_function *= "$(tab)double fbc     = 0.0;\n"  
    sb_gss_function *= "$(tab)for (int i = 0; i < gv.len_ox; i++){\n"  
    sb_gss_function *= "$(tab2)fbc += z_b.bulk_rock[i]*z_b.apo[i];\n"  
    sb_gss_function *= "$(tab)}\n\n"  

    sb_gss_function *= "$(tab)for (int i = 0; i < SS_ref_db.n_em; i++){\n"  
    sb_gss_function *= "$(tab2)SS_ref_db.ape[i] = 0.0;\n"  
    sb_gss_function *= "$(tab2)for (int j = 0; j < gv.len_ox; j++){\n"  
    sb_gss_function *= "$(tab3)SS_ref_db.ape[i] += SS_ref_db.Comp[i][j]*z_b.apo[j];\n"  
    sb_gss_function *= "$(tab2)}\n"  
    sb_gss_function *= "$(tab)}\n\n"  

    sb_gss_function *= "$(tab)SS_ref_db.fbc = z_b.fbc;\n\n"  

    sb_gss_function *= "$(tab)if (gv.verbose == 1){\n"  
    sb_gss_function *= "$(tab2)printf(\" %4s:\",name);\n"  
    sb_gss_function *= "$(tab2)for (int j = 0; j < SS_ref_db.n_em; j++){\n"  
    sb_gss_function *= "$(tab3)printf(\" %+12.5f\",SS_ref_db.gbase[j]);\n"  
    sb_gss_function *= "$(tab2)}\n"  
    sb_gss_function *= "$(tab2)printf(\"\\n\");\n"

    if sb_ver == "sb11"
        sb_gss_function *= "$(tab2)printf(\" S   C   A   F   M   N\\n\");\n"  
    else
        println("database implemented yet...")
    end
    sb_gss_function *= "$(tab2)for (int i = 0; i < SS_ref_db.n_em; i++){\n"  
    sb_gss_function *= "$(tab3)for (int j = 0; j < gv.len_ox; j++){\n"  
    sb_gss_function *= "$(tab3)$(tab)printf(\" %.1f\",SS_ref_db.Comp[i][j]);\n"  
    sb_gss_function *= "$(tab3)}\n"  
    sb_gss_function *= "$(tab3)printf(\"\\n\");\n"  
    sb_gss_function *= "$(tab2)}\n"  
    sb_gss_function *= "$(tab2)printf(\"\\n\");\n"  
    sb_gss_function *= "$(tab)}\n"  
    
    sb_gss_function *= "\n$(tab)return SS_ref_db;\n"
    sb_gss_function *= "};\n"

    return sb_gss_function
end


function get_sb_SS_xeos_PC(sb_ver,ss)
    step    = 0.1
    eps     = 1e-4
    sb_SS_xeos_PC = ""

    n_ss = length(ss)
    for ii = 1:n_ss
        mul, site_cmp    = retrieve_site_cmp(ss, ii)
        n_em             = size(site_cmp)[3]
        # Generate a simplex in n dimensions
        simplex          = generate_simplex(n_em, step)

        # Adjust the simplex points
        adjusted_simplex = adjust_simplex(simplex, eps)
        n_pc             = length(adjusted_simplex)

        sb_SS_xeos_PC   *= "struct ss_pc $(sb_ver)_$(ss[ii].abbrev)_pc_xeos[$n_pc] = {\n"
        for i=1:n_pc
            sb_SS_xeos_PC   *= "$(tab){{"
            sb_SS_xeos_PC   *= join(round.(adjusted_simplex[i],digits=6),",")
            sb_SS_xeos_PC   *= "}},\n"
        end
        sb_SS_xeos_PC   *= "};\n"
    end

    sb_SS_xeos_PC *= "\n\n"
    sb_SS_xeos_PC *= "void SB_$(sb_ver)_pc_init_function(	PC_ref 	*SS_pc_xeos,\n"
    sb_SS_xeos_PC *= "$(tab3)$(tab3)$(tab3)int 	 iss,\n"
    sb_SS_xeos_PC *= "$(tab3)$(tab3)$(tab3)char 	*name				){\n"
                            
    # sb_SS_xeos_PC *= "$(tab)for (int iss = 0; iss < gv.len_ss; iss++){\n"
    for i=1:n_ss
        if i == 1
            sb_SS_xeos_PC *= "$(tab2)if      (strcmp( name, \"$(ss[i].abbrev)\")  == 0 ){\n"
            sb_SS_xeos_PC *= "$(tab3)SS_pc_xeos[iss].ss_pc_xeos   = $(sb_ver)_$(ss[i].abbrev)_pc_xeos; 		}\n"
        else
            sb_SS_xeos_PC *= "$(tab2)else if (strcmp( name, \"$(ss[i].abbrev)\")  == 0 ){\n"
            sb_SS_xeos_PC *= "$(tab3)SS_pc_xeos[iss].ss_pc_xeos   = $(sb_ver)_$(ss[i].abbrev)_pc_xeos; 		}\n"
        end
    end
     
    sb_SS_xeos_PC *= "$(tab2)else{\n"
    sb_SS_xeos_PC *= "$(tab3)printf(\"\\nsolid solution '%s' is not in the database, cannot be initiated\\n\", name);\n"
    sb_SS_xeos_PC *= "$(tab2)}	\n"
    # sb_SS_xeos_PC *= "$(tab)};\n"	
    sb_SS_xeos_PC *= "}\n\n\n"


    return sb_SS_xeos_PC
end


function get_SB_NLopt_opt_functions(sb_ver,ss)
    SB_NLopt_opt_functions = ""

    SB_NLopt_opt_functions   *= "\n"
    SB_NLopt_opt_functions   *= "// Equality constraint function: sum(x) == 1\n"
    SB_NLopt_opt_functions   *= "double equality_constraint(unsigned n, const double *x, double *grad, void *data) {\n"
    SB_NLopt_opt_functions   *= "$(tab2)if (grad) {\n"
    SB_NLopt_opt_functions   *= "$(tab3)for (unsigned i = 0; i < n; i++) {\n"
    SB_NLopt_opt_functions   *= "$(tab3)grad[i] = 1.0;\n"
    SB_NLopt_opt_functions   *= "$(tab2)}\n"
    SB_NLopt_opt_functions   *= "$(tab)}\n"
    SB_NLopt_opt_functions   *= "$(tab)double sum = 0.0;\n"
    SB_NLopt_opt_functions   *= "$(tab)for (unsigned i = 0; i < n; i++) {\n"
    SB_NLopt_opt_functions   *= "$(tab3)sum += x[i];\n"
    SB_NLopt_opt_functions   *= "$(tab)}\n"
    SB_NLopt_opt_functions   *= "$(tab)return sum - 1.0;\n"
    SB_NLopt_opt_functions   *= "}\n"
    SB_NLopt_opt_functions   *= "\n"
    n_ss = length(ss)
    for ii = 1:n_ss
        mul, site_cmp    = retrieve_site_cmp(ss, ii)
        n_em             = size(site_cmp)[3]

        SB_NLopt_opt_functions   *= "SS_ref NLopt_opt_$(sb_ver)_$(ss[ii].abbrev)_function(global_variable gv, SS_ref SS_ref_db){\n"
        SB_NLopt_opt_functions   *= "$(tab)unsigned int    n_em     = SS_ref_db.n_em;\n"
        SB_NLopt_opt_functions   *= "$(tab)double *x  = SS_ref_db.iguess;\n"
        
        SB_NLopt_opt_functions   *= "$(tab)for (int i = 0; i < (n_em); i++){\n"
        SB_NLopt_opt_functions   *= "$(tab2)SS_ref_db.lb[i] = SS_ref_db.bounds[i][0];\n"
        SB_NLopt_opt_functions   *= "$(tab2)SS_ref_db.ub[i] = SS_ref_db.bounds[i][1];\n"
        SB_NLopt_opt_functions   *= "$(tab)}\n"
        
        SB_NLopt_opt_functions   *= "$(tab)SS_ref_db.opt = nlopt_create(NLOPT_LD_SLSQP, (n_em)); \n"
        SB_NLopt_opt_functions   *= "$(tab)nlopt_set_lower_bounds(SS_ref_db.opt, SS_ref_db.lb);\n"
        SB_NLopt_opt_functions   *= "$(tab)nlopt_set_upper_bounds(SS_ref_db.opt, SS_ref_db.ub);\n"

        SB_NLopt_opt_functions   *= "$(tab)nlopt_set_min_objective(SS_ref_db.opt, obj_$(sb_ver)_$(ss[ii].abbrev), &SS_ref_db);\n"
        SB_NLopt_opt_functions   *= "$(tab)nlopt_add_equality_constraint(SS_ref_db.opt, equality_constraint, NULL, 1e-6);\n"
        SB_NLopt_opt_functions   *= "$(tab)nlopt_set_ftol_rel(SS_ref_db.opt, gv.obj_tol);\n"
        SB_NLopt_opt_functions   *= "$(tab)nlopt_set_maxeval(SS_ref_db.opt, gv.maxeval);\n"
        
        SB_NLopt_opt_functions   *= "$(tab)double minf;\n"
        SB_NLopt_opt_functions   *= "$(tab)if (gv.maxeval==1){  \n"
        SB_NLopt_opt_functions   *= "$(tab2)minf = obj_$(sb_ver)_$(ss[ii].abbrev)(n_em, x, NULL, &SS_ref_db);\n"
        SB_NLopt_opt_functions   *= "$(tab)}\n"
        SB_NLopt_opt_functions   *= "$(tab)else{\n"
        SB_NLopt_opt_functions   *= "$(tab2)SS_ref_db.status = nlopt_optimize(SS_ref_db.opt, x, &minf);\n"
        SB_NLopt_opt_functions   *= "$(tab)}\n"
        SB_NLopt_opt_functions   *= "$(tab)/* Send back needed local solution parameters */\n"
        SB_NLopt_opt_functions   *= "$(tab)for (int i = 0; i < SS_ref_db.n_xeos; i++){\n"
        SB_NLopt_opt_functions   *= "$(tab2)SS_ref_db.xeos[i] = x[i];\n"
        SB_NLopt_opt_functions   *= "$(tab)}\n\n"
        
        SB_NLopt_opt_functions   *= "$(tab)SS_ref_db.df   = minf;\n"
        SB_NLopt_opt_functions   *= "$(tab)nlopt_destroy(SS_ref_db.opt);\n\n"
        
        SB_NLopt_opt_functions   *= "$(tab)return SS_ref_db;\n"
        SB_NLopt_opt_functions   *= "};\n"
    end
    SB_NLopt_opt_functions   *= "\n\n"

    SB_NLopt_opt_functions *= "void SB_$(sb_ver)_NLopt_opt_init(	    NLopt_type 			*NLopt_opt,\n"
    SB_NLopt_opt_functions *= "$(tab3)$(tab3)global_variable 	 gv				){\n"      
    SB_NLopt_opt_functions *= "$(tab)for (int iss = 0; iss < gv.len_ss; iss++){\n"
    for i=1:n_ss
        if i == 1
            SB_NLopt_opt_functions *= "$(tab2)if      (strcmp( gv.SS_list[iss], \"$(ss[i].abbrev)\")  == 0 ){\n"
            SB_NLopt_opt_functions *= "$(tab3)NLopt_opt[iss]  = NLopt_opt_$(sb_ver)_$(ss[i].abbrev)_function; 		}\n"
        else
            SB_NLopt_opt_functions *= "$(tab2)else if (strcmp( gv.SS_list[iss], \"$(ss[i].abbrev)\")  == 0 ){\n"
            SB_NLopt_opt_functions *= "$(tab3)NLopt_opt[iss]  = NLopt_opt_$(sb_ver)_$(ss[i].abbrev)_function; 		}\n"
        end
    end
     
    SB_NLopt_opt_functions *= "$(tab2)else{\n"
    SB_NLopt_opt_functions *= "$(tab3)printf(\"\\nsolid solution '%s' is not in the database, cannot be initiated\\n\", gv.SS_list[iss]);\n"
    SB_NLopt_opt_functions *= "$(tab2)}	\n"
    SB_NLopt_opt_functions *= "$(tab)};\n"	
    SB_NLopt_opt_functions *= "}\n\n\n"


    return SB_NLopt_opt_functions
end

function generate_C_files(sb_ver,ss,data2)

    sb_gss_init_function    = get_sb_gss_init_function(sb_ver,ss)
    sb_gss_function         = get_sb_gss_function(sb_ver,ss,data2)
    sb_objective_functions  = get_sb_objective_functions(sb_ver,ss)
    sb_SS_xeos_PC           = get_sb_SS_xeos_PC(sb_ver,ss)
    SB_NLopt_opt_functions  = get_SB_NLopt_opt_functions(sb_ver,ss)

    return sb_gss_init_function, sb_gss_function, sb_objective_functions, sb_SS_xeos_PC, SB_NLopt_opt_functions
end

