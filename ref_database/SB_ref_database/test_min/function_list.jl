# NR 12-09-24 Script to test nullspace minimization with generic solution models 

struct Phase
    id::String                  # Name id
    abbrev::String                # Abbreviation
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

struct ModelJSON
    name            :: String
    abbrev          :: String
    endmembers      :: Dict{String, Vector{String}}
    margules        :: Dict{String, Float64}
    van_laar        :: Vector{Float64}
end

function objective(x::Vector, grad::Vector, N, f_config, f_grad_config, f_mu_Gex, active)
    mu_Gex          = f_mu_Gex(x);
    S_tot           = f_config(x);
    dSdp            = f_grad_config(x);

    if length(grad) > 0
        dGdp        = dSdp+mu_Gex;                          # raw dGibbs/dp
        grad       .= (N*(dGdp'*N)') .* active;             # projected dGibbs/dp to satisfy sum(p) = 1.0
    end
    G = mu_Gex'*x + S_tot
    # println("G: $G")
    # println("sum(x): $(sum(x))")
    return  G
end

function retrieve_site_cmp(ss, i)

    em          = String.(keys(ss[i].endmembers))
    n_em        = length(em)
    elems       = ["Si", "Ca", "Al", "Fe", "Mg", "Na"]
    n_el        = length(elems)

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


function get_functions(i)
    # load stixrude informations
    @load "STIX11.jld2" data2
    ss      = JSON3.read("stx11_solution.json", Vector{ModelJSON}) 
    mul, site_cmp = retrieve_site_cmp(ss, i)
    n_sf    = size(site_cmp)[1]
    n_ox    = size(site_cmp)[2]
    n_em    = size(site_cmp)[3]
    N       = nullspace(ones(n_em)')
    em_list = String.(keys(ss[i].endmembers))

    v       = []
    if ~isempty(ss[i].van_laar)
        v   = [ ss[i].van_laar[j] for j in  keys(ss[i].van_laar)]
    end

    W       = []
    if ~isempty(ss[i].margules)
        W   = [ ss[i].margules[j] for j in  keys(ss[i].margules)]
    end

    em_comp = zeros(Float64, n_em, n_ox)
    for j=1:n_em
        em          = ss[i].endmembers[em_list[j]][1]
        id          = findfirst(data2[!,:abbrev] .== em)
        em_comp[j,:]= [ (data2[id,:oxides][k]) for k in keys(data2[id,:oxides]) ]
    end

    sym             = isempty(v) ? 1 : 0

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

    @variables p[1:n_em]
    X               = Symbolics.scalarize(p)
    R               = 8.31446261815324
    T               = 1000.0

    Xo              = C*X
    config          = R * T * (M' * Diagonal(Xo) * log.(Xo))
    grad_config     = Symbolics.gradient(config, X)

    # for g in grad_config
    #     expr = string(g)
    #     expr = replace(expr, r"(\d)(?=[\[\(a-zA-Z])" => s"\1*")
    #     println(expr)
    # end


    mu_Gex          = get_mu_Gex(W, v, n_em, sym)

    f_config        = build_function(config,        X, expression = Val{false});
    f_grad_config,  = build_function(grad_config,   X, expression = Val{false});
    f_mu_Gex,       = build_function(mu_Gex,        X, expression = Val{false});


    return n_em, N, f_config, f_grad_config, f_mu_Gex
end