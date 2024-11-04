# how to compute chear modulus:
# mu = mu0 + (P-Pr)*dmu/dP + (T-Tr)*dmu/dT
# where r is the reference pressure and temperature usually, 1bar and 298.15K

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

        println(mul, site_cmp, W, v, em)

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
    


function get_sb_gss_function(sb_ver,ss)
    sb_gss_function    = ""

    n_ss = length(ss)
    for ii = 1:n_ss

        mul, site_cmp = retrieve_site_cmp(ss, ii)

        W   = [ ss[ii].margules[j] for j in  keys(ss[ii].margules)]

        v   = [ ss[ii].van_laar[j] for j in  keys(ss[ii].van_laar)]

        em  = [ ss[ii].endmembers[j][1] for j in  keys(ss[ii].endmembers)]

        println(mul, site_cmp, W, v, em)

        sym = 1
        if ~isempty(v)
            sym = 0
        end
        
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

        for i=1:length(W)
            sb_gss_function *= "$(tab)SS_ref_db.W[$(i-1)] = $(W[i]);\n"
        end
        sb_gss_function *= "\n"
        if ~isempty(v)
            for i=1:length(v)
                sb_gss_function *= "$(tab)SS_ref_db.v[$(i-1)] = $(v[i]);\n"
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


function generate_C_files(sb_ver,ss)

    sb_gss_init_function    = get_sb_gss_init_function(sb_ver,ss)
    sb_gss_function         = get_sb_gss_function(sb_ver,ss)

    return sb_gss_init_function, sb_gss_function
end

