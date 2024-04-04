using MAGEMin_C         # load MAGEMin (needs to be loaded from main directory to pick up correct library in case it is locally compiled)

#inversion with Julia -> Optim (https://julianlsolvers.github.io/Optim.jl/stable/)

# Initialize database
# ig, igneous (HP18 -> Green et al., 2023)
# igd, igneous (T21 -> Green et al., 2023) 
# alk, alkaline (Weller et al., 2023)
# mp, metapelite (White et al 2014b) 
# mb, metabasite (Green et al., 2016)
# um, ultramafic (Evans & Frost, 2021)

db          = "mb"  
global gv, z_b, DB, splx_data      = init_MAGEMin(db);

sys_in      = "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
test        = 3         
global gv   = use_predefined_bulk_rock(gv, test, db);

P           = 4.0
T           = 1000.0
global gv.verbose  = -1        # switch off any verbose
# initialize MAGEMin up to G0 and W's point
global gv, z_b, DB, splx_data = pwm_init(P,T, gv, z_b, DB, splx_data);

# get names of the solution phases
ss_names  = unsafe_string.(unsafe_wrap(Vector{Ptr{Int8}}, gv.SS_list, gv.len_ss));

# get the solution phase structure (size gv.len_ss)
ss_struct = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);


#4: liq, 2: opx, 3: pl4tr, 13: aug, 1: sp; ilm, 6;
ss = 3; ssName = ss_names[ss];


ig = unsafe_wrap(Vector{Cdouble},ss_struct[ss].iguess, ss_struct[ss].n_xeos)
# tmp = [0.311787122    0.152332371    0.149755539    0.821682031   0.0229529747    0.391107788    0.141645662] #aug;
# tmp = [ 0.11479   0.33176   0.90264   0.10941   0.08312   0.12186   0.63121   0.07464]; #liq
# tmp = [0.114792609    0.331759108    0.902640189    0.109413663   0.0831170261    0.121861352    0.631204876   0.0746392014]; #liq2
tmp = [0.680373011  0.00109226965]; #pl4tr
# tmp = [0.338431366   0.0372603498   0.0428520434   0.0519428464    0.295506935]; #opx
# tmp = [0.621931926 0.000197970492] #ilm
# tmp = [0.837897697   0.0274371794    0.249717882]; #spn
for i =1:length(ig)
    ig[i] = tmp[i]
end
gb_lvl = unsafe_wrap(Vector{Cdouble},ss_struct[ss].gb_lvl, ss_struct[ss].n_em)
gbase = unsafe_wrap(Vector{Cdouble},ss_struct[ss].gbase, ss_struct[ss].n_em)

for i =1:length(gbase)
    gb_lvl[i] = gbase[i]
end
ss_struct[ss] = LibMAGEMin.PC_function(	gv,
                                        ss_struct[ss], 
                                        z_b,
                                        ssName 				);


# # run single point calculation
# out = pwm_run(gv, z_b, DB, splx_data);

# print(out.G_system," -> ",out.ph," ",out.ph_frac,"\n");

finalize_MAGEMin(gv,DB)
