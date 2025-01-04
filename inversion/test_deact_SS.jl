# julia> include("inversion/inversion_functions.jl")
# get_str_id (generic function with 1 method)

# julia> display_ss_infos( "ig" )

# ---------------------------------
#     Database information (ig)
# ----------------------------------
#    spl:   1;  n_W's, 28;  n_em's, 8
#    bi:    2;  n_W's, 15;  n_em's, 6
#    cd:    3;  n_W's, 3;   n_em's, 3
#    cpx:   4;  n_W's, 45;  n_em's, 10
#    ep:    5;  n_W's, 3;   n_em's, 3
#    g:     6;  n_W's, 15;  n_em's, 6
#    hb:    7;  n_W's, 55;  n_em's, 11
#    ilm:   8;  n_W's, 10;  n_em's, 5
#    liq:   9;  n_W's, 66;  n_em's, 12
#    ol:    10; n_W's, 6;   n_em's, 4
#    opx:   11; n_W's, 36;  n_em's, 9
#    pl4T:  12; n_W's, 3;   n_em's, 3
#    fl:    13; n_W's, 55;  n_em's, 11
#    fper:  14; n_W's, 1;   n_em's, 2
# ----------------------------------

# FLAGS INFORMATION
# +-------+-------+------+------+-------+
# | SS/PP |  IN	  | CSD  | HLD  |  RMV  |
# +=======+=======+======+======+=======+
# | [0]   | 0/1   |  0/1 | 0/1  | 0/1   | 
# +-------+-------+------+------+-------+
# IN -> get_active_em
# CSD -> considered
# HLD -> on hold
# RMV -> removed

using MAGEMin_C


remove_spl  = 1


db          = "ig"
gv, z_b, DB, splx_data  = init_MAGEMin(db);
sys_in      =   "mol"     #default is mol, if wt is provided conversion will be done internally (MAGEMin works on mol basis)
test        =   0         #KLB1
gv          =   use_predefined_bulk_rock(gv, test, db);
P           = 10
T           = 1100

# initialize MAGEMin up to G0 and W's point
gv, z_b, DB, splx_data = pwm_init(P, T, gv, z_b, DB, splx_data);

# get the solution phase structure (size gv.len_ss)
ss_struct   = unsafe_wrap(Vector{LibMAGEMin.SS_ref},DB.SS_ref_db,gv.len_ss);

if remove_spl == 1

    ss          = 1                                                         # this is spl see top comments on how to get the id of the solution phase
    flags       = unsafe_wrap(Vector{Cint},ss_struct[ss].ss_flags, 4)       # this get the status flags of the solution phase "ss"

    flags[1]    = 0;
    flags[2]    = 0;
    flags[3]    = 0;
    flags[4]    = 1;

end

out       = pwm_run(gv, z_b, DB, splx_data);

finalize_MAGEMin(gv,DB)

####### no spl #########
# SOL = [G: -828.507] (66 iterations, 69.27 ms)
# GAM = [-1016.348130,-1832.956600,-821.622212,-697.748899,-415.435367,-980.427740,-882.088280,-1079.665799,-282.718293,-1382.360071]

# Phase :      opx     pl4T       ol      cpx 
# Mode  :  0.21835  0.01671  0.62506  0.13988 

####### spl #########
# SOL = [G: -828.524] (74 iterations, 70.11 ms)
# GAM = [-1015.584529,-1835.452257,-822.681912,-698.129777,-415.815198,-983.490765,-880.839294,-1080.444695,-283.531656,-1404.037065]

# Phase :       ol      spl      opx      cpx 
# Mode  :  0.59845  0.01761  0.22277  0.16117 