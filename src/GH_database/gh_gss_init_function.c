/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Moccetti, N. B., Dominguez, H., Assunção J., Green E., Dolejš, D., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#include <stdio.h>
#include <string.h>

#include "../MAGEMin.h"
#include "gh_gss_init_function.h"

/**
    Metadata for the Stage-A Ghiorso/MELTS liquid ("liq"), following the
    same formulation as "sb"'s solution phases (e.g. G_SS_sb11_ol_init_function):
    n_xeos == n_em (the NLopt vector IS the endmember mole-fraction vector
    directly, p[i] = x[i], no reduced/rotated basis), with the Sigma(x)=1
    closure enforced via a single NLopt equality constraint rather than
    algebraic substitution - this matches "sb"'s LP-only solving path,
    which is what "gh" also uses (see SetupDatabase).
    - 13 real basis endmembers (n_em = n_xeos = 13; SO3 and Cl2O-1 were
      dropped - no real reference-state data exists for them anywhere in
      xMELTS - and Ni/Co/P2O5/F were already excluded, see GH_endmembers.c).
    - n_sf = n_em: the "site fractions" here are simply the 13 mole
      fractions themselves (single mixing "site"), used only for the
      generic positivity check in ss_min_function.c's SS_UPDATE_function.
    - symmetric regular solution (Margules W only, no Van Laar), n_w = C(13,2) = 78.
    - n_cat = 0: no site-mixing bookkeeping (C[]/N[]) is populated - like
      "sb", these fields are metadata-only and never read at runtime.
*/
SS_ref G_SS_gh_liq_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 1;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 0;
    if (gv.EM_database == 2){
        /* pMELTS: 12 endmembers (CO2 dropped entirely - see
           GH_init_database.c's gh_db_pmelts_dataset), n_w = C(12,2) = 66. */
        SS_ref_db.n_xeos    = 12;
        SS_ref_db.n_em      = 12;
        SS_ref_db.n_sf      = 12;
        SS_ref_db.n_w       = 66;
    }
    else {
        SS_ref_db.n_xeos    = 13;
        SS_ref_db.n_em      = 13;
        SS_ref_db.n_sf      = 13;
        SS_ref_db.n_w       = 78;
    }

    return SS_ref_db;
}

/**
    Olivine (Fo-Fa), "sb-trivial" tier: direct p=x, symmetric regular
    solution (single W term, WH1MGFE=WH2MGFE=5075 J in xMELTS'
    sources/olivine.c, i.e. the Mg-Fe binary happens to be symmetric even
    though the full 6-component model is asymmetric/van-Laar elsewhere).
    Real xMELTS olivine also carries 4 internal M1/M2 Fe-Mg-Mn-Co-Ni
    order parameters (order-disorder machinery, "needs order-param
    machinery" tier); this Fo-Fa reduction deliberately drops them and
    uses the disordered-limit treatment instead, exactly like sb's own
    "ol" (2 equivalent M-sites, Sconfig multiplicity 2, no order
    parameter) - matching the "same structure as sb" scope the user chose.
*/
SS_ref G_SS_gh_ol_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 2;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

/**
    Feldspar (Ab-An-San), Elkins & Grove (1990) ternary asymmetric
    regular solution (xMELTS sources/feldspar.c) - single framework site,
    no order parameter (Al-Si and albite/sanidine order-disorder
    corrections that xMELTS itself applies in sources/gibbs.c are not
    reproduced, same simplification already made for the SiO2 polymorphs
    and for these same endmembers' own standard states, see
    GH_PP_endmembers.h). Direct p=x, n_xeos=n_em=3.
*/
SS_ref G_SS_gh_fsp_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 7;

    return SS_ref_db;
}

/**
    Biotite (Ann-Phl), Sack & Ghiorso binary regular solution (xMELTS
    sources/biotite.c) - single octahedral site, symmetric (one W term),
    no site-multiplicity factor on the ideal entropy in xMELTS' own
    calibration (ported as-is). Direct p=x, n_xeos=n_em=2.
*/
SS_ref G_SS_gh_bi_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 3;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

/**
    Garnet (Gr-Py-Alm), Berman & Koziol (1991) ternary asymmetric
    regular solution (xMELTS sources/garnet.c) - single dodecahedral
    X-site (multiplicity 3), no order parameter. Direct p=x, n_xeos=n_em=3.
*/
SS_ref G_SS_gh_g_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 3;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 7;

    return SS_ref_db;
}

/**
    Hornblende (Parg-Fparg-Mhst), from xMELTS sources/hornblende.c - NS=0
    (no internal order parameter) despite 2 real crystallographic sites
    (M12, M3); direct p=x, n_xeos=n_em=3 (see obj_gh_hb's header comment
    for why no separate ideal Sconfig is needed on top of the site-mixing
    entropy already folded into Gex there).
*/
SS_ref G_SS_gh_hb_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 2;

    return SS_ref_db;
}

/**
    Leucite (Lc-Anl-Nlc), from xMELTS sources/leucite.c - NS=0, but a
    genuine two-site (L, S) entropy model, not a plain 3-endmember
    Sigma(p)*log(p) (see obj_gh_lc's header comment). Direct p=x,
    n_xeos=n_em=3.
*/
SS_ref G_SS_gh_lc_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 3;
    SS_ref_db.n_em      = 3;
    SS_ref_db.n_sf      = 3;
    SS_ref_db.n_w       = 3;

    return SS_ref_db;
}

/**
    Melilite (Ak-Geh-Fak-Na), from xMELTS sources/melilite.c - NS=1 (Al/Si
    ordering within the gehlenite component). Direct p=x, n_xeos=n_em=4 -
    the internal order parameter is solved by an embedded 1D Newton
    iteration inside obj_gh_mel itself, not exposed as an NLopt dimension
    (see obj_gh_mel's header comment).
*/
SS_ref G_SS_gh_mel_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 4;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_sf      = 4;
    SS_ref_db.n_w       = 10;

    return SS_ref_db;
}

SS_ref G_SS_gh_cum_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

SS_ref G_SS_gh_spn_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 5;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

/** Rhombohedral oxide (Geikielite-Hematite-Ilmenite-Pyrophanite-Corundum):
    NR=4/NS=3 reciprocal solution with a joint 3-variable order-parameter
    solve - see obj_gh_rhm in gh_objective_functions.c. */
SS_ref G_SS_gh_rhm_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 5;
    SS_ref_db.n_em      = 5;
    SS_ref_db.n_sf      = 5;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

/** Nepheline (Sack & Ghiorso 1995): NR=3/NS=1 single embedded order
    parameter - see obj_gh_nph in gh_objective_functions.c. */
SS_ref G_SS_gh_nph_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 4;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_sf      = 4;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

/** Kalsilite (Sack & Ghiorso 1995): NR=3, no order parameter (direct
    algebraic solution) - see obj_gh_kls in gh_objective_functions.c. */
SS_ref G_SS_gh_kls_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 4;
    SS_ref_db.n_em      = 4;
    SS_ref_db.n_sf      = 4;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}

/** Clinopyroxene and orthopyroxene (Di-Cen-Hed-CaTs(Al)-CaTs(Fe3+)-Ess-Jd):
    NR=6/NA=7, gh's largest phase - see obj_gh_cpx in
    gh_objective_functions.c for why both share the same init shape. */
SS_ref G_SS_gh_cpx_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 7;
    SS_ref_db.n_em      = 7;
    SS_ref_db.n_sf      = 7;
    SS_ref_db.n_w       = 1;

    return SS_ref_db;
}
SS_ref G_SS_gh_opx_init_function(SS_ref SS_ref_db, global_variable gv){
    return G_SS_gh_cpx_init_function(SS_ref_db, gv);
}

/**
    Real mixed H2O-CO2 fluid, from xMELTS sources/fluid.c's fluidPhase()
    (see obj_gh_fluid's header comment) - a genuine real-gas EOS mixture,
    not a discrete-endmember Margules solution like every other gh phase
    above. Direct p=x, n_xeos=n_em=2, n_w=0 (no Margules W is used at all -
    the mixing physics comes entirely from GH_pitzer_sterner_mix_G).
*/
SS_ref G_SS_gh_fluid_init_function(SS_ref SS_ref_db, global_variable gv){
    SS_ref_db.is_liq    = 0;
    SS_ref_db.override  = 0;
    SS_ref_db.symmetry  = 1;
    SS_ref_db.n_cat     = 1;
    SS_ref_db.n_xeos    = 2;
    SS_ref_db.n_em      = 2;
    SS_ref_db.n_sf      = 2;
    SS_ref_db.n_w       = 0;

    return SS_ref_db;
}

void GH_SS_init(            SS_init_type        *SS_init,
                            global_variable      gv              ){

    for (int iss = 0; iss < gv.len_ss; iss++){
        if (strcmp( gv.SS_list[iss], "liq") == 0 ){
            SS_init[iss]  = G_SS_gh_liq_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "ol") == 0 ){
            SS_init[iss]  = G_SS_gh_ol_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "fsp") == 0 ){
            SS_init[iss]  = G_SS_gh_fsp_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "bi") == 0 ){
            SS_init[iss]  = G_SS_gh_bi_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "g") == 0 ){
            SS_init[iss]  = G_SS_gh_g_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "hb") == 0 ){
            SS_init[iss]  = G_SS_gh_hb_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "lc") == 0 ){
            SS_init[iss]  = G_SS_gh_lc_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "mel") == 0 ){
            SS_init[iss]  = G_SS_gh_mel_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "cum") == 0 ){
            SS_init[iss]  = G_SS_gh_cum_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "spn") == 0 ){
            SS_init[iss]  = G_SS_gh_spn_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "cpx") == 0 ){
            SS_init[iss]  = G_SS_gh_cpx_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "opx") == 0 ){
            SS_init[iss]  = G_SS_gh_opx_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "fl") == 0 ){
            SS_init[iss]  = G_SS_gh_fluid_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "rhm") == 0 ){
            SS_init[iss]  = G_SS_gh_rhm_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "nph") == 0 ){
            SS_init[iss]  = G_SS_gh_nph_init_function;
        }
        else if (strcmp( gv.SS_list[iss], "kls") == 0 ){
            SS_init[iss]  = G_SS_gh_kls_init_function;
        }
        else{
            printf("\nsolid solution '%s' is not in the 'gh' database, cannot be initiated\n", gv.SS_list[iss]);
        }
    }
}
