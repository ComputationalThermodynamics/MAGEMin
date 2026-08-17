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
#include "../initialize.h"
#include "gh_gss_function.h"
#include "GH_gem_function.h"
#include "GH_fluid_eos.h"

/**
    Boil-out test for a solution-phase endmember: true if any oxide it
    needs (any nonzero entry of its own Comp[] row, already fetched via
    get_em_data) is entirely absent from the bulk-rock composition.
    Generalizes TC_database's per-endmember boil-out checks (e.g.
    G_SS_ig_bi_function's "if (z_b.bulk_rock[8]==0.){...}") to endmembers
    that span more than one oxide axis at once (common in gh - e.g.
    Fe2O3=2FeO+O, KAlSiO4=SiO2+0.5Al2O3+0.5K2O - unlike tc's biotite,
    which only ever checks one diagnostic oxide per endmember).
*/
static int GH_boiled_out(int len_ox, const double *Comp, double *bulk_rock){
    for (int j = 0; j < len_ox; j++){
        if (Comp[j] != 0.0 && bulk_rock[j] == 0.0){
            return 1;
        }
    }
    return 0;
}

/**
    Reference-state assembly for "liq": each endmember's standard-state
    Gibbs energy is fetched by name (via get_em_data, dispatched through
    GH_G_EM_function), and the 78 = C(13,2) symmetric regular-solution
    Margules W terms (enthalpy only, J - every W^S/W^V term is zero in
    this xMELTS calibration) are set directly, following the same
    self-contained, one-function-per-phase structure as TC_database's
    tc_gss_function.c (e.g. G_SS_ig_liq_function): no shared helper, no
    file-scope tables - EM_list and W live entirely inside this function.
    W upper-triangle order i=0..11,l=i+1..12:
    0=SiO2 1=TiO2 2=Al2O3 3=Fe2O3 4=MgCr2O4 5=Fe2SiO4 6=MnSi0.5O2 7=Mg2SiO4
    8=CaSiO3 9=Na2SiO3 10=KAlSiO4 11=CO2 12=H2O
    Extracted from xMELTS-master/includes/param_struct_data.h
    (originalModelParameters[], default/non-new-water-model branch).
*/
SS_ref G_SS_gh_liq_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int EM_database, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"liq_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    /* pMELTS uses real pMELTS' own rescaled-formula-unit names for 4 of
       its 12 endmembers (Si4O8/Al4O6/Ca2Si2O6/NaSi0.5O1.5, matching
       arr_em_db_gh_liq_pMELTS in GH_endmembers.c exactly) - not cosmetic,
       a mole fraction of these IS a mole fraction of the larger formula
       unit, a real physical difference from xMELTS/rMELTS' plain-oxide
       endmembers of the same name. See [[gh-multicalibration-xmelts-rmelts-pmelts]]. */
    char *EM_tmp[] = {"SiO2","TiO2","Al2O3","Fe2O3","MgCr2O4","Fe2SiO4","MnSi0.5O2","Mg2SiO4",
                       "CaSiO3","Na2SiO3","KAlSiO4","CO2","H2O"};
    /* pMELTS (EM_database==2): no CO2 entry at all - real pMELTS has no
       CO2 standard-state treatment (see GH_endmembers.c's table header
       comment) - H2O sits at index 11, not 12. */
    char *EM_tmp_pMELTS[] = {"Si4O8","TiO2","Al4O6","Fe2O3","MgCr2O4","Fe2SiO4","MnSi0.5O2","Mg2SiO4",
                       "Ca2Si2O6","NaSi0.5O1.5","KAlSiO4","H2O"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], (EM_database == 2) ? EM_tmp_pMELTS[i] : EM_tmp[i]);
    };

    /* pMELTS (EM_database==2) W-table: CORRECTED 2026-07-14, see
       [[gh-gexcess-verification]]. The ORIGINAL reasoning here (kept for
       context, now known wrong) assumed pMELTS shares xMELTS' own
       "meltsModelParameters[]" table (sources/liquid_v34.c) and just
       deletes CO2's 12 touching indices - it does NOT. Real xMELTS-master
       has a genuinely SEPARATE, independently-calibrated
       "pMeltsModelParameters[]" table in the SAME
       param_struct_data_v34.h (dispatched via
       "calculationMode==MODE_pMELTS" inside gmixLiq_v34's own WH/WS/WV
       macros - `(calculationMode==MODE__MELTS) ? meltsModelParameters[...]
       : pMeltsModelParameters[...]`), with substantially different
       calibrated values throughout (e.g. Al4O6-related terms are in the
       -100000s to -400000s J range for pMELTS vs -20000 to -180000 for
       xMELTS' Al2O3 terms - not a rescaling, a different calibration).
       The table below was ALWAYS wrong before this fix - 65 of 66 values
       (all but one, MnSi0.5O2-H2O, coincidentally 0.0 in both tables)
       didn't match real pMeltsModelParameters[] at all; the previous
       values were effectively xMELTS' own table with CO2 deleted, exactly
       matching the (wrong) reasoning above. Extracted directly and
       mechanically from pMeltsModelParameters[]'s own self-labeled "W(A,B)"
       entries (param_struct_data_v34.h) via a small Python script
       (pmelts_check2.py) - not hand-transcribed. Real pMELTS genuinely
       has no CO2 endmember at all (unlike the mistaken "delete CO2 from
       xMELTS' table" framing above implied), so this table is directly
       66 terms for pMELTS' own 12 liquid components - no deletion/
       compaction step needed or performed. */
    if (EM_database == 2){
        SS_ref_db.W[0] = 15094.7;   SS_ref_db.W[1] = -296975.2;   SS_ref_db.W[2] = -164027.4;
        SS_ref_db.W[3] = 37459.2;   SS_ref_db.W[4] = -18841.4;   SS_ref_db.W[5] = 0.0;
        SS_ref_db.W[6] = -33833.5;   SS_ref_db.W[7] = -34232.9;   SS_ref_db.W[8] = -59822.7;
        SS_ref_db.W[9] = -102706.5;   SS_ref_db.W[10] = -45181.6;   SS_ref_db.W[11] = -144804.9;
        SS_ref_db.W[12] = -212292.3;   SS_ref_db.W[13] = -22455.8;   SS_ref_db.W[14] = 9324.2;
        SS_ref_db.W[15] = 0.0;   SS_ref_db.W[16] = 16335.6;   SS_ref_db.W[17] = -9471.5;
        SS_ref_db.W[18] = 22194.2;   SS_ref_db.W[19] = -3744.0;   SS_ref_db.W[20] = 70663.0;
        SS_ref_db.W[21] = -393566.0;   SS_ref_db.W[22] = -269339.7;   SS_ref_db.W[23] = -200788.1;
        SS_ref_db.W[24] = 0.0;   SS_ref_db.W[25] = -192709.0;   SS_ref_db.W[26] = -270700.8;
        SS_ref_db.W[27] = -205068.6;   SS_ref_db.W[28] = -114506.5;   SS_ref_db.W[29] = -161944.4;
        SS_ref_db.W[30] = 201536.3;   SS_ref_db.W[31] = -211493.4;   SS_ref_db.W[32] = 0.0;
        SS_ref_db.W[33] = -196914.9;   SS_ref_db.W[34] = -146008.1;   SS_ref_db.W[35] = -123728.7;
        SS_ref_db.W[36] = -130847.5;   SS_ref_db.W[37] = -114508.6;   SS_ref_db.W[38] = -74759.0;
        SS_ref_db.W[39] = 0.0;   SS_ref_db.W[40] = -3638.5;   SS_ref_db.W[41] = 48337.5;
        SS_ref_db.W[42] = -43302.5;   SS_ref_db.W[43] = 124517.4;   SS_ref_db.W[44] = -18.9;
        SS_ref_db.W[45] = 0.0;   SS_ref_db.W[46] = -28736.4;   SS_ref_db.W[47] = -28573.8;
        SS_ref_db.W[48] = -4723.9;   SS_ref_db.W[49] = 22245.0;   SS_ref_db.W[50] = 9769.4;
        SS_ref_db.W[51] = 0.0;   SS_ref_db.W[52] = 0.0;   SS_ref_db.W[53] = 0.0;
        SS_ref_db.W[54] = 0.0;   SS_ref_db.W[55] = 0.0;   SS_ref_db.W[56] = 574.1;
        SS_ref_db.W[57] = 9272.3;   SS_ref_db.W[58] = 36512.7;   SS_ref_db.W[59] = 24630.1;
        SS_ref_db.W[60] = 7430.3;   SS_ref_db.W[61] = 19927.4;   SS_ref_db.W[62] = -1583.7;
        SS_ref_db.W[63] = -1102.3;   SS_ref_db.W[64] = 13043.1;   SS_ref_db.W[65] = 35572.8;
    }
    /* rMELTS (EM_database==1): rhyolite-MELTS v1.2.0's own liquid model
       (sources/liquid_CO2_H2O.c, includes/param_struct_data_CO2_H2O.h
       "meltsAndCO2_H2OModelParameters[]"). Real rMELTS v1.2.0 has 19
       liquid species (13 shared with xMELTS + NiSi0.5O2/CoSi0.5O2/
       Ca3(PO4)2/SO3/Cl2O-1/F2O-1) - the 6 extra track Ni/Co/P/S/Cl/F,
       none of which exist anywhere in gh's/MAGEMin's 13-oxide bulk-rock
       system, so their mole fractions are always exactly 0 and every
       Margules cross-term touching them vanishes identically - safely
       excludable rather than needing to be ported (confirmed this is a
       real simplification, not a shortcut: extending MAGEMin's oxide set
       to add Ni/Co/P/S/Cl/F would be a much larger, separate undertaking
       out of scope here). Compared programmatically (not by hand,
       same practice as the earlier W-parameter diffs) against xMELTS'
       own liquid.c calibration for all 78 relevant W-pairs AND every
       endmember's H/S/V/Cp/EOS/Vliq/Kress/Tfus/Sfus/Cpl (meltsFluidLiquid[]
       vs xMeltsLiquid[] in liq_struct_data.h): every single endmember's
       standard-state data is byte-identical (so no new endmember table
       is needed in GH_endmembers.c - EM_database==1 already correctly
       falls through to the shared arr_em_db_gh_liq array), and of the 78
       W-values, only the 9 involving H2O differ (indices 11,22,32,41,
       56,67,71,74,76) - consistent with rMELTS v1.2.0 being described as
       xMELTS plus a "new H2O model", not a wholesale liquid recalibration.
       See [[gh-multicalibration-xmelts-rmelts-pmelts]]. */
    else if (EM_database == 1){
        SS_ref_db.W[0]  =  26266.7;   SS_ref_db.W[1]  = -39120.0;   SS_ref_db.W[2]  =   8110.3;
        SS_ref_db.W[3]  =  27886.3;   SS_ref_db.W[4]  =  23660.9;   SS_ref_db.W[5]  =  18393.9;
        SS_ref_db.W[6]  =   3421.0;   SS_ref_db.W[7]  =   -863.7;   SS_ref_db.W[8]  = -99039.0;
        SS_ref_db.W[9]  = -33921.7;   SS_ref_db.W[10] =      0.0;   SS_ref_db.W[11] =  27357.0;
        SS_ref_db.W[12] = -29449.8;   SS_ref_db.W[13] = -84756.9;   SS_ref_db.W[14] = -72303.4;
        SS_ref_db.W[15] =   5209.1;   SS_ref_db.W[16] = -16123.5;   SS_ref_db.W[17] =  -4178.3;
        SS_ref_db.W[18] = -35372.5;   SS_ref_db.W[19] = -15415.6;   SS_ref_db.W[20] = -48094.6;
        SS_ref_db.W[21] = -19265.7;   SS_ref_db.W[22] =  88199.2;   SS_ref_db.W[23] = -17089.4;
        SS_ref_db.W[24] = -31770.3;   SS_ref_db.W[25] = -30509.0;   SS_ref_db.W[26] = -53874.9;
        SS_ref_db.W[27] = -32880.3;   SS_ref_db.W[28] = -57917.9;   SS_ref_db.W[29] = -130785.0;
        SS_ref_db.W[30] = -25859.2;   SS_ref_db.W[31] =      0.0;   SS_ref_db.W[32] =  11768.4;
        SS_ref_db.W[33] =  21605.9;   SS_ref_db.W[34] = -179064.9;  SS_ref_db.W[35] =   3907.9;
        SS_ref_db.W[36] = -71518.6;   SS_ref_db.W[37] =  12076.6;   SS_ref_db.W[38] = -149662.2;
        SS_ref_db.W[39] =  57555.9;   SS_ref_db.W[40] =  -3186.9;   SS_ref_db.W[41] =  50105.2;
        SS_ref_db.W[42] = -82971.8;   SS_ref_db.W[43] =    182.4;   SS_ref_db.W[44] =  46049.2;
        SS_ref_db.W[45] =  30704.7;   SS_ref_db.W[46] = 113646.0;   SS_ref_db.W[47] =  75709.1;
        SS_ref_db.W[48] =      0.0;   SS_ref_db.W[49] =      0.0;   SS_ref_db.W[50] =  -6823.9;
        SS_ref_db.W[51] = -37256.7;   SS_ref_db.W[52] = -12970.8;   SS_ref_db.W[53] = -90533.8;
        SS_ref_db.W[54] =  23649.4;   SS_ref_db.W[55] = -32464.5;   SS_ref_db.W[56] =  30936.5;
        SS_ref_db.W[57] = -13040.1;   SS_ref_db.W[58] =   2934.6;   SS_ref_db.W[59] = -15780.8;
        SS_ref_db.W[60] =  23727.4;   SS_ref_db.W[61] =      0.0;   SS_ref_db.W[62] =      0.0;
        SS_ref_db.W[63] = -31731.9;   SS_ref_db.W[64] = -41876.9;   SS_ref_db.W[65] =  22323.1;
        SS_ref_db.W[66] = -40853.6;   SS_ref_db.W[67] =  20909.9;   SS_ref_db.W[68] = -13247.1;
        SS_ref_db.W[69] =  17111.1;   SS_ref_db.W[70] =  30012.5;   SS_ref_db.W[71] =   9715.0;
        SS_ref_db.W[72] =   6522.8;   SS_ref_db.W[73] =      0.0;   SS_ref_db.W[74] = -82460.1;
        SS_ref_db.W[75] =      0.0;   SS_ref_db.W[76] =   1057.2;   SS_ref_db.W[77] =  23255.4;
    }
    /* xMELTS (EM_database==0). CO2-touching entries (indices 10,21,31,40,
       49,55,62,66,70,73,75,76,77 in this 78-term enumeration) were found
       wrong 2026-07-14 via [[gh-gexcess-verification]]'s liq H2O/CO2-
       sensitive-diff investigation: 10 of 12 CO2-touching W-values didn't
       match real xMELTS' own param_struct_data_v34.h "meltsModelParameters[]"
       table at all (not a rounding slip - e.g. SiO2-CO2 was 0.0, should be
       66600.7; Fe2SiO4-CO2 was -32464.5, should be 134680.2) - likely a
       systematic error introduced somewhere in this table's original
       transcription/generation, isolated to the CO2 column specifically
       (every non-CO2 H2O-touching entry, and every other entry, matched
       exactly). Corrected directly against the source table (all other
       78 values cross-checked clean, only the entries below changed). */
    else {
        SS_ref_db.W[0]  =  26266.7;   SS_ref_db.W[1]  = -39120.0;   SS_ref_db.W[2]  =   8110.3;
        SS_ref_db.W[3]  =  27886.3;   SS_ref_db.W[4]  =  23660.9;   SS_ref_db.W[5]  =  18393.9;
        SS_ref_db.W[6]  =   3421.0;   SS_ref_db.W[7]  =   -863.7;   SS_ref_db.W[8]  = -99039.0;
        SS_ref_db.W[9]  = -33921.7;   SS_ref_db.W[10] =  66600.7;   SS_ref_db.W[11] =  30967.3;
        SS_ref_db.W[12] = -29449.8;   SS_ref_db.W[13] = -84756.9;   SS_ref_db.W[14] = -72303.4;
        SS_ref_db.W[15] =   5209.1;   SS_ref_db.W[16] = -16123.5;   SS_ref_db.W[17] =  -4178.3;
        SS_ref_db.W[18] = -35372.5;   SS_ref_db.W[19] = -15415.6;   SS_ref_db.W[20] = -48094.6;
        SS_ref_db.W[21] =  25427.3;   SS_ref_db.W[22] =  81879.1;   SS_ref_db.W[23] = -17089.4;
        SS_ref_db.W[24] = -31770.3;   SS_ref_db.W[25] = -30509.0;   SS_ref_db.W[26] = -53874.9;
        SS_ref_db.W[27] = -32880.3;   SS_ref_db.W[28] = -57917.9;   SS_ref_db.W[29] = -130785.0;
        SS_ref_db.W[30] = -25859.2;   SS_ref_db.W[31] =  27043.7;   SS_ref_db.W[32] = -16098.1;
        SS_ref_db.W[33] =  21605.9;   SS_ref_db.W[34] = -179064.9;  SS_ref_db.W[35] =   3907.9;
        SS_ref_db.W[36] = -71518.6;   SS_ref_db.W[37] =  12076.6;   SS_ref_db.W[38] = -149662.2;
        SS_ref_db.W[39] =  57555.9;   SS_ref_db.W[40] =   1688.8;   SS_ref_db.W[41] =  31405.5;
        SS_ref_db.W[42] = -82971.8;   SS_ref_db.W[43] =    182.4;   SS_ref_db.W[44] =  46049.2;
        SS_ref_db.W[45] =  30704.7;   SS_ref_db.W[46] = 113646.0;   SS_ref_db.W[47] =  75709.1;
        SS_ref_db.W[48] =      0.0;   SS_ref_db.W[49] =      0.0;   SS_ref_db.W[50] =  -6823.9;
        SS_ref_db.W[51] = -37256.7;   SS_ref_db.W[52] = -12970.8;   SS_ref_db.W[53] = -90533.8;
        SS_ref_db.W[54] =  23649.4;   SS_ref_db.W[55] = 134680.2;   SS_ref_db.W[56] =  28873.6;
        SS_ref_db.W[57] = -13040.1;   SS_ref_db.W[58] =   2934.6;   SS_ref_db.W[59] = -15780.8;
        SS_ref_db.W[60] =  23727.4;   SS_ref_db.W[61] =      0.0;   SS_ref_db.W[62] =      0.0;
        SS_ref_db.W[63] = -31731.9;   SS_ref_db.W[64] = -41876.9;   SS_ref_db.W[65] =  22323.1;
        SS_ref_db.W[66] =  47406.1;   SS_ref_db.W[67] =  35633.7;   SS_ref_db.W[68] = -13247.1;
        SS_ref_db.W[69] =  17111.1;   SS_ref_db.W[70] =  44513.8;   SS_ref_db.W[71] =  20374.6;
        SS_ref_db.W[72] =   6522.8;   SS_ref_db.W[73] =   9846.7;   SS_ref_db.W[74] = -96937.6;
        SS_ref_db.W[75] =   1744.4;   SS_ref_db.W[76] =  10374.2;   SS_ref_db.W[77] =  43128.4;
    }

    /* pMELTS (EM_database==2) has no CO2 endmember at all (see
       GH_endmembers.c's table) - fetching "CO2" via get_em_data would hit
       find_EM_id()'s -1 "not found" return and then an out-of-bounds
       Access_GH_EM_DB(2,-1) read, since CO2 isn't in this dataset's
       endmember hashtable (gv.n_em_db==12). Branch early instead: 12
       fetches, no CO2_eq, H2O shifted from index 12 to index 11. Real
       xMELTS/rMELTS (EM_database 0/1) keep the original 13-endmember,
       CO2-inclusive path below unchanged. */
    if (EM_database == 2){
        em_data SiO2_eq    = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Si4O8",       "equilibrium");
        em_data TiO2_eq    = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "TiO2",        "equilibrium");
        em_data Al2O3_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Al4O6",       "equilibrium");
        em_data Fe2O3_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Fe2O3",       "equilibrium");
        em_data MgCr2O4_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "MgCr2O4",     "equilibrium");
        em_data Fe2SiO4_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Fe2SiO4",     "equilibrium");
        em_data MnSi_eq    = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "MnSi0.5O2",   "equilibrium");
        em_data Mg2SiO4_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Mg2SiO4",     "equilibrium");
        em_data CaSiO3_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Ca2Si2O6",    "equilibrium");
        em_data Na2SiO3_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "NaSi0.5O1.5", "equilibrium");
        em_data KAlSiO4_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "KAlSiO4",     "equilibrium");
        GH_H2O_liquid_context = 1;
        em_data H2O_eq     = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "H2O",         "equilibrium");
        GH_H2O_liquid_context = 0;

        SS_ref_db.gbase[0]  = SiO2_eq.gb;
        SS_ref_db.gbase[1]  = TiO2_eq.gb;
        SS_ref_db.gbase[2]  = Al2O3_eq.gb;
        SS_ref_db.gbase[3]  = Fe2O3_eq.gb;
        SS_ref_db.gbase[4]  = MgCr2O4_eq.gb;
        SS_ref_db.gbase[5]  = Fe2SiO4_eq.gb;
        SS_ref_db.gbase[6]  = MnSi_eq.gb;
        SS_ref_db.gbase[7]  = Mg2SiO4_eq.gb;
        SS_ref_db.gbase[8]  = CaSiO3_eq.gb;
        SS_ref_db.gbase[9]  = Na2SiO3_eq.gb;
        SS_ref_db.gbase[10] = KAlSiO4_eq.gb;
        SS_ref_db.gbase[11] = H2O_eq.gb;

        SS_ref_db.ElShearMod[0]  = SiO2_eq.ElShearMod;
        SS_ref_db.ElShearMod[1]  = TiO2_eq.ElShearMod;
        SS_ref_db.ElShearMod[2]  = Al2O3_eq.ElShearMod;
        SS_ref_db.ElShearMod[3]  = Fe2O3_eq.ElShearMod;
        SS_ref_db.ElShearMod[4]  = MgCr2O4_eq.ElShearMod;
        SS_ref_db.ElShearMod[5]  = Fe2SiO4_eq.ElShearMod;
        SS_ref_db.ElShearMod[6]  = MnSi_eq.ElShearMod;
        SS_ref_db.ElShearMod[7]  = Mg2SiO4_eq.ElShearMod;
        SS_ref_db.ElShearMod[8]  = CaSiO3_eq.ElShearMod;
        SS_ref_db.ElShearMod[9]  = Na2SiO3_eq.ElShearMod;
        SS_ref_db.ElShearMod[10] = KAlSiO4_eq.ElShearMod;
        SS_ref_db.ElShearMod[11] = H2O_eq.ElShearMod;

        SS_ref_db.ElBulkMod[0]  = SiO2_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[1]  = TiO2_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[2]  = Al2O3_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[3]  = Fe2O3_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[4]  = MgCr2O4_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[5]  = Fe2SiO4_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[6]  = MnSi_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[7]  = Mg2SiO4_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[8]  = CaSiO3_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[9]  = Na2SiO3_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[10] = KAlSiO4_eq.ElBulkMod;
        SS_ref_db.ElBulkMod[11] = H2O_eq.ElBulkMod;

        SS_ref_db.ElCp[0]  = SiO2_eq.ElCp;
        SS_ref_db.ElCp[1]  = TiO2_eq.ElCp;
        SS_ref_db.ElCp[2]  = Al2O3_eq.ElCp;
        SS_ref_db.ElCp[3]  = Fe2O3_eq.ElCp;
        SS_ref_db.ElCp[4]  = MgCr2O4_eq.ElCp;
        SS_ref_db.ElCp[5]  = Fe2SiO4_eq.ElCp;
        SS_ref_db.ElCp[6]  = MnSi_eq.ElCp;
        SS_ref_db.ElCp[7]  = Mg2SiO4_eq.ElCp;
        SS_ref_db.ElCp[8]  = CaSiO3_eq.ElCp;
        SS_ref_db.ElCp[9]  = Na2SiO3_eq.ElCp;
        SS_ref_db.ElCp[10] = KAlSiO4_eq.ElCp;
        SS_ref_db.ElCp[11] = H2O_eq.ElCp;

        SS_ref_db.ElExpansivity[0]  = SiO2_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[1]  = TiO2_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[2]  = Al2O3_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[3]  = Fe2O3_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[4]  = MgCr2O4_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[5]  = Fe2SiO4_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[6]  = MnSi_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[7]  = Mg2SiO4_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[8]  = CaSiO3_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[9]  = Na2SiO3_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[10] = KAlSiO4_eq.ElExpansivity;
        SS_ref_db.ElExpansivity[11] = H2O_eq.ElExpansivity;

        for (i = 0; i < len_ox; i++){
            SS_ref_db.Comp[0][i]  = SiO2_eq.C[i];
            SS_ref_db.Comp[1][i]  = TiO2_eq.C[i];
            SS_ref_db.Comp[2][i]  = Al2O3_eq.C[i];
            SS_ref_db.Comp[3][i]  = Fe2O3_eq.C[i];
            SS_ref_db.Comp[4][i]  = MgCr2O4_eq.C[i];
            SS_ref_db.Comp[5][i]  = Fe2SiO4_eq.C[i];
            SS_ref_db.Comp[6][i]  = MnSi_eq.C[i];
            SS_ref_db.Comp[7][i]  = Mg2SiO4_eq.C[i];
            SS_ref_db.Comp[8][i]  = CaSiO3_eq.C[i];
            SS_ref_db.Comp[9][i]  = Na2SiO3_eq.C[i];
            SS_ref_db.Comp[10][i] = KAlSiO4_eq.C[i];
            SS_ref_db.Comp[11][i] = H2O_eq.C[i];
        }

        for (i = 0; i < n_em; i++){ SS_ref_db.z_em[i] = 1.0; };

        for (i = 0; i < SS_ref_db.n_xeos; i++){
            SS_ref_db.bounds_ref[i][0] = 0.0+eps;
            SS_ref_db.bounds_ref[i][1] = 1.0-eps;
        }

        for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
        if (GH_boiled_out(len_ox, SiO2_eq.C,    z_b.bulk_rock)){ SS_ref_db.d_em[0]  = 1.0; SS_ref_db.z_em[0]  = 0.0; SS_ref_db.bounds_ref[0][0]  = 0.0; SS_ref_db.bounds_ref[0][1]  = 0.0; }
        if (GH_boiled_out(len_ox, TiO2_eq.C,    z_b.bulk_rock)){ SS_ref_db.d_em[1]  = 1.0; SS_ref_db.z_em[1]  = 0.0; SS_ref_db.bounds_ref[1][0]  = 0.0; SS_ref_db.bounds_ref[1][1]  = 0.0; }
        if (GH_boiled_out(len_ox, Al2O3_eq.C,   z_b.bulk_rock)){ SS_ref_db.d_em[2]  = 1.0; SS_ref_db.z_em[2]  = 0.0; SS_ref_db.bounds_ref[2][0]  = 0.0; SS_ref_db.bounds_ref[2][1]  = 0.0; }
        if (GH_boiled_out(len_ox, Fe2O3_eq.C,   z_b.bulk_rock)){ SS_ref_db.d_em[3]  = 1.0; SS_ref_db.z_em[3]  = 0.0; SS_ref_db.bounds_ref[3][0]  = 0.0; SS_ref_db.bounds_ref[3][1]  = 0.0; }
        if (GH_boiled_out(len_ox, MgCr2O4_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[4]  = 1.0; SS_ref_db.z_em[4]  = 0.0; SS_ref_db.bounds_ref[4][0]  = 0.0; SS_ref_db.bounds_ref[4][1]  = 0.0; }
        if (GH_boiled_out(len_ox, Fe2SiO4_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[5]  = 1.0; SS_ref_db.z_em[5]  = 0.0; SS_ref_db.bounds_ref[5][0]  = 0.0; SS_ref_db.bounds_ref[5][1]  = 0.0; }
        if (GH_boiled_out(len_ox, MnSi_eq.C,    z_b.bulk_rock)){ SS_ref_db.d_em[6]  = 1.0; SS_ref_db.z_em[6]  = 0.0; SS_ref_db.bounds_ref[6][0]  = 0.0; SS_ref_db.bounds_ref[6][1]  = 0.0; }
        if (GH_boiled_out(len_ox, Mg2SiO4_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[7]  = 1.0; SS_ref_db.z_em[7]  = 0.0; SS_ref_db.bounds_ref[7][0]  = 0.0; SS_ref_db.bounds_ref[7][1]  = 0.0; }
        if (GH_boiled_out(len_ox, CaSiO3_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[8]  = 1.0; SS_ref_db.z_em[8]  = 0.0; SS_ref_db.bounds_ref[8][0]  = 0.0; SS_ref_db.bounds_ref[8][1]  = 0.0; }
        if (GH_boiled_out(len_ox, Na2SiO3_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[9]  = 1.0; SS_ref_db.z_em[9]  = 0.0; SS_ref_db.bounds_ref[9][0]  = 0.0; SS_ref_db.bounds_ref[9][1]  = 0.0; }
        if (GH_boiled_out(len_ox, KAlSiO4_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[10] = 1.0; SS_ref_db.z_em[10] = 0.0; SS_ref_db.bounds_ref[10][0] = 0.0; SS_ref_db.bounds_ref[10][1] = 0.0; }
        if (GH_boiled_out(len_ox, H2O_eq.C,     z_b.bulk_rock)){ SS_ref_db.d_em[11] = 1.0; SS_ref_db.z_em[11] = 0.0; SS_ref_db.bounds_ref[11][0] = 0.0; SS_ref_db.bounds_ref[11][1] = 0.0; }

        return SS_ref_db;
    }

    /* xMELTS/rMELTS (EM_database 0/1): 13 endmembers including CO2. */
    em_data SiO2_eq    = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "SiO2",      "equilibrium");
    em_data TiO2_eq    = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "TiO2",      "equilibrium");
    em_data Al2O3_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Al2O3",     "equilibrium");
    em_data Fe2O3_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Fe2O3",     "equilibrium");
    em_data MgCr2O4_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "MgCr2O4",   "equilibrium");
    em_data Fe2SiO4_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Fe2SiO4",   "equilibrium");
    em_data MnSi_eq    = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "MnSi0.5O2", "equilibrium");
    em_data Mg2SiO4_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Mg2SiO4",   "equilibrium");
    em_data CaSiO3_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "CaSiO3",    "equilibrium");
    em_data Na2SiO3_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "Na2SiO3",   "equilibrium");
    em_data KAlSiO4_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "KAlSiO4",   "equilibrium");
    em_data CO2_eq     = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "CO2",       "equilibrium");
    GH_H2O_liquid_context = 1;   /* liquid's own "H2O" basis species, not the standalone "water" phase - see GH_gem_function.h */
    em_data H2O_eq     = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "H2O",       "equilibrium");
    GH_H2O_liquid_context = 0;

    SS_ref_db.gbase[0]  = SiO2_eq.gb;
    SS_ref_db.gbase[1]  = TiO2_eq.gb;
    SS_ref_db.gbase[2]  = Al2O3_eq.gb;
    SS_ref_db.gbase[3]  = Fe2O3_eq.gb;
    SS_ref_db.gbase[4]  = MgCr2O4_eq.gb;
    SS_ref_db.gbase[5]  = Fe2SiO4_eq.gb;
    SS_ref_db.gbase[6]  = MnSi_eq.gb;
    SS_ref_db.gbase[7]  = Mg2SiO4_eq.gb;
    SS_ref_db.gbase[8]  = CaSiO3_eq.gb;
    SS_ref_db.gbase[9]  = Na2SiO3_eq.gb;
    SS_ref_db.gbase[10] = KAlSiO4_eq.gb;
    SS_ref_db.gbase[11] = CO2_eq.gb;
    SS_ref_db.gbase[12] = H2O_eq.gb;

    SS_ref_db.ElShearMod[0]  = SiO2_eq.ElShearMod;
    SS_ref_db.ElShearMod[1]  = TiO2_eq.ElShearMod;
    SS_ref_db.ElShearMod[2]  = Al2O3_eq.ElShearMod;
    SS_ref_db.ElShearMod[3]  = Fe2O3_eq.ElShearMod;
    SS_ref_db.ElShearMod[4]  = MgCr2O4_eq.ElShearMod;
    SS_ref_db.ElShearMod[5]  = Fe2SiO4_eq.ElShearMod;
    SS_ref_db.ElShearMod[6]  = MnSi_eq.ElShearMod;
    SS_ref_db.ElShearMod[7]  = Mg2SiO4_eq.ElShearMod;
    SS_ref_db.ElShearMod[8]  = CaSiO3_eq.ElShearMod;
    SS_ref_db.ElShearMod[9]  = Na2SiO3_eq.ElShearMod;
    SS_ref_db.ElShearMod[10] = KAlSiO4_eq.ElShearMod;
    SS_ref_db.ElShearMod[11] = CO2_eq.ElShearMod;
    SS_ref_db.ElShearMod[12] = H2O_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0]  = SiO2_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1]  = TiO2_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2]  = Al2O3_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[3]  = Fe2O3_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[4]  = MgCr2O4_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[5]  = Fe2SiO4_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[6]  = MnSi_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[7]  = Mg2SiO4_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[8]  = CaSiO3_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[9]  = Na2SiO3_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[10] = KAlSiO4_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[11] = CO2_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[12] = H2O_eq.ElBulkMod;

    SS_ref_db.ElCp[0]  = SiO2_eq.ElCp;
    SS_ref_db.ElCp[1]  = TiO2_eq.ElCp;
    SS_ref_db.ElCp[2]  = Al2O3_eq.ElCp;
    SS_ref_db.ElCp[3]  = Fe2O3_eq.ElCp;
    SS_ref_db.ElCp[4]  = MgCr2O4_eq.ElCp;
    SS_ref_db.ElCp[5]  = Fe2SiO4_eq.ElCp;
    SS_ref_db.ElCp[6]  = MnSi_eq.ElCp;
    SS_ref_db.ElCp[7]  = Mg2SiO4_eq.ElCp;
    SS_ref_db.ElCp[8]  = CaSiO3_eq.ElCp;
    SS_ref_db.ElCp[9]  = Na2SiO3_eq.ElCp;
    SS_ref_db.ElCp[10] = KAlSiO4_eq.ElCp;
    SS_ref_db.ElCp[11] = CO2_eq.ElCp;
    SS_ref_db.ElCp[12] = H2O_eq.ElCp;

    SS_ref_db.ElExpansivity[0]  = SiO2_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1]  = TiO2_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2]  = Al2O3_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[3]  = Fe2O3_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[4]  = MgCr2O4_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[5]  = Fe2SiO4_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[6]  = MnSi_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[7]  = Mg2SiO4_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[8]  = CaSiO3_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[9]  = Na2SiO3_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[10] = KAlSiO4_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[11] = CO2_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[12] = H2O_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i]  = SiO2_eq.C[i];
        SS_ref_db.Comp[1][i]  = TiO2_eq.C[i];
        SS_ref_db.Comp[2][i]  = Al2O3_eq.C[i];
        SS_ref_db.Comp[3][i]  = Fe2O3_eq.C[i];
        SS_ref_db.Comp[4][i]  = MgCr2O4_eq.C[i];
        SS_ref_db.Comp[5][i]  = Fe2SiO4_eq.C[i];
        SS_ref_db.Comp[6][i]  = MnSi_eq.C[i];
        SS_ref_db.Comp[7][i]  = Mg2SiO4_eq.C[i];
        SS_ref_db.Comp[8][i]  = CaSiO3_eq.C[i];
        SS_ref_db.Comp[9][i]  = Na2SiO3_eq.C[i];
        SS_ref_db.Comp[10][i] = KAlSiO4_eq.C[i];
        SS_ref_db.Comp[11][i] = CO2_eq.C[i];
        SS_ref_db.Comp[12][i] = H2O_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    /* boil out any endmember whose composition needs an oxide entirely
       absent from the bulk-rock (see GH_boiled_out header comment) */
    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, SiO2_eq.C,    z_b.bulk_rock)){ SS_ref_db.d_em[0]  = 1.0; SS_ref_db.z_em[0]  = 0.0; SS_ref_db.bounds_ref[0][0]  = 0.0; SS_ref_db.bounds_ref[0][1]  = 0.0; }
    if (GH_boiled_out(len_ox, TiO2_eq.C,    z_b.bulk_rock)){ SS_ref_db.d_em[1]  = 1.0; SS_ref_db.z_em[1]  = 0.0; SS_ref_db.bounds_ref[1][0]  = 0.0; SS_ref_db.bounds_ref[1][1]  = 0.0; }
    if (GH_boiled_out(len_ox, Al2O3_eq.C,   z_b.bulk_rock)){ SS_ref_db.d_em[2]  = 1.0; SS_ref_db.z_em[2]  = 0.0; SS_ref_db.bounds_ref[2][0]  = 0.0; SS_ref_db.bounds_ref[2][1]  = 0.0; }
    if (GH_boiled_out(len_ox, Fe2O3_eq.C,   z_b.bulk_rock)){ SS_ref_db.d_em[3]  = 1.0; SS_ref_db.z_em[3]  = 0.0; SS_ref_db.bounds_ref[3][0]  = 0.0; SS_ref_db.bounds_ref[3][1]  = 0.0; }
    if (GH_boiled_out(len_ox, MgCr2O4_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[4]  = 1.0; SS_ref_db.z_em[4]  = 0.0; SS_ref_db.bounds_ref[4][0]  = 0.0; SS_ref_db.bounds_ref[4][1]  = 0.0; }
    if (GH_boiled_out(len_ox, Fe2SiO4_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[5]  = 1.0; SS_ref_db.z_em[5]  = 0.0; SS_ref_db.bounds_ref[5][0]  = 0.0; SS_ref_db.bounds_ref[5][1]  = 0.0; }
    if (GH_boiled_out(len_ox, MnSi_eq.C,    z_b.bulk_rock)){ SS_ref_db.d_em[6]  = 1.0; SS_ref_db.z_em[6]  = 0.0; SS_ref_db.bounds_ref[6][0]  = 0.0; SS_ref_db.bounds_ref[6][1]  = 0.0; }
    if (GH_boiled_out(len_ox, Mg2SiO4_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[7]  = 1.0; SS_ref_db.z_em[7]  = 0.0; SS_ref_db.bounds_ref[7][0]  = 0.0; SS_ref_db.bounds_ref[7][1]  = 0.0; }
    if (GH_boiled_out(len_ox, CaSiO3_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[8]  = 1.0; SS_ref_db.z_em[8]  = 0.0; SS_ref_db.bounds_ref[8][0]  = 0.0; SS_ref_db.bounds_ref[8][1]  = 0.0; }
    if (GH_boiled_out(len_ox, Na2SiO3_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[9]  = 1.0; SS_ref_db.z_em[9]  = 0.0; SS_ref_db.bounds_ref[9][0]  = 0.0; SS_ref_db.bounds_ref[9][1]  = 0.0; }
    if (GH_boiled_out(len_ox, KAlSiO4_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[10] = 1.0; SS_ref_db.z_em[10] = 0.0; SS_ref_db.bounds_ref[10][0] = 0.0; SS_ref_db.bounds_ref[10][1] = 0.0; }
    if (GH_boiled_out(len_ox, CO2_eq.C,     z_b.bulk_rock)){ SS_ref_db.d_em[11] = 1.0; SS_ref_db.z_em[11] = 0.0; SS_ref_db.bounds_ref[11][0] = 0.0; SS_ref_db.bounds_ref[11][1] = 0.0; }
    if (GH_boiled_out(len_ox, H2O_eq.C,     z_b.bulk_rock)){ SS_ref_db.d_em[12] = 1.0; SS_ref_db.z_em[12] = 0.0; SS_ref_db.bounds_ref[12][0] = 0.0; SS_ref_db.bounds_ref[12][1] = 0.0; }

    return SS_ref_db;
}

/**
    Olivine (Fo-Fa): W[0] is the real pressure-dependent effective
    regular-solution energy W(P) = 4*C0(P) = 20300+0.015*(Pbar-1) J,
    derived by reducing xMELTS' full Taylor-expansion olivine model to
    the Fo-Fa binary (see obj_gh_ol's header comment for the derivation).
    It is NOT the raw WH1MGFE/WH2MGFE=5075 J M1/M2 exchange-energy
    parameter from xMELTS-master/sources/olivine.c. Computed here (like
    TC_database's own P-dependent W's, e.g. G_SS_ig_liq_function's
    W[0]=9.2-0.08*SS_ref_db.P) so obj_gh_ol can just read d->W[0].
*/
SS_ref G_SS_gh_ol_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"ol_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"fo","fa"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 20300.0 + 0.015*(SS_ref_db.P*1000.0 - 1.0);

    em_data fo_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "fo", "equilibrium");
    em_data fa_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "fa", "equilibrium");

    SS_ref_db.gbase[0] = fo_eq.gb;
    SS_ref_db.gbase[1] = fa_eq.gb;

    SS_ref_db.ElShearMod[0] = fo_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = fa_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = fo_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = fa_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = fo_eq.ElCp;
    SS_ref_db.ElCp[1] = fa_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = fo_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = fa_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = fo_eq.C[i];
        SS_ref_db.Comp[1][i] = fa_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, fo_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, fa_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }

    return SS_ref_db;
}

/**
    Real mixed H2O-CO2 fluid: Duan & Zhang (2006), ported 2026-07-17 from
    real fluidPhase.c - the ACTUAL rhyolite-MELTS "fluid" solution phase
    (sol_struct_data.h's gmixFlu/actFlu function pointers + "h2oduan"/
    "co2duan" components, standard states = propertiesOfPureH2O/CO2).
    gbase[0]/gbase[1] = the DZ2006 pure-fluid standard states
    (GH_duan_pure_G, byte-exact vs real propertiesOfPure*), so the LP/NLopt
    hyperplane rotation applies to them exactly like every other gh phase;
    the mixing physics (fugacity-based, near-ideal - a genuine EOS
    cross-term, not a Margules Gex) lives in obj_gh_fluid via
    GH_duan_mix_muGex - so, unlike every other gh solution phase, W[] is
    unused here (n_w=0). The previous Pitzer-Sterner version (fluid.c's
    fluidPhase, crude linear-c mixing) was a porting mistake - real MELTS
    only uses that model for the LIQUID's H2O/CO2 standard states. The
    em_data fetches below are kept solely for Comp[]/ElMods/boil-out. */
SS_ref G_SS_gh_fluid_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"fluid_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"H2O","CO2"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    GH_H2O_liquid_context = 1;
    em_data h2o_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "H2O", "equilibrium");
    GH_H2O_liquid_context = 0;
    em_data co2_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "CO2", "equilibrium");

    SS_ref_db.gbase[0] = GH_duan_pure_G(1, SS_ref_db.T, SS_ref_db.P*1000.0)/1000.0;   /* J -> kJ */
    SS_ref_db.gbase[1] = GH_duan_pure_G(0, SS_ref_db.T, SS_ref_db.P*1000.0)/1000.0;

    SS_ref_db.ElShearMod[0] = h2o_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = co2_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = h2o_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = co2_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = h2o_eq.ElCp;
    SS_ref_db.ElCp[1] = co2_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = h2o_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = co2_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = h2o_eq.C[i];
        SS_ref_db.Comp[1][i] = co2_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, h2o_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, co2_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }

    return SS_ref_db;
}

/**
    Feldspar (Ab-An-San): Elkins & Grove (1990) asymmetric ternary
    Margules, from xMELTS-master/sources/feldspar.c. W[] here is
    display/bookkeeping only (enthalpy terms); obj_gh_fsp reads the full
    H/S/V triples directly (gh_objective_functions.c) since the
    asymmetric formula needs all three per term, not just this flat array.
    order: whaban, whanab, whabor, whorab, whanor, whoran, whabanor
*/
SS_ref G_SS_gh_fsp_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"fsp_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"ab","an","san"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] =  7924.0;
    SS_ref_db.W[1] =     0.0;
    SS_ref_db.W[2] = 18810.0;
    SS_ref_db.W[3] = 27320.0;
    SS_ref_db.W[4] = 38974.0;
    SS_ref_db.W[5] = 40317.0;
    SS_ref_db.W[6] = 12545.0;

    em_data ab_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "ab",  "equilibrium");
    em_data an_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "an",  "equilibrium");
    em_data san_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "san", "equilibrium");

    SS_ref_db.gbase[0] = ab_eq.gb;
    SS_ref_db.gbase[1] = an_eq.gb;
    SS_ref_db.gbase[2] = san_eq.gb;

    SS_ref_db.ElShearMod[0] = ab_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = an_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] = san_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = ab_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = an_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2] = san_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = ab_eq.ElCp;
    SS_ref_db.ElCp[1] = an_eq.ElCp;
    SS_ref_db.ElCp[2] = san_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = ab_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = an_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2] = san_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = ab_eq.C[i];
        SS_ref_db.Comp[1][i] = an_eq.C[i];
        SS_ref_db.Comp[2][i] = san_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, ab_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, an_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }
    if (GH_boiled_out(len_ox, san_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[2] = 1.0; SS_ref_db.z_em[2] = 0.0; SS_ref_db.bounds_ref[2][0] = 0.0; SS_ref_db.bounds_ref[2][1] = 0.0; }

    return SS_ref_db;
}

/**
    Biotite (Ann-Phl): single symmetric W term, WHANPH from
    xMELTS-master/sources/biotite.c (T,P-independent).
*/
SS_ref G_SS_gh_bi_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"bi_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"ann","phl"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 4.184*1000.0*3.21;

    em_data ann_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "ann", "equilibrium");
    em_data phl_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "phl", "equilibrium");

    SS_ref_db.gbase[0] = ann_eq.gb;
    SS_ref_db.gbase[1] = phl_eq.gb;

    SS_ref_db.ElShearMod[0] = ann_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = phl_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = ann_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = phl_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = ann_eq.ElCp;
    SS_ref_db.ElCp[1] = phl_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = ann_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = phl_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = ann_eq.C[i];
        SS_ref_db.Comp[1][i] = phl_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, ann_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, phl_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }

    return SS_ref_db;
}

/**
    Garnet (Gr-Py-Alm): Berman & Koziol (1991) asymmetric ternary
    Margules, from xMELTS-master/sources/garnet.c. Same W[] bookkeeping-
    only caveat as feldspar above; obj_gh_g reads the H/S/V triples
    directly. order: WH112, WH122, WH113, WH133, WH223, WH233, WH123
*/
SS_ref G_SS_gh_g_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"g_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"gr","py","alm"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 21560.0;
    SS_ref_db.W[1] = 69200.0;
    SS_ref_db.W[2] = 20320.0;
    SS_ref_db.W[3] =  2620.0;
    SS_ref_db.W[4] =   230.0;
    SS_ref_db.W[5] =  3720.0;
    SS_ref_db.W[6] =     0.0;

    em_data gr_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "gr",  "equilibrium");
    em_data py_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "py",  "equilibrium");
    em_data alm_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "alm", "equilibrium");

    SS_ref_db.gbase[0] = gr_eq.gb;
    SS_ref_db.gbase[1] = py_eq.gb;
    SS_ref_db.gbase[2] = alm_eq.gb;

    SS_ref_db.ElShearMod[0] = gr_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = py_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] = alm_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = gr_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = py_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2] = alm_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = gr_eq.ElCp;
    SS_ref_db.ElCp[1] = py_eq.ElCp;
    SS_ref_db.ElCp[2] = alm_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = gr_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = py_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2] = alm_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = gr_eq.C[i];
        SS_ref_db.Comp[1][i] = py_eq.C[i];
        SS_ref_db.Comp[2][i] = alm_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, gr_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, py_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }
    if (GH_boiled_out(len_ox, alm_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[2] = 1.0; SS_ref_db.z_em[2] = 0.0; SS_ref_db.bounds_ref[2][0] = 0.0; SS_ref_db.bounds_ref[2][1] = 0.0; }

    return SS_ref_db;
}

/**
    Hornblende (Parg-Fparg-Mhst): from xMELTS' sources/hornblende.c. NS=0
    (no internal order parameter), despite 2 real crystallographic sites
    (M12, M3) - the site fractions are simple affine functions of p[1],p[2]
    directly (see obj_gh_hb's header comment). W[] here is bookkeeping
    only (matches fsp/g's convention); obj_gh_hb reads its own local
    constants. order: WFEMG (M12 Fe-Mg exchange), WFEAL (M3 Fe3+-Al exchange)
*/
SS_ref G_SS_gh_hb_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"hb_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"parg","fparg","mhst"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 28116.48;
    SS_ref_db.W[1] = 16780.0;

    em_data parg_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "parg",  "equilibrium");
    em_data fparg_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "fparg", "equilibrium");
    em_data mhst_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "mhst",  "equilibrium");

    SS_ref_db.gbase[0] = parg_eq.gb;
    SS_ref_db.gbase[1] = fparg_eq.gb;
    SS_ref_db.gbase[2] = mhst_eq.gb;

    SS_ref_db.ElShearMod[0] = parg_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = fparg_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] = mhst_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = parg_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = fparg_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2] = mhst_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = parg_eq.ElCp;
    SS_ref_db.ElCp[1] = fparg_eq.ElCp;
    SS_ref_db.ElCp[2] = mhst_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = parg_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = fparg_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2] = mhst_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = parg_eq.C[i];
        SS_ref_db.Comp[1][i] = fparg_eq.C[i];
        SS_ref_db.Comp[2][i] = mhst_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, parg_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, fparg_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }
    if (GH_boiled_out(len_ox, mhst_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[2] = 1.0; SS_ref_db.z_em[2] = 0.0; SS_ref_db.bounds_ref[2][0] = 0.0; SS_ref_db.bounds_ref[2][1] = 0.0; }

    return SS_ref_db;
}

/**
    Leucite (Lc-Anl-Nlc): from xMELTS' sources/leucite.c. NS=0, but the
    ideal entropy is a genuine two-site (L: K-Na-H2O, S: K-Na-vacancy,
    multiplicity 3/2) model, not a plain 3-endmember Sigma(p)*log(p) -
    see obj_gh_lc's header comment. W[] here is bookkeeping only; obj_gh_lc
    reads its own local constants. order: WNAK, WNAH, DGR (ternary)
*/
SS_ref G_SS_gh_lc_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"lc_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"lc","anl","nlc"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 7000.0;
    SS_ref_db.W[1] = 7000.0;
    SS_ref_db.W[2] = 53000.0;

    em_data lc_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "lc",  "equilibrium");
    em_data anl_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "anl", "equilibrium");
    em_data nlc_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "nlc", "equilibrium");

    SS_ref_db.gbase[0] = lc_eq.gb;
    SS_ref_db.gbase[1] = anl_eq.gb;
    SS_ref_db.gbase[2] = nlc_eq.gb;

    SS_ref_db.ElShearMod[0] = lc_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = anl_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] = nlc_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = lc_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = anl_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2] = nlc_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = lc_eq.ElCp;
    SS_ref_db.ElCp[1] = anl_eq.ElCp;
    SS_ref_db.ElCp[2] = nlc_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = lc_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = anl_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2] = nlc_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = lc_eq.C[i];
        SS_ref_db.Comp[1][i] = anl_eq.C[i];
        SS_ref_db.Comp[2][i] = nlc_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, lc_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, anl_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }
    if (GH_boiled_out(len_ox, nlc_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[2] = 1.0; SS_ref_db.z_em[2] = 0.0; SS_ref_db.bounds_ref[2][0] = 0.0; SS_ref_db.bounds_ref[2][1] = 0.0; }

    return SS_ref_db;
}

/**
    Melilite (Ak-Geh-Fak-Na), from xMELTS' sources/melilite.c - the first
    gh phase with a genuine internal order parameter (NS=1, Al/Si
    tetrahedral ordering within the gehlenite component - see obj_gh_mel's
    header comment). Iron-akermanite and soda-melilite are reciprocal
    endmembers (see this file's header comment on "ak0"/"geh0" above):
    their gbase/Comp/ElXxx are built here at runtime as
    ak0 + 0.5*fa - 0.5*fo and geh0 + 2*ab - 2*an respectively, reusing
    fo/fa (already in this table for "ol") and ab/an (already in this
    table for "fsp") unchanged - verified byte-identical to what xMELTS'
    own gibbs.c uses for this same reaction.
    order: DG22p, W12, W12p, W13, W14, W22p, W23, W2p3, W24, W2p4
*/
SS_ref G_SS_gh_mel_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"mel_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"ak","geh","fak","na"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 12000.0;
    SS_ref_db.W[1] =  9000.0;
    SS_ref_db.W[2] = 13000.0;
    SS_ref_db.W[3] = 16000.0;
    SS_ref_db.W[4] =     0.0;
    SS_ref_db.W[5] = 38354.0;
    SS_ref_db.W[6] =  9000.0;
    SS_ref_db.W[7] = 13000.0;
    SS_ref_db.W[8] =  9000.0;
    SS_ref_db.W[9] = 13000.0;

    em_data ak_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "ak",   "equilibrium");
    em_data geh_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "geh",  "equilibrium");
    em_data ak0_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "ak0",  "equilibrium");
    em_data geh0_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "geh0", "equilibrium");
    em_data fo_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "fo",   "equilibrium");
    em_data fa_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "fa",   "equilibrium");
    em_data ab_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "ab",   "equilibrium");
    em_data an_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "an",   "equilibrium");

    SS_ref_db.gbase[0] = ak_eq.gb;
    SS_ref_db.gbase[1] = geh_eq.gb;
    SS_ref_db.gbase[2] = ak0_eq.gb  + 0.5*fa_eq.gb - 0.5*fo_eq.gb;
    SS_ref_db.gbase[3] = geh0_eq.gb + 2.0*ab_eq.gb - 2.0*an_eq.gb;

    SS_ref_db.ElShearMod[0] = ak_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = geh_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] = ak0_eq.ElShearMod  + 0.5*fa_eq.ElShearMod - 0.5*fo_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] = geh0_eq.ElShearMod + 2.0*ab_eq.ElShearMod - 2.0*an_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = ak_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = geh_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2] = ak0_eq.ElBulkMod  + 0.5*fa_eq.ElBulkMod - 0.5*fo_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[3] = geh0_eq.ElBulkMod + 2.0*ab_eq.ElBulkMod - 2.0*an_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = ak_eq.ElCp;
    SS_ref_db.ElCp[1] = geh_eq.ElCp;
    SS_ref_db.ElCp[2] = ak0_eq.ElCp  + 0.5*fa_eq.ElCp - 0.5*fo_eq.ElCp;
    SS_ref_db.ElCp[3] = geh0_eq.ElCp + 2.0*ab_eq.ElCp - 2.0*an_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = ak_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = geh_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2] = ak0_eq.ElExpansivity  + 0.5*fa_eq.ElExpansivity - 0.5*fo_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[3] = geh0_eq.ElExpansivity + 2.0*ab_eq.ElExpansivity - 2.0*an_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = ak_eq.C[i];
        SS_ref_db.Comp[1][i] = geh_eq.C[i];
        SS_ref_db.Comp[2][i] = ak0_eq.C[i]  + 0.5*fa_eq.C[i] - 0.5*fo_eq.C[i];
        SS_ref_db.Comp[3][i] = geh0_eq.C[i] + 2.0*ab_eq.C[i] - 2.0*an_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, ak_eq.C,          z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, geh_eq.C,         z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }
    if (GH_boiled_out(len_ox, SS_ref_db.Comp[2], z_b.bulk_rock)){ SS_ref_db.d_em[2] = 1.0; SS_ref_db.z_em[2] = 0.0; SS_ref_db.bounds_ref[2][0] = 0.0; SS_ref_db.bounds_ref[2][1] = 0.0; }
    if (GH_boiled_out(len_ox, SS_ref_db.Comp[3], z_b.bulk_rock)){ SS_ref_db.d_em[3] = 1.0; SS_ref_db.z_em[3] = 0.0; SS_ref_db.bounds_ref[3][0] = 0.0; SS_ref_db.bounds_ref[3][1] = 0.0; }

    return SS_ref_db;
}

/**
    Cummingtonite (Cumm-Grun, Fe-Mg amphibole): from xMELTS'
    sources/cummingtonite.c. NS=2 (M4, M1+M3, M2 site Fe-Mg ordering) -
    same "fetch by name, populate gbase/Comp/El* directly" pattern as
    every other gh phase; W[] set for bookkeeping only (the real Landau
    coefficients live in obj_gh_cum, matching how melilite's W[] is used).
*/
SS_ref G_SS_gh_cum_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"cum_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"cumm","grun"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = -25856.0;

    em_data cumm_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "cumm", "equilibrium");
    em_data grun_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "grun", "equilibrium");

    SS_ref_db.gbase[0] = cumm_eq.gb;
    SS_ref_db.gbase[1] = grun_eq.gb;

    SS_ref_db.ElShearMod[0] = cumm_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = grun_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = cumm_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = grun_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = cumm_eq.ElCp;
    SS_ref_db.ElCp[1] = grun_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = cumm_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = grun_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = cumm_eq.C[i];
        SS_ref_db.Comp[1][i] = grun_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, cumm_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, grun_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }

    return SS_ref_db;
}

/**
    Spinel (Chr-Herc-Mt-Spl-Usp): from xMELTS' sources/spinel.c (Sack &
    Ghiorso 1991). W[] set for bookkeeping only (the real Taylor
    coefficients live in obj_gh_spn, matching melilite/cummingtonite).
*/
SS_ref G_SS_gh_spn_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"spn_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"chr","herc","mt","spl","usp"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 4.5e3*4.184;

    em_data chr_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "chr",  "equilibrium");
    em_data herc_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "herc", "equilibrium");
    em_data mt_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "mt",   "equilibrium");
    em_data spl_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "spl",  "equilibrium");
    em_data usp_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "usp",  "equilibrium");

    SS_ref_db.gbase[0] = chr_eq.gb;
    SS_ref_db.gbase[1] = herc_eq.gb;
    SS_ref_db.gbase[2] = mt_eq.gb;
    SS_ref_db.gbase[3] = spl_eq.gb;
    SS_ref_db.gbase[4] = usp_eq.gb;

    SS_ref_db.ElShearMod[0] = chr_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = herc_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] = mt_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] = spl_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] = usp_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = chr_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = herc_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2] = mt_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[3] = spl_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[4] = usp_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = chr_eq.ElCp;
    SS_ref_db.ElCp[1] = herc_eq.ElCp;
    SS_ref_db.ElCp[2] = mt_eq.ElCp;
    SS_ref_db.ElCp[3] = spl_eq.ElCp;
    SS_ref_db.ElCp[4] = usp_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = chr_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = herc_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2] = mt_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[3] = spl_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[4] = usp_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = chr_eq.C[i];
        SS_ref_db.Comp[1][i] = herc_eq.C[i];
        SS_ref_db.Comp[2][i] = mt_eq.C[i];
        SS_ref_db.Comp[3][i] = spl_eq.C[i];
        SS_ref_db.Comp[4][i] = usp_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, chr_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, herc_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }
    if (GH_boiled_out(len_ox, mt_eq.C,   z_b.bulk_rock)){ SS_ref_db.d_em[2] = 1.0; SS_ref_db.z_em[2] = 0.0; SS_ref_db.bounds_ref[2][0] = 0.0; SS_ref_db.bounds_ref[2][1] = 0.0; }
    if (GH_boiled_out(len_ox, spl_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[3] = 1.0; SS_ref_db.z_em[3] = 0.0; SS_ref_db.bounds_ref[3][0] = 0.0; SS_ref_db.bounds_ref[3][1] = 0.0; }
    if (GH_boiled_out(len_ox, usp_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[4] = 1.0; SS_ref_db.z_em[4] = 0.0; SS_ref_db.bounds_ref[4][0] = 0.0; SS_ref_db.bounds_ref[4][1] = 0.0; }

    return SS_ref_db;
}

/**
    Rhombohedral oxide (Geikielite-Hematite-Ilmenite-Pyrophanite-Corundum):
    from xMELTS' sources/rhomsghiorso.c (Ghiorso & Evans 2008). W[] set for
    bookkeeping only (the real reciprocal-solution Taylor coefficients and
    the joint 3-variable order-parameter solve live in obj_gh_rhm, matching
    spinel/melilite/cummingtonite). EM_list order (gei,hem,ilm,pyr,crn)
    matches xMELTS' own FORMULAS[NA]={"MgTiO3","Fe2O3","FeTiO3","MnTiO3",
    "Al2O3"} array, i.e. r[0]=Xilm(->index2), r[1]=Xgei(->index0),
    r[2]=Xpyr(->index3), r[3]=Xcrn(->index4), hematite(->index1) dependent.
    Each endmember's own gbase already carries its intracrystalline-
    ordering (ilm/gei/pyr) or short-range-order (hem/crn) correction, added
    in GH_gem_function.c's GH_G_EM_function - see GH_rhm_pure_order_G/
    GH_rhm_SRO_G there.
*/
SS_ref G_SS_gh_rhm_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"rhm_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"gei","hem","ilm","pyr","crn"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;

    em_data gei_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "gei", "equilibrium");
    em_data hem_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "hem", "equilibrium");
    em_data ilm_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "ilm", "equilibrium");
    em_data pyr_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "pyr", "equilibrium");
    em_data crn_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "crn", "equilibrium");

    SS_ref_db.gbase[0] = gei_eq.gb;
    SS_ref_db.gbase[1] = hem_eq.gb;
    SS_ref_db.gbase[2] = ilm_eq.gb;
    SS_ref_db.gbase[3] = pyr_eq.gb;
    SS_ref_db.gbase[4] = crn_eq.gb;

    SS_ref_db.ElShearMod[0] = gei_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = hem_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] = ilm_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] = pyr_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] = crn_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = gei_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = hem_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2] = ilm_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[3] = pyr_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[4] = crn_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = gei_eq.ElCp;
    SS_ref_db.ElCp[1] = hem_eq.ElCp;
    SS_ref_db.ElCp[2] = ilm_eq.ElCp;
    SS_ref_db.ElCp[3] = pyr_eq.ElCp;
    SS_ref_db.ElCp[4] = crn_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = gei_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = hem_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2] = ilm_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[3] = pyr_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[4] = crn_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = gei_eq.C[i];
        SS_ref_db.Comp[1][i] = hem_eq.C[i];
        SS_ref_db.Comp[2][i] = ilm_eq.C[i];
        SS_ref_db.Comp[3][i] = pyr_eq.C[i];
        SS_ref_db.Comp[4][i] = crn_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, gei_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, hem_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }
    if (GH_boiled_out(len_ox, ilm_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[2] = 1.0; SS_ref_db.z_em[2] = 0.0; SS_ref_db.bounds_ref[2][0] = 0.0; SS_ref_db.bounds_ref[2][1] = 0.0; }
    if (GH_boiled_out(len_ox, pyr_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[3] = 1.0; SS_ref_db.z_em[3] = 0.0; SS_ref_db.bounds_ref[3][0] = 0.0; SS_ref_db.bounds_ref[3][1] = 0.0; }
    if (GH_boiled_out(len_ox, crn_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[4] = 1.0; SS_ref_db.z_em[4] = 0.0; SS_ref_db.bounds_ref[4][0] = 0.0; SS_ref_db.bounds_ref[4][1] = 0.0; }

    return SS_ref_db;
}

/**
    Nepheline (Sack & Ghiorso 1995): from xMELTS' sources/nepheline.c.
    NR=3/NS=1 - a single embedded order parameter, see obj_gh_nph. EM_list
    order (nane,kne,vcne,cane) matches xMELTS' own conNph component order
    m[0]=na-nepheline,m[1]=k-nepheline,m[2]=vc-nepheline,m[3]=ca-nepheline:
    r[0]=X(kne)->index1, r[1]=X(vcne)->index2, r[2]=X(cane)->index3,
    na-nepheline(->index0) dependent - a direct r-basis relabeling, same
    pattern spn/rhm already use. Uses "kne" (the nepheline-hosted
    k-nepheline row, WITH the -3572.76 offset) - kalsilite ("kls" below)
    uses "knk" instead, see GH_PP_endmembers.c's header comment.
*/
SS_ref G_SS_gh_nph_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"nph_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"nane","kne","vcne","cane"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;

    em_data nane_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "nane", "equilibrium");
    em_data kne_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "kne",  "equilibrium");
    em_data vcne_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "vcne", "equilibrium");
    em_data cane_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "cane", "equilibrium");

    SS_ref_db.gbase[0] = nane_eq.gb;
    SS_ref_db.gbase[1] = kne_eq.gb;
    SS_ref_db.gbase[2] = vcne_eq.gb;
    SS_ref_db.gbase[3] = cane_eq.gb;

    SS_ref_db.ElShearMod[0] = nane_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = kne_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] = vcne_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] = cane_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = nane_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = kne_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2] = vcne_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[3] = cane_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = nane_eq.ElCp;
    SS_ref_db.ElCp[1] = kne_eq.ElCp;
    SS_ref_db.ElCp[2] = vcne_eq.ElCp;
    SS_ref_db.ElCp[3] = cane_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = nane_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = kne_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2] = vcne_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[3] = cane_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = nane_eq.C[i];
        SS_ref_db.Comp[1][i] = kne_eq.C[i];
        SS_ref_db.Comp[2][i] = vcne_eq.C[i];
        SS_ref_db.Comp[3][i] = cane_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, nane_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, kne_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }
    if (GH_boiled_out(len_ox, vcne_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[2] = 1.0; SS_ref_db.z_em[2] = 0.0; SS_ref_db.bounds_ref[2][0] = 0.0; SS_ref_db.bounds_ref[2][1] = 0.0; }
    if (GH_boiled_out(len_ox, cane_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[3] = 1.0; SS_ref_db.z_em[3] = 0.0; SS_ref_db.bounds_ref[3][0] = 0.0; SS_ref_db.bounds_ref[3][1] = 0.0; }

    return SS_ref_db;
}

/**
    Kalsilite (Sack & Ghiorso 1995): from xMELTS' sources/kalsilite.c. NR=3,
    NO order parameter (unlike nepheline) - a direct algebraic solution,
    see obj_gh_kls. Shares 3 of nepheline's 4 endmembers verbatim (nane,
    vcne, cane) but uses its OWN k-nepheline row "knk" (WITHOUT the
    -3572.76 nepheline-only offset - see GH_PP_endmembers.c's header
    comment). Same EM_list order/r-basis convention as "nph".
*/
SS_ref G_SS_gh_kls_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"kls_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"nane","knk","vcne","cane"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;

    em_data nane_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "nane", "equilibrium");
    em_data knk_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "knk",  "equilibrium");
    em_data vcne_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "vcne", "equilibrium");
    em_data cane_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "cane", "equilibrium");

    SS_ref_db.gbase[0] = nane_eq.gb;
    SS_ref_db.gbase[1] = knk_eq.gb;
    SS_ref_db.gbase[2] = vcne_eq.gb;
    SS_ref_db.gbase[3] = cane_eq.gb;

    SS_ref_db.ElShearMod[0] = nane_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = knk_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] = vcne_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] = cane_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = nane_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = knk_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2] = vcne_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[3] = cane_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = nane_eq.ElCp;
    SS_ref_db.ElCp[1] = knk_eq.ElCp;
    SS_ref_db.ElCp[2] = vcne_eq.ElCp;
    SS_ref_db.ElCp[3] = cane_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = nane_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = knk_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2] = vcne_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[3] = cane_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = nane_eq.C[i];
        SS_ref_db.Comp[1][i] = knk_eq.C[i];
        SS_ref_db.Comp[2][i] = vcne_eq.C[i];
        SS_ref_db.Comp[3][i] = cane_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, nane_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, knk_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }
    if (GH_boiled_out(len_ox, vcne_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[2] = 1.0; SS_ref_db.z_em[2] = 0.0; SS_ref_db.bounds_ref[2][0] = 0.0; SS_ref_db.bounds_ref[2][1] = 0.0; }
    if (GH_boiled_out(len_ox, cane_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[3] = 1.0; SS_ref_db.z_em[3] = 0.0; SS_ref_db.bounds_ref[3][0] = 0.0; SS_ref_db.bounds_ref[3][1] = 0.0; }

    return SS_ref_db;
}

/**
    Clinopyroxene (Di-Cen-Hed-CaTs(Al)-CaTs(Fe3+)-Ess-Jd): from xMELTS'
    sources/clinopyroxene.c (Sack & Ghiorso 1993). W[] set for bookkeeping
    only (the real ~245 Taylor coefficients live in obj_gh_cpx, matching
    melilite/cummingtonite/spinel).
*/
SS_ref G_SS_gh_cpx_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    strcpy(SS_ref_db.fName,"cpx_MELTS");
    int i;
    int n_em = SS_ref_db.n_em;

    char *EM_tmp[] = {"di","cen","hed","cats","buff","ess","jd"};
    for (i = 0; i < n_em; i++){
        strcpy(SS_ref_db.EM_list[i], EM_tmp[i]);
    };

    SS_ref_db.W[0] = 0.0;

    em_data di_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "di",   "equilibrium");
    em_data cen_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "cen",  "equilibrium");
    em_data hed_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "hed",  "equilibrium");
    em_data cats_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "cats", "equilibrium");
    em_data buff_eq = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "buff", "equilibrium");
    em_data ess_eq  = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "ess",  "equilibrium");
    em_data jd_eq   = get_em_data(research_group, EM_dataset, len_ox, z_b, SS_ref_db.P, SS_ref_db.T, "jd",   "equilibrium");

    SS_ref_db.gbase[0] = di_eq.gb;
    SS_ref_db.gbase[1] = cen_eq.gb;
    SS_ref_db.gbase[2] = hed_eq.gb;
    SS_ref_db.gbase[3] = cats_eq.gb;
    SS_ref_db.gbase[4] = buff_eq.gb;
    SS_ref_db.gbase[5] = ess_eq.gb;
    SS_ref_db.gbase[6] = jd_eq.gb;

    SS_ref_db.ElShearMod[0] = di_eq.ElShearMod;
    SS_ref_db.ElShearMod[1] = cen_eq.ElShearMod;
    SS_ref_db.ElShearMod[2] = hed_eq.ElShearMod;
    SS_ref_db.ElShearMod[3] = cats_eq.ElShearMod;
    SS_ref_db.ElShearMod[4] = buff_eq.ElShearMod;
    SS_ref_db.ElShearMod[5] = ess_eq.ElShearMod;
    SS_ref_db.ElShearMod[6] = jd_eq.ElShearMod;

    SS_ref_db.ElBulkMod[0] = di_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[1] = cen_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[2] = hed_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[3] = cats_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[4] = buff_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[5] = ess_eq.ElBulkMod;
    SS_ref_db.ElBulkMod[6] = jd_eq.ElBulkMod;

    SS_ref_db.ElCp[0] = di_eq.ElCp;
    SS_ref_db.ElCp[1] = cen_eq.ElCp;
    SS_ref_db.ElCp[2] = hed_eq.ElCp;
    SS_ref_db.ElCp[3] = cats_eq.ElCp;
    SS_ref_db.ElCp[4] = buff_eq.ElCp;
    SS_ref_db.ElCp[5] = ess_eq.ElCp;
    SS_ref_db.ElCp[6] = jd_eq.ElCp;

    SS_ref_db.ElExpansivity[0] = di_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[1] = cen_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[2] = hed_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[3] = cats_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[4] = buff_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[5] = ess_eq.ElExpansivity;
    SS_ref_db.ElExpansivity[6] = jd_eq.ElExpansivity;

    for (i = 0; i < len_ox; i++){
        SS_ref_db.Comp[0][i] = di_eq.C[i];
        SS_ref_db.Comp[1][i] = cen_eq.C[i];
        SS_ref_db.Comp[2][i] = hed_eq.C[i];
        SS_ref_db.Comp[3][i] = cats_eq.C[i];
        SS_ref_db.Comp[4][i] = buff_eq.C[i];
        SS_ref_db.Comp[5][i] = ess_eq.C[i];
        SS_ref_db.Comp[6][i] = jd_eq.C[i];
    }

    for (i = 0; i < n_em; i++){
        SS_ref_db.z_em[i] = 1.0;
    };

    for (i = 0; i < SS_ref_db.n_xeos; i++){
        SS_ref_db.bounds_ref[i][0] = 0.0+eps;
        SS_ref_db.bounds_ref[i][1] = 1.0-eps;
    }

    for (i = 0; i < n_em; i++){ SS_ref_db.d_em[i] = 0.0; }
    if (GH_boiled_out(len_ox, di_eq.C,   z_b.bulk_rock)){ SS_ref_db.d_em[0] = 1.0; SS_ref_db.z_em[0] = 0.0; SS_ref_db.bounds_ref[0][0] = 0.0; SS_ref_db.bounds_ref[0][1] = 0.0; }
    if (GH_boiled_out(len_ox, cen_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[1] = 1.0; SS_ref_db.z_em[1] = 0.0; SS_ref_db.bounds_ref[1][0] = 0.0; SS_ref_db.bounds_ref[1][1] = 0.0; }
    if (GH_boiled_out(len_ox, hed_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[2] = 1.0; SS_ref_db.z_em[2] = 0.0; SS_ref_db.bounds_ref[2][0] = 0.0; SS_ref_db.bounds_ref[2][1] = 0.0; }
    if (GH_boiled_out(len_ox, cats_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[3] = 1.0; SS_ref_db.z_em[3] = 0.0; SS_ref_db.bounds_ref[3][0] = 0.0; SS_ref_db.bounds_ref[3][1] = 0.0; }
    if (GH_boiled_out(len_ox, buff_eq.C, z_b.bulk_rock)){ SS_ref_db.d_em[4] = 1.0; SS_ref_db.z_em[4] = 0.0; SS_ref_db.bounds_ref[4][0] = 0.0; SS_ref_db.bounds_ref[4][1] = 0.0; }
    if (GH_boiled_out(len_ox, ess_eq.C,  z_b.bulk_rock)){ SS_ref_db.d_em[5] = 1.0; SS_ref_db.z_em[5] = 0.0; SS_ref_db.bounds_ref[5][0] = 0.0; SS_ref_db.bounds_ref[5][1] = 0.0; }
    if (GH_boiled_out(len_ox, jd_eq.C,   z_b.bulk_rock)){ SS_ref_db.d_em[6] = 1.0; SS_ref_db.z_em[6] = 0.0; SS_ref_db.bounds_ref[6][0] = 0.0; SS_ref_db.bounds_ref[6][1] = 0.0; }

    return SS_ref_db;
}

/**
    Orthopyroxene: same 7 endmembers/energetics as clinopyroxene (see
    obj_gh_cpx's header comment) - registered separately so MAGEMin's own
    phase competition can pick cpx vs opx by stability.
*/
SS_ref G_SS_gh_opx_function(SS_ref SS_ref_db, char* research_group, int EM_dataset, int len_ox, bulk_info z_b, double eps){
    SS_ref_db = G_SS_gh_cpx_function(SS_ref_db, research_group, EM_dataset, len_ox, z_b, eps);
    strcpy(SS_ref_db.fName,"opx_MELTS");
    return SS_ref_db;
}

/**
    Per-dataset dispatch (mirrors G_SS_sb11_EM_function): runs the
    phase-specific reference-state function gv.n_Diff times at the P/T
    finite-difference stencil, storing each pass's gbase[] for later
    FD-derived volume/Cp/expansivity of the solution phase (toolkit.c).
*/
SS_ref G_SS_gh_EM_function(    global_variable  gv,
                                SS_ref           SS_ref_db,
                                int              EM_dataset,
                                bulk_info        z_b,
                                char            *name                  ){

    double eps  = gv.bnd_val;
    double P    = SS_ref_db.P;
    double T    = SS_ref_db.T;

    SS_ref_db.ss_flags[0] = 1;

    for (int FD = 0; FD < gv.n_Diff; FD++){
        SS_ref_db.P = P + gv.gb_P_eps*gv.pdev[0][FD];
        SS_ref_db.T = T + gv.gb_T_eps*gv.pdev[1][FD];

        if (strcmp( name, "liq") == 0 ){
            SS_ref_db = G_SS_gh_liq_function(SS_ref_db, gv.research_group, EM_dataset, gv.EM_database, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "ol") == 0 ){
            SS_ref_db = G_SS_gh_ol_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "fsp") == 0 ){
            SS_ref_db = G_SS_gh_fsp_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "bi") == 0 ){
            SS_ref_db = G_SS_gh_bi_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "g") == 0 ){
            SS_ref_db = G_SS_gh_g_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "hb") == 0 ){
            SS_ref_db = G_SS_gh_hb_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "lc") == 0 ){
            SS_ref_db.ss_flags[0]  = 0;
            SS_ref_db = G_SS_gh_lc_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "mel") == 0 ){
            SS_ref_db = G_SS_gh_mel_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "cum") == 0 ){
            SS_ref_db = G_SS_gh_cum_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "spn") == 0 ){
            SS_ref_db = G_SS_gh_spn_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "cpx") == 0 ){
            SS_ref_db = G_SS_gh_cpx_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "opx") == 0 ){
            SS_ref_db = G_SS_gh_opx_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "fl") == 0 ){
            SS_ref_db = G_SS_gh_fluid_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "rhm") == 0 ){
            SS_ref_db = G_SS_gh_rhm_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "nph") == 0 ){
            SS_ref_db = G_SS_gh_nph_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else if (strcmp( name, "kls") == 0 ){
            SS_ref_db = G_SS_gh_kls_function(SS_ref_db, gv.research_group, EM_dataset, gv.len_ox, z_b, eps);
        }
        else{
            printf("\nsolid solution '%s' is not in the 'gh' database\n", name);
        }

        for (int j = 0; j < SS_ref_db.n_em; j++){
            SS_ref_db.mu_array[FD][j] = SS_ref_db.gbase[j];
        }
    }

    for (int j = 0; j < SS_ref_db.n_xeos; j++){
        SS_ref_db.bounds[j][0] = SS_ref_db.bounds_ref[j][0];
        SS_ref_db.bounds[j][1] = SS_ref_db.bounds_ref[j][1];
    }

    double fbc = 0.0;
    for (int i = 0; i < gv.len_ox; i++){
        fbc += z_b.bulk_rock[i]*z_b.apo[i];
    }

    for (int i = 0; i < SS_ref_db.n_em; i++){
        SS_ref_db.ape[i] = 0.0;
        for (int j = 0; j < gv.len_ox; j++){
            SS_ref_db.ape[i] += SS_ref_db.Comp[i][j]*z_b.apo[j];
        }
    }

    SS_ref_db.fbc = z_b.fbc;

    /* mirrors the verbose endmember-name/gbase/composition dump that TC's
       per-dataset EM_function dispatchers print (e.g. G_SS_ig_EM_function
       in TC_database/tc_gss_function.c) - using the actual gv.ox[] oxide
       names rather than TC's hardcoded single-letter abbreviations, since
       "gh"'s oxide list/order differs from any TC dataset's.             */
    if (gv.verbose == 1){
        printf(" %4s:\n", name);
        printf("----\n");
        for (int j = 0; j < SS_ref_db.n_em; j++){
            printf(" %12s", SS_ref_db.EM_list[j]);
        }
        printf("\n");
        for (int j = 0; j < SS_ref_db.n_em; j++){
            printf(" %+12.5f", SS_ref_db.gbase[j]);
        }
        printf("\n");

        for (int j = 0; j < gv.len_ox; j++){
            printf(" %8s", gv.ox[j]);
        }
        printf("\n");
        for (int i = 0; i < SS_ref_db.n_em; i++){
            for (int j = 0; j < gv.len_ox; j++){
                printf(" %+8.2f", SS_ref_db.Comp[i][j]);
            }
            printf("\n");
        }
        printf("\n");
    }

    return SS_ref_db;
};
