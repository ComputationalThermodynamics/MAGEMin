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
#include "GH_endmembers.h"

static const EM_db_gh arr_em_db_gh_liq[13] = {
    /* SiO2 */
    { "SiO2", { 1.0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -906377.0, 46.029, 2.730,
        { 83.51, -374.7, -2.4554e6, 2.8007e8, 0.0, 0.0, 0.0, 0.0 },
        2.690, { 0.0, -1.89e-5, 1.3e-8, 3.6e-10 }, 1999.0, 4.46, 82.6 },
    /* TiO2 */
    { "TiO2", { 0,0,0,0,0,0,0, 1.0, 0,0,0,0,0,0,0,0 },
        -944750.0, 50.460, 1.882,
        { 77.84, 0.0, -3.3678e6, 4.0294e8, 0.0, 0.0, 0.0, 0.0 },
        2.316, { 7.246e-4, -2.310e-5, 0.0, 5.0e-10 }, 1870.0, 35.824, 109.2 },
    /* Al2O3 */
    { "Al2O3", { 0, 1.0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -1675700.0, 50.82, 0.0,
        { 155.02, -828.4, -3.8614e6, 4.0908e8, 0.0, 0.0, 0.0, 0.0 },
        3.711, { 2.62e-4, -2.26e-5, 2.7e-8, 4.0e-10 }, 2319.65, 48.61, 170.3 },
    /* Fe2O3 = 2 FeO + O */
    { "Fe2O3", { 0,0,0,0, 2.0, 0,0,0, 1.0, 0,0,0,0,0,0,0 },
        -822000.0, 87.40, 3.027,
        { 146.86, 0.0, -5.5768e6, 5.2563e8, 955.0, 1287.0, -7.403e-2, 27.921e-5 },
        4.213, { 9.09e-4, -2.53e-5, 3.1e-8, 4.4e-10 }, 1895.0, 60.41, 240.9 },
    /* MgCr2O4 = MgO + Cr2O3 */
    { "MgCr2O4", { 0,0,0, 1.0, 0,0,0,0,0,0, 1.0, 0,0,0,0,0 },
        -1783640.0, 106.02, 4.3560,
        { 201.981, -551.9, -5.7844e6, 5.7729e8, 0.0, 0.0, 0.0, 0.0 },
        5.358, { 11.71e-4, -2.26e-5, 1.8e-8, 4.67e-10 }, 2673.15, 73.22, 335.1 },
    /* Fe2SiO4 = 2 FeO + SiO2 */
    { "Fe2SiO4", { 1.0, 0,0,0, 2.0, 0,0,0,0,0,0,0,0,0,0,0 },
        -1479360.0, 150.930, 0.0,
        { 248.93, -1923.9, 0.0, -1.391e8, 0.0, 0.0, 0.0, 0.0 },
        5.420, { 5.84e-4, -2.79e-5, -2.3e-8, 14.6e-10 }, 1490.0, 59.9, 240.2 },
    /* MnSi0.5O2 = MnO + 0.5 SiO2 */
    { "MnSi0.5O2", { 0.5, 0,0,0,0,0,0,0,0, 1.0, 0,0,0,0,0,0 },
        -866000.0, 77.950, 2.4445,
        { 109.945, -635.5, -1.0248e6, 8.826e7, 0.0, 0.0, 0.0, 0.0 },
        2.840, { 2.92e-4, -1.395e-5, -1.15e-8, 7.3e-10 }, 1620.0, 27.6, 121.6 },
    /* Mg2SiO4 = 2 MgO + SiO2 */
    { "Mg2SiO4", { 1.0, 0,0, 2.0, 0,0,0,0,0,0,0,0,0,0,0,0 },
        -2174420.0, 94.010, 0.0,
        { 238.64, -2001.3, 0.0, -1.1624e8, 0.0, 0.0, 0.0, 0.0 },
        4.980, { 5.24e-4, -1.35e-5, -1.3e-8, 4.14e-10 }, 2163.0, 57.2, 271.0 },
    /* CaSiO3 = CaO + SiO2 */
    { "CaSiO3", { 1.0, 0, 1.0, 0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -1627427.0, 85.279, 4.016,
        { 141.16, -417.2, -5.8576e6, 9.4074e8, 0.0, 0.0, 0.0, 0.0 },
        4.347, { 2.92e-4, -1.55e-5, -1.6e-8, 3.89e-10 }, 1817.0, 31.5, 172.4 },
    /* Na2SiO3 = Na2O + SiO2 */
    { "Na2SiO3", { 1.0, 0,0,0,0,0, 1.0, 0,0,0,0,0,0,0,0,0 },
        -1561426.96, 113.84664, 0.0,
        { 234.77, -2218.9, 0.0, 1.353e8, 0.0, 0.0, 0.0, 0.0 },
        5.568, { 7.41e-4, -4.29e-5, -5.3e-8, 8.4e-10 }, 1361.0, 38.34, 180.2 },
    /* KAlSiO4 = 0.5 K2O + 0.5 Al2O3 + SiO2 */
    { "KAlSiO4", { 1.0, 0.5, 0,0,0, 0.5, 0,0,0,0,0,0,0,0,0,0 },
        -2111813.55, 133.9653, 5.989,
        { 186.0, 0.0, -1.31067e7, 2.13893e9, 800.15, 1154.0, -7.096454e-2, 21.682e-5 },
        6.8375, { 7.265e-4, -6.395e-5, -4.6e-8, 12.1e-10 }, 2023.15, 24.5, 217.0 },
    /* CO2 - G computed via GH_pitzer_sterner_G(), not this placeholder data */
    { "CO2", { 0,0,0,0,0,0,0,0,0,0,0,0, 1.0, 0,0,0 },
        0.0, 0.0, 0.0, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        0.0, { 0.0, 0.0, 0.0, 0.0 }, 1000.0, 0.0, 0.0 },
    /* H2O - G computed via GH_pitzer_sterner_G(), not this placeholder data */
    { "H2O", { 0,0,0,0,0,0,0,0,0,0,0, 1.0, 0,0,0,0 },
        0.0, 0.0, 0.0, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        0.0, { 0.0, 0.0, 0.0, 0.0 }, 1000.0, 0.0, 0.0 },
};

/** pMELTS liquid endmembers, ported from xMELTS-master's
    includes/liq_struct_data.h "pMeltsLiquid[]" table (sources/liquid_v34.c,
    calculationMode==MODE_pMELTS - shared with rhyolite-MELTS v1.0.2, but
    NOT what gh's rMELTS uses, since rMELTS replicates v1.2.0/liquid_CO2_H2O.c
    instead, per [[gh-multicalibration-xmelts-rmelts-pmelts]]).

    Differs from the xMELTS table above in 2 ways, both confirmed directly
    against source (H0/S0/V0/Cp/Vliq/kress numbers recomputed from scratch,
    not hand-adjusted from the xMELTS table):
    1. 4 of the 11 real oxide endmembers use a DIFFERENT, rescaled formula
       unit as pMELTS' own "one mole" convention (a real physical choice,
       not a bug - larger units for SiO2/Al2O3/CaSiO3 model silica/alumina/
       wollastonite polymerization in the melt, smaller for Na2SiO3):
       SiO2->Si4O8 (x4), Al2O3->Al4O6 (x2), CaSiO3->Ca2Si2O6 (x2),
       Na2SiO3->NaSi0.5O1.5 (x0.5); the other 9 stay x1. H0/S0/V0/Cp/Vliq/
       kress/Sfus/Cpl all scale by this same multiplier (Tfus does not -
       it's an intensive per-mole quantity, not per-formula-unit); gh keeps
       the plain oxide names (SiO2 not Si4O8) as EM_list labels for
       consistency with the xMELTS table, only the underlying VALUES are
       rescaled.
    2. A small additive "corrHxx"/"corrSxx" calibration correction on top
       of the (possibly rescaled) xMELTS base value for most endmembers -
       all corrVxx are 0 in the real table (confirmed), so V0 needs no
       correction, only H0/S0 do.

    CO2 has NO entry here at all: real pMELTS has no CO2 standard-state
    treatment whatsoever (its own gibbs.c/liquid_v34.c CO2 branch is a
    literal zero placeholder - not a gh gap, pMELTS' calibration focus is
    anhydrous peridotite melting), so gh drops CO2 entirely from this
    dataset's oxide/endmember list (12 oxides, no "fl", no CO2-bearing
    carbonates - see GH_init_database.c's gh_db_pmelts_dataset) rather than
    carry a permanently-inert placeholder around. H2O's entry carries real
    (non-dummy) corrH19/corrS19/V data, but real pMELTS still special-cases
    H2O via its own fluidPhase()-based branch in gibbs.c rather than using
    this liquid-endmember table directly for water's true standard state -
    see the H2O investigation note in [[gh-multicalibration-xmelts-rmelts-pmelts]].
*/
static const EM_db_gh arr_em_db_gh_liq_pMELTS[12] = {
    /* SiO2 (pMELTS: Si4O8, x4 mult) + corrH01/corrS01 - real pMELTS' own
       name for this rescaled formula unit, matching liq_struct_data.h's
       "pMeltsLiquid[]" table exactly (not left as "SiO2" like the other
       9 unrescaled endmembers) - lets GH_G_EM_function's name-based
       dispatch distinguish this from xMELTS/rMELTS' plain "SiO2" without
       an extra EM_database check, and avoids reporting misleading xeos
       labels (a mole fraction of THIS endmember is a mole fraction of
       Si4O8-as-one-particle, not of SiO2 - a real, not cosmetic,
       difference - see this table's own header comment). */
    { "Si4O8", { 4.0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -3555755.403, 229.696, 10.920,
        { 334.04, -1498.8, -9.8216e6, 1.12028e9, 0.0, 0.0, 0.0, 0.0 },
        10.760, { 0.0, -7.56e-5, 5.2e-8, 1.44e-9 }, 1999.0, 17.84, 330.4 },
    /* TiO2 (x1 mult) + corrH02/corrS02(=0) */
    { "TiO2", { 0,0,0,0,0,0,0, 1.0, 0,0,0,0,0,0,0,0 },
        -961690.793, 50.460, 1.882,
        { 77.84, 0.0, -3.3678e6, 4.0294e8, 0.0, 0.0, 0.0, 0.0 },
        2.316, { 7.246e-4, -2.310e-5, 0.0, 5.0e-10 }, 1870.0, 35.824, 109.2 },
    /* Al2O3 (pMELTS: Al4O6, x2 mult) + corrH03/corrS03 - real pMELTS name */
    { "Al4O6", { 0, 2.0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -3307936.987, 76.474, 0.0,
        { 310.04, -1656.8, -7.7228e6, 8.1816e8, 0.0, 0.0, 0.0, 0.0 },
        7.422, { 5.24e-4, -4.52e-5, 5.4e-8, 8.0e-10 }, 2319.65, 97.22, 340.6 },
    /* Fe2O3 = 2 FeO + O (x1 mult) + corrH04/corrS04 */
    { "Fe2O3", { 0,0,0,0, 2.0, 0,0,0, 1.0, 0,0,0,0,0,0,0 },
        -658347.093, 112.813, 3.027,
        { 146.86, 0.0, -5.5768e6, 5.2563e8, 955.0, 1287.0, -7.403e-2, 27.921e-5 },
        4.213, { 9.09e-4, -2.53e-5, 3.1e-8, 4.4e-10 }, 1895.0, 60.41, 240.9 },
    /* MgCr2O4 = MgO + Cr2O3 (x1 mult) + corrH05/corrS05(=0) */
    { "MgCr2O4", { 0,0,0, 1.0, 0,0,0,0,0,0, 1.0, 0,0,0,0,0 },
        -1776401.419, 106.02, 4.3560,
        { 201.981, -551.9, -5.7844e6, 5.7729e8, 0.0, 0.0, 0.0, 0.0 },
        5.358, { 11.71e-4, -2.26e-5, 1.8e-8, 4.67e-10 }, 2673.15, 73.22, 335.1 },
    /* Fe2SiO4 = 2 FeO + SiO2 (x1 mult) + corrH06/corrS06 */
    { "Fe2SiO4", { 1.0, 0,0,0, 2.0, 0,0,0,0,0,0,0,0,0,0,0 },
        -1465620.986, 157.354, 0.0,
        { 248.93, -1923.9, 0.0, -1.391e8, 0.0, 0.0, 0.0, 0.0 },
        5.420, { 5.84e-4, -2.79e-5, -2.3e-8, 14.6e-10 }, 1490.0, 59.9, 240.2 },
    /* MnSi0.5O2 = MnO + 0.5 SiO2 (x1 mult) + corrH07(=0)/corrS07(=0) */
    { "MnSi0.5O2", { 0.5, 0,0,0,0,0,0,0,0, 1.0, 0,0,0,0,0,0 },
        -866000.0, 77.950, 2.4445,
        { 109.945, -635.5, -1.0248e6, 8.826e7, 0.0, 0.0, 0.0, 0.0 },
        2.840, { 2.92e-4, -1.395e-5, -1.15e-8, 7.3e-10 }, 1620.0, 27.6, 121.6 },
    /* Mg2SiO4 = 2 MgO + SiO2 (x1 mult) + corrH08/corrS08 */
    { "Mg2SiO4", { 1.0, 0,0, 2.0, 0,0,0,0,0,0,0,0,0,0,0,0 },
        -2166508.063, 106.787, 0.0,
        { 238.64, -2001.3, 0.0, -1.1624e8, 0.0, 0.0, 0.0, 0.0 },
        4.980, { 5.24e-4, -1.35e-5, -1.3e-8, 4.14e-10 }, 2163.0, 57.2, 271.0 },
    /* CaSiO3 (pMELTS: Ca2Si2O6, x2 mult) + corrH11/corrS11 - real pMELTS name
       (pseudowollastonite) */
    { "Ca2Si2O6", { 2.0, 0, 2.0, 0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -3276716.658, 168.292, 8.032,
        { 282.32, -834.4, -1.17152e7, 1.88148e9, 0.0, 0.0, 0.0, 0.0 },
        8.694, { 5.84e-4, -3.1e-5, -3.2e-8, 7.78e-10 }, 1817.0, 63.0, 344.8 },
    /* Na2SiO3 (pMELTS: NaSi0.5O1.5, x0.5 mult) + corrH12/corrS12(=0) - real
       pMELTS name */
    { "NaSi0.5O1.5", { 0.5, 0,0,0,0,0, 0.5, 0,0,0,0,0,0,0,0,0 },
        -800742.230, 56.9233, 0.0,
        { 117.385, -1109.45, 0.0, 6.765e7, 0.0, 0.0, 0.0, 0.0 },
        2.784, { 3.705e-4, -2.145e-5, -2.65e-8, 4.2e-10 }, 1361.0, 19.17, 90.1 },
    /* KAlSiO4 = 0.5 K2O + 0.5 Al2O3 + SiO2 (x1 mult) + corrH13/corrS13 */
    { "KAlSiO4", { 1.0, 0.5, 0,0,0, 0.5, 0,0,0,0,0,0,0,0,0,0 },
        -2127828.836, 126.6543, 5.989,
        { 186.0, 0.0, -1.31067e7, 2.13893e9, 800.15, 1154.0, -7.096454e-2, 21.682e-5 },
        6.8375, { 7.265e-4, -6.395e-5, -4.6e-8, 12.1e-10 }, 2023.15, 24.5, 217.0 },
    /* H2O (index 11, not 12 - CO2's row is gone entirely) - real pMELTS
       special-cases this via fluidPhase(), not this table's own (non-
       dummy but unused for the true standard state) corrH19/corrS19/V
       data - kept as dummy placeholder pending the fluidPhase()-
       equivalent investigation noted above. */
    { "H2O", { 0,0,0,0,0,0,0,0,0,0,0, 1.0, 0,0,0,0 },
        0.0, 0.0, 0.0, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
        0.0, { 0.0, 0.0, 0.0, 0.0 }, 1000.0, 0.0, 0.0 },
};

EM_db_gh Access_GH_EM_DB(int EM_database, int id){
    if (EM_database == 2){ return arr_em_db_gh_liq_pMELTS[id]; }
    return arr_em_db_gh_liq[id];
}
