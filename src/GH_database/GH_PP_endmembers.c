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
#include <string.h>

#include "GH_PP_endmembers.h"

/**
    Comp[] axis order (identical to GH_endmembers.h / GH_init_database.c):
    0 SiO2, 1 Al2O3, 2 CaO, 3 MgO, 4 FeO, 5 K2O, 6 Na2O, 7 TiO2, 8 O,
    9 MnO, 10 Cr2O3, 11 H2O, 12 CO2

    O2 and H2O pure phases are handled separately in GH_gem_function.c
    (a dedicated ideal-gas formula for O2, ported from xMELTS'
    sources/gibbs.c; the same Pitzer & Sterner (1994) EOS already used for
    the liquid's H2O basis species) and so are not part of this table.
**/
static const PP_db_gh arr_pp_db_gh[GH_N_PP] = {
    /* quartz - SiO2 */
    { "q", { 1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -910700.0, 41.460, 2.269,
        { 80.01, -2.403E2, -35.467E5, 49.157E7, 848.0, 0.0, -9.187E-2, 24.607E-5 },
        { -2.434E-6, 10.137E-12, 23.895E-6, 0.0 } },
    /* cristobalite - SiO2 */
    { "crst", { 1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -907753.0, 43.394, 2.587,
        { 83.51, -3.747E2, -24.554E5, 28.007E7, 535.0, 0.0, -14.216E-2, 44.142E-5 },
        { -2.515E-6, 0.0, 20.824E-6, 0.0 } },
    /* tridymite - SiO2 */
    { "trd", { 1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -907750.0, 43.770, 2.675,
        { 75.37, 0.0, -59.581E5, 95.825E7, 383.0, 130.0, 42.670E-2, -144.575E-5 },
        { -2.508E-6, 0.0, 19.339E-6, 0.0 } },
    /* corundum - Al2O3 */
    { "cor", { 0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -1675700.0, 50.820, 2.558,
        { 155.02, -8.284E2, -38.614E5, 40.908E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.385E-6, 0.375E-12, 21.342E-6, 47.180E-10 } },
    /* sillimanite - Al2SiO5 = Al2O3 + SiO2 */
    { "sill", { 1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -2586091.0, 95.930, 4.983,
        { 256.73, -18.872E2, -29.774E5, 25.096E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.753E-6, 0.0, 13.431E-6, 0.0 } },
        
    /* andalusite - Al2SiO5 = Al2O3 + SiO2. xMELTS has no data for this
       (not modeled by MELTS itself); H/S/V/Cp ported from Theriak-Domino's
       JUN92d.bs Berman (1988) database (../theriak-domino/src/JUN92d.bs),
       cross-validated exactly against xMELTS' own sillimanite entry above
       (same H/S/V/Cp to 5+ sig figs - confirms same underlying Berman
       calibration). Theriak's volume-EOS line uses a different internal
       parameterization than xMELTS' (v1,v2,v3,v4) polynomial that could
       not be safely reverse-mapped (checked against corundum: the
       per-term scale factors were inconsistent) - reuses sillimanite's
       own real EOS terms as a physically-reasonable stand-in, since all
       three Al2SiO5 polymorphs have similar compressibility and this only
       affects a small pressure correction to G. */
    { "and", { 1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -2589972.17, 91.4337, 5.147,
        { 236.47818, -1102.941, -7526810.0, 936442368.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.753E-6, 0.0, 13.431E-6, 0.0 } },
    /* kyanite - Al2SiO5 = Al2O3 + SiO2. Same source/caveat as andalusite
       above (Theriak-Domino JUN92d.bs, Berman 1988; EOS borrowed from
       sillimanite). */
    { "ky", { 1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -2594220.46, 82.4300, 4.412,
        { 262.68478, -2001.407, -1999740.0, -63181880.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.753E-6, 0.0, 13.431E-6, 0.0 } },
    /* rutile - TiO2 */
    { "ru", { 0,0,0,0,0,0,0,1.0,0,0,0,0,0,0,0,0 },
        -944750.0, 50.460, 1.882,
        { 77.84, 0.0, -33.678E5, 40.294E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.454E-6, 0.584E-12, 25.716E-6, 15.409E-10 } },
    /* sphene (titanite) - CaTiSiO5 = CaO + TiO2 + SiO2 */
    { "sph", { 1.0,0,1.0,0,0,0,0,1.0,0,0,0,0,0,0,0,0 },
        -2596652.0, 129.290, 5.565,
        { 234.62, -10.403E2, -51.183E5, 59.146E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.590E-6, 0.0, 25.200E-6, 0.0 } },

    /* --- olivine endmembers (Fo-Fa, "sb-trivial" tier) --- */
    /* forsterite - Mg2SiO4 */
    { "fo", { 1.0,0,0,2.0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -2174420.0, 94.010, 4.366,
        { 238.64, -20.013E2, 0.0, -11.624E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.791E-6, 1.351E-12, 29.464E-6, 88.633E-10 } },
    /* fayalite - Fe2SiO4 */
    { "fa", { 1.0,0,0,0,2.0,0,0,0,0,0,0,0,0,0,0,0 },
        -1479360.0, 150.930, 4.630,
        { 248.93, -19.239E2, 0.0, -13.910E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.730E-6, 0.0, 26.546E-6, 79.482E-10 } },

    /* --- feldspar endmembers ("sb-trivial" tier) --- */
    /* albite - NaAlSi3O8 = 0.5 Na2O + 0.5 Al2O3 + 3 SiO2 (monalbite reference state) */
    { "ab", { 3.0,0.5,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0 },
        -3921618.0, 224.412, 10.083,
        { 393.64, -24.155E2, -78.928E5, 107.064E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.945E-6, 4.861E-12, 26.307E-6, 32.407E-10 } },
    /* anorthite - CaAl2Si2O8 = CaO + Al2O3 + 2 SiO2 */
    { "an", { 2.0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -4228730.0+3.7*4184.0, 200.186+3.7*4184.0/2200.0, 10.075,
        { 439.37, -37.341E2, 0.0, -31.702E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.272E-6, 3.176E-12, 10.918E-6, 41.985E-10 } },
    /* sanidine - KAlSi3O8 = 0.5 K2O + 0.5 Al2O3 + 3 SiO2 (fully-disordered reference state) */
    { "san", { 3.0,0.5,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0 },
        -3970791.0, 214.145, 10.869,
        { 381.37, -19.410E2, -120.373E5, 183.643E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.805E-6, 5.112E-12, 15.145E-6, 54.850E-10 } },

    /* --- biotite endmembers ("sb-trivial" tier) --- */
    /* annite - KFe3AlSi3O10(OH)2 = 0.5 K2O + 0.5 Al2O3 + 3 FeO + 3 SiO2 + H2O */
    { "ann", { 3.0,0.5,0,0,3.0,0.5,0,0,0,0,0,1.0,0,0,0,0 },
        -5142800.0, 420.0, 15.408,
        { 727.208, -47.75040E2, -138.319E5, 211.906E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.6969784E-6, 0.0, 34.4473262E-6, 0.0 } },
    /* phlogopite - KMg3AlSi3O10(OH)2 = 0.5 K2O + 0.5 Al2O3 + 3 MgO + 3 SiO2 + H2O */
    { "phl", { 3.0,0.5,0,3.0,0,0.5,0,0,0,0,0,1.0,0,0,0,0 },
        -6210391.0, 334.346, 14.977,
        { 610.37988, -20.83781E2, -215.33008E5, 284.1040896E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.6969784E-6, 0.0, 34.4473262E-6, 0.0 } },

    /* --- garnet endmembers ("sb-trivial" tier) --- */
    /* almandine - Fe3Al2Si3O12 = 3 FeO + Al2O3 + 3 SiO2 */
    { "alm", { 3.0,1.0,0,0,3.0,0,0,0,0,0,0,0,0,0,0,0 },
        -5267216.0, 340.007, 11.511,
        { 573.96, -14.831E2, -292.920E5, 502.208E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.558E-6, 0.321E-12, 18.613E-6, 74.539E-10 } },
    /* grossular - Ca3Al2Si3O12 = 3 CaO + Al2O3 + 3 SiO2 */
    { "gr", { 3.0,1.0,3.0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -6632859.0, 255.150, 12.538,
        { 573.43, -20.394E2, -188.872E5, 231.931E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.654E-6, 1.635E-12, 18.994E-6, 79.756E-10 } },
    /* pyrope - Mg3Al2Si3O12 = 3 MgO + Al2O3 + 3 SiO2 */
    { "py", { 3.0,1.0,0,3.0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -6286548.0, 266.359, 11.316,
        { 640.72, -45.421E2, -47.019E5, 0.0E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.576E-6, 0.442E-12, 22.519E-6, 37.044E-10 } },

    /* --- additional common MELTS pure phases --- */
    /* perovskite - CaTiO3 = CaO + TiO2 */
    { "perov", { 0,0,1.0,0,0,0,0,1.0,0,0,0,0,0,0,0,0 },
        -1660630.0, 93.64, 3.3626,
        { 150.49, -6.213E2, 0.0, -43.010E7, 1530.0, 550.0*4.184, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0 } },
    /* calcite - CaCO3 = CaO + CO2 */
    { "cc", { 0,0,1.0,0,0,0,0,0,0,0,0,0,1.0,0,0,0 },
        -1206819.0, 91.725, 3.690,
        { 178.19, -16.577E2, -4.827E5, 16.660E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.4E-6, 0.0, 8.907E-6, 227.402E-10 } },
    /* aragonite - CaCO3 = CaO + CO2 (polymorph of calcite) */
    { "arag", { 0,0,1.0,0,0,0,0,0,0,0,0,0,1.0,0,0,0 },
        -1206819.0+1100.0, 88.0, 3.415,
        { 166.62, -14.994E2, 0.0, 5.449E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.4E-6, 0.0, 8.907E-6, 227.402E-10 } },
    /* magnesite - MgCO3 = MgO + CO2 */
    { "mgs", { 0,0,0,1.0,0,0,0,0,0,0,0,0,1.0,0,0,0 },
        -1113636.0, 65.210, 2.803,
        { 162.30, -11.093E2, -48.826E5, 87.466E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.890E-6, 2.212E-12, 18.436E-6, 415.968E-10 } },
    /* siderite - FeCO3 = FeO + CO2 */
    { "sid", { 0,0,0,0,1.0,0,0,0,0,0,0,0,1.0,0,0,0 },
        -755900.0, 95.5, 2.938,
        { 177.36, -16.694E2, -3.551E5, 15.078E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.890E-6, 2.212E-12, 18.436E-6, 415.968E-10 } },
    /* dolomite - CaMg(CO3)2 = CaO + MgO + 2 CO2 */
    { "dol", { 0,0,1.0,1.0,0,0,0,0,0,0,0,0,2.0,0,0,0 },
        -2324500.0+1100.0, 155.2, 6.434,
        { 368.02, -37.508E2, 0.0, 18.079E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.4E-6, 0.0, 8.907E-6, 227.402E-10 } },
    /* spurrite - Ca5Si2O8CO3 = 2 SiO2 + 5 CaO + CO2 */
    { "spu", { 2.0,0,5.0,0,0,0,0,0,0,0,0,0,1.0,0,0,0 },
        -5840200.0, 331.0, 14.712,
        { 597.163, -36.929E2, -50.5712E5, 43.382E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.890E-6, 2.212E-12, 18.436E-6, 415.968E-10 } },
    /* tilleyite - Ca5Si2O7(CO3)2 = 2 SiO2 + 5 CaO + 2 CO2, companion to spurrite
       (one more CO2, one less "framework" O) - EOS left to parallel calcite/
       spurrite (Berman 1988), same convention xMELTS' own entry uses. */
    { "til", { 2.0,0,5.0,0,0,0,0,0,0,0,0,0,2.0,0,0,0 },
        -6372200.0, 394.0, 17.43,
        { 716.789, -51.992E2, -50.5712E5, 60.769E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.890E-6, 2.212E-12, 18.436E-6, 415.968E-10 } },
    /* muscovite - KAl2(AlSi3O10)(OH)2 = 3 SiO2 + 1.5 Al2O3 + 0.5 K2O + H2O */
    { "mu", { 3.0,1.5,0,0,0,0.5,0,0,0,0,0,1.0,0,0,0,0 },
        -5976740.0, 293.157, 14.087,
        { 651.49, -38.732E2, -185.232E5, 274.247E7, 0.0, 0.0, 0.0, 0.0 },
        { -1.717E-6, 4.295E-12, 33.527E-6, 0.0 } },
    /* aegirine (acmite) - NaFeSi2O6 (Fe3+) = 2 SiO2 + 0.5 Na2O + FeO + 0.5 O,
       Sack & Ghiorso (1994). Cp is xMELTS' own composite (jadeite + hedenbergite
       - diopside reaction sum, "311.29+303.909-297.499" etc in the source) -
       ported as the already-summed numeric result, not re-derived. */
    { "aeg", { 2.0,0,0,0,1.0,0,0.5,0,0.5,0,0,0,0,0,0,0 },
        -2576800.0, 170.57, 63.997,
        { 317.700, -2066.81, -3004650.0, 256011000.0, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0 } },
    /* aenigmatite - Na2Fe5TiSi6O20 (all Fe2+, charge-balances exactly with no
       Fe3+/O-excess needed) = Na2O + 5 FeO + TiO2 + 6 SiO2. H is xMELTS' own
       "Guess" (no calorimetric measurement exists), S from sum-of-oxides. */
    { "aen", { 6.0,0,0,0,5.0,0,1.0,1.0,0,0,0,0,0,0,0,0 },
        -8472805.0, 740.23, 22.8546,
        { 1092.073, -6159.81, -22526790.0, 3268659000.0, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0 } },

    /* --- hornblende endmembers (pargasite-ferropargasite-magnesiohastingsite) --- */
    /* pargasite - NaCa2Mg4AlAl2Si6O22(OH)2 = 6 SiO2 + 1.5 Al2O3 + 2 CaO + 4 MgO + 0.5 Na2O + H2O */
    { "parg", { 6.0,1.5,2.0,4.0,0,0,0.5,0,0,0,0,1.0,0,0,0,0 },
        -12621554.82, 669.44, 27.35,
        { 1267.25, -6654.34, -30378700.0, 3913530000.0, 0.0, 0.0, 0.0, 0.0 },
        { -1.392E-6, 3.481E-12, 24.374E-6, 98.338E-10 } },
    /* ferropargasite - NaCa2Fe4AlAl2Si6O22(OH)2 = 6 SiO2 + 1.5 Al2O3 + 2 CaO + 4 FeO + 0.5 Na2O + H2O */
    { "fparg", { 6.0,1.5,2.0,0,4.0,0,0.5,0,0,0,0,1.0,0,0,0,0 },
        -11224859.39, 776.132, 27.989,
        { 1342.61, -8348.62, -24760400.0, 3485070000.0, 0.0, 0.0, 0.0, 0.0 },
        { -1.392E-6, 3.481E-12, 24.374E-6, 98.338E-10 } },
    /* magnesiohastingsite - NaCa2Mg4FeAl2Si6O22(OH)2 (Fe3+) = 6 SiO2 + Al2O3 + 2 CaO + 4 MgO + FeO + 0.5 O + 0.5 Na2O + H2O */
    { "mhst", { 6.0,1.0,2.0,4.0,1.0,0,0.5,0,0.5,0,0,1.0,0,0,0,0 },
        -12166501.59, 685.3392, 27.38,
        { 1273.66, -6716.06, -28033100.0, 3506970000.0, 0.0, 0.0, 0.0, 0.0 },
        { -1.392E-6, 3.481E-12, 24.374E-6, 98.338E-10 } },

    /* --- leucite endmembers (leucite-analcime-na-leucite) --- */
    /* leucite - KAlSi2O6 = 2 SiO2 + 0.5 Al2O3 + 0.5 K2O */
    { "lc", { 2.0,0.5,0,0,0,0.5,0,0,0,0,0,0,0,0,0,0 },
        -3010279.63, 205.846, 8.8390,
        { 271.14, -944.1, -7857200.0, 959200000.0, 955.0, 256.0, -0.09731, 0.0003373 },
        { -1.5568254E-6, 1.2739E-12, 12.5167877E-6, 0.0 } },
    /* analcime - NaAlSi2O5(OH)2 = 2 SiO2 + 0.5 Al2O3 + 0.5 Na2O + H2O */
    { "anl", { 2.0,0.5,0,0,0,0,0.5,0,0,0,0,1.0,0,0,0,0 },
        -3325900.0, 228.10, 9.71,
        { 571.83, -7188.7, 0.0, 1493060000.0, 0.0, 0.0, 0.0, 0.0 },
        { -1.5568254E-6, 1.2739E-12, 12.5167877E-6, 0.0 } },
    /* na-leucite - NaAlSi2O6 = 2 SiO2 + 0.5 Al2O3 + 0.5 Na2O */
    { "nlc", { 2.0,0.5,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0 },
        -3001975.55, 165.74, 8.91,
        { 401.27, -4248.0, 0.0, 216300000.0, 0.0, 0.0, 0.0, 0.0 },
        { -1.5568254E-6, 1.2739E-12, 12.5167877E-6, 0.0 } },

    /* --- melilite endmembers (akermanite-gehlenite-iron_akermanite-soda_melilite) ---
       iron-akermanite and soda-melilite are reciprocal endmembers in real xMELTS
       (gibbs.c): iron-akermanite = "ak0" (own-row Ca2MgSi2O7-shaped energy, NOT
       the same H298 as "ak" below) + 0.5*fa - 0.5*fo; soda-melilite = "geh0"
       (own-row Ca2Al2SiO7-shaped energy, NOT the same H298 as "geh" below) +
       2*ab - 2*an. fo/fa/ab/an already exist in this table (ol/fsp) with H/S/V/Cp/EOS
       verified byte-identical to what gibbs.c uses for this same reaction, so the
       reaction is applied at RUNTIME in G_SS_gh_mel_function (gh_gss_function.c),
       reusing those 4 existing get_em_data lookups rather than duplicating them. */
    /* akermanite - Ca2MgSi2O7 = 2 SiO2 + 2 CaO + MgO */
    { "ak", { 2.0,0,2.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -3855541.0, 212.000, 9.252,
        { 387.06, -2938.8, 0.0, -40790000.0, 358.0, 452.0, 0.0, 0.0 },
        { -0.785E-6, 0.0, 25.011E-6, 67.224E-10 } },
    /* gehlenite - Ca2Al2SiO7 = 1 SiO2 + 1 Al2O3 + 2 CaO */
    { "geh", { 1.0,1.0,2.0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -3969258.0, 198.600, 9.033,
        { 373.09, -2276.8, -4778500.0, 477910000.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.996E-6, 2.488E-12, 24.926E-6, 5.664E-10 } },
    /* iron-akermanite "own row" (Ca2MgSi2O7-shaped, pre-reaction) - same Comp/Cp/EOS as
       "ak" above, distinct H298 (xMELTS' own two separate calibration constants) */
    { "ak0", { 2.0,0,2.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -3856441.0, 212.000, 9.252,
        { 387.06, -2938.8, 0.0, -40790000.0, 358.0, 452.0, 0.0, 0.0 },
        { -0.785E-6, 0.0, 25.011E-6, 67.224E-10 } },
    /* soda-melilite "own row" (Ca2Al2SiO7-shaped, pre-reaction) - same Comp/Cp/EOS as
       "geh" above, distinct H298 */
    { "geh0", { 1.0,1.0,2.0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -3981458.0, 198.600, 9.033,
        { 373.09, -2276.8, -4778500.0, 477910000.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.996E-6, 2.488E-12, 24.926E-6, 5.664E-10 } },

    /* --- cummingtonite endmembers (cummingtonite-grunerite, Fe-Mg amphibole) --- */
    /* cummingtonite - Mg7Si8O22(OH)2 = 8 SiO2 + 7 MgO + H2O */
    { "cumm", { 8.0,0,0,7.0,0,0,0,0,0,0,0,1.0,0,0,0,0 },
        -12067517.0, 540.2587, 26.33,
        { 1233.792, -7133.980, -22163800.0, 2333937490.0, 0.0, 0.0, 0.0, 0.0 },
        { -1.1393847E-6, 0.0, 28.1048234E-6, 62.894E-10 } },
    /* grunerite - Fe7Si8O22(OH)2 = 8 SiO2 + 7 FeO + H2O */
    { "grun", { 8.0,0,0,0,7.0,0,0,0,0,0,0,1.0,0,0,0,0 },
        -9623300.0, 725.0, 27.840,
        { 1347.83, -9356.91, -20228480.0, 3039190000.0, 0.0, 0.0, 0.0, 0.0 },
        { -1.670299E-6, 8.68919E-12, 28.400E-6, 0.0 } },

    /* --- spinel endmembers (chromite-hercynite-magnetite-spinel-ulvospinel) ---
       Fe3+ tracked via gh's "1 FeO + 0.5 O per mole Fe3+" convention (matching
       hb's mhst): magnetite Fe3O4 = 1 Fe2+ + 2 Fe3+ -> 3 FeO + 1.0 O. */
    /* chromite - FeCr2O4 = FeO + Cr2O3 */
    { "chr", { 0,0,0,0,1.0,0,0,0,0,0,1.0,0,0,0,0,0 },
        -1445490.0, 142.676, 4.4010,
        { 236.874, -1679.6, 0.0, -167650000.0, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0 } },
    /* hercynite - FeAl2O4 = FeO + Al2O3 */
    { "herc", { 0,1.0,0,0,1.0,0,0,0,0,0,0,0,0,0,0,0 },
        -1947681.0, 115.362, 4.075018,
        { 235.190, -1437.0, -4691300.0, 645640000.0, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0 } },
    /* magnetite - Fe3O4 = 3 FeO + 1.0 O (2 of the 3 Fe are Fe3+) */
    { "mt", { 0,0,0,0,3.0,0,0,0,1.0,0,0,0,0,0,0,0 },
        -1117403.0, 146.114, 4.452,
        { 207.93, 0.0, -7243300.0, 664360000.0, 848.0, 1565.0, -0.19502, 0.00061037 },
        { -0.582E-6, 1.751E-12, 30.291E-6, 138.470E-10 } },
    /* spinel (sensu stricto) - MgAl2O4 = MgO + Al2O3 */
    { "spl", { 0,1.0,0,1.0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -2300313.0, 84.535, 3.977,
        { 235.90, -1766.6, -1710400.0, 40620000.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.489E-6, 0.0, 21.691E-6, 50.528E-10 } },
    /* ulvospinel - Fe2TiO4 = 2 FeO + TiO2 */
    { "usp", { 0,0,0,0,2.0,0,0,1.0,0,0,0,0,0,0,0,0 },
        -1488500.0, 185.447, 4.682,
        { 249.63, -1817.4, 0.0, -54530000.0, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0 } },

    /* --- clinopyroxene/orthopyroxene endmembers (Di-Cen-Hed-CaTs(Al)-
       CaTs(Fe3+)-Ess-Jd, Sack & Ghiorso 1993) - from xMELTS'
       includes/sol_struct_data.h, same 7 endmembers shared by both the
       "clinopyroxene" and "orthopyroxene" registry entries (see
       obj_gh_cpx's header comment in gh_objective_functions.c for why
       both phases use the identical energetics and endmember set).
       Fe3+ tracked via gh's usual "1 FeO + 0.5 O per mole Fe3+" convention
       (buffonite, essenite). alumino-buffonite/buffonite standard states
       include xMELTS' own "Guess" H/S offset (no calorimetric data exists
       for these rare accessory components) - a normal, accepted practice
       for minor endmembers in an otherwise real, calibrated (Berman 1988/
       Cosca & Peacor 1987) solution model, not a red flag like ortho-
       oxide's fully-guessed W-parameters were. */
    /* diopside - CaMgSi2O6 = 2 SiO2 + CaO + MgO */
    { "di", { 2.0,0,1.0,1.0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -3200583.0, 142.5, 6.620,
        { 305.41, -1604.9, -7166000.0, 921840000.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.872E-6, 1.707E-12, 27.795E-6, 83.082E-10 } },
    /* clinoenstatite - Mg2Si2O6 = 2 SiO2 + 2 MgO */
    { "cen", { 2.0,0,0,2.0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -3086083.0, 135.164, 6.3279,
        { 333.16, -2401.2, -4541200.0, 558300000.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.749E-6, 0.447E-12, 24.656E-6, 74.670E-10 } },
    /* hedenbergite - CaFeSi2O6 = 2 SiO2 + CaO + FeO */
    { "hed", { 2.0,0,1.0,0,1.0,0,0,0,0,0,0,0,0,0,0,0 },
        -2842221.0, 174.2, 6.7894,
        { 307.89, -1597.3, -6992500.0, 935220000.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.9925E-6, 1.4835E-12, 31.371E-6, 83.672E-10 } },
    /* alumino-buffonite (Ca-Ti-Al Tschermak) - CaTi0.5Mg0.5AlSiO6 =
       SiO2 + 0.5 Al2O3 + CaO + 0.5 MgO + 0.5 TiO2. H/S ref include
       xMELTS' own "Guess" offset (H31=-8565.18 J, S31=0) - no
       calorimetric data for this rare component (see header comment). */
    { "cats", { 1.0,0.5,1.0,0.5,0,0,0,0.5,0,0,0,0,0,0,0,0 },
        -3283830.18, 143.745, 6.356,
        { 297.499, -1355.96, -6702190.0, 759082000.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.870E-6, 2.171E-12, 22.250E-6, 52.863E-10 } },
    /* buffonite (Ca-Ti-Fe3+ Tschermak) - CaTi0.5Mg0.5FeSiO6 = SiO2 +
       CaO + 0.5 MgO + 0.5 TiO2 + 1 Fe3+ (-> FeO + 0.5 O). H/S ref
       include xMELTS' own "Guess" offset (H41=+7932.05 J,
       S41=+18.756 J/K) - no calorimetric data for this rare component. */
    { "buff", { 1.0,0,1.0,0.5,1.0,0,0,0.5,0.5,0,0,0,0,0,0,0 },
        -2828776.95, 180.741, 6.7217,
        { 303.909, -1417.67, -4356540.0, 352523000.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.870E-6, 2.171E-12, 22.250E-6, 52.863E-10 } },
    /* essenite - CaFeAlSiO6 = SiO2 + 0.5 Al2O3 + CaO + 1 Fe3+
       (-> FeO + 0.5 O). H/S ref include xMELTS' own offset relative to
       CaTs (Berman 1988) - real Cosca & Peacor (1987) volume, not
       guessed. */
    { "ess", { 1.0,0.5,1.0,0,1.0,0,0,0,0.5,0,0,0,0,0,0,0 },
        -2867478.01, 177.747, 6.7217,
        { 317.11, -1733.3, -5109650.0, 542220000.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.870E-6, 2.171E-12, 22.250E-6, 52.863E-10 } },
    /* jadeite - NaAlSi2O6 = 2 SiO2 + 0.5 Al2O3 + 0.5 Na2O */
    { "jd", { 2.0,0.5,0,0,0,0,0.5,0,0,0,0,0,0,0,0,0 },
        -3025118.0, 133.574, 6.034,
        { 311.29, -2005.1, -5350300.0, 662570000.0, 0.0, 0.0, 0.0, 0.0 },
        { -0.860E-6, 2.149E-12, 23.118E-6, 25.785E-10 } },

    /* --- rhombohedral oxide endmembers (geikielite-hematite-ilmenite-
       pyrophanite-corundum, Ghiorso & Evans 2008) - from xMELTS'
       includes/sol_struct_data.h (#ifdef RHYOLITE_ADJUSTMENTS block).
       Each of these 5 gets an additional intracrystalline-ordering (or,
       for hematite/corundum, short-range-order) correction added at
       runtime in GH_gem_function.c (GH_ilm_gei_pyr_ordering_G /
       GH_rhm_SRO_G), mirroring albite/sanidine's ordering corrections
       above - NOT baked into H/S here. "crn" is a distinct entry from
       "cor" above: xMELTS gives this Al2O3 reference state used inside
       the rhm-oxide solution a +20000 J offset relative to the standalone
       corundum phase ("cor"), so it must not be merged with it. */
    /* geikielite - MgTiO3 = MgO + TiO2 */
    { "gei", { 0,0,0,1.0,0,0,0,1.0,0,0,0,0,0,0,0,0 },
        -1572560.0, 74.56, 3.086,
        { 146.20, -4.160E2, -39.998E5, 40.233E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.584e-6, 1.230e-12, 27.248e-6, 29.968e-10 } },
    /* hematite (rhm-oxide's own entry) - Fe2O3 = 2 FeO + 1.0 O (2 Fe3+) */
    { "hem", { 0,0,0,0,2.0,0,0,0,1.0,0,0,0,0,0,0,0 },
        -825627.0, 87.437, 3.027,
        { 146.86, 0.0, -55.768E5, 52.563E7, 955.0, 1287.0, -7.403E-2, 27.921E-5 },
        { -0.479e-6, 0.304e-12, 38.310e-6, 1.650e-10 } },
    /* ilmenite - FeTiO3 = FeO + TiO2 */
    { "ilm", { 0,0,0,0,1.0,0,0,1.0,0,0,0,0,0,0,0,0 },
        -1231947.0, 108.628, 3.170,
        { 150.00, -4.416E2, -33.237E5, 34.815E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.584e-6, 1.230e-12, 27.248e-6, 29.968e-10 } },
    /* pyrophanite - MnTiO3 = MnO + TiO2 */
    { "pyr", { 0,0,0,0,0,0,0,1.0,0,1.0,0,0,0,0,0,0 },
        -1350707.0, 104.935, 2.8859,
        { 150.00, -4.416E2, -33.237E5, 34.815E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.584e-6, 1.230e-12, 27.248e-6, 29.968e-10 } },
    /* corundum (rhm-oxide's own entry, +20000 J vs standalone "cor") - Al2O3 */
    { "crn", { 0,1.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 },
        -1675700.0+20000.0, 50.820, 2.558,
        { 155.02, -8.284E2, -38.614E5, 40.908E7, 0.0, 0.0, 0.0, 0.0 },
        { -0.385E-6, 0.375E-12, 21.342E-6, 47.180E-10 } },

    /* na-nepheline - Na4Al4Si4O16 = 2 Na2O + 2 Al2O3 + 4 SiO2 */
    { "nane", { 4.0,2.0,0,0,0,0,2.0,0,0,0,0,0,0,0,0,0 },
        -2093004.0*4.0, 124.641*4.0, 5.4131*4.0,
        { 205.24*4.0, -7.599E2*4.0, -108.383E5*4.0, 208.182E7*4.0, 467.0, 241.0*4.0, -50.249E-2*2.0, 165.95E-5*2.0 },
        { -2.0500e-6, 5.2000e-12, 31.802e-6, 213.0e-10 } },
    /* k-nepheline (nepheline-hosted) - K4Al4Si4O16 = 2 K2O + 2 Al2O3 + 4 SiO2.
       H ref includes the real xMELTS -3572.76 J nepheline-only offset (on
       top of the shared -(DH2) term, DH2=-1.35*1000.0*4.184=-5648.4 J). */
    { "kne", { 4.0,2.0,0,0,0,2.0,0,0,0,0,0,0,0,0,0,0 },
        -2109563.55*4.0-(-5648.4)-3572.76, 133.9653*4.0-(0.0), 6.043478*4.0-(-0.00001*1000.0*4.184),
        { 186.0*4.0, 0.0, -131.067E5*4.0, 213.893E7*4.0, 800.15, 1154.0*4.0, -7.096454E-2*2.0, 21.682E-5*2.0 },
        { -2.0500e-6, 5.2000e-12, 31.802e-6, 213.0e-10 } },
    /* k-nepheline (kalsilite-hosted) - same as "kne" but WITHOUT the
       -3572.76 nepheline-only offset (real xMELTS: two distinct
       registrations of the "same" endmember, see header comment above). */
    { "knk", { 4.0,2.0,0,0,0,2.0,0,0,0,0,0,0,0,0,0,0 },
        -2109563.55*4.0-(-5648.4), 133.9653*4.0-(0.0), 6.043478*4.0-(-0.00001*1000.0*4.184),
        { 186.0*4.0, 0.0, -131.067E5*4.0, 213.893E7*4.0, 800.15, 1154.0*4.0, -7.096454E-2*2.0, 21.682E-5*2.0 },
        { -2.0500e-6, 5.2000e-12, 31.802e-6, 213.0e-10 } },
    /* vc-nepheline - Na3Al3Si5O16 = 1.5 Na2O + 1.5 Al2O3 + 5 SiO2. H/S/Cp
       here are real xMELTS' own table values (= high-albite's Berman 1988
       calibration, "ab"'s H/S/Cp above); V=0/EOS=0 deliberately - the
       Vinet-EOS phantom-reaction correction is added on top by name in
       GH_gem_function.c, not folded into this table row. */
    { "vcne", { 5.0,1.5,0,0,0,0,1.5,0,0,0,0,0,0,0,0,0 },
        -3921618.0, 224.412, 0.0,
        { 393.64, -24.155E2, -78.928E5, 107.064E7, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0 } },
    /* ca-nepheline - CaNa2Al4Si4O16 = 1 CaO + 1 Na2O + 2 Al2O3 + 4 SiO2.
       H/S/Cp here are high-anorthite's Berman 1988 calibration ("an"'s H/S/
       Cp above); V=0/EOS=0 deliberately, same phantom-reaction convention
       as "vcne" (GH_caneph_G in GH_gem_function.c adds the +23096 J / zero-
       point-entropy / Vinet-EOS na-nepheline correction on top). */
    { "cane", { 4.0,2.0,1.0,0,0,0,1.0,0,0,0,0,0,0,0,0,0 },
        -4228730.0+3.7*4184.0, 200.186+3.7*4184.0/2200.0, 0.0,
        { 439.37, -37.341E2, 0.0, -31.702E7, 0.0, 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0, 0.0 } },
};

/**
    Beta (high-T) polymorph reference states for the 3 SiO2 phases whose
    alpha-form data is in arr_pp_db_gh[] above (Berman 1988), used by
    GH_SiO2_polymorph_G() (GH_gem_function.c) when T exceeds the
    pressure-shifted transition temperature Tt+dTdP*(P-1).
**/
static const PP_db_gh_beta arr_pp_db_gh_beta[3] = {
    /* beta-quartz */
    { "q",    -908627.0, 44.207, 2.370, { -1.238E-6,  7.087E-13, 0.0,     0.0 }, 0.0237 },
    /* beta-cristobalite */
    { "crst", -906377.0, 46.029, 2.730, { -1.100E-6,   5.535E-12, 3.189E-6, 0.0 }, 0.0480 },
    /* beta-tridymite (no pressure shift on Tt: dTdP=0) */
    { "trd",  -907045.0, 45.524, 2.737, { -0.740E-6,   3.735E-12, 4.829E-6, 0.0 }, 0.0 },
};

PP_db_gh Access_GH_PP_DB(int id){
    return arr_pp_db_gh[id];
}

int GH_find_PP_id(char *name){
    for (int i = 0; i < GH_N_PP; i++){
        if (strcmp(arr_pp_db_gh[i].Name, name) == 0){
            return i;
        }
    }
    return -1;
}

PP_db_gh_beta Access_GH_SiO2_beta_DB(int id){
    return arr_pp_db_gh_beta[id];
}

int GH_find_SiO2_beta_id(char *name){
    for (int i = 0; i < 3; i++){
        if (strcmp(arr_pp_db_gh_beta[i].Name, name) == 0){
            return i;
        }
    }
    return -1;
}
