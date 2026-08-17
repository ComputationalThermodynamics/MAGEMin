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
#ifndef __ENDMEMBER_DATABASE_GH_H_
#define __ENDMEMBER_DATABASE_GH_H_

    /**
        Store standard-state thermodynamic data for a Ghiorso/MELTS liquid
        basis species. Ported from xMELTS includes/liq_struct_data.h
        (xMeltsLiquid[]). The endmember's liquid standard-state Gibbs energy
        is built from a *solid* reference state (H, S at Tr=298.15K/Pr=1bar,
        4-term Berman Cp) integrated up to a fusion temperature, plus a
        fusion correction (Tfus, Sfus), plus the liquid's own constant Cp
        above Tfus, plus a pressure correction from a linear-in-(T,P)
        "Kress" liquid-volume polynomial. See GH_gem_function.c.
    **/
    typedef struct EM_db_gh_ {
        char   Name[16];
        double Comp[16];        /** oxide composition, same 16-oxide axis order/convention as
                                     TC_database/TC_init_database.c's oxide_info table              */
        double H;                /** solid reference enthalpy at Tr,Pr (J)                          */
        double S;                /** solid reference entropy at Tr,Pr (J/K)                         */
        double Vs;               /** solid reference volume at Tr,Pr (J/bar), informational only    */
        double cp_berman[8];     /** Berman Cp coefficients: k0,k1,k2,k3,Tt,deltah,l1,l2 (solid)     */
        double Vl;               /** liquid volume at Tr (J/bar)                                     */
        double kress[4];         /** liquid volume EOS: dV/dT, dV/dP, d2V/dTdP, d2V/dP2 (J/bar units)*/
        double Tfus;             /** fusion temperature (K)                                          */
        double Sfus;             /** entropy of fusion (J/K)                                         */
        double Cpl;              /** liquid heat capacity, constant (J/K)                            */
    } EM_db_gh;

    EM_db_gh Access_GH_EM_DB(int EM_database, int id);

#endif
