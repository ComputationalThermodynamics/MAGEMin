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
#ifndef __PP_ENDMEMBER_DATABASE_GH_H_
#define __PP_ENDMEMBER_DATABASE_GH_H_


    typedef struct PP_db_gh_ {
        char   Name[16];
        double Comp[16];            /** oxide composition, same axis order as GH_endmembers.h   */
        double H;                   /** reference enthalpy at Tr,Pr (J)                         */
        double S;                   /** reference entropy at Tr,Pr (J/K)                        */
        double V;                   /** reference volume at Tr,Pr (J/bar)                       */
        double cp_berman[8];        /** Berman Cp coefficients: k0,k1,k2,k3,Tt,deltah,l1,l2     */
        double eos_berman[4];       /** Berman volume EOS coefficients: v1,v2,v3,v4             */
    } PP_db_gh;


    typedef struct PP_db_gh_beta_ {
        char   Name[16];
        double H;
        double S;
        double V;
        double eos_berman[4];
        double dTdP;
    } PP_db_gh_beta;

    PP_db_gh_beta Access_GH_SiO2_beta_DB(int id);
    int           GH_find_SiO2_beta_id(char *name);

    #define GH_N_PP 64

    PP_db_gh Access_GH_PP_DB(int id);
    int      GH_find_PP_id(char *name);

#endif
