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
#ifndef __PP_ENDMEMBER_DATABASE_BR_H_
#define __PP_ENDMEMBER_DATABASE_BR_H_

    typedef struct PP_db_br_ {
        char   Name[16];
        double Comp[16];
        double H;
        double S;
        double V;
        double cp_berman[11];
        double eos_berman[4];
    } PP_db_br;

    #define BR_N_PP 87

    PP_db_br Access_BR_PP_DB(int id);
    int      BR_find_PP_id(char *name);

#endif
