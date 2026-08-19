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
#ifndef __ENDMEMBER_DATABASE_TC_H_
#define __ENDMEMBER_DATABASE_TC_H_


    /** store endmember database **/
    typedef struct EM_db_ {
        char   Name[20];			/** pure species name 														*/
        double Comp[16];       	 	/** pure species composition [0-10] + number of atom [11] 					*/
        double input_1[3];          /** first line of the thermodynamics datable 								*/
        double input_2[4];          /** second line of the thermodynamics datable 								*/
        double input_3[11];         /** third line of the thermodynamics datable 								*/
        double input_4[3];         	/** third line of the thermodynamics datable 								*/
    } EM_db;

    /**
        Store DEW2019 aqueous species database (Sverjensky et al., Deep Earth Water model).
        Universal solvent-EOS constants shared by every species (T_r=298.15K, P_r=1bar,
        Psi=2600bar, eta=694657 J*Angstrom/mol, theta=228K) are NOT stored here - they are
        hardcoded in G_DEW_function (TC_gem_function.c).
    **/
    typedef struct DEW_db_ {
        char   Name[20];			/** aqueous species name (DEW2019 notation)								*/
        double Comp[16];       	 	/** oxide-basis composition (O-corrected), canonical oxide order [0-14]; [15] unused	*/
        double MuComp[16];          /** formation-reaction stoichiometry vs oxide components, same canonical order [0-14]; [15] unused */
        double input_1[4];          /** G_ref, S_ref, a1, a2													*/
        double input_2[4];          /** a3, a4, c1, c2															*/
        double input_3[3];          /** rH, omega0, z (charge)													*/
        double input_4[1];          /** mu_Hp - H+ stoichiometric coefficient of the formation reaction		*/
    } DEW_db;

    EM_db Access_EM_DB(int id, int EM_dataset);

    DEW_db Access_DEW_DB(int id);

    // extern EM_db arr_em_db_tc_ds62[257];
    // extern EM_db arr_em_db_tc_ds633[289];
    // extern EM_db arr_em_db_tc_ds634[291];
    // extern EM_db arr_em_db_tc_ds635[291];
    // extern EM_db arr_em_db_tc_ds636[291];
    EM_db* get_arr_em_db_tc( int EM_dataset);
#endif