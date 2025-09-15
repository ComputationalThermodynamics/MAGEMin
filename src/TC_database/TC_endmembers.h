/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
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

    /** store endmember database **/
    typedef struct FS_db_ {
        char   Name[20];			/** pure species name 														*/
        double Comp[17];       	 	/** pure species composition [0-10] + number of atom [11] + charge[12]	    */
        double input_1[4];          /** first line of the thermodynamics datable 								*/
        double input_2[7];          /** second line of the thermodynamics datable 								*/
        double input_3[1];         	/** third line of the thermodynamics datable 								*/
    } FS_db;

    EM_db Access_EM_DB(int id, int EM_dataset);

    FS_db Access_FS_DB(int id);
    
    // extern EM_db arr_em_db_tc_ds62[257];
    // extern EM_db arr_em_db_tc_ds633[289];
    // extern EM_db arr_em_db_tc_ds634[291];
    // extern EM_db arr_em_db_tc_ds635[291];
    // extern EM_db arr_em_db_tc_ds636[291];
    EM_db* get_arr_em_db_tc( int EM_dataset);
#endif