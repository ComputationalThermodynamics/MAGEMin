/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus, Jamison Assunção
 **   Contributors : Dominguez, H., Green E., Berlie N., and Rummel L.
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#ifndef __ENDMEMBER_DATABASE_SB_H_
#define __ENDMEMBER_DATABASE_SB_H_

    /** store endmember database **/
    typedef struct EM_db_sb_ {
        char   Name[20];			/** pure species name 														*/
        double Comp[16];       	 	/** pure species composition [0-10] + number of atom [11] 					*/
        double input_1[3];          /** first line of the thermodynamics datable 								*/
        double input_2[4];          /** second line of the thermodynamics datable 								*/
        double input_3[11];         /** third line of the thermodynamics datable 								*/
        double input_4[3];         	/** third line of the thermodynamics datable 								*/
    } EM_db_sb;

    EM_db_sb Access_SB_EM_DB(int id, int EM_dataset);

#endif