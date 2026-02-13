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
#ifndef __ENDMEMBER_DATABASE_SB_H_
#define __ENDMEMBER_DATABASE_SB_H_

    /** store endmember database **/
    typedef struct EM_db_sb_ {
        char   Name[50];			/** pure species name 														*/
        char   FullName[80];		/** pure species name 														*/
        char   Equation[90];		/** pure species name 														*/
        double Comp[17];       	 	/** pure species composition [0-16] + number of atom [17] 					*/
        double input_1[10];         /** second line of the thermodynamics datable 								*/
        double input_2[3];          /** second line of the thermodynamics datable 								*/
    } EM_db_sb;

    EM_db_sb Access_SB_EM_DB(int id, int EM_dataset);

#endif