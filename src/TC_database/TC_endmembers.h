#ifndef __ENDMEMBER_DATABASE_TC_H_
#define __ENDMEMBER_DATABASE_TC_H_

    /** store endmember database **/
    struct EM_db {
        char   Name[20];			/** pure species name 														*/
        double Comp[16];       	 	/** pure species composition [0-10] + number of atom [11] 					*/
        double input_1[3];          /** first line of the thermodynamics datable 								*/
        double input_2[4];          /** second line of the thermodynamics datable 								*/
        double input_3[11];         /** third line of the thermodynamics datable 								*/
        double input_4[3];         	/** third line of the thermodynamics datable 								*/
    };

    /** store endmember database **/
    struct FS_db {
        char   Name[20];			/** pure species name 														*/
        double Comp[16];       	 	/** pure species composition [0-10] + number of atom [11] 					*/
        double input_1[4];          /** first line of the thermodynamics datable 								*/
        double input_2[7];          /** second line of the thermodynamics datable 								*/
        double input_3[1];         	/** third line of the thermodynamics datable 								*/
    };

    struct EM_db Access_EM_DB(int id, int EM_database);

    struct FS_db Access_FS_DB(int id);

#endif