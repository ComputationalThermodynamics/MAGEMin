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
#ifndef __HASH_INIT_H_
#define __HASH_INIT_H_

    #include "uthash.h"

    /*---------------------------------------------------------------------------*/ 
    /*  Hashtable for endmember in thermodynamic database                        */
    typedef struct EM2id_{
        char EM_tag[20];           /* key (string is WITHIN the structure)       */
        int id;                    /* id of the key (array index)                */
        UT_hash_handle hh;         /* makes this structure hashable              */
    } EM2id;
    EM2id *EM = NULL;

    int find_EM_id(char *EM_tag) {
        EM2id *p_s;
        HASH_FIND_STR ( EM, EM_tag, p_s );
        return(p_s->id);	
    }

    /*  Hashtable for fluid species in thermodynamic database                    */
    typedef struct FS2id_{
        char FS_tag[20];           /* key (string is WITHIN the structure)       */
        int id;                    /* id of the key (array index)                */
        UT_hash_handle hh;         /* makes this structure hashable              */
    } FS2id;
    FS2id *FS = NULL;

    int find_FS_id(char *FS_tag) {
        FS2id *fs_s;
        HASH_FIND_STR ( FS, FS_tag, fs_s );
        return(fs_s->id);	
    }

    /**
        Hashtable for Pure-phases in the pure-phases list                        
    */
    typedef struct PP2id_{
        char PP_tag[20];           /* key (string is WITHIN the structure)       */
        int id;                    /* id of the key (array index)                */
        UT_hash_handle hh;         /* makes this structure hashable              */
    } PP2id;
    PP2id *PP = NULL;

    int find_PP_id(char *PP_tag) {
        PP2id *pp_s;
        HASH_FIND_STR ( PP, PP_tag, pp_s );
        return(pp_s->id);	
    }

#endif