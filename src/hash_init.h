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
#ifndef __HASH_INIT_H_
#define __HASH_INIT_H_

    #include <stdatomic.h>
    #include <stdio.h>
    #include <stdlib.h>
    #include <string.h>
    #include "uthash.h"

    /*---------------------------------------------------------------------------*/
    /*  Hashtable for endmember in thermodynamic database                        */
    typedef struct EM2id_{
        char EM_tag[20];           /* key (string is WITHIN the structure)       */
        int id;                    /* id of the key (array index)                */
        UT_hash_handle hh;         /* makes this structure hashable              */
    } EM2id;

    #define n_EM_tables_max 64
    typedef struct EM_tables_{
        char    research_group[20];
        int     EM_dataset;
        int     n_names;
        EM2id  *table;
    } EM_table;
    static EM_table     EM_tables[n_EM_tables_max];
    static atomic_int   n_EM_tables     = 0;
    static atomic_int   EM_tables_lock  = 0;

    static int get_EM_index(char *research_group, int EM_dataset, int n_names) {
        int n = atomic_load_explicit(&n_EM_tables, memory_order_acquire);
        for (int i = 0; i < n; i++){
            if (EM_tables[i].EM_dataset == EM_dataset && (n_names < 0 || EM_tables[i].n_names == n_names) && strcmp(EM_tables[i].research_group, research_group) == 0){
                return i;
            }
        }
        return -1;
    }

    static EM2id *get_EM_table(char *research_group, int EM_dataset, int n_names) {
        int i = get_EM_index(research_group, EM_dataset, n_names);
        return (i < 0) ? NULL : EM_tables[i].table;
    }

    static void warn_tag_truncation(const char *kind, const char *name, size_t cap) {
        if (strlen(name) >= cap){
            fprintf(stderr, "MAGEMin warning: %s name '%s' is %zu characters, longer than the %zu available; it will be truncated to '%.*s' and may collide with another name sharing that prefix\n", kind, name, strlen(name), cap - 1, (int)(cap - 1), name);
        }
    }

    void register_EM_table(char *research_group, int EM_dataset, char **names, int n_names) {
        while (atomic_exchange_explicit(&EM_tables_lock, 1, memory_order_acquire)){}
        int n = atomic_load_explicit(&n_EM_tables, memory_order_relaxed);
        if (get_EM_index(research_group, EM_dataset, n_names) < 0){
            if (n >= n_EM_tables_max){
                fprintf(stderr, "MAGEMin: endmember table registry full (%d entries), cannot register research group '%s', EM_dataset %d, n_em_db %d. Increase n_EM_tables_max in hash_init.h\n", n_EM_tables_max, research_group, EM_dataset, n_names);
                abort();
            }
            EM2id *table = NULL;
            for (int i = 0; i < n_names; i++){
                warn_tag_truncation("endmember", names[i], sizeof(((EM2id *)0)->EM_tag));
                EM2id *p_s = (EM2id *)malloc(sizeof *p_s);
                strncpy(p_s->EM_tag, names[i], sizeof(p_s->EM_tag) - 1);
                p_s->EM_tag[sizeof(p_s->EM_tag) - 1] = '\0';
                p_s->id = i;
                HASH_ADD_STR( table, EM_tag, p_s );
            }
            strncpy(EM_tables[n].research_group, research_group, sizeof(EM_tables[n].research_group) - 1);
            EM_tables[n].research_group[sizeof(EM_tables[n].research_group) - 1] = '\0';
            EM_tables[n].EM_dataset = EM_dataset;
            EM_tables[n].n_names    = n_names;
            EM_tables[n].table      = table;
            atomic_store_explicit(&n_EM_tables, n + 1, memory_order_release);
        }
        atomic_store_explicit(&EM_tables_lock, 0, memory_order_release);
    }

    int find_EM_id(char *research_group, int EM_dataset, char *EM_tag) {
        EM2id *table = get_EM_table(research_group, EM_dataset, -1);
        EM2id *p_s   = NULL;
        HASH_FIND_STR ( table, EM_tag, p_s );
        return (p_s == NULL) ? -1 : p_s->id;
    }

    void check_lookup_id(int id, const char *kind, const char *name) {
        if (id < 0){
            fprintf(stderr, "\nMAGEMin fatal: %s '%s' was not found in the registered name table (lookup returned %d).\n", kind, name ? name : "(null)", id);
            fprintf(stderr, "Refusing to index the thermodynamic database with an invalid id.\n");
            abort();
        }
    }

    /*  Hashtable for DEW2019 aqueous species in thermodynamic database          */
    typedef struct DEW2id_{
        char DEW_tag[20];          /* key (string is WITHIN the structure)       */
        int id;                    /* id of the key (array index)                */
        UT_hash_handle hh;         /* makes this structure hashable              */
    } DEW2id;
    static DEW2id      *DEW             = NULL;
    static atomic_int   DEW_ready       = 0;

    void register_DEW_table(char **names, int n_names) {
        while (atomic_exchange_explicit(&EM_tables_lock, 1, memory_order_acquire)){}
        if (atomic_load_explicit(&DEW_ready, memory_order_relaxed) == 0){
            DEW2id *table = NULL;
            for (int i = 0; i < n_names; i++){
                warn_tag_truncation("DEW species", names[i], sizeof(((DEW2id *)0)->DEW_tag));
                DEW2id *dew_s = (DEW2id *)malloc(sizeof *dew_s);
                strncpy(dew_s->DEW_tag, names[i], sizeof(dew_s->DEW_tag) - 1);
                dew_s->DEW_tag[sizeof(dew_s->DEW_tag) - 1] = '\0';
                dew_s->id = i;
                HASH_ADD_STR( table, DEW_tag, dew_s );
            }
            DEW = table;
            atomic_store_explicit(&DEW_ready, 1, memory_order_release);
        }
        atomic_store_explicit(&EM_tables_lock, 0, memory_order_release);
    }

    int find_DEW_id(char *DEW_tag) {
        DEW2id *dew_s = NULL;
        if (atomic_load_explicit(&DEW_ready, memory_order_acquire) == 1){
            HASH_FIND_STR ( DEW, DEW_tag, dew_s );
        }
        return (dew_s == NULL) ? -1 : dew_s->id;
    }

#endif
