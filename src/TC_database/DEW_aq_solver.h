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
#ifndef __DEW_AQ_SOLVER_H_
#define __DEW_AQ_SOLVER_H_

#include "TC_endmembers.h"
#include "TC_gem_function.h"

#define DEW_N_SPECIES_DB 120


typedef struct AQ_data_t {
    char     name[20];
    int      n_sp;                  /** number of active aqueous species (excludes water)   */
    int      len_ox;                /** number of active oxide components                   */
    int      id_H2O;                /** position of H2O within the active oxide list, or -1 */

    double **em_comp;               /** [n_sp+1][len_ox] oxide-basis composition            */
    char   **em_names;               /** [n_sp]                                              */
    double  *gbase;                  /** [n_sp] kJ/mol, HSC convention (AQ_ref.gbase)        */

    double **mu_comp;                /** [n_sp+1][len_ox+1] formation-reaction stoichiometry vs oxide components + H+ */
    double  *z;                      /** [n_sp] charge                                       */

    double   rho_w;                  /** g/cm3                                               */
    double   gb_w;                   /** kJ/mol, water standard-state G                      */
} AQ_data;

/**
    Iterative-solver scratch + output state, sized to one AQ_data's n_sp. Mirrors
    MAGEMin.jl's local_optimizer_aq_ + the relevant fields of AQ_act_. Every field that
    is a bare scalar in Julia (wrapped in a length-1 Vector there, a mutation idiom not
    needed in C) is a plain double here.
**/
typedef struct AQ_solver_t {
    double  *m;                     /** [n_sp] molality, current iterate                    */
    double   sum_m;
    double   sum_charged_p;
    double   sum_charged_m;

    double   Gamma;                  /** Setschenow-type water term                          */
    double   gamma_e[5];             /** Debye-Hückel activity coefficient, indexed by |z| (0..4) */
    double   I_str;                  /** ionic strength                                      */
    double   Lambda;
    double   sigma;

    double   coef;
    double   a_coef;                 /** water activity coefficient                          */
    double   mu_w;                   /** water chemical potential, kJ/mol                    */
    double   mu_Hp;                  /** H+ chemical potential, kJ/mol                       */

    /* outputs, valid after DEW_aq_min_iterative returns */
    double  *m1;                     /** [n_sp] converged molality                           */
    double  *x;                      /** [n_sp+1] mole fraction (index n_sp = water)         */
    double  *a;                      /** [n_sp+1] activity (index n_sp = water)              */
    double  *mu;                     /** [n_sp+1] chemical potential, kJ/mol (index n_sp = water) */
    double   G;                      /** total Gibbs energy of the fluid, kJ/mol             */
    double   z_res;                  /** |sum(z_i*m_i)| at exit - charge-balance residual actually reached, whether or not z_res_tol was met */
    int      converged;              /** 1 if z_res <= z_res_tol at exit, 0 if the stall detector or max_iter cut this start off first */
} AQ_solver;

/* TEMPORARY investigation instrumentation - see DEW_aq_solver.c */
extern _Thread_local long DEW_stat_picard_passes;
extern _Thread_local long DEW_stat_residual_evals;
void DEW_stat_reset(void);

int  DEW_species_active(   DEW_db   species,
                            int      len_ox,
                            int     *id                );

int  DEW_find_id_H2O(      int      len_ox,
                            int     *id                );

int  DEW_canonical_oxide_index(    char    *oxide_name       );

void DEW_build_id_map(     int      len_ox,
                            char   **ox_names,
                            int     *id_out            );

AQ_data  init_DEW_aqueous_model(   int      n_dew_db,
                                    int      len_ox,
                                    int     *id,
                                    double  *ElEntropy,
                                    double   P,
                                    double   T,
                                    double   rho_w,
                                    double   gb_w              );

void DEW_get_water_reference(      int      EM_dataset,
                                    int      len_ox,
                                    int     *id,
                                    double  *bulk_rock,
                                    double  *apo,
                                    double   P,
                                    double   T,
                                    double  *rho_w,
                                    double  *gb_w              );

AQ_data  init_DEW_aqueous_model_at_point(  int      n_dew_db,
                                            int      EM_dataset,
                                            int      len_ox,
                                            int     *id,
                                            double  *bulk_rock,
                                            double  *apo,
                                            double  *ElEntropy,
                                            double   P,
                                            double   T              );

void free_DEW_aqueous_model(       AQ_data *AQ                );

AQ_solver init_DEW_solver(         AQ_data *AQ                );

void free_DEW_solver(              AQ_solver *S               );

void DEW_aq_min_iterative(         AQ_data     *AQ,
                                    AQ_solver   *S,
                                    double      *Gamma_ox,
                                    double       R,
                                    double       T,
                                    double       P,
                                    int          max_iter,
                                    double       z_res_tol,
                                    double       Hp_eps,
                                    const double *x_warm,
                                    int          algorithm      );

void DEW_aq_min_iterative_mixed(   AQ_data     *AQ,
                                    AQ_solver   *S,
                                    double      *Gamma_ox,
                                    double       R,
                                    double       T,
                                    double       P,
                                    int          max_iter,
                                    double       z_res_tol,
                                    double       Hp_eps,
                                    const double *x_warm       );

void DEW_aq_min_multistart(        AQ_data     *AQ,
                                    AQ_solver   *S,
                                    double      *Gamma_ox,
                                    double       R,
                                    double       T,
                                    double       P,
                                    int          max_iter,
                                    double       z_res_tol,
                                    int          algorithm      );


int  DEW_aq_min_warmstart(         AQ_data      *AQ,
                                    AQ_solver    *S,
                                    double       *Gamma_ox,
                                    double        R,
                                    double        T,
                                    double        P,
                                    int           max_iter,
                                    double        z_res_tol,
                                    const double *x_warm,
                                    int           algorithm      );

void DEW_aq_evaluate(              AQ_data      *AQ,
                                    const double *x,
                                    double        R,
                                    double        T,
                                    double        P,
                                    double       *mu_out,
                                    double       *G_out        );

void DEW_aq_activity_state(        AQ_data      *AQ,
                                    const double *x,
                                    double        T,
                                    double        P,
                                    double       *m,
                                    double        gamma_e[5],
                                    double       *a_coef_out       );

double DEW_aq_pH(                  AQ_data      *AQ,
                                    const double *x,
                                    double        T,
                                    double        P               );

#endif
