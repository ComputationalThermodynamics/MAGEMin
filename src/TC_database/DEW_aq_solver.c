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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "DEW_aq_solver.h"

/* TEMPORARY investigation instrumentation (2026-08-14, low-T DEW convergence-cost
   investigation) - counts Picard outer-loop passes (DEW_aq_min_iterative's own while
   loop) and charge-residual evaluations (DEW_charge_residual, called from both the
   bracket-expansion and bisection phases of the mu_Hp root-find, each internally O(n_sp))
   accumulated across an entire DEW_aq_min_multistart call (all Hp_eps starts summed).
   Not meant to be permanent - reset via DEW_stat_reset() by the caller before each
   multistart call it wants to measure.

   _Thread_local: MAGEMin_wrappers.jl's point_wise_minimization_series runs point batches
   across multiple Julia threads (Threads.@threads :static), each with its own private
   DB/gv/SS_ref_db - but these two counters were plain global C variables, so every thread
   was incrementing the SAME two `long`s concurrently with no synchronization (a genuine
   data race, found 2026-08-14 while investigating a reported segfault in
   swap_PGE_pseudocompounds during a 6876-point --db=mpe --DEW_activate=1 batch run).
   Thread-local storage gives every thread its own independent counters instead, matching
   the per-call, per-thread instrumentation semantics DEW_stat_reset()/the [DEW stat]
   printf were always meant to have. */
_Thread_local long DEW_stat_picard_passes  = 0;
_Thread_local long DEW_stat_residual_evals = 0;
void DEW_stat_reset(void){ DEW_stat_picard_passes = 0; DEW_stat_residual_evals = 0; }

/**
    Whether a DEW species is usable given the current active-oxide set: its
    (O-corrected) oxide-basis composition must be exactly zero on every canonical
    oxide slot the system doesn't track. Mirrors DEW19_import_database.jl's
    `active = issubset(els, gv.cstPointVariables.elements)` check (there performed on
    the species' raw element_formula - equivalent here since Comp[] only carries a
    nonzero coefficient on slots the raw formula actually touches, see
    tools/DEW_implementation_plan.md Phase 1).
*/
int DEW_species_active( DEW_db species, int len_ox, int *id ){
    int active_mask[16] = {0};
    for (int i = 0; i < len_ox; i++){
        if (id[i] >= 0 && id[i] < 16){ active_mask[id[i]] = 1; }   /* id[i]<0: oxide outside DEW's 15-oxide canonical set (e.g. from DEW_build_id_map) - contributes nothing to the mask, correctly excluding species that would need it */
    }
    for (int k = 0; k < 15; k++){
        if (!active_mask[k] && species.Comp[k] != 0.0){
            return 0;
        }
    }
    return 1;
}

/** Position (within the active oxide list) of the H2O component, or -1 if absent. */
int DEW_find_id_H2O( int len_ox, int *id ){
    for (int i = 0; i < len_ox; i++){
        if (id[i] == 11){ return i; }   /* canonical oxide_info index of "H2O" */
    }
    return -1;
}

/* Canonical oxide_info order (0-14) that DEW_db.Comp[]/MuComp[] are indexed against -
   same order as oxide_info in TC_init_database.c and z_b.id[]'s target indices. */
static const char *DEW_CANON_OXIDES[15] = {
    "SiO2","Al2O3","CaO","MgO","FeO","K2O","Na2O","TiO2","O","MnO","Cr2O3","H2O","CO2","S","Cl"
};

int DEW_canonical_oxide_index( char *oxide_name ){
    for (int k = 0; k < 15; k++){
        if (strcmp(oxide_name, DEW_CANON_OXIDES[k]) == 0){ return k; }
    }
    return -1;
}

void DEW_build_id_map( int len_ox, char **ox_names, int *id_out ){
    for (int i = 0; i < len_ox; i++){
        id_out[i] = DEW_canonical_oxide_index(ox_names[i]);
    }
}

AQ_data init_DEW_aqueous_model(    int      n_dew_db,
                                    int      len_ox,
                                    int     *id,
                                    double  *ElEntropy,
                                    double   P,
                                    double   T,
                                    double   rho_w,
                                    double   gb_w              ){

    int i, j, s, k;

    int n_active = 0;
    for (s = 0; s < n_dew_db; s++){
        DEW_db sp = Access_DEW_DB(s);
        if (DEW_species_active(sp, len_ox, id)){ n_active++; }
    }

    AQ_data AQ;
    strcpy(AQ.name, "DEW");
    AQ.n_sp   = n_active;
    AQ.len_ox = len_ox;
    AQ.id_H2O = DEW_find_id_H2O(len_ox, id);
    AQ.rho_w  = rho_w;
    AQ.gb_w   = gb_w;

    AQ.em_names = malloc(n_active * sizeof(char*));
    AQ.gbase    = malloc(n_active * sizeof(double));
    AQ.z        = malloc(n_active * sizeof(double));
    AQ.em_comp  = malloc((n_active+1) * sizeof(double*));
    AQ.mu_comp  = malloc((n_active+1) * sizeof(double*));
    for (i = 0; i <= n_active; i++){
        AQ.em_comp[i] = malloc(len_ox * sizeof(double));
        AQ.mu_comp[i] = malloc((len_ox+1) * sizeof(double));
    }

    k = 0;
    for (s = 0; s < n_dew_db; s++){
        DEW_db sp = Access_DEW_DB(s);
        if (!DEW_species_active(sp, len_ox, id)){ continue; }

        AQ_ref ref = G_DEW_function(len_ox, id, ElEntropy, P, T, rho_w, sp.Name);

        AQ.em_names[k] = malloc(20 * sizeof(char));
        strcpy(AQ.em_names[k], ref.Name);
        AQ.gbase[k] = ref.gbase;
        AQ.z[k]     = ref.z;
        for (j = 0; j < len_ox; j++){
            AQ.em_comp[k][j] = ref.Comp[j];
            AQ.mu_comp[k][j] = ref.MuComp[j];
        }
        AQ.mu_comp[k][len_ox] = ref.MuComp[len_ox];
        k++;
    }

    /* water's own row: pure H2O composition, zero formation-reaction stoichiometry
       (its chemical potential is self-consistently solved, not levelled against Gamma) */
    for (j = 0; j < len_ox; j++){
        AQ.em_comp[n_active][j] = (id[j] == 11) ? 1.0 : 0.0;
        AQ.mu_comp[n_active][j] = 0.0;
    }
    AQ.mu_comp[n_active][len_ox] = 0.0;

    return AQ;
}

/**
    Water's reference-state density/G from TC_G_EM_function's own H2O standard state
    (Pitzer & Sterner 1994), consistent with the rest of the TC-group database's H2O
    endmember (see tools/DEW_implementation_plan.md, decision 3).
*/
void DEW_get_water_reference(  int      EM_dataset,
                                int      len_ox,
                                int     *id,
                                double  *bulk_rock,
                                double  *apo,
                                double   P,
                                double   T,
                                double  *rho_w,
                                double  *gb_w              ){

    PP_ref water = TC_G_EM_function(EM_dataset, len_ox, id, bulk_rock, apo, P, T, "H2O", "equilibrium");
    *rho_w = water.phase_density;
    *gb_w  = water.gbase;
}

AQ_data init_DEW_aqueous_model_at_point(   int      n_dew_db,
                                            int      EM_dataset,
                                            int      len_ox,
                                            int     *id,
                                            double  *bulk_rock,
                                            double  *apo,
                                            double  *ElEntropy,
                                            double   P,
                                            double   T              ){

    double rho_w, gb_w;
    DEW_get_water_reference(EM_dataset, len_ox, id, bulk_rock, apo, P, T, &rho_w, &gb_w);
    return init_DEW_aqueous_model(n_dew_db, len_ox, id, ElEntropy, P, T, rho_w, gb_w);
}

void free_DEW_aqueous_model( AQ_data *AQ ){
    for (int i = 0; i <= AQ->n_sp; i++){
        free(AQ->em_comp[i]);
        free(AQ->mu_comp[i]);
    }
    free(AQ->em_comp);
    free(AQ->mu_comp);
    for (int i = 0; i < AQ->n_sp; i++){
        free(AQ->em_names[i]);
    }
    free(AQ->em_names);
    free(AQ->gbase);
    free(AQ->z);
}

AQ_solver init_DEW_solver( AQ_data *AQ ){
    AQ_solver S;
    int n = AQ->n_sp;

    S.m  = calloc(n, sizeof(double));
    S.m1 = calloc(n, sizeof(double));
    S.x  = calloc(n+1, sizeof(double));
    S.a  = calloc(n+1, sizeof(double));
    S.mu = calloc(n+1, sizeof(double));

    for (int i = 0; i < 5; i++){ S.gamma_e[i] = 1.0; }

    S.sum_m         = 0.0;
    S.sum_charged_p = 0.0;
    S.sum_charged_m = 0.0;
    S.Gamma         = 0.0;
    S.I_str         = 0.0;
    S.Lambda        = 0.0;
    S.sigma         = 0.0;
    S.coef          = 0.0;
    S.a_coef        = 1.0;
    S.mu_w          = 0.0;
    S.mu_Hp         = 0.0;
    S.G             = 0.0;
    S.z_res         = 1.0;
    S.converged     = 0;

    return S;
}

void free_DEW_solver( AQ_solver *S ){
    free(S->m);
    free(S->m1);
    free(S->x);
    free(S->a);
    free(S->mu);
}

/** Set to 0 to roll back to the original single-step Newton-style mu_Hp charge-balance
    correction (still present below, unmodified, behind this flag) instead of the
    bracketed-bisection solve. The Newton step can overshoot by an unbounded amount in
    one iteration when the crude initial guess is badly charge-imbalanced (confirmed
    root cause of a real "DEW spuriously absent" bug, see tools/DEW_implementation_plan.md
    2026-08-12 entries) - charge balance vs mu_Hp is a monotonic 1D root-find (same
    problem structure as pH vs charge balance in any speciation code), so bisection can't
    overshoot into the same overflow territory: every trial mu_Hp it evaluates uses the
    same clamped exp() as everywhere else in this file, so the search itself is bounded. */
#define DEW_MU_HP_ROOT_FIND 1

#if DEW_MU_HP_ROOT_FIND
/**
    Recomputes every species' molality at a trial mu_Hp, holding mu_w/gamma_e fixed
    (matches the outer Picard loop's own convention - those are only refreshed once per
    outer pass, same as the original single-step correction did). Returns the signed
    charge residual sum(m*z); m_out is filled with the resulting molalities. Used only
    as the residual function for DEW_solve_mu_Hp_bisection below.
*/
static double DEW_charge_residual( AQ_data *AQ, double *Gamma_ox, double R, double T,
                                    int id_H2O, int len_ox, int n_sp,
                                    double mu_w, double mu_Hp, double *gamma_e,
                                    double *m_out ){
    DEW_stat_residual_evals++;
    double z_res = 0.0;
    for (int m = 0; m < n_sp; m++){
        int z_id = (int)fabs(AQ->z[m]);
        double mu_exp = -AQ->gbase[m];
        for (int j = 0; j < len_ox; j++){
            if (j == id_H2O){ mu_exp += AQ->mu_comp[m][j]*mu_w; }
            else             { mu_exp += AQ->mu_comp[m][j]*Gamma_ox[j]; }
        }
        mu_exp += AQ->mu_comp[m][len_ox]*mu_Hp;
        m_out[m] = exp(fmin(mu_exp/(R*T), 50.0)) / gamma_e[z_id];
        z_res += m_out[m]*AQ->z[m];
    }
    return z_res;
}

/**
    Bracketed bisection for mu_Hp that drives the charge residual to ~0, holding
    mu_w/gamma_e fixed for this outer iteration. Expands outward from the current mu_Hp
    guess (doubling step, capped at 60 expansions) until a sign change is bracketed,
    then bisects to tol. If no bracket is found (a genuinely degenerate case, e.g. every
    active species carries the same-sign charge - charge balance is then physically
    impossible and no mu_Hp fixes it), returns mu_Hp0 unchanged and lets the outer loop's
    own z_res/max_iter check surface the non-convergence, same as it always could.
*/
static double DEW_solve_mu_Hp_bisection( AQ_data *AQ, double *Gamma_ox, double R, double T,
                                          int id_H2O, int len_ox, int n_sp,
                                          double mu_w, double mu_Hp0, double *gamma_e,
                                          double tol, int max_iter ){
    double *m_scratch = malloc(n_sp * sizeof(double));

    double f0 = DEW_charge_residual(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp, mu_w, mu_Hp0, gamma_e, m_scratch);
    if (fabs(f0) < tol){ free(m_scratch); return mu_Hp0; }

    double lo = mu_Hp0, hi = mu_Hp0, f_lo = f0, f_hi = f0;
    double step = fmax(fabs(R*T), 1.0);
    int bracketed = 0;
    for (int k = 0; k < 60; k++){
        lo = mu_Hp0 - step;
        hi = mu_Hp0 + step;
        f_lo = DEW_charge_residual(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp, mu_w, lo, gamma_e, m_scratch);
        f_hi = DEW_charge_residual(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp, mu_w, hi, gamma_e, m_scratch);
        if ((f_lo <= 0.0 && f_hi >= 0.0) || (f_lo >= 0.0 && f_hi <= 0.0)){ bracketed = 1; break; }
        step *= 2.0;
    }
    if (!bracketed){ free(m_scratch); return mu_Hp0; }

    double result = 0.5*(lo+hi);
    for (int k = 0; k < max_iter; k++){
        double mid = 0.5*(lo+hi);
        double f_mid = DEW_charge_residual(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp, mu_w, mid, gamma_e, m_scratch);
        result = mid;
        if (fabs(f_mid) < tol || 0.5*(hi-lo) < 1e-10){ break; }
        if ((f_lo < 0.0) == (f_mid < 0.0)){ lo = mid; f_lo = f_mid; }
        else                               { hi = mid; }
    }

    free(m_scratch);
    return result;
}

/**
    Same residual as DEW_charge_residual, but also returns d(z_res)/d(mu_Hp) analytically:
    each species' molality is m[i] = exp(mu_exp[i]/RT)/gamma_e, and mu_exp[i] is LINEAR in
    mu_Hp with slope AQ->mu_comp[i][len_ox] - so dm[i]/dmu_Hp = m[i]*mu_comp[i][len_ox]/(RT),
    a free byproduct of the same loop that computes m[i] itself (no extra residual
    evaluations needed, unlike a finite-difference derivative would require). Used only by
    DEW_solve_mu_Hp_safeguarded below (gv.DEW_solve_algorithm==2).
*/
static double DEW_charge_residual_deriv( AQ_data *AQ, double *Gamma_ox, double R, double T,
                                          int id_H2O, int len_ox, int n_sp,
                                          double mu_w, double mu_Hp, double *gamma_e,
                                          double *m_out, double *dres_out ){
    DEW_stat_residual_evals++;
    double z_res = 0.0;
    double dres  = 0.0;
    double RT = R*T;
    for (int m = 0; m < n_sp; m++){
        int z_id = (int)fabs(AQ->z[m]);
        double mu_exp = -AQ->gbase[m];
        for (int j = 0; j < len_ox; j++){
            if (j == id_H2O){ mu_exp += AQ->mu_comp[m][j]*mu_w; }
            else             { mu_exp += AQ->mu_comp[m][j]*Gamma_ox[j]; }
        }
        double slope = AQ->mu_comp[m][len_ox];
        mu_exp += slope*mu_Hp;
        double raw = mu_exp/RT;
        if (raw >= 50.0){
            m_out[m] = exp(50.0) / gamma_e[z_id];
        }
        else{
            m_out[m] = exp(raw) / gamma_e[z_id];
            dres += m_out[m]*AQ->z[m]*slope/RT;
        }
        z_res += m_out[m]*AQ->z[m];
    }
    *dres_out = dres;
    return z_res;
}

/**
    Newton-Raphson safeguarded by bisection ("rtsafe"-style) for mu_Hp - same bracketing
    guarantee as DEW_solve_mu_Hp_bisection (a trial mu_Hp is only ever accepted if it stays
    strictly within a verified sign-changing bracket), but converges quadratically instead
    of linearly once close to the root, since the analytic derivative above is essentially
    free. Falls back to a plain bisection step whenever the Newton step would leave the
    current bracket (or the derivative is too small to trust) - this is exactly the
    robustness property the codebase's own switch away from raw Newton (see the comment on
    DEW_MU_HP_ROOT_FIND above) was missing, so it cannot reintroduce that class of bug:
    worst case, every iteration degrades to plain bisection. gv.DEW_solve_algorithm==2
    selects this over DEW_solve_mu_Hp_bisection at both call sites below - prototyped and
    benchmarked in isolation first (40-1000x fewer residual evaluations per solve; total
    wall-clock win is largest exactly where bisection was struggling, e.g. deep/hot points
    needing hundreds of bisection passes - roughly neutral where bisection was already
    cheap), then verified for 0 regressions across an 84-point P-T sweep over mpe/ume/mbe
    before being wired in as a runtime option here.
*/
static double DEW_solve_mu_Hp_safeguarded( AQ_data *AQ, double *Gamma_ox, double R, double T,
                                            int id_H2O, int len_ox, int n_sp,
                                            double mu_w, double mu_Hp0, double *gamma_e,
                                            double tol, int max_iter ){
    double *m_scratch = malloc(n_sp * sizeof(double));

    double f0 = DEW_charge_residual(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp, mu_w, mu_Hp0, gamma_e, m_scratch);
    if (fabs(f0) < tol){ free(m_scratch); return mu_Hp0; }

    double lo = mu_Hp0, hi = mu_Hp0, f_lo = f0, f_hi = f0;
    double step = fmax(fabs(R*T), 1.0);
    int bracketed = 0;
    for (int k = 0; k < 60; k++){
        lo = mu_Hp0 - step;
        hi = mu_Hp0 + step;
        f_lo = DEW_charge_residual(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp, mu_w, lo, gamma_e, m_scratch);
        f_hi = DEW_charge_residual(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp, mu_w, hi, gamma_e, m_scratch);
        if ((f_lo <= 0.0 && f_hi >= 0.0) || (f_lo >= 0.0 && f_hi <= 0.0)){ bracketed = 1; break; }
        step *= 2.0;
    }
    if (!bracketed){ free(m_scratch); return mu_Hp0; }
    if (f_lo > f_hi){
        double tmp = lo; lo = hi; hi = tmp;
        tmp = f_lo; f_lo = f_hi; f_hi = tmp;
    }

    double x    = 0.5*(lo+hi);
    double dres;
    double f_x  = DEW_charge_residual_deriv(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp, mu_w, x, gamma_e, m_scratch, &dres);
    double dx_prev = fabs(hi-lo);
    double dx_cur  = dx_prev;

    for (int k = 0; k < max_iter; k++){
        if (fabs(f_x) < tol || 0.5*(hi-lo) < 1e-10){ break; }

        if (f_x < 0.0){ lo = x; } else { hi = x; }

        int use_bisection = 0;
        double x_new;
        if (fabs(dres) < 1e-300){
            use_bisection = 1;
        }
        else{
            double dx_newton = -f_x/dres;
            x_new = x + dx_newton;
            if (x_new <= lo || x_new >= hi){ use_bisection = 1; }
            else if (fabs(2.0*f_x) > fabs(dx_prev*dres)){ use_bisection = 1; }
        }

        dx_prev = dx_cur;
        if (use_bisection){
            x_new  = 0.5*(lo+hi);
            dx_cur = fabs(x_new - x);
        }
        else{
            dx_cur = fabs(x_new - x);
        }

        x   = x_new;
        f_x = DEW_charge_residual_deriv(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp, mu_w, x, gamma_e, m_scratch, &dres);
    }

    free(m_scratch);
    return x;
}
#endif

/**
    Fixed-point aqueous speciation solve, given the system's current oxide-component
    chemical potentials (Gamma_ox, kJ/mol). Direct port of aq_min_iterative (see file
    header for the one deliberate deviation, P unit handling).

    R must be in kJ/mol/K (MAGEMin's convention, 0.0083144), T in K, P in kbar.

    x_warm: NULL for the normal Hp_eps-derived cold start (matches original behaviour
    exactly). Non-NULL (see DEW_aq_min_warmstart below): a previously converged
    mole-fraction vector [n_sp+1], index n_sp=water - seeds S->m[] directly from it
    (skipping the Hp_eps-derived initial guess) so pass 1 already starts from something
    close to the true fixed point instead of a generic starting guess. Hp_eps is ignored
    when x_warm is given.

    algorithm: gv.DEW_solve_algorithm as forwarded by DEW_aq_min_iterative_dispatch - only
    the mu_Hp inner-solve choice is read here (2 = DEW_solve_mu_Hp_safeguarded, anything
    else = DEW_solve_mu_Hp_bisection); the outer Picard loop below is always the plain
    (undamped) variant regardless of this value.
*/
void DEW_aq_min_iterative(     AQ_data     *AQ,
                                AQ_solver   *S,
                                double      *Gamma_ox,
                                double       R,
                                double       T,
                                double       P,
                                int          max_iter,
                                double       z_res_tol,
                                double       Hp_eps,
                                const double *x_warm,
                                int          algorithm      ){

    int n_sp   = AQ->n_sp;
    int len_ox = AQ->len_ox;
    int id_H2O = AQ->id_H2O;

    int i, j, m;
    double Pbar = P * 1000.0;

    double Mw  = 18.01528/1000.0;      /* kg/mol */
    double a_i = 3.72;                  /* Angstrom, fixed ion-size parameter (Truesdell-Jones), shared by every species */
    double b   = 0.03;                  /* B-dot parameter */

    double Ag = DEW_Agamma(T, Pbar, AQ->rho_w);
    double Bg = DEW_Bgamma(T, Pbar, AQ->rho_w);

    for (i = 0; i < n_sp; i++){ S->m[i] = 0.0; }
    S->I_str      = 0.0;
    S->a_coef = 1.0;
    S->coef   = 1.0;
    for (i = 0; i < 5; i++){ S->gamma_e[i] = 1.0; }

    S->mu_w  = AQ->gb_w + R*T*log(S->a_coef);

    if (x_warm != NULL){
        double Omega = 1.0/Mw;
        double x_w   = x_warm[n_sp];
        if (x_w > 0.0 && isfinite(x_w)){
            for (m = 0; m < n_sp; m++){
                double m_val = x_warm[m]*Omega/x_w;
                S->m[m] = (isfinite(m_val) && m_val > 0.0) ? m_val : 0.0;
            }
        }
        S->mu_Hp = R*T*log(1e-6);
    }
    else{
        S->mu_Hp = R*T*log(Hp_eps);

        for (m = 0; m < n_sp; m++){
            int z_id = (int)fabs(AQ->z[m]);
            double mu_exp = -AQ->gbase[m];
            for (j = 0; j < len_ox; j++){
                if (j == id_H2O){ mu_exp += AQ->mu_comp[m][j]*S->mu_w; }
                else             { mu_exp += AQ->mu_comp[m][j]*Gamma_ox[j]; }
            }
            mu_exp += AQ->mu_comp[m][len_ox]*S->mu_Hp;
            S->m[m] = exp(fmin(mu_exp/(R*T), 50.0)) / S->gamma_e[z_id];
        }
    }

    double z_res = 1.0;
    int ite = 1;
    double best_z_res     = z_res;
    int    best_z_res_ite = ite;
    const int stall_window = 100;
    while (z_res > z_res_tol && ite <= max_iter){
        DEW_stat_picard_passes++;

        double sum_m_charged = 0.0;
        for (i = 0; i < n_sp; i++){ sum_m_charged += S->m[i]*AQ->z[i]*AQ->z[i]; }

        S->sum_m = 0.0;
        for (i = 0; i < n_sp; i++){ S->sum_m += S->m[i]; }

        S->Gamma  = -log10(1.0 + Mw*S->sum_m);
        S->I_str      = 0.5*sum_m_charged;
        S->Lambda = 1.0 + a_i*Bg*sqrt(S->I_str);

        if (S->sum_m > 0.0){

            S->sigma  = (S->I_str > 1e-8) ? 3.0/(pow(a_i,3.0)*pow(Bg,3.0)*pow(S->I_str,1.5)) * (S->Lambda - 1.0/S->Lambda - 2.0*log(S->Lambda)) : 1.0;
            double Gamma_over_Mw_sum_m = (Mw*S->sum_m > 1e-8) ? S->Gamma/(Mw*S->sum_m) : -1.0/log(10.0);
            S->coef   = -log(10.0) * ( Gamma_over_Mw_sum_m - b*S->I_str/2.0 + sum_m_charged*Ag*(sqrt(S->I_str*S->sigma))/(3.0*S->sum_m) );
            S->a_coef = exp( fmax(fmin(-S->sum_m*S->coef / (1.0/Mw), 700.0), -700.0) );
        }
        else{

            S->sigma  = 1.0;
            S->coef   = 0.0;
            S->a_coef = 1.0;
        }

        for (i = 0; i < 5; i++){
            double z_i = (double)i;
            double term1 = -Ag*z_i*z_i*sqrt(S->I_str)/S->Lambda;
            S->gamma_e[i] = pow(10.0, fmax(fmin(term1 + S->Gamma + b*S->I_str, 300.0), -300.0));
        }

        S->sum_charged_p = 0.0;
        S->sum_charged_m = 0.0;
        for (i = 0; i < n_sp; i++){
            if (AQ->z[i] > 0.0){ S->sum_charged_p += fabs(AQ->z[i])*S->m[i]; }
            else if (AQ->z[i] < 0.0){ S->sum_charged_m += fabs(AQ->z[i])*S->m[i]; }
        }

        S->mu_w = AQ->gb_w + R*T*log(S->a_coef);

#if DEW_MU_HP_ROOT_FIND
        if (algorithm == 2){
            S->mu_Hp = DEW_solve_mu_Hp_safeguarded(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp,
                                                    S->mu_w, S->mu_Hp, S->gamma_e,
                                                    z_res_tol, 80);
        }
        else{
            S->mu_Hp = DEW_solve_mu_Hp_bisection(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp,
                                                  S->mu_w, S->mu_Hp, S->gamma_e,
                                                  z_res_tol, 80);
        }
#else
        double mu_Hp_cor = R*T*log( sqrt(S->sum_charged_p*S->sum_charged_m) / S->sum_charged_p );
        if (isfinite(mu_Hp_cor)){
            mu_Hp_cor = fmax(fmin(mu_Hp_cor, 700.0*R*T), -700.0*R*T);
            S->mu_Hp += mu_Hp_cor;
        }
#endif

        for (m = 0; m < n_sp; m++){
            int z_id = (int)fabs(AQ->z[m]);
            double mu_exp = -AQ->gbase[m];
            for (j = 0; j < len_ox; j++){
                if (j == id_H2O){ mu_exp += AQ->mu_comp[m][j]*S->mu_w; }
                else             { mu_exp += AQ->mu_comp[m][j]*Gamma_ox[j]; }
            }
            mu_exp += AQ->mu_comp[m][len_ox]*S->mu_Hp;
            S->m[m] = exp(fmin(mu_exp/(R*T), 50.0)) / S->gamma_e[z_id];
        }

        z_res = 0.0;
        for (i = 0; i < n_sp; i++){ z_res += S->m[i]*AQ->z[i]; }
        z_res = fabs(z_res);

        if (z_res < best_z_res){
            best_z_res     = z_res;
            best_z_res_ite = ite;
        }
        else if (ite - best_z_res_ite > stall_window){
            ite += 1;
            break;
        }

        ite += 1;
    }

    S->z_res     = z_res;
    S->converged = (z_res <= z_res_tol) ? 1 : 0;

    double Omega = 1.0/Mw;
    double sum_m = 0.0;
    for (i = 0; i < n_sp; i++){ sum_m += S->m[i]; }
    double x_w = Omega/(Omega + sum_m);

    for (i = 0; i < n_sp; i++){
        S->m1[i] = S->m[i];
        S->x[i]  = S->m[i]/(Omega + sum_m);
    }
    S->x[n_sp] = x_w;

    for (i = 0; i < n_sp; i++){
        int z_id = (int)fabs(AQ->z[i]);
        S->a[i]  = S->m[i]*S->gamma_e[z_id];
        S->mu[i] = AQ->gbase[i] + R*T*log(S->m[i]*S->gamma_e[z_id]);
    }
    S->a[n_sp]  = S->a_coef;
    S->mu[n_sp] = S->mu_w;

    S->G = 0.0;
    for (i = 0; i <= n_sp; i++){
        if (S->x[i] > 0.0){ S->G += S->mu[i]*S->x[i]; }
    }
}

/**
    Algorithm 1 (gv.DEW_solve_algorithm==1): damped/mixed Picard variant of
    DEW_aq_min_iterative above. Same equilibrium condition, same per-pass computation of
    Gamma/I_str/sigma/coef/a_coef/gamma_e/mu_w/mu_Hp - deliberately kept as a full,
    independent copy rather than threading a branch through the original (numerically
    delicate, heavily-clamped) function, so a bug here cannot touch the default algorithm
    0 path at all.

    Root cause this targets (see NLopt_opt_DEW_function's 2026-08-14 diagnostic: I_str/
    sum_m spiking to unphysical values, a_coef collapsing to 0, before the stall detector
    cuts the pass off - confirmed to happen ONLY in the concentrated/high-ionic-strength
    regime, e.g. 20kbar/500C, never in the dilute one, e.g. 2kbar/200C): the composition
    (S->m) -> activity-coefficient (a_coef/gamma_e via Gamma/I_str/sigma/coef) -> mu_w/
    mu_Hp -> composition feedback loop is only weakly nonlinear at low ionic strength
    (plain Picard converges in a few passes) but becomes strongly, tightly coupled once
    ionic strength reaches ~1 mol/kg+ - the classic regime where an unmixed fixed-point
    iteration oscillates/overshoots instead of converging (the same issue self-consistent-
    field electronic-structure codes hit, solved the same way there: damped/mixed updates).

    Fix: the newly computed molality vector is linearly mixed with the previous pass's
    molality before being used as next pass's input, instead of being used directly:
        m_(k+1) = beta*m_new + (1-beta)*m_k
    beta starts at 1.0 (identical to plain Picard - a well-behaved start pays no damping
    penalty). Whenever a pass's mixed z_res fails to improve on the previous pass's z_res
    (a direct oscillation/divergence signal), the SAME already-computed m_new is re-mixed
    at half beta and re-checked (cheap - no mu_exp/bisection recompute needed) until it
    improves or beta hits a floor; beta is then slowly relaxed back up on sustained
    improvement, so it only stays small while actually needed.

    Verified 2026-08-14 (8-point P-T sweep, mpe + mbe/ume sanity): same G/phase assemblage
    as algorithm 0 to within normal outer-LP path-dependent noise, no I_str/a_coef
    pathology at any point tested, warm-start hit rate at 20kbar/500C rose from 44% to 98%
    (the damped solve tracks the outer PGE loop's small Gamma_ox steps far more reliably),
    ~34% faster overall on the sweep (up to 61% at the hardest point).
*/
void DEW_aq_min_iterative_mixed(   AQ_data     *AQ,
                                    AQ_solver   *S,
                                    double      *Gamma_ox,
                                    double       R,
                                    double       T,
                                    double       P,
                                    int          max_iter,
                                    double       z_res_tol,
                                    double       Hp_eps,
                                    const double *x_warm       ){

    int n_sp   = AQ->n_sp;
    int len_ox = AQ->len_ox;
    int id_H2O = AQ->id_H2O;

    int i, j, m;
    double Pbar = P * 1000.0;

    double Mw  = 18.01528/1000.0;      /* kg/mol */
    double a_i = 3.72;                  /* Angstrom, fixed ion-size parameter (Truesdell-Jones), shared by every species */
    double b   = 0.03;                  /* B-dot parameter */

    double Ag = DEW_Agamma(T, Pbar, AQ->rho_w);
    double Bg = DEW_Bgamma(T, Pbar, AQ->rho_w);

    double *m_prev = malloc(n_sp * sizeof(double));
    double *m_raw  = malloc(n_sp * sizeof(double));

    for (i = 0; i < n_sp; i++){ S->m[i] = 0.0; }
    S->I_str      = 0.0;
    S->a_coef = 1.0;
    S->coef   = 1.0;
    for (i = 0; i < 5; i++){ S->gamma_e[i] = 1.0; }

    S->mu_w  = AQ->gb_w + R*T*log(S->a_coef);

    if (x_warm != NULL){
        double Omega = 1.0/Mw;
        double x_w   = x_warm[n_sp];
        if (x_w > 0.0 && isfinite(x_w)){
            for (m = 0; m < n_sp; m++){
                double m_val = x_warm[m]*Omega/x_w;
                S->m[m] = (isfinite(m_val) && m_val > 0.0) ? m_val : 0.0;
            }
        }
        S->mu_Hp = R*T*log(1e-6);
    }
    else{
        S->mu_Hp = R*T*log(Hp_eps);

        for (m = 0; m < n_sp; m++){
            int z_id = (int)fabs(AQ->z[m]);
            double mu_exp = -AQ->gbase[m];
            for (j = 0; j < len_ox; j++){
                if (j == id_H2O){ mu_exp += AQ->mu_comp[m][j]*S->mu_w; }
                else             { mu_exp += AQ->mu_comp[m][j]*Gamma_ox[j]; }
            }
            mu_exp += AQ->mu_comp[m][len_ox]*S->mu_Hp;
            S->m[m] = exp(fmin(mu_exp/(R*T), 50.0)) / S->gamma_e[z_id];
        }
    }

    double z_res = 1.0;
    double z_res_prev = 1.0;
    int ite = 1;
    double best_z_res     = z_res;
    int    best_z_res_ite = ite;
    const int stall_window = 100;

    double beta = 1.0;
    const double beta_floor    = 0.03125;   /* 1/32 - five halvings from 1.0 */
    const double beta_relax_up = 1.10;      /* slow creep back toward 1.0 on sustained improvement */
    const int    max_mix_tries = 6;

    while (z_res > z_res_tol && ite <= max_iter){
        DEW_stat_picard_passes++;

        for (i = 0; i < n_sp; i++){ m_prev[i] = S->m[i]; }

        double sum_m_charged = 0.0;
        for (i = 0; i < n_sp; i++){ sum_m_charged += S->m[i]*AQ->z[i]*AQ->z[i]; }

        S->sum_m = 0.0;
        for (i = 0; i < n_sp; i++){ S->sum_m += S->m[i]; }

        S->Gamma  = -log10(1.0 + Mw*S->sum_m);
        S->I_str      = 0.5*sum_m_charged;
        S->Lambda = 1.0 + a_i*Bg*sqrt(S->I_str);

        if (S->sum_m > 0.0){
            S->sigma  = (S->I_str > 1e-8) ? 3.0/(pow(a_i,3.0)*pow(Bg,3.0)*pow(S->I_str,1.5)) * (S->Lambda - 1.0/S->Lambda - 2.0*log(S->Lambda)) : 1.0;
            double Gamma_over_Mw_sum_m = (Mw*S->sum_m > 1e-8) ? S->Gamma/(Mw*S->sum_m) : -1.0/log(10.0);
            S->coef   = -log(10.0) * ( Gamma_over_Mw_sum_m - b*S->I_str/2.0 + sum_m_charged*Ag*(sqrt(S->I_str*S->sigma))/(3.0*S->sum_m) );
            S->a_coef = exp( fmax(fmin(-S->sum_m*S->coef / (1.0/Mw), 700.0), -700.0) );
        }
        else{
            S->sigma  = 1.0;
            S->coef   = 0.0;
            S->a_coef = 1.0;
        }

        for (i = 0; i < 5; i++){
            double z_i = (double)i;
            double term1 = -Ag*z_i*z_i*sqrt(S->I_str)/S->Lambda;
            S->gamma_e[i] = pow(10.0, fmax(fmin(term1 + S->Gamma + b*S->I_str, 300.0), -300.0));
        }

        S->sum_charged_p = 0.0;
        S->sum_charged_m = 0.0;
        for (i = 0; i < n_sp; i++){
            if (AQ->z[i] > 0.0){ S->sum_charged_p += fabs(AQ->z[i])*S->m[i]; }
            else if (AQ->z[i] < 0.0){ S->sum_charged_m += fabs(AQ->z[i])*S->m[i]; }
        }

        S->mu_w = AQ->gb_w + R*T*log(S->a_coef);

        S->mu_Hp = DEW_solve_mu_Hp_bisection(AQ, Gamma_ox, R, T, id_H2O, len_ox, n_sp,
                                              S->mu_w, S->mu_Hp, S->gamma_e,
                                              z_res_tol, 80);

        for (m = 0; m < n_sp; m++){
            int z_id = (int)fabs(AQ->z[m]);
            double mu_exp = -AQ->gbase[m];
            for (j = 0; j < len_ox; j++){
                if (j == id_H2O){ mu_exp += AQ->mu_comp[m][j]*S->mu_w; }
                else             { mu_exp += AQ->mu_comp[m][j]*Gamma_ox[j]; }
            }
            mu_exp += AQ->mu_comp[m][len_ox]*S->mu_Hp;
            m_raw[m] = exp(fmin(mu_exp/(R*T), 50.0)) / S->gamma_e[z_id];
        }

        /* Adaptive linear mixing: retry at half beta (cheap - no mu_exp/bisection redo)
           until the mixed result improves on the previous pass's z_res, or beta floors out. */
        double z_res_trial;
        int mix_tries = 0;
        while (1){
            for (m = 0; m < n_sp; m++){ S->m[m] = beta*m_raw[m] + (1.0-beta)*m_prev[m]; }
            z_res_trial = 0.0;
            for (i = 0; i < n_sp; i++){ z_res_trial += S->m[i]*AQ->z[i]; }
            z_res_trial = fabs(z_res_trial);
            if (z_res_trial <= z_res_prev || beta <= beta_floor || mix_tries >= max_mix_tries){ break; }
            beta = fmax(beta*0.5, beta_floor);
            mix_tries++;
        }
        z_res = z_res_trial;
        if (z_res < z_res_prev){ beta = fmin(1.0, beta*beta_relax_up); }
        z_res_prev = z_res;

        if (z_res < best_z_res){
            best_z_res     = z_res;
            best_z_res_ite = ite;
        }
        else if (ite - best_z_res_ite > stall_window){
            ite += 1;
            break;
        }

        ite += 1;
    }

    free(m_prev);
    free(m_raw);

    S->z_res     = z_res;
    S->converged = (z_res <= z_res_tol) ? 1 : 0;

    double Omega = 1.0/Mw;
    double sum_m = 0.0;
    for (i = 0; i < n_sp; i++){ sum_m += S->m[i]; }
    double x_w = Omega/(Omega + sum_m);

    for (i = 0; i < n_sp; i++){
        S->m1[i] = S->m[i];
        S->x[i]  = S->m[i]/(Omega + sum_m);
    }
    S->x[n_sp] = x_w;

    for (i = 0; i < n_sp; i++){
        int z_id = (int)fabs(AQ->z[i]);
        S->a[i]  = S->m[i]*S->gamma_e[z_id];
        S->mu[i] = AQ->gbase[i] + R*T*log(S->m[i]*S->gamma_e[z_id]);
    }
    S->a[n_sp]  = S->a_coef;
    S->mu[n_sp] = S->mu_w;

    S->G = 0.0;
    for (i = 0; i <= n_sp; i++){
        if (S->x[i] > 0.0){ S->G += S->mu[i]*S->x[i]; }
    }
}

/**
    Shared core of DEW_aq_evaluate/DEW_aq_pH (and, externally, dump_function.c's
    DEW output reporting): derives per-species molality m[] (caller-allocated,
    size AQ->n_sp) and the Debye-Hückel activity coefficient gamma_e[5] (indexed by
    |z|, 0..4) from a converged mole-fraction vector x[n_sp+1] (index n_sp =
    water). Also returns a_coef (water activity coefficient), needed by
    DEW_aq_evaluate for mu_w and by dump_function.c for the reported water
    activity, but not by DEW_aq_pH. Not `static` so it can be called directly
    from outside this file (see DEW_aq_solver.h) - a fully-populated AQ_data
    (em_names set) is required by callers that need to identify species by name.
*/
void DEW_aq_activity_state(        AQ_data      *AQ,
                                    const double *x,
                                    double        T,
                                    double        P,
                                    double       *m,
                                    double        gamma_e[5],
                                    double       *a_coef_out       ){

    int n_sp = AQ->n_sp;
    double Pbar = P * 1000.0;

    double Mw  = 18.01528/1000.0;
    double a_i = 3.72;
    double b   = 0.03;
    double Omega = 1.0/Mw;

    double Ag = DEW_Agamma(T, Pbar, AQ->rho_w);
    double Bg = DEW_Bgamma(T, Pbar, AQ->rho_w);

    double x_w = x[n_sp];
    for (int i = 0; i < n_sp; i++){ m[i] = x[i]*Omega/x_w; }

    double sum_m_charged = 0.0;
    double sum_m = 0.0;
    for (int i = 0; i < n_sp; i++){
        sum_m_charged += m[i]*AQ->z[i]*AQ->z[i];
        sum_m         += m[i];
    }

    double Gamma  = -log10(1.0 + Mw*sum_m);
    double I_str  = 0.5*sum_m_charged;
    double Lambda = 1.0 + a_i*Bg*sqrt(I_str);
    double sigma, coef, a_coef;
    if (sum_m > 0.0){
        sigma  = (I_str > 1e-8) ? 3.0/(pow(a_i,3.0)*pow(Bg,3.0)*pow(I_str,1.5)) * (Lambda - 1.0/Lambda - 2.0*log(Lambda)) : 1.0;
        double Gamma_over_Mw_sum_m = (Mw*sum_m > 1e-8) ? Gamma/(Mw*sum_m) : -1.0/log(10.0);
        coef   = -log(10.0) * ( Gamma_over_Mw_sum_m - b*I_str/2.0 + sum_m_charged*Ag*(sqrt(I_str*sigma))/(3.0*sum_m) );
        a_coef = exp( fmax(fmin(-sum_m*coef / (1.0/Mw), 700.0), -700.0) );
    }
    else{
        sigma  = 1.0;
        coef   = 0.0;
        a_coef = 1.0;
    }

    for (int i = 0; i < 5; i++){
        double z_i = (double)i;
        double term1 = -Ag*z_i*z_i*sqrt(I_str)/Lambda;
        gamma_e[i] = pow(10.0, fmax(fmin(term1 + Gamma + b*I_str, 300.0), -300.0));
    }

    *a_coef_out = a_coef;
}

void DEW_aq_evaluate(      AQ_data      *AQ,
                            const double *x,
                            double        R,
                            double        T,
                            double        P,
                            double       *mu_out,
                            double       *G_out        ){

    int n_sp = AQ->n_sp;

    double *m = malloc(n_sp * sizeof(double));
    double gamma_e[5];
    double a_coef;
    DEW_aq_activity_state(AQ, x, T, P, m, gamma_e, &a_coef);

    for (int i = 0; i < n_sp; i++){
        int z_id  = (int)fabs(AQ->z[i]);
        mu_out[i] = AQ->gbase[i] + R*T*log(m[i]*gamma_e[z_id]);
    }
    double mu_w = AQ->gb_w + R*T*log(a_coef);
    mu_out[n_sp] = mu_w;

    *G_out = 0.0;
    for (int i = 0; i <= n_sp; i++){
        if (x[i] > 0.0){ *G_out += mu_out[i]*x[i]; }
    }

    free(m);
}

/** Set to 0 to roll back to a single DEW_aq_min_iterative call (Hp_eps fixed at 1e-6)
    instead of the multi-start scan below. Confirmed via direct testing (not assumed) that
    DEW_aq_min_iterative's fixed-point iteration has multiple distinct self-consistent
    fixed points depending on the Hp_eps starting guess - not just different convergence
    paths to the same answer, but genuinely different converged G (e.g. one starting guess
    finding DEW clearly favoured at 600-800 C, another converging to a DEW-absent
    state with measurably HIGHER G at the exact same conditions). A Picard-type iteration
    like this has no general guarantee of a unique fixed point, so trusting any single
    arbitrary starting guess is not safe - see tools/DEW_implementation_plan.md 2026-08-12. */
#define DEW_MULTISTART 1

static void DEW_aq_min_iterative_dispatch( int           algorithm,
                                            AQ_data      *AQ,
                                            AQ_solver    *S,
                                            double       *Gamma_ox,
                                            double        R,
                                            double        T,
                                            double        P,
                                            int           max_iter,
                                            double        z_res_tol,
                                            double        Hp_eps,
                                            const double *x_warm       ){
    if (algorithm == 1){
        DEW_aq_min_iterative_mixed(AQ, S, Gamma_ox, R, T, P, max_iter, z_res_tol, Hp_eps, x_warm);
    }
    else{
        DEW_aq_min_iterative(AQ, S, Gamma_ox, R, T, P, max_iter, z_res_tol, Hp_eps, x_warm, algorithm);
    }
}

/**
    Multi-start wrapper around DEW_aq_min_iterative: tries a spread of Hp_eps starting
    guesses and keeps whichever converges to the lowest G, evaluated self-consistently via
    DEW_aq_evaluate using AQ's own (absolute) gbase for every candidate - the same
    quantity DEW_aq_min_iterative's own equilibrium condition targets, so the ranking is
    valid regardless of whatever hyperplane-rotation convention the caller applies
    afterward. Re-solves once more at the winning Hp_eps at the end so every field of S
    (not just x/m1) is self-consistently populated by that final run, not stitched
    together from different candidates. DEW_MULTISTART=0 rolls back to a single
    Hp_eps=1e-6 call (the pre-2026-08-12 behavior) - the symbol always exists either way
    so no caller needs its own #if.
*/
void DEW_aq_min_multistart(    AQ_data     *AQ,
                                AQ_solver   *S,
                                double      *Gamma_ox,
                                double       R,
                                double       T,
                                double       P,
                                int          max_iter,
                                double       z_res_tol,
                                int          algorithm    ){

#if DEW_MULTISTART
    static const double Hp_eps_starts[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-10, 1e-12};
    const int n_starts = (int)(sizeof(Hp_eps_starts)/sizeof(Hp_eps_starts[0]));

    int n_sp = AQ->n_sp;
    double *mu_tmp = malloc((n_sp+1) * sizeof(double));

    double best_G = 0.0;
    double best_Hp_eps = Hp_eps_starts[0];
    int have_best = 0;
    int best_converged = 0;

    for (int s = 0; s < n_starts; s++){
        DEW_aq_min_iterative_dispatch(algorithm, AQ, S, Gamma_ox, R, T, P, max_iter, z_res_tol, Hp_eps_starts[s], NULL);

        double G_tmp;
        DEW_aq_evaluate(AQ, S->x, R, T, P, mu_tmp, &G_tmp);

        int take = !have_best
                   || (S->converged && !best_converged)
                   || (S->converged == best_converged && G_tmp < best_G);
        if (take){
            best_G         = G_tmp;
            best_Hp_eps    = Hp_eps_starts[s];
            have_best      = 1;
            best_converged = S->converged;
        }
    }

    free(mu_tmp);

    /* leave S fully, self-consistently populated by the winning starting guess */
    DEW_aq_min_iterative_dispatch(algorithm, AQ, S, Gamma_ox, R, T, P, max_iter, z_res_tol, best_Hp_eps, NULL);
#else
    DEW_aq_min_iterative_dispatch(algorithm, AQ, S, Gamma_ox, R, T, P, max_iter, z_res_tol, 1e-6, NULL);
#endif
}

/**
    Warm-start accelerator for the outer PGE loop. Attempts a single Picard solve seeded
    directly from a previous outer-iteration's converged mole-fraction vector (x_warm),
    instead of re-running the full 8-point Hp_eps multistart grid DEW_aq_min_multistart does
    above. Only worth trusting once the caller has already verified x_warm came from a real
    multistart-grid solve for THIS point (see NLopt_opt_DEW_function / SS_ref.dew_warm_ok)
    - relies on Gamma_ox changing only incrementally between consecutive outer PGE
    iterations of the same point (the premise the whole PGE architecture already runs on for
    every other phase's per-iteration NLopt warm start), NOT on x_warm being a
    globally-safe guess in general - a single Picard start from an arbitrary/unverified
    guess is exactly what motivated multistart in the first place (see the DEW_MULTISTART
    comment above).

    Returns S->converged (1/0). Caller should fall back to DEW_aq_min_multistart's full grid
    whenever this returns 0 - a stalled/capped warm-started solve carries no more special
    trust than any other non-converged candidate.
*/
int DEW_aq_min_warmstart(      AQ_data      *AQ,
                                AQ_solver    *S,
                                double       *Gamma_ox,
                                double        R,
                                double        T,
                                double        P,
                                int           max_iter,
                                double        z_res_tol,
                                const double *x_warm,
                                int           algorithm      ){
    DEW_aq_min_iterative_dispatch(algorithm, AQ, S, Gamma_ox, R, T, P, max_iter, z_res_tol, 0.0, x_warm);
    return S->converged;
}

/**
    pH = -log10(activity of H+) = -log10(m_H+ * gamma_H+), the standard definition -
    NOT a port of MAGEMin.jl's print_local_minimum_AQ (local_optimizer.jl:888), which
    computes -log10(act.m1[id_Hp] * act.a[id_Hp]) where act.a[id_Hp] is ALREADY
    m[id_Hp]*gamma_e[id_Hp] (local_optimizer.jl:1053) - i.e. it multiplies molality in
    twice (m^2*gamma instead of m*gamma). That looks like a bug unique to that one debug
    print statement (pH isn't computed anywhere else in MAGEMin.jl to cross-check it
    against), so this is a from-scratch correct implementation, not a faithful port.

    Only ever called with a fully-populated AQ_data (never obj_DEW's lightweight
    AQ_shim, which leaves em_names unset) since it needs to find "H+" by name.
*/
double DEW_aq_pH(          AQ_data      *AQ,
                            const double *x,
                            double        T,
                            double        P               ){

    int n_sp = AQ->n_sp;
    int id_Hp = -1;
    for (int i = 0; i < n_sp; i++){
        if (strcmp(AQ->em_names[i], "H+") == 0){ id_Hp = i; break; }
    }
    if (id_Hp == -1){ return NAN; }

    double *m = malloc(n_sp * sizeof(double));
    double gamma_e[5];
    double a_coef;
    DEW_aq_activity_state(AQ, x, T, P, m, gamma_e, &a_coef);

    int z_id = (int)fabs(AQ->z[id_Hp]);
    double pH = -log10(m[id_Hp]*gamma_e[z_id]);

    free(m);
    return pH;
}
