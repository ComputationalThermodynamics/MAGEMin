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
/*@
 **   Minimal C API to call MAGEMin as a library from external C/C++ code.
 @*/

#ifndef __MAGEMIN_API_H_
#define __MAGEMIN_API_H_

#include "MAGEMin.h"

#ifdef __cplusplus
extern "C" {
#endif

/** opaque-ish handle bundling the persistent state needed across calls */
typedef struct MAGEMin_Handles {
	global_variable gv;
	bulk_info       z_b;
	Databases       DB;
	simplex_data    splx_data;

	/* phases excluded from minimization, see MAGEMin_SetSuppressedPhases */
	int            *ss_suppressed;		/** [gv.len_ss] 1 = suppressed		*/
	int            *pp_suppressed;		/** [gv.len_pp] 1 = suppressed		*/
	int             n_suppressed;
} MAGEMin_Handle;

/**
 * Allocate and initialize a MAGEMin instance for the given thermodynamic
 * database (same acronyms as the CLI --db flag, e.g. "ig","mp","mb","mbe",
 * "um","ume","mtl","igd","igad","all","sb11","sb21","sb24",...). The research
 * group (CLI --rg) is inferred from the acronym: "sb11"/"sb21"/"sb24" -> "sb",
 * "xMELTS"/"pMELTS"/"rMELTS" -> "gh", "po" -> "br", anything else -> "tc".
 * verbose follows the CLI semantics (0 = silent minimization internals,
 * 1 = verbose); -1 has no special meaning here (only the CLI banner/progress
 * prints use it).
 * Returns NULL on failure.
 */
MAGEMin_Handle *MAGEMin_Init(	const char *database,
								int         verbose			);

/**
 * Like MAGEMin_Init, but with the chemical potential of the n oxides named in
 * mu_fix_oxides (MAGEMin_OxideNames, e.g. {"MgO","K2O"}) fixed directly instead
 * of solved for from the bulk (native mu-mu mechanism, CLI --n_mu_fix/
 * --mu_fix_idx). Each fixed oxide gets a fictive pure phase "m<oxide>",
 * appended as the last n entries of MAGEMin_PurePhaseNames, which appears in
 * the stable assemblage when its target is active. The targets are set with
 * MAGEMin_SetMuFix. The set of fixed oxides is structural (it sizes the pure-
 * phase list), so it is fixed for the handle's lifetime.
 *
 * Each fixed oxide's bulk content must be set generously in excess of what its
 * target implies, or the fictive phase may not become stable and the
 * oxide's chemical potential silently misses the target. The final Newton
 * step (MAGEMin_SetFinalNewtonStep) is skipped while any oxide is fixed.
 *
 * Returns NULL on failure, including an unknown or repeated oxide name.
 */
MAGEMin_Handle *MAGEMin_InitMuFix(	const char  *database,
									int          verbose,
									const char **mu_fix_oxides,
									int          n				);

/** number of oxides with a fixed chemical potential (0 unless MAGEMin_InitMuFix) */
int MAGEMin_NMuFix(				MAGEMin_Handle *h				);

/**
 * Set the target chemical potentials of the oxides fixed by MAGEMin_InitMuFix
 * (CLI --mu_fix_val): n = MAGEMin_NMuFix(h) values, same order and same units
 * as the chemical potentials reported in stb_system.gamma. Persists across
 * MAGEMin_ComputeEquilibrium calls until changed.
 *
 * Returns 0 on success, -1 if n does not match MAGEMin_NMuFix(h).
 */
int MAGEMin_SetMuFix(			MAGEMin_Handle *h,
								const double   *mu_fix_val,
								int             n				);

/**
 * Select how the elastic moduli of the phases are averaged into the system's
 * seismic properties: scheme 0 = Voigt-Reuss-Hill (default), 1 = Hashin-
 * Shtrikman; weight (0..1) positions the result between the soft and the stiff
 * bound of that scheme (VRH: 0 = Reuss, 0.5 = classic VRH (default), 1 =
 * Voigt; HS: 0 = HS-, 0.5 = HS average, 1 = HS+). Stays in effect for the
 * handle's lifetime.
 *
 * Returns 0 on success, -1 on an invalid scheme or weight.
 */
int MAGEMin_SetSeismicScheme(	MAGEMin_Handle *h,
								int             scheme,
								double          weight			);

/**
 * Select the minimization solver (CLI --solver): 0 = legacy, 1 = PGE + legacy
 * hybrid, 2 = hybrid PGE/LP (default), 3 = metastable calculation without
 * minimization (diagnostic). Stays in effect for the handle's lifetime.
 *
 * Returns 0 on success, 1 if the request was ignored because the database's
 * research group ("sb", "gh", "br") only supports the legacy solver (0 stays
 * active), -1 on an invalid solver value.
 */
int MAGEMin_SetSolver(			MAGEMin_Handle *h,
								int             solver			);

int MAGEMin_SetSSSolver(		MAGEMin_Handle *h,
								int             ss_solver		);

/**
 * Newton polish of the converged assemblage (CLI --final_Newton_step): 1 = active
 * (default), 0 = inactive. Applies to the "tc" and "sb" research groups.
 *
 * Returns 0 on success, 1 if the request has no effect for the database's
 * research group, -1 on an invalid value.
 */
int MAGEMin_SetFinalNewtonStep(	MAGEMin_Handle *h,
								int             final_Newton_step	);

/**
 * Store the final Newton system at the converged state (1) or not (0, default).
 * When stored, MAGEMin_GetFinalNewtonSystem returns the Jacobian (column-major,
 * nz x nz), the residual, the layout and the converged state of the last call.
 *
 * Returns 0 on success, -1 on an invalid value.
 */
int MAGEMin_SetStoreKKT(		MAGEMin_Handle *h,
								int             fn_store_kkt		);

/**
 * Final Newton system of the last MAGEMin_ComputeEquilibrium call, or NULL when
 * the polish was not accepted or storage is off. Owned by the handle, valid
 * until the next call.
 */
const fn_kkt *MAGEMin_GetFinalNewtonSystem(	MAGEMin_Handle *h			);

/**
 * Set (or clear) the oxygen buffer / fixed-activity constraint (CLI --buffer,
 * --buffer_n). Persists across MAGEMin_ComputeEquilibrium calls until changed.
 *
 * buffer   : a redox buffer "O2","qfm","mw","qif","nno","hm","iw","cco" (buffer_n
 *            is an offset in log units), a fixed activity "aH2O","aO2",
 *            "aMgO","aFeO","aAl2O3","aTiO2" (buffer_n is the activity), or
 *            "none"/NULL to clear. Must be available in the handle's database.
 *
 * Returns 0 on success, -1 on a name that is not one of the above or is
 * not available in the database.
 */
int MAGEMin_SetBuffer(			MAGEMin_Handle *h,
								const char     *buffer,
								double          buffer_n		);

/**
 * Exclude phases from minimization. names holds n solution-phase model names
 * (MAGEMin_SolutionPhaseNames) and/or pure-phase names (MAGEMin_PurePhaseNames);
 * the set replaces any previous one and persists across
 * MAGEMin_ComputeEquilibrium calls. Pass NULL/0 to clear it. The fictive
 * "m<oxide>" phases of MAGEMin_InitMuFix cannot be suppressed.
 *
 * While any phase is suppressed, both variants of the near-degenerate model
 * pairs gated by gv.mbCpx/mbIlm/mpSp/mpIlm are generated (set to 2), so that
 * suppressing one variant (e.g. "ilm" in "mp") lets the other ("ilmm") replace
 * it; otherwise they keep their default (0).
 *
 * Returns 0 on success, -1 if any name is unknown (the previous set is kept).
 */
int MAGEMin_SetSuppressedPhases(	MAGEMin_Handle  *h,
									const char     **names,
									int              n				);

/** number of solution-phase models of the handle's database */
int MAGEMin_NSolutionPhases(	MAGEMin_Handle *h				);

/** solution-phase model names. Owned by the handle, do not free. */
char **MAGEMin_SolutionPhaseNames(	MAGEMin_Handle *h			);

/** number of pure phases of the handle's database */
int MAGEMin_NPurePhases(		MAGEMin_Handle *h				);

/** pure-phase names. Owned by the handle, do not free. */
char **MAGEMin_PurePhaseNames(	MAGEMin_Handle *h				);

/** number of oxide/system components expected in the bulk[] array */
int MAGEMin_NOxides(			MAGEMin_Handle *h				);

/** names of the oxide/system components, in the order expected by bulk[] in
 *  MAGEMin_ComputeEquilibrium. Owned by the handle, do not free. */
char **MAGEMin_OxideNames(		MAGEMin_Handle *h				);

/**
 * Compute the stable equilibrium phase assemblage at one (P,T,bulk) point.
 *
 * P      : pressure [kbar]
 * T      : temperature [Celsius]
 * bulk   : MAGEMin_NOxides(h) values, in MAGEMin_OxideNames(h) order
 * sys_in : "mol" or "wt", composition unit of bulk[]
 *
 * Returns a pointer to the handle's internal stb_system. The pointer is
 * owned by the handle: it stays valid until the next call to
 * MAGEMin_ComputeEquilibrium or until MAGEMin_Free, so copy out whatever
 * fields you need before calling again.
 */
stb_system *MAGEMin_ComputeEquilibrium(	MAGEMin_Handle *h,
											double          P,
											double          T,
											const double   *bulk,
											const char     *sys_in		);

/** free everything allocated by MAGEMin_Init and the handle itself */
void MAGEMin_Free(				MAGEMin_Handle *h				);

#ifdef __cplusplus
}
#endif

#endif
