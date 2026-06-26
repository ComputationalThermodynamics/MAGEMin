/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Contributors : Nickolas B. Moccetti, Dominguez, H., Assunção J., Green E., Berlie N., and Rummel L.
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
} MAGEMin_Handle;

/**
 * Allocate and initialize a MAGEMin instance for the given thermodynamic
 * database (same acronyms as the CLI --db flag, e.g. "ig","mp","mb","mbe",
 * "um","ume","mtl","igd","igad",...). verbose follows the CLI semantics
 * (0 = silent minimization internals, 1 = verbose); -1 has no special
 * meaning here (only the CLI banner/progress prints use it).
 * Returns NULL on failure.
 */
MAGEMin_Handle *MAGEMin_Init(	const char *database,
								int         verbose			);

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
