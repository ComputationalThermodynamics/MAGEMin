/*
   Minimal example showing how to call MAGEMin as a library from external
   C/C++ code using the MAGEMin_api.h wrapper (src/MAGEMin_api.h/.c).

   Build (after "make lib" in the repo root):
     gcc -Isrc examples/call_as_library.c -L. -lMAGEMin \
         -lm -framework Accelerate /opt/homebrew/lib/libnlopt.dylib \
         -o call_as_library                      (macOS, adapt LIBS to your Makefile)

   Run:
     ./call_as_library
*/
#include <stdio.h>
#include "MAGEMin_api.h"

int main(void){

	MAGEMin_Handle *h = MAGEMin_Init("ig", 0);
	if (h == NULL){
		fprintf(stderr,"MAGEMin_Init failed\n");
		return 1;
	}

	int     n_ox  = MAGEMin_NOxides(h);
	char  **names = MAGEMin_OxideNames(h);

	printf("Oxide order (%d): ", n_ox);
	for (int i = 0; i < n_ox; i++){ printf("%s ",names[i]); }
	printf("\n\n");

	/* KLB1 peridotite bulk composition (mol%), Holland et al. 2018 / Green */
	double bulk[] = { 38.494, 1.776, 2.824, 50.566, 5.886, 0.01, 0.250, 0.10, 0.096, 0.109, 0.0 };

	double P_points[] = { 10.0, 20.0, 30.0 };
	double T_points[] = { 1100.0, 1200.0, 1300.0 };

	for (int i = 0; i < 3; i++){
		stb_system *sp = MAGEMin_ComputeEquilibrium(	h,
														P_points[i],
														T_points[i],
														bulk,
														"mol"			);
		if (sp == NULL){
			fprintf(stderr,"MAGEMin_ComputeEquilibrium failed\n");
			continue;
		}

		printf("P = %5.1f kbar, T = %6.1f C -> G = %+f, status = %d, %d phase(s):\n",
				P_points[i], T_points[i], sp->G, sp->status, sp->n_ph);
		for (int j = 0; j < sp->n_ph; j++){
			printf("    %-10s frac (mol) = %.5f\n", sp->ph[j], sp->ph_frac[j]);
		}
		printf("\n");
	}

	MAGEMin_Free(h);
	return 0;
}
