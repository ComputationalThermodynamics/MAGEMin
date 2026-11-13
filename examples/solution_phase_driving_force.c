/*
   Example: evaluating a single solution phase with the MAGEMin_api.h calls
   MAGEMin_SetModelVariant, MAGEMin_EvaluateSolutionPhase and
   MAGEMin_MinimizeSolutionPhase.

   Use case: driving force of garnet nucleation in an overstepped metapelite.
   The garnet-free matrix is computed with garnet suppressed; garnet is then
   minimised against the matrix chemical potentials (gamma). -df is the
   maximum driving force (parallel-tangent construction) and p the nucleus
   composition. The model-variant switches are fixed first, so that the
   suppressed and unsuppressed calculations use the same ilm/ilmm and sp/mt
   models.

   Build (after "make lib" in the repo root; Linux shown, see
   call_as_library.c for macOS):
     gcc -Isrc examples/solution_phase_driving_force.c -L. -lMAGEMin \
         -lm -llapacke -lnlopt -o solution_phase_driving_force

   Run:
     LD_LIBRARY_PATH=. ./solution_phase_driving_force
*/
#include <stdio.h>
#include "MAGEMin_api.h"

/* minimum of df over a small grid of starting points (several local minima
   are possible); returns 0 on success */
static int minimize_multistart(	MAGEMin_Handle *h, const char *name, const double *gamma, int n_ox,
								double *df_best, double *p_best, int n_em, int n_xeos			){
	double starts[3] = { 0.15, 0.50, 0.85 };
	int    n_start   = 1;
	for (int k = 0; k < n_xeos; k++){ n_start *= 3; }

	*df_best = 1e30;
	for (int s = 0; s < n_start; s++){
		double x0[n_xeos], df, p[n_em];
		int    status, code = s;
		for (int k = 0; k < n_xeos; k++){ x0[k] = starts[code % 3]; code /= 3; }
		if (MAGEMin_MinimizeSolutionPhase(h, name, gamma, n_ox, x0, n_xeos,
										  &df, NULL, p, NULL, &status) != 0){ return -1; }
		if (df < *df_best){
			*df_best = df;
			for (int k = 0; k < n_em; k++){ p_best[k] = p[k]; }
		}
	}
	return 0;
}

int main(void){

	MAGEMin_Handle *h = MAGEMin_Init("mp", 0);
	if (h == NULL){ fprintf(stderr,"MAGEMin_Init failed\n"); return 1; }

	int n_ox = MAGEMin_NOxides(h);

	/* FPWorldMedian pelite (mol), Forshaw & Pattison (2023), MAGEMin "mp" test 0:
	   SiO2 Al2O3 CaO MgO FeO K2O Na2O TiO2 O MnO H2O */
	double bulk[] = { 70.999, 12.8065, 0.771, 3.978, 6.342, 2.7895, 1.481, 0.758, 0.72933, 0.075, 30.0 };

	/* same model variants in every call below, with or without suppression */
	MAGEMin_SetModelVariant(h, NULL, 2);

	int n_em, n_xeos;
	MAGEMin_SolutionPhaseSize(h, "g", &n_em, &n_xeos);

	/* 1. consistency: against the equilibrium gamma of a garnet-bearing point
	      the minimum is df ~ 0 at the stable garnet composition */
	stb_system *sp = MAGEMin_ComputeEquilibrium(h, 8.0, 600.0, bulk, "mol");
	if (sp == NULL){ return 1; }
	char **em = MAGEMin_SolutionPhaseEndmemberNames(h, "g");	/* set up by the first equilibrium */
	double gamma_eq[n_ox], df, p[n_em];
	for (int k = 0; k < n_ox; k++){ gamma_eq[k] = sp->gamma[k]; }
	minimize_multistart(h, "g", gamma_eq, n_ox, &df, p, n_em, n_xeos);
	printf("8 kbar, 600 C, equilibrium gamma: df = %+.6f kJ/mol\n   p:", df);
	for (int k = 0; k < n_em; k++){ printf(" %s %.4f", em[k], p[k]); }
	printf("\n\n");

	/* absolute G of that composition, and its distance from the hyperplane */
	double G, comp[n_ox];
	int    sf_ok;
	MAGEMin_EvaluateSolutionPhase(h, "g", p, n_em, &G, NULL, comp, &sf_ok);
	double plane = 0.0;
	for (int k = 0; k < n_ox; k++){ plane += comp[k] * gamma_eq[k]; }
	printf("   G = %.4f kJ/mol, sum(comp*gamma) = %.4f, difference = %+.6f, sf_ok = %d\n\n",
			G, plane, G - plane, sf_ok);

	/* 2. overstep: garnet suppressed, matrix gamma, maximum driving force */
	const char *suppress[] = { "g" };
	MAGEMin_SetSuppressedPhases(h, suppress, 1);
	double T_points[] = { 550.0, 560.0, 570.0, 580.0 };
	for (int i = 0; i < 4; i++){
		sp = MAGEMin_ComputeEquilibrium(h, 6.0, T_points[i], bulk, "mol");
		if (sp == NULL){ continue; }
		double gamma_m[n_ox];
		for (int k = 0; k < n_ox; k++){ gamma_m[k] = sp->gamma[k]; }
		minimize_multistart(h, "g", gamma_m, n_ox, &df, p, n_em, n_xeos);
		printf("6 kbar, %5.1f C, garnet-free matrix: driving force -df = %+8.4f kJ/mol (%+7.1f J/mol O)\n   p:",
				T_points[i], -df, -df * 1000.0 / 12.0);
		for (int k = 0; k < n_em; k++){ printf(" %s %.4f", em[k], p[k]); }
		printf("\n");
	}

	MAGEMin_SetSuppressedPhases(h, NULL, 0);
	MAGEMin_SetModelVariant(h, NULL, -1);			/* back to automatic */
	MAGEMin_Free(h);
	return 0;
}
