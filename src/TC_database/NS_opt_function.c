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
#include <complex.h>

#if __APPLE__
	extern void dgesvd( char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info );
#else
	#include <lapacke.h>
#endif

#include "nlopt.h"
#include "../MAGEMin.h"
#include "objective_functions.h"
#include "NLopt_opt_function.h"
#include "NS_opt_function.h"
#include "../SB_database/sb_objective_functions.h"

#define NS_N_SAMP_ADD 		4
#define NS_RES_TOL 			1e-9
#define NS_GN_TOL 			1e-12
#define NS_RANK_TOL 		1e-10
#define NS_FIX_TOL 			1e-12
#define NS_ABSENT_EPS 		1e-12

static int ns_n_samp(	SS_ref *d	){
	return 2*d->n_em + NS_N_SAMP_ADD;
}

static int ns_rows_max(	SS_ref *d	){
	int n_samp = ns_n_samp(d);
	return (n_samp > d->ns_nc) ? n_samp : d->ns_nc;
}

static int ns_cols_max(	SS_ref *d	){
	return (d->ns_nc + 1 > d->n_em) ? d->ns_nc + 1 : d->n_em;
}

static int ns_svd(		SS_ref 		*d,
						int 		 m,
						int 		 n				){
	char 	jobu  = 'A';
	char 	jobvt = 'A';
	int 	lda   = (m > 1) ? m : 1;
	int 	ldvt  = (n > 1) ? n : 1;
	int 	lwork = d->ns_lwork;
	int 	info  = 0;

	#if __APPLE__
		dgesvd(&jobu, &jobvt, &m, &n, d->ns_svdA, &lda, d->ns_svdS, d->ns_svdU, &lda, d->ns_svdVt, &ldvt, d->ns_work, &lwork, &info);
	#else
		info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, jobu, jobvt, m, n, d->ns_svdA, lda, d->ns_svdS, d->ns_svdU, lda, d->ns_svdVt, ldvt, d->ns_work, lwork);
	#endif

	return info;
}

static int ns_rank(		SS_ref 		*d,
						int 		 m,
						int 		 n				){
	int 	mn = (m < n) ? m : n;
	int 	r  = 0;
	if (mn == 0){ return 0; }
	double 	tol = NS_RANK_TOL*((d->ns_svdS[0] > 1.0) ? d->ns_svdS[0] : 1.0);
	for (int k = 0; k < mn; k++){
		if (d->ns_svdS[k] > tol){ r += 1; }
	}
	return r;
}

static int ns_lstsq(	SS_ref 		*d,
						int 		 m,
						int 		 n,
						int 		 k,
						const double *Y,
						double 	   **B				){
	if (ns_svd(d, m, n) != 0){ return 0; }

	int 	mn  = (m < n) ? m : n;
	double 	tol = 1e-12*d->ns_svdS[0];
	double 	z[n];

	for (int c = 0; c < k; c++){
		for (int r = 0; r < n; r++){
			z[r] = 0.0;
			if (r < mn && d->ns_svdS[r] > tol){
				double a = 0.0;
				for (int i = 0; i < m; i++){ a += d->ns_svdU[i + r*m]*Y[i*k + c]; }
				z[r] = a/d->ns_svdS[r];
			}
		}
		for (int j = 0; j < n; j++){
			double a = 0.0;
			for (int r = 0; r < n; r++){ a += d->ns_svdVt[r + j*n]*z[r]; }
			B[j][c] = a;
		}
	}
	return 1;
}

double ns_x_of_sf(			global_variable 	 gv,
							SS_ref 				*d,
							const double 		*sf				){
	int 	n_sf = d->ns_nc;
	int 	n_x  = d->n_xeos;
	int 	n_em = d->n_em;

	if (d->ns_mode == 4){
		for (int j = 0; j < n_x; j++){ d->ns_x[j] = sf[j]; }
		return 0.0;
	}
	if (d->ns_mode == 2){
		for (int i = 0; i < n_sf; i++){ d->ns_sft[i] = sf[i]; }
		for (int j = 0; j < n_x; j++){
			double a = d->ns_M[n_sf][j];
			for (int i = 0; i < n_sf; i++){ a += sf[i]*d->ns_M[i][j]; }
			d->ns_x[j] = a;
		}
		return 0.0;
	}

	for (int e = 0; e < n_em; e++){
		double a = d->ns_M[n_sf][e];
		for (int i = 0; i < n_sf; i++){ a += sf[i]*d->ns_M[i][e]; }
		d->p[e]      = a;
		d->ns_pt[e]  = a;
	}
	d->ns_p2x(d, gv.bnd_val);
	for (int j = 0; j < n_x; j++){ d->ns_x[j] = d->iguess[j]; }

	return 0.0;
}

static double ns_p_residual(	SS_ref 		*d				){
	double res = 0.0;
	if (d->ns_mode == 4){ return 0.0; }
	if (d->ns_mode == 2){
		for (int i = 0; i < d->ns_nc; i++){
			double r = fabs(d->sf[i] - d->ns_sft[i]);
			if (!(r <= res)){ res = r; }
		}
		return res;
	}
	for (int e = 0; e < d->n_em; e++){
		double r = fabs(d->p[e] - d->ns_pt[e]);
		if (!(r <= res)){ res = r; }
	}
	return res;
}

static double ns_eval_fd(	global_variable 	 gv,
							SS_ref 				*d,
							const double 		*sf,
							double 				*grad			);

double ns_eval(				global_variable 	 gv,
							SS_ref 				*d,
							const double 		*sf,
							double 				*grad			){
	int 	n_sf = d->ns_nc;
	int 	n_x  = d->n_xeos;
	int 	n_em = d->n_em;

	if (grad != NULL && d->ns_fd == 1){ return ns_eval_fd(gv, d, sf, grad); }

	ns_x_of_sf(gv, d, sf);
	double G = d->ns_obj(n_x, d->ns_x, (grad != NULL && (d->ns_mode == 2 || d->ns_mode == 4)) ? d->ns_w : NULL, d);
	if (!(ns_p_residual(d) <= NS_RES_TOL)){ return NAN; }
	if (d->ns_mode == 3){
		for (int i = 0; i < n_sf; i++){ if (d->sf[i] < -1e-12){ return NAN; } }
	}

	if (grad != NULL && d->ns_mode == 4){
		for (int i = 0; i < n_sf; i++){ grad[i] = d->ns_w[i]; }
	}
	else if (grad != NULL && d->ns_mode == 2){
		for (int i = 0; i < n_sf; i++){
			double a = 0.0;
			for (int j = 0; j < n_x; j++){ a += d->ns_M[i][j]*d->ns_w[j]; }
			grad[i] = a;
		}
	}
	else if (grad != NULL){
		double c = d->df_raw/d->sum_apep;
		for (int e = 0; e < n_em; e++){ d->ns_w[e] = d->factor*(d->mu[e] - d->ape[e]*c); }
		for (int i = 0; i < n_sf; i++){
			double a = 0.0;
			for (int e = 0; e < n_em; e++){ a += d->ns_M[i][e]*d->ns_w[e]; }
			grad[i] = a;
		}
	}
	return G;
}

static double ns_eval_fd(	global_variable 	 gv,
							SS_ref 				*d,
							const double 		*sf,
							double 				*grad			){
	int 	n_sf  = d->ns_nc;
	int 	n_dir = d->ns_n_dir0;
	double 	sp[n_sf];

	for (int i = 0; i < n_sf; i++){ grad[i] = 0.0; }
	for (int k = 0; k < n_dir; k++){
		double h  = 1e-6;
		double hp = h, hm = h;
		for (int i = 0; i < n_sf; i++){
			double v = d->ns_N0[i][k];
			if (d->ns_sf_state[i] != 0 || v == 0.0){ continue; }
			double room_p = (v > 0.0) ? (d->ns_ub[i] - sf[i])/v : -sf[i]/v;
			double room_m = (v > 0.0) ? sf[i]/v : (sf[i] - d->ns_ub[i])/v;
			if (0.5*room_p < hp){ hp = 0.5*room_p; }
			if (0.5*room_m < hm){ hm = 0.5*room_m; }
		}
		if (!(hp > 0.0) && !(hm > 0.0)){ continue; }
		double Gp = 0.0, Gm = 0.0;
		if (hp > 0.0){
			for (int i = 0; i < n_sf; i++){ sp[i] = sf[i] + hp*d->ns_N0[i][k]; }
			Gp = ns_eval(gv, d, sp, NULL);
		}
		if (hm > 0.0){
			for (int i = 0; i < n_sf; i++){ sp[i] = sf[i] - hm*d->ns_N0[i][k]; }
			Gm = ns_eval(gv, d, sp, NULL);
		}
		double dk;
		if (hp > 0.0 && hm > 0.0){ dk = (Gp - Gm)/(hp + hm); }
		else if (hp > 0.0){
			double G0 = ns_eval(gv, d, sf, NULL);
			dk = (Gp - G0)/hp;
		}
		else{
			double G0 = ns_eval(gv, d, sf, NULL);
			dk = (G0 - Gm)/hm;
		}
		if (isnan(dk) || isinf(dk)){ return NAN; }
		for (int i = 0; i < n_sf; i++){ grad[i] += dk*d->ns_N0[i][k]; }
	}
	return ns_eval(gv, d, sf, NULL);
}

static void ns_proj(		SS_ref 		 *d,
							double 		**N,
							int 		  n_dir,
							const double *v,
							double 		 *out			){
	int 	n_sf = d->ns_nc;

	for (int k = 0; k < n_dir; k++){
		double a = 0.0;
		for (int i = 0; i < n_sf; i++){ a += N[i][k]*v[i]; }
		d->ns_z[k] = a;
	}
	for (int i = 0; i < n_sf; i++){
		double a = 0.0;
		for (int k = 0; k < n_dir; k++){ a += N[i][k]*d->ns_z[k]; }
		out[i] = a;
	}
}

static int ns_active_basis(	SS_ref 		*d				){
	int 	n_sf  = d->ns_nc;
	int 	n_dir = d->ns_n_dir0;
	int 	nC    = d->ns_n_C;

	if (nC == 0){
		for (int i = 0; i < n_sf; i++){
			for (int k = 0; k < n_dir; k++){ d->ns_N[i][k] = d->ns_N0[i][k]; }
		}
		return n_dir;
	}

	for (int k = 0; k < n_dir; k++){
		for (int q = 0; q < nC; q++){ d->ns_svdA[q + k*nC] = d->ns_N0[d->ns_C[q]][k]; }
	}
	if (ns_svd(d, nC, n_dir) != 0){ return -1; }

	int 	r  = ns_rank(d, nC, n_dir);
	int 	nn = n_dir - r;
	for (int k = 0; k < nn; k++){
		for (int i = 0; i < n_sf; i++){
			double a = 0.0;
			for (int l = 0; l < n_dir; l++){ a += d->ns_N0[i][l]*d->ns_svdVt[(r + k) + l*n_dir]; }
			d->ns_N[i][k] = a;
		}
	}
	return nn;
}

static int ns_release_index(	global_variable 	 gv,
								SS_ref 				*d				){
	int 	n_sf  = d->ns_nc;
	int 	n_dir = d->ns_n_dir0;
	int 	nC    = d->ns_n_C;
	double 	b[n_dir], lam[nC];

	if (nC == 0 || n_dir == 0){ return -1; }

	for (int k = 0; k < n_dir; k++){
		double a = 0.0;
		for (int i = 0; i < n_sf; i++){ a += d->ns_N0[i][k]*d->ns_g[i]; }
		b[k] = a;
		for (int q = 0; q < nC; q++){ d->ns_svdA[k + q*n_dir] = d->ns_N0[d->ns_C[q]][k]; }
	}
	if (ns_svd(d, n_dir, nC) != 0){ return -1; }

	int 	r = ns_rank(d, n_dir, nC);
	for (int q = 0; q < nC; q++){ lam[q] = 0.0; }
	for (int k = 0; k < r; k++){
		double c = 0.0;
		for (int i = 0; i < n_dir; i++){ c += d->ns_svdU[i + k*n_dir]*b[i]; }
		c /= d->ns_svdS[k];
		for (int q = 0; q < nC; q++){ lam[q] += d->ns_svdVt[k + q*nC]*c; }
	}

	int 	best = -1;
	double 	lmin = -100.0*gv.ns_tol_pg;
	for (int q = 0; q < nC; q++){
		double l = (d->ns_sf[d->ns_C[q]] > 0.5*d->ns_ub[d->ns_C[q]]) ? -lam[q] : lam[q];
		if (l < lmin){ lmin = l; best = q; }
	}
	return best;
}

static int ns_init_phase(	global_variable 	 gv,
							SS_ref 				*d				);

static int ns_absent_comp(	global_variable 	 gv,
							bulk_info 	 		 z_b,
							SS_ref 				*d,
							const double 		*x,
							double 				*c				){
	int 	n_c = 0;
	d->ns_obj(d->n_xeos, x, NULL, d);
	for (int j = 0; j < gv.len_ox; j++){
		if (z_b.bulk_rock[j] != 0.0){ continue; }
		int 	in = 0;
		double 	a  = 0.0;
		for (int k = 0; k < d->n_em; k++){
			if (fabs(d->Comp[k][j]) > 1e-12){ in = 1; }
			a += d->p[k]*d->Comp[k][j];
		}
		if (in){ c[n_c] = a; n_c += 1; }
	}
	return n_c;
}

static void ns_reduce_x(	global_variable 	 gv,
							bulk_info 	 		 z_b,
							SS_ref 				*d				){
	int 	n_sf   = d->ns_nc;
	int 	n_x    = d->n_xeos;
	int 	n_free = 0;
	int 	n_dir  = 0;
	int 	fr[n_x];
	double 	xc[n_x];
	double 	xi[n_x];
	double 	Z[n_x][n_x];

	for (int j = 0; j < n_x; j++){
		if (d->bounds_ref[j][1] - d->bounds_ref[j][0] > 1e-14){ fr[n_free] = j; n_free += 1; }
	}
	for (int j = 0; j < n_x; j++){
		double lo = d->bounds_ref[j][0], hi = d->bounds_ref[j][1];
		xc[j] = (hi - lo > 1e-14) ? lo + (hi - lo)/(double)(n_free + 1) : lo;
		xi[j] = (hi - lo > 1e-14) ? 0.5*(lo + hi) : lo;
	}
	for (int l = 0; l < n_free; l++){
		for (int k = 0; k < n_free; k++){ Z[l][k] = (l == k) ? 1.0 : 0.0; }
	}
	n_dir = n_free;

	if (d->ns_absent == 1 && n_free > 0){
		int 	n_ox = gv.len_ox;
		double 	c0[n_ox], ch[n_ox], K[n_ox][n_x], xh[n_x];
		int 	n_c  = ns_absent_comp(gv, z_b, d, xc, c0);
		for (int l = 0; l < n_free; l++){
			double h = 0.25*(d->bounds_ref[fr[l]][1] - d->bounds_ref[fr[l]][0]);
			for (int j = 0; j < n_x; j++){ xh[j] = xc[j]; }
			xh[fr[l]] += h;
			ns_absent_comp(gv, z_b, d, xh, ch);
			for (int q = 0; q < n_c; q++){ K[q][l] = (ch[q] - c0[q])/h; }
		}
		for (int l = 0; l < n_free; l++){
			for (int k = 0; k < n_free; k++){
				double a = 0.0;
				for (int q = 0; q < n_c; q++){ a += K[q][l]*K[q][k]; }
				d->ns_svdA[l + k*n_free] = a;
			}
		}
		if (n_c > 0){
			if (ns_svd(d, n_free, n_free) != 0){ d->ns_ok = 0; return; }
			int 	rk = ns_rank(d, n_free, n_free);
			double 	g[n_free];
			for (int l = 0; l < n_free; l++){
				double a = 0.0;
				for (int q = 0; q < n_c; q++){ a += K[q][l]*c0[q]; }
				g[l] = a;
			}
			for (int r = 0; r < rk; r++){
				double a = 0.0;
				for (int l = 0; l < n_free; l++){ a += d->ns_svdVt[r + l*n_free]*g[l]; }
				a /= d->ns_svdS[r];
				for (int l = 0; l < n_free; l++){ xc[fr[l]] -= a*d->ns_svdVt[r + l*n_free]; }
			}
			n_dir = n_free - rk;
			for (int k = 0; k < n_dir; k++){
				for (int l = 0; l < n_free; l++){ Z[l][k] = d->ns_svdVt[(rk + k) + l*n_free]; }
			}
			n_c = ns_absent_comp(gv, z_b, d, xc, c0);
			for (int q = 0; q < n_c; q++){
				if (fabs(c0[q]) > 1e-9){ d->ns_ok = 0; return; }
			}
			for (int j = 0; j < n_x; j++){ xc[j] = (1.0 - NS_ABSENT_EPS)*xc[j] + NS_ABSENT_EPS*xi[j]; }
		}
	}

	for (int i = 0; i < n_sf; i++){
		double a = d->ns_V[n_x][i];
		for (int j = 0; j < n_x; j++){ a += d->ns_V[j][i]*xc[j]; }
		if (d->ns_absent == 1 && a < -1e-9){ d->ns_ok = 0; return; }
		d->ns_sf0[i] = a;
		d->ns_sfc[i] = a;
	}
	d->ns_n_em_act = -1;
	if (n_dir == 0){
		for (int i = 0; i < n_sf; i++){ d->ns_sf_state[i] = 1; }
		d->ns_n_dir0 = 0;
		d->ns_ok     = 1;
		return;
	}
	for (int i = 0; i < n_sf; i++){
		double cmax = 0.0;
		for (int k = 0; k < n_dir; k++){
			double a = 0.0;
			for (int l = 0; l < n_free; l++){ a += d->ns_V[fr[l]][i]*Z[l][k]; }
			d->ns_svdA[k + i*n_dir] = a;
			if (fabs(a) > cmax){ cmax = fabs(a); }
		}
		d->ns_sf_state[i] = (cmax < NS_FIX_TOL) ? 1 : 0;
	}
	if (ns_svd(d, n_dir, n_sf) != 0){ d->ns_ok = 0; return; }
	int 	rk = ns_rank(d, n_dir, n_sf);
	for (int i = 0; i < n_sf; i++){
		for (int k = 0; k < rk; k++){ d->ns_N0[i][k] = d->ns_svdVt[k + i*n_sf]; }
	}
	d->ns_n_dir0 = rk;
	d->ns_ok     = 1;
}

static int ns_has_absent(	global_variable 	 gv,
							bulk_info 	 		 z_b,
							SS_ref 				*d				){
	for (int k = 0; k < d->n_em; k++){
		if (d->z_em[k] != 1.0){ continue; }
		for (int j = 0; j < gv.len_ox; j++){
			if (z_b.bulk_rock[j] == 0.0 && fabs(d->Comp[k][j]) > 1e-12){ return 1; }
		}
	}
	return 0;
}

static void ns_reduce_phase(	global_variable 	 gv,
								bulk_info 	 		 z_b,
								SS_ref 				*d				){
	if (d->ns_ok_init == 0){
		d->ns_ok_init = ns_init_phase(gv, d) ? 1 : -1;
		if (gv.verbose == 1 && d->ns_ok_init != 1){
			printf(" NS: solution phase model not eligible for nullspace minimization (NLopt only)\n");
		}
	}
	d->ns_absent = (gv.ss_solver > 0 && gv.solver == 0) ? ns_has_absent(gv, z_b, d) : 0;
	if (d->ns_mode == 4 && gv.ss_solver > 0 && gv.solver == 0 && d->ns_ok_init == 1){
		for (int k = 0; k < d->n_em; k++){
			if (d->z_em[k] != 1.0 || d->bounds_ref[k][1] - d->bounds_ref[k][0] <= 1e-14){ d->ns_absent = 1; }
		}
	}
	if (d->ns_absent == 1){
		for (int j = 0; j < d->n_xeos; j++){
			if (d->bounds_ref[j][1] - d->bounds_ref[j][0] <= 1e-14){ continue; }
			for (int b = 0; b < 2; b++){
				double v = d->bounds_ref[j][b];
				double r = round(v);
				if (v != r && fabs(v - r) <= 2.0*gv.bnd_val){ d->bounds_ref[j][b] = r; }
			}
		}
	}
	if (d->ns_ok_init != 1){ d->ns_ok = 0; return; }
	if (d->ss_flags[0] == 0){ d->ns_ok = 0; return; }
	if (d->ns_mode == 2){ ns_reduce_x(gv, z_b, d); return; }

	int 	n_sf  = d->ns_nc;
	int 	n_em  = d->n_em;
	int 	R[n_em];
	int 	nR    = 0;

	for (int k = 0; k < n_em; k++){
		int ok = (d->z_em[k] == 1.0);
		if (d->ns_mode == 4 && d->bounds_ref[k][1] - d->bounds_ref[k][0] <= 1e-14){ ok = 0; }
		for (int j = 0; j < gv.len_ox && ok; j++){
			if (z_b.bulk_rock[j] == 0.0 && fabs(d->Comp[k][j]) > 1e-12){ ok = 0; }
		}
		if (ok){ R[nR] = k; nR += 1; }
	}

	int 	same = (nR == d->ns_n_em_act);
	for (int r = 0; r < nR && same; r++){
		if (d->ns_em_act[r] != R[r]){ same = 0; }
	}
	if (same && d->ns_absent == 0){
		d->ns_ok = (nR > 0);
		return;
	}

	d->ns_n_em_act = nR;
	for (int r = 0; r < nR; r++){ d->ns_em_act[r] = R[r]; }

	if (nR == 0){
		d->ns_ok     = 0;
		d->ns_n_dir0 = 0;
		return;
	}

	for (int i = 0; i < n_sf; i++){
		d->ns_sf0[i] = d->ns_V[R[0]][i];
		d->ns_sfc[i] = 0.0;
		for (int r = 0; r < nR; r++){ d->ns_sfc[i] += d->ns_V[R[r]][i]/(double)nR; }
		double dmax = 0.0;
		for (int r = 1; r < nR; r++){
			double dv = fabs(d->ns_V[R[r]][i] - d->ns_V[R[0]][i]);
			if (dv > dmax){ dmax = dv; }
		}
		d->ns_sf_state[i] = (dmax < NS_FIX_TOL) ? 1 : 0;
	}
	if (d->ns_absent == 1){
		int 	nz = 0;
		double 	cf[n_sf];
		for (int i = 0; i < n_sf; i++){ cf[i] = 0.0; }
		for (int k = 0; k < n_em; k++){
			if (d->z_em[k] != 1.0 && d->ns_mode != 4){ continue; }
			for (int i = 0; i < n_sf; i++){ cf[i] += d->ns_V[k][i]; }
			nz += 1;
		}
		for (int i = 0; i < n_sf; i++){
			double c = cf[i]/(double)nz;
			d->ns_sf0[i] = (1.0 - NS_ABSENT_EPS)*d->ns_sf0[i] + NS_ABSENT_EPS*c;
			d->ns_sfc[i] = (1.0 - NS_ABSENT_EPS)*d->ns_sfc[i] + NS_ABSENT_EPS*c;
		}
	}

	if (nR == 1){
		d->ns_n_dir0 = 0;
		d->ns_ok     = 1;
		return;
	}

	int 	m = nR - 1;
	for (int i = 0; i < n_sf; i++){
		for (int r = 0; r < m; r++){
			d->ns_svdA[r + i*m] = d->ns_V[R[r + 1]][i] - d->ns_V[R[0]][i];
		}
	}
	if (ns_svd(d, m, n_sf) != 0){
		d->ns_ok = 0;
		d->ns_n_em_act = -1;
		return;
	}

	int 	rk = ns_rank(d, m, n_sf);
	for (int i = 0; i < n_sf; i++){
		for (int k = 0; k < rk; k++){ d->ns_N0[i][k] = d->ns_svdVt[k + i*n_sf]; }
	}
	d->ns_n_dir0 = rk;
	d->ns_ok     = 1;
}

void ns_reduce_system(		global_variable 	 gv,
							bulk_info 	 		 z_b,
							SS_ref 				*d				){
	ns_reduce_phase(gv, z_b, d);
	if (d->ns_absent != 1 || d->ns_ok_init != 1){ return; }
	for (int k = 0; k < d->n_em; k++){
		for (int j = 0; j < gv.len_ox; j++){
			if (z_b.bulk_rock[j] == 0.0 && fabs(d->Comp[k][j]) > 1e-12){ d->z_em[k] = 0.0; }
		}
	}
}

static void ns_zgrad(		SS_ref 		 *d,
							int 		  n_dir,
							const double *g,
							double 		 *gz			){
	int 	n_sf = d->ns_nc;
	for (int k = 0; k < n_dir; k++){
		double a = 0.0;
		for (int i = 0; i < n_sf; i++){ a += d->ns_N[i][k]*g[i]; }
		gz[k] = a;
	}
}

static void ns_restart_H(	SS_ref 		 *d,
							int 		  n_dir,
							double 		  h0				){
	int 	 n_sf = d->ns_nc;
	double **H    = d->ns_BkI;
	for (int a = 0; a < n_dir; a++){
		for (int b = 0; b <= a; b++){
			double v = 0.0;
			for (int i = 0; i < n_sf; i++){
				if (d->ns_sf_state[i] == 0){ v += d->ns_N[i][a]*d->ns_sf[i]*d->ns_N[i][b]; }
			}
			H[a][b] = h0*v;
			H[b][a] = h0*v;
		}
	}
}

static void ns_start_sf(	global_variable 	 gv,
							SS_ref 				*d				){
	int 	 n_sf   = d->ns_nc;
	double 	*sf     = d->ns_sf;
	double 	 eps_b  = gv.ns_eps_bnd;
	double 	 sfi[n_sf];
	double 	 r[n_sf];

	d->ns_obj(d->n_xeos, d->iguess, NULL, d);
	for (int i = 0; i < n_sf; i++){ sfi[i] = (d->ns_mode == 3 || d->ns_mode == 4) ? ((i < d->n_em) ? d->p[i] : 0.0) : d->sf[i]; }
	for (int i = 0; i < n_sf; i++){ r[i] = sfi[i] - d->ns_sf0[i]; }
	ns_proj(d, d->ns_N0, d->ns_n_dir0, r, sf);
	for (int i = 0; i < n_sf; i++){
		sf[i] += d->ns_sf0[i];
		if (d->ns_sf_state[i] == 1){ sf[i] = d->ns_sf0[i]; }
	}

	double theta = 0.0;
	double lo    = 10.0*eps_b;
	for (int i = 0; i < n_sf; i++){
		if (d->ns_sf_state[i] == 1){ continue; }
		if (sf[i] < lo && d->ns_sfc[i] > sf[i]){
			double t = (lo - sf[i])/(d->ns_sfc[i] - sf[i]);
			if (t > theta){ theta = t; }
		}
		double hi = d->ns_ub[i] - 10.0*eps_b;
		if (sf[i] > hi && d->ns_sfc[i] < sf[i]){
			double t = (sf[i] - hi)/(sf[i] - d->ns_sfc[i]);
			if (t > theta){ theta = t; }
		}
	}
	if (theta > 1.0){ theta = 1.0; }
	for (int i = 0; i < n_sf; i++){
		if (d->ns_sf_state[i] == 0){ sf[i] = (1.0 - theta)*sf[i] + theta*d->ns_sfc[i]; }
	}
}

static int ns_minimize(		global_variable 	 gv,
							SS_ref 				*d				){
	int 	 n_sf   = d->ns_nc;
	double 	*sf     = d->ns_sf;
	double 	**BkI   = d->ns_BkI;
	double 	*gz0    = d->ns_pg0;
	double 	*gz1    = d->ns_pg1;
	double 	*pz     = d->ns_z;
	double 	*dsf    = d->ns_pk;
	double 	*sz     = d->ns_sk;
	double 	*yz     = d->ns_yk;
	double 	*vt     = d->ns_vt;
	double 	*ny     = d->ns_r;
	double 	 eps    = gv.ns_eps;
	double 	 eps_b  = gv.ns_eps_bnd;

	ns_start_sf(gv, d);

	d->ns_ite = 0;
	if (d->ns_n_dir0 == 0){
		for (int i = 0; i < n_sf; i++){ sf[i] = d->ns_sf0[i]; }
		double G0 = ns_eval(gv, d, sf, d->ns_g);
		return (isnan(G0) || isinf(G0)) ? 103 : 3;
	}

	d->ns_n_C   = 0;
	d->ns_n_dir = ns_active_basis(d);
	int 	n_dir = d->ns_n_dir;

	double G = ns_eval(gv, d, sf, d->ns_g);
	if (isnan(G) || isinf(G)){ return 103; }
	ns_zgrad(d, n_dir, d->ns_g, gz0);

	double h0   = 1.0;
	double gmax = 0.0;
	for (int i = 0; i < n_sf; i++){
		double a = 0.0;
		for (int k = 0; k < n_dir; k++){ a += d->ns_N[i][k]*gz0[k]; }
		if (fabs(sf[i]*a) > gmax){ gmax = fabs(sf[i]*a); }
	}
	if (gmax > 0.0){ h0 = 0.1/gmax; }

	int 	o         = 1;
	int 	ite       = 0;
	int 	red       = 1;
	int 	converged = 0;
	int 	n_rel     = 0;
	int 	reset_h   = 1;

	while (red){
		double G_out = G;
		double dxo   = 1.0;
		converged    = 0;

		while (dxo > gv.ns_tol && o < gv.ns_omax && ite < gv.ns_max_ite){
			double beta = (o == 1) ? gv.ns_beta_ini : gv.ns_beta;
			if (reset_h){
				ns_restart_H(d, n_dir, h0);
				reset_h = 0;
			}
			double dxi = 1.0;
			int    it  = 1;

			while (dxi > gv.ns_tol && it < gv.ns_imax && ite < gv.ns_max_ite){
				for (int k = 0; k < n_dir; k++){
					double a = 0.0;
					for (int l = 0; l < n_dir; l++){ a += BkI[k][l]*gz0[l]; }
					pz[k] = a;
				}
				for (int i = 0; i < n_sf; i++){
					double a = 0.0;
					for (int k = 0; k < n_dir; k++){ a += d->ns_N[i][k]*pz[k]; }
					dsf[i] = a;
				}

				double dl_min = 1.0;
				for (int i = 0; i < n_sf; i++){
					if (d->ns_sf_state[i] != 0){ continue; }
					if (dsf[i] > 0.0){
						double dl = (sf[i] - eps_b)/dsf[i];
						if (dl < dl_min){ dl_min = dl; }
					}
					else if (dsf[i] < 0.0){
						double dl = (sf[i] - (d->ns_ub[i] - eps_b))/dsf[i];
						if (dl < dl_min){ dl_min = dl; }
					}
				}
				double scale = (dl_min <= 0.0) ? 0.0 : ((gv.ns_frac_bnd*dl_min < 1.0) ? gv.ns_frac_bnd*dl_min : 1.0);
				double stp   = 0.0;
				for (int k = 0; k < n_dir; k++){ stp -= pz[k]*scale*gz0[k]; }
				if (stp < -1.0){ stp = -1.0; }

				double t  = 1.0;
				for (int i = 0; i < n_sf; i++){ d->ns_trial[i] = sf[i] - t*scale*dsf[i]; }
				double G1 = ns_eval(gv, d, d->ns_trial, NULL);
				int    ls = 0;
				while ((G1 > G + gv.ns_alpha*t*stp || isnan(G1)) && ls < 60){
					t  *= beta;
					ls += 1;
					for (int i = 0; i < n_sf; i++){ d->ns_trial[i] = sf[i] - t*scale*dsf[i]; }
					G1 = ns_eval(gv, d, d->ns_trial, NULL);
				}
				if (isnan(G1)){ return 103; }
				for (int i = 0; i < n_sf; i++){ sf[i] = d->ns_trial[i]; }

				double Gn = ns_eval(gv, d, sf, d->ns_g);
				if (isnan(Gn) || isinf(Gn)){ return 103; }
				ns_zgrad(d, n_dir, d->ns_g, gz1);

				double sy = 0.0, ss = 0.0, yvt = 0.0;
				for (int k = 0; k < n_dir; k++){
					sz[k] = -t*scale*pz[k];
					yz[k] = gz1[k] - gz0[k];
					ss   += sz[k]*sz[k];
					sy   += sz[k]*yz[k];
				}
				double yy = 0.0;
				double pgm = 0.0;
				for (int i = 0; i < n_sf; i++){
					double a = 0.0, b = 0.0;
					for (int k = 0; k < n_dir; k++){ a += d->ns_N[i][k]*yz[k]; b += d->ns_N[i][k]*gz1[k]; }
					ny[i] = a;
					if (d->ns_sf_state[i] == 0){ yy += a*sf[i]*a; }
					if (fabs(b) > pgm){ pgm = fabs(b); }
				}
				if (sy > 0.0 && yy > 0.0){ h0 = sy/yy; }
				if (ss > 0.0 && sy > 1e-300){
					for (int k = 0; k < n_dir; k++){
						double a = 0.0;
						for (int l = 0; l < n_dir; l++){ a += BkI[k][l]*yz[l]; }
						vt[k] = a;
						yvt  += yz[k]*a;
					}
					double c1 = (sy + yvt)/(sy*sy);
					for (int k = 0; k < n_dir; k++){
						for (int l = 0; l < n_dir; l++){
							BkI[k][l] += c1*sz[k]*sz[l] - (vt[k]*sz[l] + sz[k]*vt[l])/sy;
						}
					}
				}

				double gs  = (fabs(G) > 1.0) ? fabs(G) : 1.0;
				dxi = (pgm < gv.ns_tol_pg) ? 0.0 : fabs(Gn - G)/gs;
				if (dxi < gv.ns_tol && pgm > 100.0*gv.ns_tol_pg){
					dxi = 1.0;
					if (it > 3 && fabs(Gn - G) < 1e-14*gs){ dxi = 0.0; }
				}

				G = Gn;
				for (int k = 0; k < n_dir; k++){ gz0[k] = gz1[k]; }
				it  += 1;
				ite += 1;
			}
			double gs = (fabs(G_out) > 1.0) ? fabs(G_out) : 1.0;
			dxo       = fabs(G - G_out)/gs;
			G_out     = G;
			o        += 1;
			if (dxi <= gv.ns_tol){ converged = 1; }
		}

		int 	added = 0;
		for (int i = 0; i < n_sf; i++){
			if (d->ns_sf_state[i] != 0){ continue; }
			if ((sf[i] < eps && sf[i] != 0.0) || (sf[i] > d->ns_ub[i] - eps && sf[i] != d->ns_ub[i])){
				d->ns_sf_state[i]       = 2;
				d->ns_C[d->ns_n_C]      = i;
				d->ns_n_C              += 1;
				added                  += 1;
			}
		}

		red = 0;
		if (added > 0){
			int nn = ns_active_basis(d);
			if (nn > 0 && nn == d->ns_n_dir0 - d->ns_n_C){
				n_dir       = nn;
				d->ns_n_dir = nn;
				red         = 1;
				reset_h     = 1;
				ns_zgrad(d, n_dir, d->ns_g, gz0);
			}
		}
		else if (converged && d->ns_n_C > 0 && n_rel < 2*n_sf && o < gv.ns_omax && ite < gv.ns_max_ite){
			int q = ns_release_index(gv, d);
			if (q >= 0){
				d->ns_sf_state[d->ns_C[q]] = 0;
				for (int k = q; k < d->ns_n_C - 1; k++){ d->ns_C[k] = d->ns_C[k + 1]; }
				d->ns_n_C  -= 1;
				n_rel      += 1;
				int nn      = ns_active_basis(d);
				n_dir       = nn;
				d->ns_n_dir = nn;
				red         = 1;
				reset_h     = 1;
				ns_zgrad(d, n_dir, d->ns_g, gz0);
			}
		}
	}

	for (int q = 0; q < d->ns_n_C; q++){ d->ns_sf_state[d->ns_C[q]] = 0; }
	d->ns_ite = ite;

	return converged ? 3 : 102;
}

static int ns_out_valid(		global_variable 	 gv,
							SS_ref 				*d,
							double 				 theta			){
	int 	n_sf = d->ns_nc;
	int 	n_x  = d->n_xeos;
	double 	sft[n_sf];

	for (int i = 0; i < n_sf; i++){
		sft[i] = (d->ns_sf_state[i] == 0) ? (1.0 - theta)*d->ns_sf[i] + theta*d->ns_sfc[i] : d->ns_sf[i];
	}
	ns_x_of_sf(gv, d, sft);
	for (int j = 0; j < n_x; j++){
		double tol = (d->ns_absent == 1 && d->bounds_ref[j][1] - d->bounds_ref[j][0] <= 1e-14) ? 2.0*gv.bnd_val : 1e-9;
		if (d->ns_x[j] < d->bounds_ref[j][0] - tol || d->ns_x[j] > d->bounds_ref[j][1] + tol){ return 0; }
	}
	for (int j = 0; j < n_x; j++){
		double v = d->ns_x[j];
		int    pin = (d->ns_mode == 4 && d->bounds_ref[j][1] - d->bounds_ref[j][0] <= 1e-14);
		if (!pin && v < d->bounds_ref[j][0]){ v = d->bounds_ref[j][0]; }
		if (!pin && v > d->bounds_ref[j][1]){ v = d->bounds_ref[j][1]; }
		d->xeos[j] = v;
	}
	d->df = d->ns_obj(n_x, d->xeos, NULL, d);
	return (!isnan(d->df) && !isinf(d->df) && ns_p_residual(d) <= ((d->ns_absent == 1) ? 2.0*gv.bnd_val : 1e-8));
}

static int ns_out_fit(		global_variable 	 gv,
							SS_ref 				*d				){
	double 	lo    = 0.0;
	double 	hi    = -1.0;

	if (ns_out_valid(gv, d, 0.0)){ return 1; }

	double theta = 1e-9;
	for (int k = 0; k < 8 && hi < 0.0; k++){
		if (ns_out_valid(gv, d, theta)){ hi = theta; }
		else{ lo = theta; theta *= 10.0; }
	}
	if (hi < 0.0){ return 0; }
	for (int k = 0; k < 30; k++){
		double mid = 0.5*(lo + hi);
		if (ns_out_valid(gv, d, mid)){ hi = mid; } else { lo = mid; }
	}
	return ns_out_valid(gv, d, hi);
}

static int ns_out_clamp(		global_variable 	 gv,
							SS_ref 				*d				){
	int 	n_x = d->n_xeos;

	ns_x_of_sf(gv, d, d->ns_sf);
	for (int j = 0; j < n_x; j++){
		double v = d->ns_x[j];
		int    pin = (d->ns_mode == 4 && d->bounds_ref[j][1] - d->bounds_ref[j][0] <= 1e-14);
		if (!pin && v < d->bounds_ref[j][0]){ v = d->bounds_ref[j][0]; }
		if (!pin && v > d->bounds_ref[j][1]){ v = d->bounds_ref[j][1]; }
		d->xeos[j] = v;
	}
	d->df = d->ns_obj(n_x, d->xeos, NULL, d);
	if (isnan(d->df) || isinf(d->df)){ return 0; }
	for (int i = 0; i < d->ns_nc && d->ns_mode != 4; i++){
		if (!(d->sf[i] >= -1e-12)){ return 0; }
	}
	return 1;
}

int ns_project_x(				global_variable 	 gv,
								SS_ref 				*d,
								double 				*x				){
	int 	n_x = d->n_xeos;
	double 	xs[n_x];

	if (d->ns_ok != 1){ return 0; }
	for (int j = 0; j < n_x; j++){
		xs[j]           = d->iguess[j];
		d->iguess[j]    = x[j];
		d->bounds[j][0] = d->bounds_ref[j][0] - 1.0;
		d->bounds[j][1] = d->bounds_ref[j][1] + 1.0;
	}
	ns_start_sf(gv, d);
	int 	ok = ns_out_fit(gv, d);
	for (int j = 0; j < n_x; j++){
		if (ok){ x[j] = d->xeos[j]; }
		d->iguess[j]    = xs[j];
		d->bounds[j][0] = d->bounds_ref[j][0];
		d->bounds[j][1] = d->bounds_ref[j][1];
	}
	return ok;
}

int ns_pc_mode(					global_variable 	 gv,
								SS_ref 				*d				){
	if (!(gv.ss_solver > 0 && gv.solver == 0) || d->ns_absent != 1 || d->ns_ok_init != 1){ return 0; }
	return (d->ns_ok == 1) ? 1 : 2;
}

SS_ref NS_opt_function(		global_variable 	 gv,
							SS_ref 				 SS_ref_db		){
	SS_ref *d = &SS_ref_db;

	d->ns_ite = 0;
	if (d->ns_ok != 1){
		d->status    = 100;
		d->ns_status = 100;
		return SS_ref_db;
	}

	int 	n_x = d->n_xeos;
	for (int j = 0; j < n_x; j++){
		d->ns_x0[j]     = d->iguess[j];
		d->bounds[j][0] = d->bounds_ref[j][0] - 1.0;
		d->bounds[j][1] = d->bounds_ref[j][1] + 1.0;
	}

	int 	status = ns_minimize(gv, d);

	if (status == 3 && !ns_out_fit(gv, d) && !(d->ns_absent == 1 && ns_out_clamp(gv, d))){ status = 103; }
	if (status != 3 && d->ns_absent == 1){
		for (int j = 0; j < n_x; j++){ d->xeos[j] = d->ns_x0[j]; }
		d->df = d->ns_obj(n_x, d->xeos, NULL, d);
	}

	for (int j = 0; j < n_x; j++){
		d->iguess[j]    = (status == 3) ? d->xeos[j] : d->ns_x0[j];
		d->bounds[j][0] = d->bounds_ref[j][0];
		d->bounds[j][1] = d->bounds_ref[j][1];
	}

	d->status    = status;
	d->ns_status = status;

	return SS_ref_db;
}

static double **ns_alloc_mat(	int 	nr,
								int 	nc			){
	double **A = malloc(nr*sizeof(double*));
	for (int i = 0; i < nr; i++){ A[i] = calloc(nc, sizeof(double)); }
	return A;
}

static void ns_free_mat(		double **A,
								int 	 nr			){
	if (A == NULL){ return; }
	for (int i = 0; i < nr; i++){ free(A[i]); }
	free(A);
}

static void ns_alloc(		global_variable 	 gv,
							SS_ref 				*d				){
	int 	n_sf   = d->ns_nc;
	int 	n_x    = d->n_xeos;
	int 	n_em   = d->n_em;
	int 	n_samp = ns_n_samp(d);
	int 	mr     = ns_rows_max(d);
	int 	mc     = ns_cols_max(d);

	d->ns_V        = ns_alloc_mat(n_em, n_sf);
	d->ns_M        = ns_alloc_mat(n_sf + 1, n_em);
	d->ns_N0       = ns_alloc_mat(n_sf, n_sf);
	d->ns_N        = ns_alloc_mat(n_sf, n_sf);
	d->ns_BkI      = ns_alloc_mat(n_sf, n_sf);

	d->ns_sf0      = calloc(n_sf, sizeof(double));
	d->ns_sfc      = calloc(n_sf, sizeof(double));
	d->ns_sf       = calloc(n_sf, sizeof(double));
	d->ns_sf_prev  = calloc(n_sf, sizeof(double));
	d->ns_trial    = calloc(n_sf, sizeof(double));
	d->ns_g        = calloc(n_sf, sizeof(double));
	d->ns_pg0      = calloc(n_sf, sizeof(double));
	d->ns_pg1      = calloc(n_sf, sizeof(double));
	d->ns_pk       = calloc(n_sf, sizeof(double));
	d->ns_sk       = calloc(n_sf, sizeof(double));
	d->ns_yk       = calloc(n_sf, sizeof(double));
	d->ns_vt       = calloc(n_sf, sizeof(double));
	d->ns_r        = calloc(n_sf, sizeof(double));
	d->ns_z        = calloc(n_sf, sizeof(double));

	d->ns_x        = calloc(n_x, sizeof(double));
	d->ns_x0       = calloc(n_x, sizeof(double));
	d->ns_w        = calloc((n_em > n_x) ? n_em : n_x, sizeof(double));
	d->ns_pt       = calloc(n_em, sizeof(double));
	d->ns_sft      = calloc(n_sf, sizeof(double));
	d->ns_gb       = calloc(n_em, sizeof(double));
	d->ns_bsave    = calloc(2*n_x, sizeof(double));
	d->ns_ub       = calloc(n_sf, sizeof(double));
	d->ns_S        = calloc(n_samp*n_sf, sizeof(double));
	d->ns_P        = calloc(n_samp*n_em, sizeof(double));

	d->ns_svdA     = calloc(mr*mc, sizeof(double));
	d->ns_svdS     = calloc(mc, sizeof(double));
	d->ns_svdU     = calloc(mr*mr, sizeof(double));
	d->ns_svdVt    = calloc(mc*mc, sizeof(double));

	d->ns_sf_state = calloc(n_sf, sizeof(int));
	d->ns_C        = calloc(n_sf, sizeof(int));
	d->ns_em_act   = calloc(n_em, sizeof(int));
	d->ns_n_em_act = -1;

	double 	wq   = 0.0;
	char 	ja   = 'A';
	int 	lda  = mr;
	int 	ldvt = mc;
	int 	lw   = -1;
	int 	info = 0;

	#if __APPLE__
		dgesvd(&ja, &ja, &mr, &mc, d->ns_svdA, &lda, d->ns_svdS, d->ns_svdU, &lda, d->ns_svdVt, &ldvt, &wq, &lw, &info);
	#else
		info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, ja, ja, mr, mc, d->ns_svdA, lda, d->ns_svdS, d->ns_svdU, lda, d->ns_svdVt, ldvt, &wq, lw);
	#endif

	d->ns_lwork = (int)wq + 1;
	d->ns_work  = calloc(d->ns_lwork, sizeof(double));
}

void NS_free(				SS_ref 				*d				){
	if (d->ns_V == NULL){ return; }

	int 	n_sf = d->ns_nc;
	int 	n_em = d->n_em;

	ns_free_mat(d->ns_V,   n_em);
	ns_free_mat(d->ns_M,   n_sf + 1);
	ns_free_mat(d->ns_N0,  n_sf);
	ns_free_mat(d->ns_N,   n_sf);
	ns_free_mat(d->ns_BkI, n_sf);

	free(d->ns_sf0);
	free(d->ns_sfc);
	free(d->ns_sf);
	free(d->ns_sf_prev);
	free(d->ns_trial);
	free(d->ns_g);
	free(d->ns_pg0);
	free(d->ns_pg1);
	free(d->ns_pk);
	free(d->ns_sk);
	free(d->ns_yk);
	free(d->ns_vt);
	free(d->ns_r);
	free(d->ns_z);
	free(d->ns_x);
	free(d->ns_x0);
	free(d->ns_w);
	free(d->ns_pt);
	free(d->ns_sft);
	free(d->ns_gb);
	free(d->ns_bsave);
	free(d->ns_ub);
	free(d->ns_S);
	free(d->ns_P);
	free(d->ns_svdA);
	free(d->ns_svdS);
	free(d->ns_svdU);
	free(d->ns_svdVt);
	free(d->ns_work);
	free(d->ns_sf_state);
	free(d->ns_C);
	free(d->ns_em_act);

	d->ns_V       = NULL;
	d->ns_ok      = 0;
	d->ns_ok_init = -1;
}

static int ns_check_gradient(	global_variable 	 gv,
								SS_ref 				*d				){
	int 	n_sf  = d->ns_nc;
	int 	n_em  = d->n_em;
	int 	n_dir = d->ns_n_dir0;
	double 	sfc[n_sf], sp[n_sf], sm[n_sf], g[n_sf];

	int 	n_samp = ns_n_samp(d);
	int 	qbest  = -1;
	double 	mbest  = 0.0;
	for (int q = 0; q < n_samp; q++){
		double m = 1.0;
		for (int i = 0; i < n_sf; i++){ if (d->ns_S[q*n_sf + i] < m){ m = d->ns_S[q*n_sf + i]; } }
		if (m > mbest){ mbest = m; qbest = q; }
	}
	if (qbest < 0){
		for (int i = 0; i < n_sf; i++){
			sfc[i] = 0.0;
			for (int e = 0; e < n_em; e++){ sfc[i] += d->ns_V[e][i]/(double)n_em; }
		}
	}
	else{
		for (int i = 0; i < n_sf; i++){ sfc[i] = d->ns_S[qbest*n_sf + i]; }
	}
	double G0 = ns_eval(gv, d, sfc, g);
	if (isnan(G0) || isinf(G0)){ return 0; }

	int 	n_test = (n_dir < 3) ? n_dir : 3;
	for (int k = 0; k < n_test; k++){
		double dmax = 0.0;
		for (int i = 0; i < n_sf; i++){ if (fabs(d->ns_N0[i][k]) > dmax){ dmax = fabs(d->ns_N0[i][k]); } }
		double h = 1e-6;
		for (int i = 0; i < n_sf; i++){
			if (fabs(d->ns_N0[i][k]) > 0.0 && sfc[i] - h*fabs(d->ns_N0[i][k]) < 0.0){ h = 0.5*sfc[i]/fabs(d->ns_N0[i][k]); }
		}
		if (!(h > 1e-12)){ continue; }
		for (int i = 0; i < n_sf; i++){ sp[i] = sfc[i] + h*d->ns_N0[i][k]; sm[i] = sfc[i] - h*d->ns_N0[i][k]; }
		double Gp = ns_eval(gv, d, sp, NULL);
		double Gm = ns_eval(gv, d, sm, NULL);
		if (isnan(Gp) || isnan(Gm)){ return 0; }
		double fd = (Gp - Gm)/(2.0*h);
		double an = 0.0;
		for (int i = 0; i < n_sf; i++){ an += g[i]*d->ns_N0[i][k]; }
		double sc = (fabs(fd) > 1.0) ? fabs(fd) : 1.0;
		if (fabs(fd - an) > 1e-4*sc){ return 0; }
	}
	return 1;
}

static int ns_pcoord_ok(	SS_ref 				*d				){
	int 	n_sf  = d->ns_nc;
	int 	n_x   = d->n_xeos;
	int 	n_em  = d->n_em;
	int 	n_in  = 0;
	int 	n_out = 0;
	unsigned seed = 2654435761u;
	double 	x[n_x];

	for (int q = 0; q < 600; q++){
		for (int j = 0; j < n_x; j++){
			seed = seed*1103515245u + 12345u;
			double u = ((seed >> 8) + 0.5)/16777216.0;
			x[j] = (q < 300) ? u : -log(u);
		}
		if (q >= 300){
			double se = 1.0;
			for (int j = 0; j < n_x; j++){ se += x[j]; }
			for (int j = 0; j < n_x; j++){ x[j] /= se; }
		}
		d->ns_obj(n_x, x, NULL, d);
		int sf_ok = 1, p_ok = 1;
		for (int i = 0; i < n_sf; i++){ if (!(d->sf[i] >= -1e-12)){ sf_ok = 0; } }
		for (int e = 0; e < n_em; e++){ if (!(d->p[e] >= -1e-12)){ p_ok = 0; } }
		if (sf_ok != p_ok){ return 0; }
		if (sf_ok){ n_in += 1; } else { n_out += 1; }
	}
	return (n_in >= 20 && n_out >= 20);
}

static int ns_init_phase_p(	global_variable 	 gv,
							SS_ref 				*d				){
	int 	nc     = d->ns_nc;
	int 	n_em   = d->n_em;
	int 	n_samp = ns_n_samp(d);
	int 	m      = n_em - 1;
	int 	ok     = 1;
	double 	primes[16] = {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0, 53.0};

	if (d->n_xeos != n_em || n_em < 2 || n_em > 16){ return 0; }

	for (int e = 0; e < n_em; e++){
		d->ns_gb[e]  = d->gb_lvl[e];
		d->gb_lvl[e] = d->gbase[e];
		for (int i = 0; i < nc; i++){ d->ns_V[e][i] = (i == e) ? 1.0 : 0.0; }
	}
	for (int i = 0; i < nc; i++){
		for (int r = 0; r < m; r++){ d->ns_svdA[r + i*m] = d->ns_V[r + 1][i] - d->ns_V[0][i]; }
	}
	if (ns_svd(d, m, nc) != 0 || ns_rank(d, m, nc) != m){ ok = 0; }
	else{
		for (int i = 0; i < nc; i++){
			for (int k = 0; k < m; k++){ d->ns_N0[i][k] = d->ns_svdVt[k + i*nc]; }
		}
		d->ns_n_dir0 = m;
	}
	for (int q = 0; q < n_samp && ok; q++){
		double se = 0.0;
		for (int i = 0; i < nc; i++){
			d->ns_S[q*nc + i] = 0.1 + 0.8*fmod((double)(q + 1)*sqrt(primes[i]), 1.0);
			se += d->ns_S[q*nc + i];
		}
		for (int i = 0; i < nc; i++){ d->ns_S[q*nc + i] /= se; }
	}
	if (ok){
		d->ns_fd = 0;
		d->ns_absent = 0;
		if (!ns_check_gradient(gv, d)){ d->ns_fd = 1; }
		for (int i = 0; i < nc; i++){ d->ns_ub[i] = 1.0; }
	}
	for (int e = 0; e < n_em; e++){ d->gb_lvl[e] = d->ns_gb[e]; }
	d->ns_n_em_act = -1;

	return ok;
}

static int ns_init_phase(	global_variable 	 gv,
							SS_ref 				*d				){
	if (d->ns_mode == 4){ return ns_init_phase_p(gv, d); }
	int 	n_sf   = d->ns_nc;
	int 	n_x    = d->n_xeos;
	int 	n_em   = d->n_em;
	int 	n_samp = ns_n_samp(d);
	int 	ok     = 1;
	double 	*S     = d->ns_S;
	double 	*P     = d->ns_P;
	double 	primes[16] = {2.0, 3.0, 5.0, 7.0, 11.0, 13.0, 17.0, 19.0, 23.0, 29.0, 31.0, 37.0, 41.0, 43.0, 47.0, 53.0};
	double 	x[n_x];

	if (n_x > 16 || n_em < 2 || n_sf < 1){ return 0; }

	for (int j = 0; j < n_x; j++){
		d->ns_bsave[2*j]     = d->bounds[j][0];
		d->ns_bsave[2*j + 1] = d->bounds[j][1];
		d->bounds[j][0]      = -1.0e3;
		d->bounds[j][1]      =  1.0e3;
	}
	for (int e = 0; e < n_em; e++){
		d->ns_gb[e]     = d->gb_lvl[e];
		d->gb_lvl[e]    = d->gbase[e];
	}

	for (int q = 0; q < n_samp; q++){
		for (int j = 0; j < n_x; j++){
			x[j] = 0.1 + 0.8*fmod((double)(q + 1)*sqrt(primes[j]), 1.0);
		}
		d->ns_obj(n_x, x, NULL, d);
		for (int i = 0; i < n_sf; i++){ S[q*n_sf + i] = d->sf[i]; }
		for (int e = 0; e < n_em; e++){ P[q*n_em + e] = d->p[e]; }
	}

	for (int e = 0; e < n_em; e++){
		for (int q = 0; q < n_samp; q++){ d->ns_svdA[q + e*n_samp] = P[q*n_em + e]; }
	}
	if (!ns_lstsq(d, n_samp, n_em, n_sf, S, d->ns_V)){ ok = 0; }

	for (int i = 0; i <= n_sf && ok; i++){
		for (int q = 0; q < n_samp; q++){ d->ns_svdA[q + i*n_samp] = (i < n_sf) ? S[q*n_sf + i] : 1.0; }
	}
	if (ok && !ns_lstsq(d, n_samp, n_sf + 1, n_em, P, d->ns_M)){ ok = 0; }

	double 	err = 0.0;
	for (int q = 0; q < n_samp && ok; q++){
		for (int i = 0; i < n_sf; i++){
			double a = 0.0;
			for (int e = 0; e < n_em; e++){ a += P[q*n_em + e]*d->ns_V[e][i]; }
			if (!(fabs(a - S[q*n_sf + i]) <= err)){ err = fabs(a - S[q*n_sf + i]); }
		}
		for (int e = 0; e < n_em; e++){
			double a = d->ns_M[n_sf][e];
			for (int i = 0; i < n_sf; i++){ a += S[q*n_sf + i]*d->ns_M[i][e]; }
			if (!(fabs(a - P[q*n_em + e]) <= err)){ err = fabs(a - P[q*n_em + e]); }
		}
	}
	if (!(err < NS_RES_TOL)){ ok = 0; }
	d->ns_mode = 1;

	if (!ok && n_x + 1 <= n_em){
		ok = 1;
		double *X = P;
		for (int q = 0; q < n_samp; q++){
			for (int j = 0; j < n_x; j++){ X[q*n_x + j] = 0.1 + 0.8*fmod((double)(q + 1)*sqrt(primes[j]), 1.0); }
		}
		for (int j = 0; j <= n_x; j++){
			for (int q = 0; q < n_samp; q++){ d->ns_svdA[q + j*n_samp] = (j < n_x) ? X[q*n_x + j] : 1.0; }
		}
		if (!ns_lstsq(d, n_samp, n_x + 1, n_sf, S, d->ns_V)){ ok = 0; }
		for (int i = 0; i <= n_sf && ok; i++){
			for (int q = 0; q < n_samp; q++){ d->ns_svdA[q + i*n_samp] = (i < n_sf) ? S[q*n_sf + i] : 1.0; }
		}
		if (ok && !ns_lstsq(d, n_samp, n_sf + 1, n_x, X, d->ns_M)){ ok = 0; }
		err = 0.0;
		for (int q = 0; q < n_samp && ok; q++){
			for (int i = 0; i < n_sf; i++){
				double a = d->ns_V[n_x][i];
				for (int j = 0; j < n_x; j++){ a += X[q*n_x + j]*d->ns_V[j][i]; }
				if (!(fabs(a - S[q*n_sf + i]) <= err)){ err = fabs(a - S[q*n_sf + i]); }
			}
			for (int j = 0; j < n_x; j++){
				double a = d->ns_M[n_sf][j];
				for (int i = 0; i < n_sf; i++){ a += S[q*n_sf + i]*d->ns_M[i][j]; }
				if (!(fabs(a - X[q*n_x + j]) <= err)){ err = fabs(a - X[q*n_x + j]); }
			}
		}
		if (!(err < NS_RES_TOL)){ ok = 0; }
		if (ok){
			for (int i = 0; i < n_sf; i++){
				for (int j = 0; j < n_x; j++){ d->ns_svdA[j + i*n_x] = d->ns_V[j][i]; }
			}
			if (ns_svd(d, n_x, n_sf) != 0 || ns_rank(d, n_x, n_sf) != n_x){ ok = 0; }
			else{
				for (int i = 0; i < n_sf; i++){
					for (int k = 0; k < n_x; k++){ d->ns_N0[i][k] = d->ns_svdVt[k + i*n_sf]; }
				}
				d->ns_n_dir0 = n_x;
				d->ns_mode   = 2;
			}
		}
	}
	if (!ok && n_em <= n_sf && n_em - 1 == n_x && ns_pcoord_ok(d)){
		ok = 1;
		for (int e = 0; e < n_em; e++){
			for (int i = 0; i < n_sf; i++){ d->ns_V[e][i] = (i == e) ? 1.0 : 0.0; }
		}
		for (int i = 0; i <= n_sf; i++){
			for (int e = 0; e < n_em; e++){ d->ns_M[i][e] = (i == e) ? 1.0 : 0.0; }
		}
		for (int q = 0; q < n_samp; q++){
			for (int j = 0; j < n_x; j++){ x[j] = 0.1 + 0.8*fmod((double)(q + 1)*sqrt(primes[j]), 1.0); }
			d->ns_obj(n_x, x, NULL, d);
			for (int i = 0; i < n_sf; i++){ S[q*n_sf + i] = (i < n_em) ? d->p[i] : 0.0; }
		}
		int m = n_em - 1;
		for (int i = 0; i < n_sf; i++){
			for (int r = 0; r < m; r++){ d->ns_svdA[r + i*m] = d->ns_V[r + 1][i] - d->ns_V[0][i]; }
		}
		if (ns_svd(d, m, n_sf) != 0 || ns_rank(d, m, n_sf) != n_x){ ok = 0; }
		else{
			for (int i = 0; i < n_sf; i++){
				for (int k = 0; k < n_x; k++){ d->ns_N0[i][k] = d->ns_svdVt[k + i*n_sf]; }
			}
			d->ns_n_dir0 = n_x;
			d->ns_mode   = 3;
		}
	}
	else if (ok && d->ns_mode == 1){
		int m = n_em - 1;
		for (int i = 0; i < n_sf; i++){
			for (int r = 0; r < m; r++){ d->ns_svdA[r + i*m] = d->ns_V[r + 1][i] - d->ns_V[0][i]; }
		}
		if (ns_svd(d, m, n_sf) != 0 || ns_rank(d, m, n_sf) != n_x){ ok = 0; }
		else{
			for (int i = 0; i < n_sf; i++){
				for (int k = 0; k < n_x; k++){ d->ns_N0[i][k] = d->ns_svdVt[k + i*n_sf]; }
			}
			d->ns_n_dir0 = n_x;
		}
	}

	for (int q = 0; q < n_samp && ok; q++){
		ns_x_of_sf(gv, d, &S[q*n_sf]);
		d->ns_obj(n_x, d->ns_x, NULL, d);
		if (!(ns_p_residual(d) <= NS_RES_TOL)){ ok = 0; }
	}

	for (int i = 0; i < n_sf && ok; i++){
		d->ns_ub[i] = 1.0;
		for (int e = 0; e < n_em && d->ns_mode == 1; e++){
			if (d->ns_V[e][i] > d->ns_ub[i]){ d->ns_ub[i] = d->ns_V[e][i]; }
		}
	}

	if (ok){
		for (int e = 0; e < n_em; e++){
			P[e]           = d->d_em[e];
			P[n_em + e]    = d->ox_penalty[e];
			d->d_em[e]     = 0.0;
			d->ox_penalty[e] = 0.0;
		}
		d->ns_fd = 0;
		d->ns_absent = 0;
		if (!ns_check_gradient(gv, d)){ d->ns_fd = 1; }
		for (int e = 0; e < n_em; e++){
			d->d_em[e]       = P[e];
			d->ox_penalty[e] = P[n_em + e];
		}
	}

	for (int j = 0; j < n_x; j++){
		d->bounds[j][0] = d->ns_bsave[2*j];
		d->bounds[j][1] = d->ns_bsave[2*j + 1];
	}
	for (int e = 0; e < n_em; e++){ d->gb_lvl[e] = d->ns_gb[e]; }
	d->ns_n_em_act = -1;

	return ok;
}

void TC_NS_init(			global_variable 	 gv,
							SS_ref 				*SS_ref_db		){
	obj_type 	SS_objective[gv.len_ss];
	P2X_type 	P2X_read[gv.len_ss];

	for (int iss = 0; iss < gv.len_ss; iss++){
		SS_objective[iss] = NULL;
		P2X_read[iss]     = NULL;
	}
	TC_SS_objective_init_function(	SS_objective,
									gv							);
	TC_P2X_init(					P2X_read,
									gv							);

	for (int iss = 0; iss < gv.len_ss; iss++){
		SS_ref *d     = &SS_ref_db[iss];
		d->ns_obj     = NULL;
		d->ns_p2x     = NULL;
		d->ns_ok      = 0;
		d->ns_ok_init = -1;
		d->ns_mode    = 0;
		d->ns_fd      = 0;
		d->ns_absent  = 0;

		if (strcmp(gv.SS_list[iss], "DEW") == 0 || d->n_em < 2 || d->n_sf < 1 || d->n_xeos > 16){ continue; }
		if (SS_objective[iss] == NULL || P2X_read[iss] == NULL){ continue; }

		d->ns_obj     = SS_objective[iss];
		d->ns_p2x     = P2X_read[iss];
		d->ns_nc      = d->n_sf;
		ns_alloc(gv, d);
		d->ns_ok_init = 0;
	}
}

void SB_NS_init(			global_variable 	 gv,
							SS_ref 				*SS_ref_db		){
	obj_type 	SS_objective[gv.len_ss];

	for (int iss = 0; iss < gv.len_ss; iss++){ SS_objective[iss] = NULL; }
	SB_SS_objective_init_function(	SS_objective,
									gv							);

	for (int iss = 0; iss < gv.len_ss; iss++){
		SS_ref *d     = &SS_ref_db[iss];
		d->ns_obj     = NULL;
		d->ns_p2x     = NULL;
		d->ns_ok      = 0;
		d->ns_ok_init = -1;
		d->ns_mode    = 0;
		d->ns_fd      = 0;
		d->ns_absent  = 0;

		if (d->n_em < 2 || d->n_xeos != d->n_em || d->n_em > 16 || SS_objective[iss] == NULL){ continue; }

		d->ns_obj     = SS_objective[iss];
		d->ns_mode    = 4;
		d->ns_nc      = d->n_em;
		ns_alloc(gv, d);
		d->ns_ok_init = 0;
	}
}
