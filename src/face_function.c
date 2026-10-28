/*@ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 **
 **   Project      : MAGEMin
 **   License      : GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 **   Developers   : Nicolas Riel, Boris Kaus
 **   Organization : Institute of Geosciences, Johannes-Gutenberg University, Mainz
 **   Contact      : nriel[at]uni-mainz.de, kaus[at]uni-mainz.de
 **
 ** ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ @*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if __APPLE__
	extern void dgesvd( char* jobu, char* jobvt, int* m, int* n, double* a, int* lda, double* s, double* u, int* ldu, double* vt, int* ldvt, double* work, int* lwork, int* info );
#else
	#include <lapacke.h>
#endif

#include "MAGEMin.h"
#include "face_function.h"

#define FACE_PIV_TOL 	1e-11
#define FACE_RC_TOL 	1e-12
#define FACE_SVD_TOL 	1e-9
#define FACE_ON_TOL 	1e-9

typedef struct {
	int 	 nr, nc, ncol;
	double 	*T;
	int 	*basis;
} face_tab;

static double *tab(face_tab *t, int r, int c){ return &t->T[r*t->ncol + c]; }

static void face_pivot(face_tab *t, int pr, int pc){
	double p = *tab(t, pr, pc);
	for (int c = 0; c < t->ncol; c++){ *tab(t, pr, c) /= p; }
	for (int r = 0; r <= t->nr; r++){
		if (r == pr){ continue; }
		double f = *tab(t, r, pc);
		if (f == 0.0){ continue; }
		for (int c = 0; c < t->ncol; c++){ *tab(t, r, c) -= f*(*tab(t, pr, c)); }
	}
	t->basis[pr] = pc;
}

static int face_simplex(face_tab *t, int n_enter){
	int rhs = t->ncol - 1;
	for (int it = 0; it < 50*(t->nr + n_enter) + 100; it++){
		int pc = -1;
		for (int c = 0; c < n_enter; c++){
			if (*tab(t, t->nr, c) < -FACE_RC_TOL){ pc = c; break; }
		}
		if (pc < 0){ return 0; }
		int 	pr   = -1;
		double 	best = 0.0;
		for (int r = 0; r < t->nr; r++){
			double a = *tab(t, r, pc);
			if (a > FACE_PIV_TOL){
				double ratio = *tab(t, r, rhs)/a;
				if (pr < 0 || ratio < best - 1e-14 || (fabs(ratio - best) <= 1e-14 && t->basis[r] < t->basis[pr])){ pr = r; best = ratio; }
			}
		}
		if (pr < 0){ return 2; }
		face_pivot(t, pr, pc);
	}
	return 3;
}

/*
	max cx s.t. Ax = b, x >= 0 (A row-major nr x nc). Returns 0 optimal, 1 infeasible, 2 unbounded, 3 iteration limit.
	With c == NULL only feasibility is tested (phase I). x (length nc) may be NULL.
*/
static int face_lp(int nr, int nc, const double *A, const double *b, const double *c, double *x, double feas_tol){
	face_tab t;
	t.nr 	= nr;
	t.nc 	= nc;
	t.ncol 	= nc + nr + 1;
	t.T 	= calloc((size_t)(nr + 1)*t.ncol, sizeof(double));
	t.basis = malloc(nr*sizeof(int));
	int rhs = t.ncol - 1;

	for (int r = 0; r < nr; r++){
		double s = (b[r] < 0.0) ? -1.0 : 1.0;
		for (int j = 0; j < nc; j++){ *tab(&t, r, j) = s*A[r*nc + j]; }
		*tab(&t, r, nc + r) = 1.0;
		*tab(&t, r, rhs) 	= s*b[r];
		t.basis[r] 			= nc + r;
	}
	for (int j = 0; j < t.ncol; j++){
		double a = 0.0;
		for (int r = 0; r < nr; r++){ a += *tab(&t, r, j); }
		*tab(&t, nr, j) = (j >= nc && j < nc + nr) ? 0.0 : -a;
	}
	int st = face_simplex(&t, nc + nr);
	if (st == 3 || -(*tab(&t, nr, rhs)) > feas_tol){
		free(t.T); free(t.basis);
		return (st == 3) ? 3 : 1;
	}
	for (int r = 0; r < nr; r++){
		if (t.basis[r] < nc){ continue; }
		for (int j = 0; j < nc; j++){
			if (fabs(*tab(&t, r, j)) > FACE_PIV_TOL){ face_pivot(&t, r, j); break; }
		}
	}
	if (c != NULL){
		for (int j = 0; j < t.ncol; j++){ *tab(&t, nr, j) = (j < nc) ? -c[j] : 0.0; }
		for (int r = 0; r < nr; r++){
			int bj = t.basis[r];
			if (bj < nc && c[bj] != 0.0){
				double f = c[bj];
				for (int j = 0; j < t.ncol; j++){ *tab(&t, nr, j) += f*(*tab(&t, r, j)); }
			}
		}
		st = face_simplex(&t, nc);
	}
	if (x != NULL){
		for (int j = 0; j < nc; j++){ x[j] = 0.0; }
		for (int r = 0; r < nr; r++){ if (t.basis[r] < nc){ x[t.basis[r]] = *tab(&t, r, rhs); } }
	}
	free(t.T); free(t.basis);
	return st;
}

static int face_svd_basis(int m, int n, const double *Acm, double *U, double *S){
	if (n == 0){ return 0; }
	char 	jobu = 'A', jobvt = 'N';
	int 	lda = m, ldu = m, ldvt = 1, info = 0;
	int 	mn = (m < n) ? m : n;
	double 	*a = malloc((size_t)m*n*sizeof(double));
	memcpy(a, Acm, (size_t)m*n*sizeof(double));
	double 	vt[1];
	int 	lwork = 5*(m + n) + 64;
	double 	*work = malloc(lwork*sizeof(double));
	#if __APPLE__
		dgesvd(&jobu, &jobvt, &m, &n, a, &lda, S, U, &ldu, vt, &ldvt, work, &lwork, &info);
	#else
		info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, jobu, jobvt, m, n, a, lda, S, U, ldu, vt, ldvt, work, lwork);
	#endif
	free(a); free(work);
	if (info != 0){ return -1; }
	int 	r = 0;
	double 	tol = FACE_SVD_TOL*((S[0] > 1.0) ? S[0] : 1.0);
	for (int k = 0; k < mn; k++){ if (S[k] > tol){ r += 1; } }
	return r;
}

double face_deviation(					global_variable 	 gv,
										const double 		*comp			){
	double dev = 0.0;
	for (int k = 0; k < gv.face_n; k++){
		double a = 0.0;
		for (int j = 0; j < gv.len_ox; j++){ a += gv.face_l[k*gv.len_ox + j]*comp[j]; }
		dev += fabs(a);
	}
	return dev;
}

void face_project_gamma(				global_variable 	 gv,
										double 				*gam			){
	for (int k = 0; k < gv.face_n; k++){
		double a = 0.0;
		for (int j = 0; j < gv.len_ox; j++){ a += gv.face_l[k*gv.len_ox + j]*gam[j]; }
		for (int j = 0; j < gv.len_ox; j++){ gam[j] -= a*gv.face_l[k*gv.len_ox + j]; }
	}
}

global_variable face_detect(			bulk_info 	 		 z_b,
										global_variable 	 gv,
										PP_ref  			*PP_ref_db,
										SS_ref  			*SS_ref_db		){
	gv.face_n = 0;
	if (gv.face_l == NULL){ return gv; }
	if (!(gv.solver == 0 && gv.ss_solver > 0)){ return gv; }

	int 	nox = gv.len_ox;
	int 	act[nox], na = 0;
	double 	bsum = 0.0;
	for (int j = 0; j < nox; j++){ if (z_b.bulk_rock[j] > 0.0){ act[na] = j; na += 1; bsum += z_b.bulk_rock[j]; } }
	if (na < 2 || !(bsum > 0.0)){ return gv; }

	int 	m_max = gv.len_pp;
	for (int i = 0; i < gv.len_ss; i++){ m_max += SS_ref_db[i].n_em; }
	int 	*g_ph = malloc(m_max*sizeof(int));
	int 	*g_em = malloc(m_max*sizeof(int));
	double 	*C 	  = malloc((size_t)na*m_max*sizeof(double));
	int 	 m 	  = 0;

	for (int i = 0; i < gv.len_pp; i++){
		if (gv.pp_flags[i][0] != 1){ continue; }
		if (gv.n_mu_fix > 0 && i >= gv.len_pp - gv.n_mu_fix){ continue; }
		int ok = 1;
		for (int j = 0; j < nox && ok; j++){ if (z_b.bulk_rock[j] == 0.0 && fabs(PP_ref_db[i].Comp[j]) > 1e-12){ ok = 0; } }
		if (!ok){ continue; }
		for (int a = 0; a < na; a++){ C[a*m_max + m] = PP_ref_db[i].Comp[act[a]]; }
		g_ph[m] = i; g_em[m] = -1; m += 1;
	}
	for (int i = 0; i < gv.len_ss; i++){
		if (SS_ref_db[i].ss_flags[0] != 1){ continue; }
		for (int k = 0; k < SS_ref_db[i].n_em; k++){
			if (SS_ref_db[i].z_em[k] != 1.0){ continue; }
			int ok = 1;
			for (int j = 0; j < nox && ok; j++){ if (z_b.bulk_rock[j] == 0.0 && fabs(SS_ref_db[i].Comp[k][j]) > 1e-12){ ok = 0; } }
			if (!ok){ continue; }
			for (int a = 0; a < na; a++){ C[a*m_max + m] = SS_ref_db[i].Comp[k][act[a]]; }
			g_ph[m] = gv.len_pp + i; g_em[m] = k; m += 1;
		}
	}
	if (m == 0){ free(g_ph); free(g_em); free(C); return gv; }

	double 	*A = malloc((size_t)na*m*sizeof(double));
	double 	 b[na], bs[na];
	double 	*bud = malloc(m*sizeof(double));
	for (int a = 0; a < na; a++){
		b[a] = z_b.bulk_rock[act[a]]/bsum;
		for (int j = 0; j < m; j++){ A[a*m + j] = C[a*m_max + j]; }
	}
	for (int j = 0; j < m; j++){
		double u = INFINITY;
		for (int a = 0; a < na; a++){ if (A[a*m + j] > 1e-12){ double v = b[a]/A[a*m + j]; if (v < u){ u = v; } } }
		bud[j] = isfinite(u) ? u : 1.0;
	}

	/* fast path: a representation with every generator strictly positive exists -> interior */
	for (int a = 0; a < na; a++){
		double s = 0.0;
		for (int j = 0; j < m; j++){ s += 1e-4*bud[j]/(double)m*A[a*m + j]; }
		bs[a] = b[a] - s;
	}
	int st_fast = face_lp(na, m, A, bs, NULL, NULL, 1e-12);
	if (st_fast == 0){
		free(g_ph); free(g_em); free(C); free(A); free(bud);
		return gv;
	}

	/* exact path */
	double 	*x  = malloc(m*sizeof(double));
	double 	*cj = calloc(m, sizeof(double));
	int 	*on = calloc(m, sizeof(int));
	int st = face_lp(na, m, A, b, NULL, x, 1e-10);
	if (st == 1){
		gv.face_n = -1;
		if (gv.verbose == 1){ printf(" FACE: bulk outside the composition space of the selected phases\n"); }
		free(g_ph); free(g_em); free(C); free(A); free(bud); free(x); free(cj); free(on);
		return gv;
	}
	for (int j = 0; j < m; j++){ if (x[j] > FACE_ON_TOL*bud[j]){ on[j] = 1; } }
	for (int j = 0; j < m; j++){
		if (on[j]){ continue; }
		cj[j] = 1.0;
		st 	  = face_lp(na, m, A, b, cj, x, 1e-10);
		cj[j] = 0.0;
		if (st == 2){ on[j] = 1; continue; }
		if (st != 0){ continue; }
		for (int k = 0; k < m; k++){ if (x[k] > FACE_ON_TOL*bud[k]){ on[k] = 1; } }
	}

	int n_off = 0;
	for (int j = 0; j < m; j++){ if (!on[j]){ n_off += 1; } }
	if (n_off > 0){
		double 	*Acm = malloc((size_t)na*m*sizeof(double));
		double 	*Ua  = malloc((size_t)na*na*sizeof(double));
		double 	*Uo  = malloc((size_t)na*na*sizeof(double));
		double 	 S[na];
		for (int j = 0; j < m; j++){ for (int a = 0; a < na; a++){ Acm[a + j*na] = A[a*m + j]; } }
		int r_all = face_svd_basis(na, m, Acm, Ua, S);
		int n_on  = 0;
		for (int j = 0; j < m; j++){
			if (!on[j]){ continue; }
			for (int a = 0; a < na; a++){ Acm[a + n_on*na] = A[a*m + j]; }
			n_on += 1;
		}
		int r_on = (n_on > 0) ? face_svd_basis(na, n_on, Acm, Uo, S) : 0;
		int nf   = 0;
		if (r_all > 0 && r_on >= 0 && r_all > r_on){
			double v[na];
			for (int q = 0; q < r_all && nf < r_all - r_on; q++){
				for (int a = 0; a < na; a++){ v[a] = Ua[a + q*na]; }
				for (int k = 0; k < r_on; k++){
					double d = 0.0;
					for (int a = 0; a < na; a++){ d += Uo[a + k*na]*v[a]; }
					for (int a = 0; a < na; a++){ v[a] -= d*Uo[a + k*na]; }
				}
				for (int k = 0; k < nf; k++){
					double d = 0.0;
					for (int a = 0; a < na; a++){ d += gv.face_l[k*nox + act[a]]*v[a]; }
					for (int a = 0; a < na; a++){ v[a] -= d*gv.face_l[k*nox + act[a]]; }
				}
				double nv = 0.0;
				for (int a = 0; a < na; a++){ nv += v[a]*v[a]; }
				nv = sqrt(nv);
				if (nv < 1e-6){ continue; }
				for (int j = 0; j < nox; j++){ gv.face_l[nf*nox + j] = 0.0; }
				for (int a = 0; a < na; a++){ gv.face_l[nf*nox + act[a]] = v[a]/nv; }
				nf += 1;
			}
		}
		gv.face_n = nf;
		free(Acm); free(Ua); free(Uo);
	}

	if (gv.face_n > 0){
		for (int j = 0; j < m; j++){
			if (on[j] || g_em[j] >= 0){ continue; }
			for (int f = 0; f < gv.n_flags; f++){ gv.pp_flags[g_ph[j]][f] = 0; }
		}
		for (int i = 0; i < gv.len_ss; i++){
			int has = 0, any = 0;
			for (int j = 0; j < m; j++){
				if (g_ph[j] != gv.len_pp + i){ continue; }
				any = 1;
				if (on[j]){ has = 1; }
			}
			if (any && !has){
				for (int f = 0; f < gv.n_flags; f++){ SS_ref_db[i].ss_flags[f] = 0; }
			}
		}
		if (gv.verbose == 1){
			printf(" FACE: bulk on a face of codimension %d; %d of %d generators off the face\n", gv.face_n, n_off, m);
		}
	}

	free(g_ph); free(g_em); free(C); free(A); free(bud); free(x); free(cj); free(on);
	return gv;
}
