/**
	List of objective functions used for non-linear minimization and to generate pseudocompounds
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <complex.h> 
#include "nlopt.h"
#include "MAGEMin.h"


/**************************************************************************************/
/**************************************************************************************/
/**********************METABASITE DATABASE (Gree et al., 2016)*************************/
/**************************************************************************************/
/**************************************************************************************/

/**
    Endmember to xeos for L
*/
void p2x_mb_liq(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[7]  = d->p[8];
    d->iguess[2]   = d->p[1]/(d->p[1] + d->p[2]);
    d->iguess[3]   = (d->p[3] + d->iguess[7])/(d->iguess[7] + 1.0);
    d->iguess[4]  = (d->p[4] + d->iguess[7])/(d->iguess[7] + 1.0);
    d->iguess[0]    = d->p[0]/(d->iguess[7] + 1.0);
    d->iguess[6]   = d->p[5]/(d->p[5] + d->p[6]);
    d->iguess[1]  = (d->p[1] + d->p[2])/(d->iguess[7] + 1.0);
    d->iguess[5]  = (d->p[5] + d->p[6])/(d->iguess[7] + 1.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for hb
*/
void p2x_mb_hb(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[7]   = d->p[10];
    d->iguess[6]   = d->p[8];
    d->iguess[2]   = d->iguess[6] + d->p[3];
    d->iguess[3]   = d->p[2] + d->p[9];
    d->iguess[4]   = d->p[9]/(d->p[2] + d->p[9]);
    d->iguess[5]   = d->iguess[3] + d->p[0] + d->p[1] + d->iguess[7];
    d->iguess[1]   = -0.5*d->iguess[3] + d->iguess[5] - 1.0*d->iguess[6] - 1.0*d->p[0] - 1.0*d->iguess[7] + d->iguess[2];
    d->iguess[0]  = (5.0*d->iguess[5] + 5.0*d->p[4] - 2.0*d->p[5] + d->p[6] + 5.0*d->iguess[2] - 5.0)/(2.0*d->iguess[5] + 2.0*d->iguess[6] + 2.0*d->iguess[7] + 2.0*d->iguess[1] + 2.0*d->iguess[2] - 7.0);
    d->iguess[8]  = -0.4*d->iguess[5]*d->iguess[0] + 2.0*d->iguess[5] - 0.4*d->iguess[6]*d->iguess[0] + 2.0*d->p[4] - 0.4*d->p[5] + 1.2*d->p[6] - 0.4*d->iguess[7]*d->iguess[0] - 0.4*d->iguess[0]*d->iguess[1] - 0.4*d->iguess[0]*d->iguess[2] + 2.4*d->iguess[0] + 2.0*d->iguess[2] - 2.0;
    d->iguess[9]  = (-2.0*d->iguess[5]*d->iguess[0] + 5.0*d->iguess[5] + 5.0*d->p[4] + 3.0*d->p[6] - 2.0*d->iguess[0]*d->iguess[2] + 5.0*d->iguess[0] + 5.0*d->iguess[2] - 5.0)/(2.0*d->iguess[6] + 2.0*d->iguess[7] + 2.0*d->iguess[1] - 2.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for aug
*/
void p2x_mb_aug(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[6]  = d->p[5];
    d->iguess[1]    = d->p[6] + d->iguess[6];
    d->iguess[2]    = d->p[4];
    d->iguess[4]    = d->iguess[2] + d->p[3];
    d->iguess[3]    = d->p[0] + d->iguess[1];
    d->iguess[0]   = (2.0*d->iguess[4] + 2.0*d->p[1] + d->p[7] + 2.0*d->iguess[3] - 2.0)/(2.0*d->iguess[4] + d->iguess[1] + d->iguess[3] - 2.0);
    d->iguess[5]   = (4.0*d->iguess[4]*d->p[1] + 4.0*d->iguess[4]*d->p[2] + 2.0*d->iguess[4]*d->p[7] + 4.0*d->iguess[4]*d->iguess[1] + 4.0*d->iguess[4]*d->iguess[3] - 8.0*d->iguess[4] + 4.0*d->iguess[4]*d->iguess[4] + 4.0*d->p[1]*d->iguess[1] - 4.0*d->p[1] + 2.0*d->p[2]*d->iguess[1] + 2.0*d->p[2]*d->iguess[3] - 4.0*d->p[2] + 2.0*d->p[7]*d->iguess[1] - 2.0*d->p[7] + 4.0*d->iguess[1]*d->iguess[3] - 4.0*d->iguess[1] - 4.0*d->iguess[3] + 4.0)/(d->iguess[4]*d->iguess[1] + 3.0*d->iguess[4]*d->iguess[3] - 4.0*d->iguess[4] + 2.0*d->iguess[4]*d->iguess[4] + d->iguess[1]*d->iguess[3] - 1.0*d->iguess[1] - 3.0*d->iguess[3] + d->iguess[3]*d->iguess[3] + 2.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for dio
*/
void p2x_mb_dio(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[4]  = 0.5*d->p[6];
    d->iguess[3]  = 0.5*d->p[4];
    d->iguess[1]    = d->iguess[3] + d->p[0] + d->p[3] + 2.0*d->iguess[4];
    d->iguess[5]   = (d->iguess[3]*d->p[2] + 0.5*d->iguess[3]*d->p[5] + 0.5*d->iguess[1]*d->p[5] - 0.5*d->p[5])/(d->iguess[3]*d->iguess[1] - 1.0*d->iguess[3] - 2.0*d->iguess[1] + d->iguess[1]*d->iguess[1] + 1.0);
    d->iguess[0]   = (-1.0*d->iguess[3]*d->iguess[5] - 1.0*d->iguess[1]*d->iguess[5] + 0.5*d->p[5] + d->iguess[5])/d->iguess[3];
    d->iguess[2]    = (-1.0*d->iguess[3] + d->iguess[1] - 1.0*d->p[0] - 1.0*d->iguess[4])/d->iguess[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for opx
*/
void p2x_mb_opx(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[1]  = d->p[3];
    d->iguess[2]  = d->p[4];
    d->iguess[3]  = d->p[5];
    d->iguess[0]  = (d->iguess[3] + d->iguess[2] + d->p[0] - 1.0*d->p[1] + d->iguess[1] - 1.0)/(d->iguess[3] + d->iguess[2] + d->iguess[1] - 2.0);
    d->iguess[4]  = (-1.0*d->iguess[3]*d->iguess[0] + d->iguess[2]*d->iguess[0] - 1.0*d->p[2] + d->iguess[0]*d->iguess[1])/(d->iguess[3] - 1.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for g
*/
void p2x_mb_g(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]  = d->p[3];
    d->iguess[1]  = d->p[2];
    d->iguess[0]  = d->p[1]/(1.0 - 1.0*d->iguess[1]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for ol
*/
void p2x_mb_ol(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = d->p[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for pl4tr
*/
void p2x_mb_pl4tr(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[1]  = d->p[2];
    d->iguess[0]  = d->p[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for abc
*/
void p2x_mb_abc(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = d->p[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for k4tr
*/
void p2x_mb_k4tr(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = d->p[0];
    d->iguess[1]  = d->p[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for sp
*/
void p2x_mb_sp(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]   = d->p[3];
    d->iguess[1]  = -1.0*d->p[2] - 1.0*d->iguess[2] + 1.0;
    d->iguess[0]  = (-1.0*d->p[1] + d->iguess[2] + 1.0)/(d->iguess[2] + 1.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for ilm
*/
void p2x_mb_ilm(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[1]   = d->p[0];
    d->iguess[0]  = d->p[1] + d->iguess[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for ilmm
*/
void p2x_mb_ilmm(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]   = d->p[0];
    d->iguess[0]  = 1.0 - 1.0*d->p[2];
    d->iguess[1]   = d->p[3];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for ep
*/
void p2x_mb_ep(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[1]  = 0.5*d->p[1];
    d->iguess[0]  = d->p[2] + d->iguess[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for bi
*/
void p2x_mb_bi(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]  = d->p[5];
    d->iguess[3]  = d->p[4];
    d->iguess[1]  = d->p[3];
    d->iguess[0]  = (-3.0*d->p[1] - 1.0*d->p[2])/(d->iguess[2] + d->iguess[3] + d->iguess[1] - 3.0);
    d->iguess[4]  = 1.5*d->iguess[2]*d->iguess[0] - 1.5*d->iguess[2] + 1.5*d->p[0] + 1.5*d->iguess[3]*d->iguess[0] - 1.5*d->iguess[3] + 1.5*d->iguess[0]*d->iguess[1] - 1.5*d->iguess[0] - 1.5*d->iguess[1] + 1.5;
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for mu
*/
void p2x_mb_mu(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]  = d->p[5];
    d->iguess[4]  = d->p[4];
    d->iguess[3]  = d->p[3];
    d->iguess[1]  = d->iguess[4] + d->iguess[2] + d->iguess[3] + d->p[0];
    d->iguess[0]  = (d->p[1] + d->iguess[1] - 1.0)/(d->iguess[1] - 1.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for chl
*/
void p2x_mb_chl(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]  = d->p[6];
    d->iguess[3]  = 0.5*d->p[0] + 0.5*d->p[3];
    d->iguess[1]  = d->p[2] + d->iguess[3];
    d->iguess[5]  = (5.0*d->iguess[2]*d->p[0] + d->iguess[2]*d->p[1] - 3.0*d->iguess[2]*d->p[4] - 8.0*d->iguess[2]*d->iguess[3] + 2.0*d->iguess[2]*d->iguess[1] - 2.0*d->iguess[2] + d->iguess[2]*d->iguess[2] + 5.0*d->p[0]*d->iguess[3] + 5.0*d->p[0]*d->iguess[1] - 5.0*d->p[0] + d->p[1]*d->iguess[3] + d->p[1]*d->iguess[1] - 1.0*d->p[1] - 4.0*d->p[4]*d->iguess[3] - 2.0*d->p[4]*d->iguess[1] - 2.0*d->p[4] - 8.0*d->iguess[3]*d->iguess[1] + 8.0*d->iguess[3] - 9.0*d->iguess[3]*d->iguess[3] - 2.0*d->iguess[1] + d->iguess[1]*d->iguess[1] + 1.0)/(d->iguess[2]*d->iguess[3] + 3.0*d->iguess[2]*d->iguess[1] - 7.0*d->iguess[2] + d->iguess[2]*d->iguess[2] + 2.0*d->iguess[3]*d->iguess[1] - 6.0*d->iguess[3] - 8.0*d->iguess[1] + 2.0*d->iguess[1]*d->iguess[1] + 6.0);
    d->iguess[0]  = (d->iguess[2]*d->iguess[5] - 1.0*d->p[4] + d->iguess[5]*d->iguess[3] + d->iguess[5]*d->iguess[1] - 1.0*d->iguess[5])/(d->iguess[2] + d->iguess[3] + d->iguess[1] - 1.0);
    d->iguess[4]  = (d->iguess[2]*d->iguess[5] - 0.8*d->iguess[2]*d->iguess[0] - 0.8*d->p[5] + d->iguess[5]*d->iguess[3] + d->iguess[5]*d->iguess[1] - 1.0*d->iguess[5] - 1.6*d->iguess[0]*d->iguess[1] + 0.8*d->iguess[0])/(d->iguess[3] - 1.0*d->iguess[1] + 1.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Update dpdx matrix of L
*/
void dpdx_mb_liq(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[7] + 1.0;      dp_dx[0][1] = 0.0;      dp_dx[0][2] = 0.0;      dp_dx[0][3] = 0.0;      dp_dx[0][4] = 0.0;      dp_dx[0][5] = 0.0;      dp_dx[0][6] = 0.0;      dp_dx[0][7] = x[0];      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = x[2]*(x[7] + 1.0);      dp_dx[1][2] = x[1]*(x[7] + 1.0);      dp_dx[1][3] = 0.0;      dp_dx[1][4] = 0.0;      dp_dx[1][5] = 0.0;      dp_dx[1][6] = 0.0;      dp_dx[1][7] = x[1]*x[2];      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = (1.0 -x[2])*(x[7] + 1.0);      dp_dx[2][2] = -x[1]*(x[7] + 1.0);      dp_dx[2][3] = 0.0;      dp_dx[2][4] = 0.0;      dp_dx[2][5] = 0.0;      dp_dx[2][6] = 0.0;      dp_dx[2][7] = x[1]*(1.0 -x[2]);      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = x[7] + 1.0;      dp_dx[3][4] = 0.0;      dp_dx[3][5] = 0.0;      dp_dx[3][6] = 0.0;      dp_dx[3][7] = x[3] - 1.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = x[7] + 1.0;      dp_dx[4][5] = 0.0;      dp_dx[4][6] = 0.0;      dp_dx[4][7] = x[4] - 1.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      dp_dx[5][5] = x[6]*(x[7] + 1.0);      dp_dx[5][6] = x[5]*(x[7] + 1.0);      dp_dx[5][7] = x[5]*x[6];      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 0.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = (1.0 -x[6])*(x[7] + 1.0);      dp_dx[6][6] = -x[5]*(x[7] + 1.0);      dp_dx[6][7] = x[5]*(1.0 -x[6]);      
    dp_dx[7][0] = -x[7] - 1.0;      dp_dx[7][1] = -x[7] - 1.0;      dp_dx[7][2] = 0.0;      dp_dx[7][3] = -x[7] - 1.0;      dp_dx[7][4] = -x[7] - 1.0;      dp_dx[7][5] = -x[7] - 1.0;      dp_dx[7][6] = 0.0;      dp_dx[7][7] = -x[0] -x[1] -x[3] -x[4] -x[5] + 1.0;      
    dp_dx[8][0] = 0.0;      dp_dx[8][1] = 0.0;      dp_dx[8][2] = 0.0;      dp_dx[8][3] = 0.0;      dp_dx[8][4] = 0.0;      dp_dx[8][5] = 0.0;      dp_dx[8][6] = 0.0;      dp_dx[8][7] = 1.00;      
}


/**
    Update dpdx matrix of hb
*/
void dpdx_mb_hb(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = -1.00;      dp_dx[0][2] = 1.00;      dp_dx[0][3] = -0.500;      dp_dx[0][4] = 0.0;      dp_dx[0][5] = 1.00;      dp_dx[0][6] = -1.00;      dp_dx[0][7] = -1.00;      dp_dx[0][8] = 0.0;      dp_dx[0][9] = 0.0;      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = 1.00;      dp_dx[1][2] = -1.00;      dp_dx[1][3] = -0.500;      dp_dx[1][4] = 0.0;      dp_dx[1][5] = 0.0;      dp_dx[1][6] = 1.00;      dp_dx[1][7] = 0.0;      dp_dx[1][8] = 0.0;      dp_dx[1][9] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 0.0;      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 1.0 -x[4];      dp_dx[2][4] = -x[3];      dp_dx[2][5] = 0.0;      dp_dx[2][6] = 0.0;      dp_dx[2][7] = 0.0;      dp_dx[2][8] = 0.0;      dp_dx[2][9] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.00;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      dp_dx[3][5] = 0.0;      dp_dx[3][6] = -1.00;      dp_dx[3][7] = 0.0;      dp_dx[3][8] = 0.0;      dp_dx[3][9] = 0.0;      
    dp_dx[4][0] = x[2] + x[5] - 1.0;      dp_dx[4][1] = x[9];      dp_dx[4][2] = x[0] - 1.0;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = 0.0;      dp_dx[4][5] = x[0] - 1.0;      dp_dx[4][6] = x[9];      dp_dx[4][7] = x[9];      dp_dx[4][8] = -1.50;      dp_dx[4][9] = x[1] + x[6] + x[7] - 1.0;      
    dp_dx[5][0] = -x[1] + x[2] + x[5] -x[6] -x[7] + 1.0;      dp_dx[5][1] = -x[0] + 2.0*x[9];      dp_dx[5][2] = x[0];      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      dp_dx[5][5] = x[0];      dp_dx[5][6] = -x[0] + 2.0*x[9];      dp_dx[5][7] = -x[0] + 2.0*x[9];      dp_dx[5][8] = -2.50;      dp_dx[5][9] = 2.0*x[1] + 2.0*x[6] + 2.0*x[7] - 2.0;      
    dp_dx[6][0] = -x[2] -x[5];      dp_dx[6][1] = -x[9];      dp_dx[6][2] = -x[0];      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = -x[0];      dp_dx[6][6] = -x[9];      dp_dx[6][7] = -x[9];      dp_dx[6][8] = 2.50;      dp_dx[6][9] = -x[1] -x[6] -x[7] + 1.0;      
    dp_dx[7][0] = x[1] -x[2] -x[5] + x[6] + x[7];      dp_dx[7][1] = x[0] - 2.0*x[9];      dp_dx[7][2] = -x[0];      dp_dx[7][3] = 0.0;      dp_dx[7][4] = 0.0;      dp_dx[7][5] = -x[0];      dp_dx[7][6] = x[0] - 2.0*x[9];      dp_dx[7][7] = x[0] - 2.0*x[9];      dp_dx[7][8] = 1.50;      dp_dx[7][9] = -2.0*x[1] - 2.0*x[6] - 2.0*x[7] + 2.0;      
    dp_dx[8][0] = 0.0;      dp_dx[8][1] = 0.0;      dp_dx[8][2] = 0.0;      dp_dx[8][3] = 0.0;      dp_dx[8][4] = 0.0;      dp_dx[8][5] = 0.0;      dp_dx[8][6] = 1.00;      dp_dx[8][7] = 0.0;      dp_dx[8][8] = 0.0;      dp_dx[8][9] = 0.0;      
    dp_dx[9][0] = 0.0;      dp_dx[9][1] = 0.0;      dp_dx[9][2] = 0.0;      dp_dx[9][3] = x[4];      dp_dx[9][4] = x[3];      dp_dx[9][5] = 0.0;      dp_dx[9][6] = 0.0;      dp_dx[9][7] = 0.0;      dp_dx[9][8] = 0.0;      dp_dx[9][9] = 0.0;      
    dp_dx[10][0] = 0.0;      dp_dx[10][1] = 0.0;      dp_dx[10][2] = 0.0;      dp_dx[10][3] = 0.0;      dp_dx[10][4] = 0.0;      dp_dx[10][5] = 0.0;      dp_dx[10][6] = 0.0;      dp_dx[10][7] = 1.00;      dp_dx[10][8] = 0.0;      dp_dx[10][9] = 0.0;      
}


/**
    Update dpdx matrix of aug
*/
void dpdx_mb_aug(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = -1.00;      dp_dx[0][2] = 0.0;      dp_dx[0][3] = 1.00;      dp_dx[0][4] = 0.0;      dp_dx[0][5] = 0.0;      dp_dx[0][6] = 0.0;      
    dp_dx[1][0] = x[3] + x[4] - 1.0;      dp_dx[1][1] = 0.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = x[0] + 0.5*x[5] - 1.0;      dp_dx[1][4] = x[0] + 0.5*x[5] - 1.0;      dp_dx[1][5] = 0.5*x[3] + 0.5*x[4] - 0.5;      dp_dx[1][6] = 0.0;      
    dp_dx[2][0] = -x[1] -x[4] + 1.0;      dp_dx[2][1] = -x[0];      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 0.5*x[5];      dp_dx[2][4] = -x[0] + 0.5*x[5];      dp_dx[2][5] = 0.5*x[3] + 0.5*x[4] - 0.5;      dp_dx[2][6] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = -1.00;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 1.00;      dp_dx[3][5] = 0.0;      dp_dx[3][6] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 1.00;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = 0.0;      dp_dx[4][5] = 0.0;      dp_dx[4][6] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      dp_dx[5][5] = 0.0;      dp_dx[5][6] = 1.00;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 1.00;      dp_dx[6][2] = 0.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = 0.0;      dp_dx[6][6] = -1.00;      
    dp_dx[7][0] = x[1] -x[3];      dp_dx[7][1] = x[0];      dp_dx[7][2] = 0.0;      dp_dx[7][3] = -x[0] -x[5];      dp_dx[7][4] = -x[5];      dp_dx[7][5] = -x[3] -x[4] + 1.0;      dp_dx[7][6] = 0.0;      
}


/**
    Update dpdx matrix of dio
*/
void dpdx_mb_dio(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = 1.0 -x[2];      dp_dx[0][2] = -x[1];      dp_dx[0][3] = -1.00;      dp_dx[0][4] = -1.00;      dp_dx[0][5] = 0.0;      
    dp_dx[1][0] = x[1] -x[3] - 1.0;      dp_dx[1][1] = x[0] -x[5] - 1.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = -x[0] -x[5] - 1.0;      dp_dx[1][4] = 0.0;      dp_dx[1][5] = -x[1] -x[3] + 1.0;      
    dp_dx[2][0] = -x[1] -x[3] + 1.0;      dp_dx[2][1] = -x[0] -x[5];      dp_dx[2][2] = 0.0;      dp_dx[2][3] = -x[0] -x[5];      dp_dx[2][4] = 0.0;      dp_dx[2][5] = -x[1] -x[3] + 1.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = x[2];      dp_dx[3][2] = x[1];      dp_dx[3][3] = 0.0;      dp_dx[3][4] = -1.00;      dp_dx[3][5] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 2.00;      dp_dx[4][4] = 0.0;      dp_dx[4][5] = 0.0;      
    dp_dx[5][0] = 2.0*x[3];      dp_dx[5][1] = 2.0*x[5];      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 2.0*x[0] + 2.0*x[5];      dp_dx[5][4] = 0.0;      dp_dx[5][5] = 2.0*x[1] + 2.0*x[3] - 2.0;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 0.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 2.00;      dp_dx[6][5] = 0.0;      
}


/**
    Update dpdx matrix of opx
*/
void dpdx_mb_opx(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[3] - 1.0;      dp_dx[0][1] = -1.00;      dp_dx[0][2] = -1.00;      dp_dx[0][3] = x[0] + 0.5*x[4] - 1.0;      dp_dx[0][4] = 0.5*x[3] - 0.5;      
    dp_dx[1][0] = -x[1] -x[2] + 1.0;      dp_dx[1][1] = -x[0];      dp_dx[1][2] = -x[0];      dp_dx[1][3] = 0.5*x[4];      dp_dx[1][4] = 0.5*x[3] - 0.5;      
    dp_dx[2][0] = x[1] + x[2] -x[3];      dp_dx[2][1] = x[0];      dp_dx[2][2] = x[0];      dp_dx[2][3] = -x[0] -x[4];      dp_dx[2][4] = 1.0 -x[3];      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 1.00;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 1.00;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 1.00;      dp_dx[5][4] = 0.0;      
}


/**
    Update dpdx matrix of g
*/
void dpdx_mb_g(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[1] - 1.0;      dp_dx[0][1] = x[0] - 1.0;      dp_dx[0][2] = -1.00;      
    dp_dx[1][0] = 1.0 -x[1];      dp_dx[1][1] = -x[0];      dp_dx[1][2] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.00;      dp_dx[2][2] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.00;      
}


/**
    Update dpdx matrix of ol
*/
void dpdx_mb_ol(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.00;      
    dp_dx[1][0] = 1.00;      
}


/**
    Update dpdx matrix of pl4tr
*/
void dpdx_mb_pl4tr(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.00;      dp_dx[0][1] = -1.00;      
    dp_dx[1][0] = 1.00;      dp_dx[1][1] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.00;      
}


/**
    Update dpdx matrix of abc
*/
void dpdx_mb_abc(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.00;      
    dp_dx[1][0] = 1.00;      
}


/**
    Update dpdx matrix of k4tr
*/
void dpdx_mb_k4tr(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 1.00;      dp_dx[0][1] = 0.0;      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = 1.00;      
    dp_dx[2][0] = -1.00;      dp_dx[2][1] = -1.00;      
}


/**
    Update dpdx matrix of sp
*/
void dpdx_mb_sp(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[2] + 1.0;      dp_dx[0][1] = 1.00;      dp_dx[0][2] = x[0] - 1.0;      
    dp_dx[1][0] = -x[2] - 1.0;      dp_dx[1][1] = 0.0;      dp_dx[1][2] = 1.0 -x[0];      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = -1.00;      dp_dx[2][2] = -1.00;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.00;      
}


/**
    Update dpdx matrix of ilm
*/
void dpdx_mb_ilm(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = 1.00;      
    dp_dx[1][0] = 1.00;      dp_dx[1][1] = -1.00;      
    dp_dx[2][0] = -1.00;      dp_dx[2][1] = 0.0;      
}


/**
    Update dpdx matrix of ilmm
*/
void dpdx_mb_ilmm(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = 0.0;      dp_dx[0][2] = 1.00;      
    dp_dx[1][0] = 1.00;      dp_dx[1][1] = -1.00;      dp_dx[1][2] = -1.00;      
    dp_dx[2][0] = -1.00;      dp_dx[2][1] = 0.0;      dp_dx[2][2] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 1.00;      dp_dx[3][2] = 0.0;      
}


/**
    Update dpdx matrix of ep
*/
void dpdx_mb_ep(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.00;      dp_dx[0][1] = -1.00;      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = 2.00;      
    dp_dx[2][0] = 1.00;      dp_dx[2][1] = -1.00;      
}


/**
    Update dpdx matrix of bi
*/
void dpdx_mb_bi(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[1] + x[2] + x[3] - 1.0;      dp_dx[0][1] = x[0] - 1.0;      dp_dx[0][2] = x[0] - 1.0;      dp_dx[0][3] = x[0] - 1.0;      dp_dx[0][4] = -2.0/3.0;      
    dp_dx[1][0] = 1.00;      dp_dx[1][1] = 0.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      dp_dx[1][4] = -1.0/3.0;      
    dp_dx[2][0] = -x[1] -x[2] -x[3];      dp_dx[2][1] = -x[0];      dp_dx[2][2] = -x[0];      dp_dx[2][3] = -x[0];      dp_dx[2][4] = 1.00;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 1.00;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 1.00;      dp_dx[4][4] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 1.00;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      
}


/**
    Update dpdx matrix of mu
*/
void dpdx_mb_mu(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = 1.00;      dp_dx[0][2] = -1.00;      dp_dx[0][3] = -1.00;      dp_dx[0][4] = -1.00;      
    dp_dx[1][0] = x[1] - 1.0;      dp_dx[1][1] = x[0] - 1.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      dp_dx[1][4] = 0.0;      
    dp_dx[2][0] = 1.0 -x[1];      dp_dx[2][1] = -x[0];      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 0.0;      dp_dx[2][4] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 1.00;      dp_dx[3][4] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = 1.00;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 1.00;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      
}


/**
    Update dpdx matrix of chl
*/
void dpdx_mb_chl(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -x[1] -x[2] -x[3];      dp_dx[0][1] = -x[0] + 0.25*x[4] + 1.25*x[5];      dp_dx[0][2] = -x[0] + 1.25*x[5];      dp_dx[0][3] = -x[0] - 0.25*x[4] + 1.25*x[5] + 2.0;      dp_dx[0][4] = 0.25*x[1] - 0.25*x[3] - 0.25;      dp_dx[0][5] = 1.25*x[1] + 1.25*x[2] + 1.25*x[3] - 1.25;      
    dp_dx[1][0] = 3.0*x[1] + 2.0*x[2] + x[3] - 2.0;      dp_dx[1][1] = 3.0*x[0] - 1.25*x[4] - 2.25*x[5] - 1.0;      dp_dx[1][2] = 2.0*x[0] - 2.25*x[5] - 1.0;      dp_dx[1][3] = x[0] + 1.25*x[4] - 2.25*x[5] - 1.0;      dp_dx[1][4] = -1.25*x[1] + 1.25*x[3] + 1.25;      dp_dx[1][5] = -2.25*x[1] - 2.25*x[2] - 2.25*x[3] + 2.25;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.00;      dp_dx[2][2] = 0.0;      dp_dx[2][3] = -1.00;      dp_dx[2][4] = 0.0;      dp_dx[2][5] = 0.0;      
    dp_dx[3][0] = x[1] + x[2] + x[3];      dp_dx[3][1] = x[0] - 0.25*x[4] - 1.25*x[5];      dp_dx[3][2] = x[0] - 1.25*x[5];      dp_dx[3][3] = x[0] + 0.25*x[4] - 1.25*x[5];      dp_dx[3][4] = -0.25*x[1] + 0.25*x[3] + 0.25;      dp_dx[3][5] = -1.25*x[1] - 1.25*x[2] - 1.25*x[3] + 1.25;      
    dp_dx[4][0] = -x[1] -x[2] -x[3] + 1.0;      dp_dx[4][1] = -x[0] + x[5];      dp_dx[4][2] = -x[0] + x[5];      dp_dx[4][3] = -x[0] + x[5];      dp_dx[4][4] = 0.0;      dp_dx[4][5] = x[1] + x[2] + x[3] - 1.0;      
    dp_dx[5][0] = -2.0*x[1] -x[2] + 1.0;      dp_dx[5][1] = -2.0*x[0] + 1.25*x[4] + 1.25*x[5];      dp_dx[5][2] = -x[0] + 1.25*x[5];      dp_dx[5][3] = -1.25*x[4] + 1.25*x[5];      dp_dx[5][4] = 1.25*x[1] - 1.25*x[3] - 1.25;      dp_dx[5][5] = 1.25*x[1] + 1.25*x[2] + 1.25*x[3] - 1.25;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 1.00;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = 0.0;      
}

    
/**
    Endmember fraction of L
*/
void px_mb_liq(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[0]*(x[7] + 1.0);
        p[1]           = x[1]*x[2]*(x[7] + 1.0);
        p[2]           = x[1]*(1.0 -x[2])*(x[7] + 1.0);
        p[3]           = x[3]*(x[7] + 1.0) -x[7];
        p[4]           = x[4]*(x[7] + 1.0) -x[7];
        p[5]           = x[5]*x[6]*(x[7] + 1.0);
        p[6]           = x[5]*(1.0 -x[6])*(x[7] + 1.0);
        p[7]           = x[7] + (x[7] + 1.0)*(-1.0*x[0] -x[1] -x[3] -x[4] -x[5]) + 1.0;
        p[8]           = x[7];
}

    
/**
    Endmember fraction of hb
*/
void px_mb_hb(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[1] + x[2] - 0.5*x[3] + x[5] -x[6] -x[7];
        p[1]           = x[1] -x[2] - 0.5*x[3] + x[6];
        p[2]           = -x[3]*x[4] + x[3];
        p[3]           = x[2] -x[6];
        p[4]           = x[0]*x[2] + x[0]*x[5] -x[0] + x[1]*x[9] -x[2] -x[5] + x[6]*x[9] + x[7]*x[9] - 1.5*x[8] -x[9] + 1.0;
        p[5]           = -x[0]*x[1] + x[0]*x[2] + x[0]*x[5] -x[0]*x[6] -x[0]*x[7] + x[0] + 2.0*x[1]*x[9] + 2.0*x[6]*x[9] + 2.0*x[7]*x[9] - 2.5*x[8] - 2.0*x[9];
        p[6]           = -x[0]*x[2] -x[0]*x[5] -x[1]*x[9] -x[6]*x[9] -x[7]*x[9] + 2.5*x[8] + x[9];
        p[7]           = x[0]*x[1] -x[0]*x[2] -x[0]*x[5] + x[0]*x[6] + x[0]*x[7] - 2.0*x[1]*x[9] - 2.0*x[6]*x[9] - 2.0*x[7]*x[9] + 1.5*x[8] + 2.0*x[9];
        p[8]           = x[6];
        p[9]           = x[3]*x[4];
        p[10]           = x[7];
}

    
/**
    Endmember fraction of aug
*/
void px_mb_aug(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[1] + x[3];
        p[1]           = x[0]*x[3] + x[0]*x[4] -x[0] + 0.5*x[3]*x[5] -x[3] + 0.5*x[4]*x[5] -x[4] - 0.5*x[5] + 1.0;
        p[2]           = -x[0]*x[1] -x[0]*x[4] + x[0] + 0.5*x[3]*x[5] + 0.5*x[4]*x[5] - 0.5*x[5];
        p[3]           = -x[2] + x[4];
        p[4]           = x[2];
        p[5]           = x[6];
        p[6]           = x[1] -x[6];
        p[7]           = x[0]*x[1] -x[0]*x[3] -x[3]*x[5] -x[4]*x[5] + x[5];
}

    
/**
    Endmember fraction of dio
*/
void px_mb_dio(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[1]*x[2] + x[1] -x[3] -x[4];
        p[1]           = x[0]*x[1] -x[0]*x[3] -x[0] -x[1]*x[5] -x[1] -x[3]*x[5] -x[3] + x[5] + 1.0;
        p[2]           = -x[0]*x[1] -x[0]*x[3] + x[0] -x[1]*x[5] -x[3]*x[5] + x[5];
        p[3]           = x[1]*x[2] -x[4];
        p[4]           = 2.0*x[3];
        p[5]           = 2.0*x[0]*x[3] + 2.0*x[1]*x[5] + 2.0*x[3]*x[5] - 2.0*x[5];
        p[6]           = 2.0*x[4];
}

    
/**
    Endmember fraction of opx
*/
void px_mb_opx(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[0]*x[3] -x[0] -x[1] -x[2] + 0.5*x[3]*x[4] -x[3] - 0.5*x[4] + 1.0;
        p[1]           = -x[0]*x[1] -x[0]*x[2] + x[0] + 0.5*x[3]*x[4] - 0.5*x[4];
        p[2]           = x[0]*x[1] + x[0]*x[2] -x[0]*x[3] -x[3]*x[4] + x[4];
        p[3]           = x[1];
        p[4]           = x[2];
        p[5]           = x[3];
}

    
/**
    Endmember fraction of g
*/
void px_mb_g(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[0]*x[1] -x[0] -x[1] -x[2] + 1.0;
        p[1]           = -x[0]*x[1] + x[0];
        p[2]           = x[1];
        p[3]           = x[2];
}

    
/**
    Endmember fraction of ol
*/
void px_mb_ol(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 1.0 -x[0];
        p[1]           = x[0];
}

    
/**
    Endmember fraction of pl4tr
*/
void px_mb_pl4tr(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[0] -x[1] + 1.0;
        p[1]           = x[0];
        p[2]           = x[1];
}

    
/**
    Endmember fraction of abc
*/
void px_mb_abc(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 1.0 -x[0];
        p[1]           = x[0];
}

    
/**
    Endmember fraction of k4tr
*/
void px_mb_k4tr(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[0];
        p[1]           = x[1];
        p[2]           = -x[0] -x[1] + 1.0;
}

    
/**
    Endmember fraction of sp
*/
void px_mb_sp(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[1] + (x[0] - 1.0)*(x[2] + 1.0);
        p[1]           = (1.0 -x[0])*(x[2] + 1.0);
        p[2]           = -x[1] -x[2] + 1.0;
        p[3]           = x[2];
}

    
/**
    Endmember fraction of ilm
*/
void px_mb_ilm(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[1];
        p[1]           = x[0] -x[1];
        p[2]           = 1.0 -x[0];
}

    
/**
    Endmember fraction of ilmm
*/
void px_mb_ilmm(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[2];
        p[1]           = x[0] -x[1] -x[2];
        p[2]           = 1.0 -x[0];
        p[3]           = x[1];
}

    
/**
    Endmember fraction of ep
*/
void px_mb_ep(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[0] -x[1] + 1.0;
        p[1]           = 2.0*x[1];
        p[2]           = x[0] -x[1];
}

    
/**
    Endmember fraction of bi
*/
void px_mb_bi(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[0]*x[1] + x[0]*x[2] + x[0]*x[3] -x[0] -x[1] -x[2] -x[3] - 2.0/3.0*x[4] + 1.0;
        p[1]           = x[0] - 1.0/3.0*x[4];
        p[2]           = -x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + x[4];
        p[3]           = x[1];
        p[4]           = x[3];
        p[5]           = x[2];
}

    
/**
    Endmember fraction of mu
*/
void px_mb_mu(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[1] -x[2] -x[3] -x[4];
        p[1]           = x[0]*x[1] -x[0] -x[1] + 1.0;
        p[2]           = -x[0]*x[1] + x[0];
        p[3]           = x[3];
        p[4]           = x[4];
        p[5]           = x[2];
}

    
/**
    Endmember fraction of chl
*/
void px_mb_chl(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + 0.25*x[1]*x[4] + 1.25*x[1]*x[5] + 1.25*x[2]*x[5] - 0.25*x[3]*x[4] + 1.25*x[3]*x[5] + 2.0*x[3] - 0.25*x[4] - 1.25*x[5];
        p[1]           = 3.0*x[0]*x[1] + 2.0*x[0]*x[2] + x[0]*x[3] - 2.0*x[0] - 1.25*x[1]*x[4] - 2.25*x[1]*x[5] -x[1] - 2.25*x[2]*x[5] -x[2] + 1.25*x[3]*x[4] - 2.25*x[3]*x[5] -x[3] + 1.25*x[4] + 2.25*x[5] + 1.0;
        p[2]           = x[1] -x[3];
        p[3]           = x[0]*x[1] + x[0]*x[2] + x[0]*x[3] - 0.25*x[1]*x[4] - 1.25*x[1]*x[5] - 1.25*x[2]*x[5] + 0.25*x[3]*x[4] - 1.25*x[3]*x[5] + 0.25*x[4] + 1.25*x[5];
        p[4]           = -x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + x[0] + x[1]*x[5] + x[2]*x[5] + x[3]*x[5] -x[5];
        p[5]           = -2.0*x[0]*x[1] -x[0]*x[2] + x[0] + 1.25*x[1]*x[4] + 1.25*x[1]*x[5] + 1.25*x[2]*x[5] - 1.25*x[3]*x[4] + 1.25*x[3]*x[5] - 1.25*x[4] - 1.25*x[5];
        p[6]           = x[2];
}


    
/**
    Objective function of L
*/
double obj_mb_liq(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_liq(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = -x[7] +(x[7] + 1.0)*(x[0] + x[1] + x[3] + x[4] + x[5]);
    sf[1]          = x[0]*(x[7] + 1.0);
    sf[2]          = x[1]*x[2]*(x[7] + 1.0);
    sf[3]          = x[1]*(1.0 -x[2])*(x[7] + 1.0);
    sf[4]          = x[3]*(x[7] + 1.0) - x[7];
    sf[5]          = x[4]*(x[7] + 1.0) - x[7];
    sf[6]          = x[7] +(x[7] + 1.0)*(-1.0*x[0] -x[1] -x[3] -x[4] -x[5]) + 1.0;
    sf[7]          = x[7];
    sf[8]          = x[5]*(x[7] + 1.0);
    sf[9]          = x[6];
    sf[10]         = 1.0 - x[6];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*sf[1])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[0]*sf[2])) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[0]*sf[3])) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(sf[0]*sf[4])) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(sf[0]*sf[5])) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(sf[0]*sf[8]*cpow(sf[9], 5.0))) + mu_Gex[5];
    mu[6]          = gb[6] + R*T*creal(clog(sf[0]*cpow(sf[10], 5.0)*sf[8])) + mu_Gex[6];
    mu[7]          = gb[7] + R*T*creal(clog(cpow(sf[6], 2.0))) + mu_Gex[7];
    mu[8]          = gb[8] + R*T*creal(clog(sf[0]*sf[7])) + mu_Gex[8];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_liq(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of hb
*/
double obj_mb_hb(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_hb(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0 - x[3];
    sf[1]          = -x[3]*x[4] +x[3];
    sf[2]          =x[3]*x[4];
    sf[3]          = -x[0] +x[8] + 1.0;
    sf[4]          =x[0] - x[8];
    sf[5]          =x[0]*x[1] +x[0]*x[6] +x[0]*x[7] - x[0] - x[1]*x[9] - x[1] - x[6]*x[9] - x[6] - x[7]*x[9] - x[7] +x[9] + 1.0;
    sf[6]          = -x[0]*x[1] - x[0]*x[6] - x[0]*x[7] +x[0] +x[1]*x[9] +x[6]*x[9] +x[7]*x[9] - x[9];
    sf[7]          =x[1];
    sf[8]          =x[6];
    sf[9]          =x[7];
    sf[10]          =x[5];
    sf[11]          =x[0]*x[2] +x[0]*x[5] - x[0] +x[1]*x[9] - x[2] - x[5] +x[6]*x[9] +x[7]*x[9] - 1.5*x[8] - x[9] + 1.0;
    sf[12]          = -x[0]*x[2] - x[0]*x[5] +x[0] - x[1]*x[9] - x[6]*x[9] - x[7]*x[9] + 1.5*x[8] +x[9];
    sf[13]          =x[2];
    sf[14]          = -0.5*x[1] + 0.5*x[2] - 0.25*x[3] - 0.5*x[6] - 0.5*x[7] + 1.0;
    sf[15]          = 0.5*x[1] - 0.5*x[2] + 0.25*x[3] + 0.5*x[6] + 0.5*x[7];
    sf[16]          = 1.0 - x[7];
    sf[17]          =x[7];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*cpow(sf[10], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[5], 2.0))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(2.0*sf[0]*cpow(sf[10], 2.0)*sqrt(sf[14])*sqrt(sf[15])*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[7], 2.0))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(8.0*cpow(sf[10], 2.0)*sqrt(sf[14])*sqrt(sf[15])*cpow(sf[16], 2.0)*sf[1]*cpow(sf[3], 3.0)*sf[5]*sf[7])) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(sf[0]*cpow(sf[13], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[7], 2.0))) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(sf[0]*cpow(sf[11], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[5], 2.0))) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(sf[0]*cpow(sf[12], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[4], 3.0)*cpow(sf[6], 2.0))) + mu_Gex[5];
    mu[6]          = gb[6] + R*T*creal(clog(sf[0]*cpow(sf[12], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[6], 2.0))) + mu_Gex[6];
    mu[7]          = gb[7] + R*T*creal(clog(sf[0]*cpow(sf[12], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[4], 3.0)*cpow(sf[5], 2.0))) + mu_Gex[7];
    mu[8]          = gb[8] + R*T*creal(clog(sf[0]*cpow(sf[13], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[8], 2.0))) + mu_Gex[8];
    mu[9]          = gb[9] + R*T*creal(clog(8.0*cpow(sf[10], 2.0)*sqrt(sf[14])*sqrt(sf[15])*cpow(sf[16], 2.0)*sf[2]*cpow(sf[3], 3.0)*sf[5]*sf[7])) + mu_Gex[9];
    mu[10]          = gb[10] + R*T*creal(clog(2.0*sf[0]*cpow(sf[10], 2.0)*sqrt(sf[14])*sqrt(sf[15])*cpow(sf[17], 2.0)*cpow(sf[3], 3.0)*cpow(sf[9], 2.0))) + mu_Gex[10];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_hb(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of aug
*/
double obj_mb_aug(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_aug(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          =x[0]*x[1] +x[0]*x[4] - x[0] - x[1] - 0.5*x[3]*x[5] - 0.5*x[4]*x[5] - x[4] + 0.5*x[5] + 1.0;
    sf[1]          = -x[0]*x[1] - x[0]*x[4] +x[0] + 0.5*x[3]*x[5] + 0.5*x[4]*x[5] - 0.5*x[5];
    sf[2]          =x[1] - x[2] +x[4];
    sf[3]          =x[2];
    sf[4]          =x[0]*x[3] +x[0]*x[4] - x[0] + 0.5*x[3]*x[5] - x[3] + 0.5*x[4]*x[5] - x[4] - 0.5*x[5] + 1.0;
    sf[5]          = -x[0]*x[3] - x[0]*x[4] +x[0] - 0.5*x[3]*x[5] - 0.5*x[4]*x[5] + 0.5*x[5];
    sf[6]          =x[3];
    sf[7]          =x[4];
    sf[8]          = -0.5*x[1] + 0.5*x[6] + 1.0;
    sf[9]          = 0.5*x[1] - 0.5*x[6];
    sf[10]          = -0.5*x[1] - 0.5*x[6] + 1.0;
    sf[11]          = 0.5*x[1] + 0.5*x[6];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*cpow(sf[10], 0.25)*sf[6]*cpow(sf[8], 0.25))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[0]*cpow(sf[10], 0.25)*sf[4]*cpow(sf[8], 0.25))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(cpow(sf[10], 0.25)*sf[1]*sf[5]*cpow(sf[8], 0.25))) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(cpow(sf[10], 0.25)*sf[2]*sf[7]*cpow(sf[8], 0.25))) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(cpow(sf[10], 0.25)*sf[3]*sf[7]*cpow(sf[8], 0.25))) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(cpow(sf[11], 0.25)*sf[2]*sf[6]*cpow(sf[8], 0.25))) + mu_Gex[5];
    mu[6]          = gb[6] + R*T*creal(clog(1.4142*cpow(sf[10], 0.125)*cpow(sf[11], 0.125)*sf[2]*sf[6]*cpow(sf[8], 0.125)*cpow(sf[9], 0.125))) + mu_Gex[6];
    mu[7]          = gb[7] + R*T*creal(clog(sf[0]*cpow(sf[10], 0.25)*sf[5]*cpow(sf[8], 0.25))) + mu_Gex[7];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_aug(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of dio
*/
double obj_mb_dio(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_dio(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          =x[0]*x[1] - x[0]*x[3] - x[0] - x[1]*x[5] - x[1] - x[3]*x[5] +x[3] +x[5] + 1.0;
    sf[1]          = -x[0]*x[1] +x[0]*x[3] +x[0] +x[1]*x[5] +x[3]*x[5] - x[5];
    sf[2]          =x[1]*x[2] - x[4];
    sf[3]          = -x[1]*x[2] +x[1] - x[3] +x[4];
    sf[4]          =x[0]*x[1] +x[0]*x[3] - x[0] +x[1]*x[5] - x[1] +x[3]*x[5] - x[3] - x[5] + 1.0;
    sf[5]          = -x[0]*x[1] - x[0]*x[3] +x[0] - x[1]*x[5] - x[3]*x[5] +x[5];
    sf[6]          =x[1]*x[2] +x[4];
    sf[7]          = -x[1]*x[2] +x[1] +x[3] - x[4];
    sf[8]          =x[1] - x[3];
    sf[9]          = -x[1] +x[3] + 1.0;
    sf[10]          =x[1] +x[3];
    sf[11]          = -x[1] - x[3] + 1.0;
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sqrt(sf[10])*sqrt(sf[3])*sqrt(sf[7])*sqrt(sf[8]))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sqrt(sf[0])*sqrt(sf[11])*sqrt(sf[4])*sqrt(sf[9]))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sqrt(sf[11])*sqrt(sf[1])*sqrt(sf[5])*sqrt(sf[9]))) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(sqrt(sf[10])*sqrt(sf[2])*sqrt(sf[6])*sqrt(sf[8]))) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(sqrt(sf[0])*sqrt(sf[10])*sqrt(sf[7])*sqrt(sf[9]))) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(sqrt(sf[11])*sqrt(sf[1])*sqrt(sf[4])*sqrt(sf[9]))) + mu_Gex[5];
    mu[6]          = gb[6] + R*T*creal(clog(sqrt(sf[10])*sqrt(sf[3])*sqrt(sf[6])*sqrt(sf[8]))) + mu_Gex[6];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_dio(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of opx
*/
double obj_mb_opx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_opx(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          =x[0]*x[1] +x[0]*x[2] - x[0] - x[1] - x[2] - 0.5*x[3]*x[4] + 0.5*x[4] + 1.0;
    sf[1]          = -x[0]*x[1] - x[0]*x[2] +x[0] + 0.5*x[3]*x[4] - 0.5*x[4];
    sf[2]          =x[2];
    sf[3]          =x[1];
    sf[4]          =x[0]*x[3] - x[0] + 0.5*x[3]*x[4] - x[3] - 0.5*x[4] + 1.0;
    sf[5]          = -x[0]*x[3] +x[0] - 0.5*x[3]*x[4] + 0.5*x[4];
    sf[6]          =x[3];
    sf[7]          = 0.5*x[1] + 0.5*x[2];
    sf[8]          = -0.5*x[1] - 0.5*x[2] + 1.0;
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*sf[4]*sqrt(sf[8]))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[1]*sf[5]*sqrt(sf[8]))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[0]*sf[5]*sqrt(sf[8]))) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(1.4142*sf[3]*sf[4]*cpow(sf[7], 0.25)*cpow(sf[8], 0.25))) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(1.4142*sf[2]*sf[4]*cpow(sf[7], 0.25)*cpow(sf[8], 0.25))) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(sf[0]*sf[6]*sqrt(sf[8]))) + mu_Gex[5];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_opx(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of g
*/
double obj_mb_g(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_g(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          =x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[1]          = -x[0]*x[1] +x[0];
    sf[2]          =x[1];
    sf[3]          = 1.0 - x[2];
    sf[4]          =x[2];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(cpow(sf[0], 3.0)*cpow(sf[3], 2.0))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(cpow(sf[1], 3.0)*cpow(sf[3], 2.0))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(cpow(sf[2], 3.0)*cpow(sf[3], 2.0))) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(cpow(sf[0], 3.0)*cpow(sf[4], 2.0))) + mu_Gex[3];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_g(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of ol
*/
double obj_mb_ol(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_ol(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0 - x[0];
    sf[1]          =x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(cpow(sf[0], 2.0))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(cpow(sf[1], 2.0))) + mu_Gex[1];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_ol(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of pl4tr
*/
double obj_mb_pl4tr(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_pl4tr(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          = -x[0] - x[1] + 1.0;
    sf[1]          =x[0];
    sf[2]          =x[1];
    sf[3]          = 0.25*x[0] + 0.25;
    sf[4]          = 0.75 - 0.25*x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(1.7548*sf[0]*cpow(sf[3], 0.25)*cpow(sf[4], 0.75))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(2.0*sf[1]*sqrt(sf[3])*sqrt(sf[4]))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(1.7548*sf[2]*cpow(sf[3], 0.25)*cpow(sf[4], 0.75))) + mu_Gex[2];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_pl4tr(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of abc
*/
double obj_mb_abc(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_abc(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0 - x[0];
    sf[1]          =x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[1])) + mu_Gex[1];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_abc(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of k4tr
*/
double obj_mb_k4tr(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_k4tr(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          =x[0];
    sf[1]          =x[1];
    sf[2]          = -x[0] - x[1] + 1.0;
    sf[3]          = 0.25*x[1] + 0.25;
    sf[4]          = 0.75 - 0.25*x[1];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(1.7548*sf[0]*cpow(sf[3], 0.25)*cpow(sf[4], 0.75))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(2.0*sf[1]*sqrt(sf[3])*sqrt(sf[4]))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(1.7548*sf[2]*cpow(sf[3], 0.25)*cpow(sf[4], 0.75))) + mu_Gex[2];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_k4tr(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of sp
*/
double obj_mb_sp(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_sp(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          =x[1];
    sf[1]          = -x[1] - x[2] + 1.0;
    sf[2]          =x[2];
    sf[3]          = 1.0 - x[0];
    sf[4]          =x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*sf[4])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[0]*sf[3])) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[1]*sf[4])) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(sf[2]*sf[4])) + mu_Gex[3];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_sp(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of ilm
*/
double obj_mb_ilm(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_ilm(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 0.5*x[0] + 0.5*x[1];
    sf[1]          = 0.5*x[0] - 0.5*x[1];
    sf[2]          = 1.0 - x[0];
    sf[3]          = 0.5*x[0] - 0.5*x[1];
    sf[4]          = 0.5*x[0] + 0.5*x[1];
    sf[5]          = 1.0 - x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*sf[4])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(4.0*sqrt(sf[0])*sqrt(sf[1])*sqrt(sf[3])*sqrt(sf[4]))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[2]*sf[5])) + mu_Gex[2];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_ilm(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of ilmm
*/
double obj_mb_ilmm(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_ilmm(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 0.5*x[0] - 0.5*x[1] + 0.5*x[2];
    sf[1]          = 0.5*x[0] - 0.5*x[1] - 0.5*x[2];
    sf[2]          =x[1];
    sf[3]          = 1.0 - x[0];
    sf[4]          = 0.5*x[0] - 0.5*x[1] - 0.5*x[2];
    sf[5]          = 0.5*x[0] + 0.5*x[1] + 0.5*x[2];
    sf[6]          = 1.0 - x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*sf[5])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(4.0*sqrt(sf[0])*sqrt(sf[1])*sqrt(sf[4])*sqrt(sf[5]))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[3]*sf[6])) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(sf[2]*sf[5])) + mu_Gex[3];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_ilmm(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of ep
*/
double obj_mb_ep(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_ep(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          =x[0] - x[1];
    sf[1]          = -x[0] +x[1] + 1.0;
    sf[2]          =x[0] +x[1];
    sf[3]          = -x[0] - x[1] + 1.0;
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[1]*sf[3])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[1]*sf[2])) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[0]*sf[2])) + mu_Gex[2];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_ep(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of bi
*/
double obj_mb_bi(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_bi(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          =x[0]*x[1] +x[0]*x[2] +x[0]*x[3] - x[0] - x[1] - x[2] - x[3] - 2.0/3.0*x[4] + 1.0;
    sf[1]          = -x[0]*x[1] - x[0]*x[2] - x[0]*x[3] +x[0] + 2.0/3.0*x[4];
    sf[2]          =x[2];
    sf[3]          =x[3];
    sf[4]          =x[1];
    sf[5]          = -x[0] + 1.0/3.0*x[4] + 1.0;
    sf[6]          =x[0] - 1.0/3.0*x[4];
    sf[7]          = -0.5*x[1] - 0.5*x[2] + 0.5;
    sf[8]          = 0.5*x[1] + 0.5*x[2] + 0.5;
    sf[9]          = 1.0 - x[3];
    sf[10]          =x[3];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(4.0*sf[0]*cpow(sf[5], 2.0)*sf[7]*sf[8]*cpow(sf[9], 2.0))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(4.0*sf[1]*cpow(sf[6], 2.0)*sf[7]*sf[8]*cpow(sf[9], 2.0))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(4.0*sf[1]*cpow(sf[5], 2.0)*sf[7]*sf[8]*cpow(sf[9], 2.0))) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(sf[4]*cpow(sf[5], 2.0)*cpow(sf[8], 2.0)*cpow(sf[9], 2.0))) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(4.0*cpow(sf[10], 2.0)*sf[3]*cpow(sf[5], 2.0)*sf[7]*sf[8])) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(sf[2]*cpow(sf[5], 2.0)*cpow(sf[8], 2.0)*cpow(sf[9], 2.0))) + mu_Gex[5];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_bi(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of mu
*/
double obj_mb_mu(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_mu(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          = -x[3] - x[4] + 1.0;
    sf[1]          =x[3];
    sf[2]          =x[4];
    sf[3]          =x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[4]          = -x[0]*x[1] +x[0];
    sf[5]          =x[1];
    sf[6]          = 1.0 - x[2];
    sf[7]          =x[2];
    sf[8]          = -0.5*x[1] - 0.5*x[4] + 1.0;
    sf[9]          = 0.5*x[1] + 0.5*x[4];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(4.0*sf[0]*sf[5]*sf[6]*sf[8]*sf[9])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[0]*sf[3]*sf[6]*cpow(sf[8], 2.0))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[0]*sf[4]*sf[6]*cpow(sf[8], 2.0))) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(4.0*sf[1]*sf[5]*sf[6]*sf[8]*sf[9])) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(sf[2]*sf[5]*sf[6]*cpow(sf[9], 2.0))) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(4.0*sf[0]*sf[5]*sf[7]*sf[8]*sf[9])) + mu_Gex[5];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_mu(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of chl
*/
double obj_mb_chl(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mb_chl(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          =x[0]*x[1] - x[0]*x[3] - x[0] - x[1]*x[4] - x[1] +x[3]*x[4] +x[3] +x[4] + 1.0;
    sf[1]          = -x[0]*x[1] +x[0]*x[3] +x[0] +x[1]*x[4] - x[3]*x[4] - x[4];
    sf[2]          =x[1] - x[3];
    sf[3]          = -x[0] + 0.25*x[1]*x[4] + 0.25*x[1]*x[5] + 0.25*x[2]*x[5] - 0.25*x[3]*x[4] + 0.25*x[3]*x[5] - 0.25*x[4] - 0.25*x[5] + 1.0;
    sf[4]          =x[0] - 0.25*x[1]*x[4] - 0.25*x[1]*x[5] - 0.25*x[2]*x[5] + 0.25*x[3]*x[4] - 0.25*x[3]*x[5] + 0.25*x[4] + 0.25*x[5];
    sf[5]          =x[0]*x[1] +x[0]*x[2] +x[0]*x[3] - x[0] - x[1]*x[5] - x[1] - x[2]*x[5] - x[2] - x[3]*x[5] - x[3] +x[5] + 1.0;
    sf[6]          = -x[0]*x[1] - x[0]*x[2] - x[0]*x[3] +x[0] +x[1]*x[5] +x[2]*x[5] +x[3]*x[5] - x[5];
    sf[7]          =x[2];
    sf[8]          =x[1] +x[3];
    sf[9]          = -x[1] - 0.5*x[2] + 1.0;
    sf[10]          =x[1] + 0.5*x[2];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(4.0*sf[0]*sf[10]*cpow(sf[3], 4.0)*sf[8]*sf[9])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[0]*cpow(sf[3], 4.0)*sf[5]*cpow(sf[9], 2.0))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(cpow(sf[10], 2.0)*sf[2]*cpow(sf[3], 4.0)*sf[8])) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(4.0*sf[10]*sf[1]*cpow(sf[4], 4.0)*sf[8]*sf[9])) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(sf[0]*cpow(sf[4], 4.0)*sf[6]*cpow(sf[9], 2.0))) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(sf[1]*cpow(sf[3], 4.0)*sf[5]*cpow(sf[9], 2.0))) + mu_Gex[5];
    mu[6]          = gb[6] + R*T*creal(clog(4.0*sf[0]*sf[10]*cpow(sf[3], 4.0)*sf[7]*sf[9])) + mu_Gex[6];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mb_chl(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}



/**************************************************************************************/
/**************************************************************************************/
/**********************ULTRAMAFIC DATABASE (Evans & Frost., 2021)**********************/
/**************************************************************************************/
/**************************************************************************************/

/**
    Update dpdx matrix of fluid
*/
void dpdx_um_fluid(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 1.0;      
    dp_dx[1][0] = -1.0;      
}


/**
    Update dpdx matrix of ol
*/
void dpdx_um_ol(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.0;      
    dp_dx[1][0] = 1.0;      
}


/**
    Update dpdx matrix of br
*/
void dpdx_um_br(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.0;      
    dp_dx[1][0] = 1.0;      
}


/**
    Update dpdx matrix of ch
*/
void dpdx_um_ch(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.0;      
    dp_dx[1][0] = 1.0;      
}


/**
    Update dpdx matrix of atg
*/
void dpdx_um_atg(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[1] + x[2] - 1.0;      dp_dx[0][1] = x[0] -x[3] - 1.0;      dp_dx[0][2] = x[0] -x[3] - 1.0;      dp_dx[0][3] = -x[1] -x[2] + 1.0;      
    dp_dx[1][0] = 1.0;      dp_dx[1][1] = -0.5*x[3];      dp_dx[1][2] = -0.5*x[3];      dp_dx[1][3] = -0.5*x[1] - 0.5*x[2] + 0.5;      
    dp_dx[2][0] = -x[1] -x[2];      dp_dx[2][1] = -x[0] + 1.5*x[3];      dp_dx[2][2] = -x[0] + 1.5*x[3];      dp_dx[2][3] = 1.5*x[1] + 1.5*x[2] - 1.5;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 1.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 1.0;      dp_dx[4][3] = 0.0;      
}


/**
    Update dpdx matrix of g
*/
void dpdx_um_g(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.0;      
    dp_dx[1][0] = 1.0;      
}


/**
    Update dpdx matrix of ta
*/
void dpdx_um_ta(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[1] + x[2] - 1.0;      dp_dx[0][1] = x[0] - 2.0;      dp_dx[0][2] = x[0] - 2.0;      dp_dx[0][3] = 0.5*x[4] + 1.0;      dp_dx[0][4] = 0.5*x[3] - 0.5;      
    dp_dx[1][0] = 1.0 -x[3];      dp_dx[1][1] = 0.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = -x[0] + x[4];      dp_dx[1][4] = x[3] - 1.0;      
    dp_dx[2][0] = -x[1] -x[2] + x[3];      dp_dx[2][1] = -x[0];      dp_dx[2][2] = -x[0];      dp_dx[2][3] = x[0] - 1.5*x[4];      dp_dx[2][4] = 1.5 - 1.5*x[3];      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 2.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = -2.0;      dp_dx[3][4] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 2.0;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 1.0;      dp_dx[5][4] = 0.0;      
}


/**
    Update dpdx matrix of chl
*/
void dpdx_um_chl(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -x[1] -x[2] -x[3];      dp_dx[0][1] = -x[0] + 0.25*x[4] + 1.25*x[5];      dp_dx[0][2] = -x[0] + 1.25*x[5];      dp_dx[0][3] = -x[0] - 0.25*x[4] + 1.25*x[5] + 2.0;      dp_dx[0][4] = 0.25*x[1] - 0.25*x[3] - 0.25;      dp_dx[0][5] = 1.25*x[1] + 1.25*x[2] + 1.25*x[3] - 1.25;      
    dp_dx[1][0] = 3.0*x[1] + 2.0*x[2] + x[3] - 2.0;      dp_dx[1][1] = 3.0*x[0] - 1.25*x[4] - 2.25*x[5] - 1.0;      dp_dx[1][2] = 2.0*x[0] - 2.25*x[5] - 1.0;      dp_dx[1][3] = x[0] + 1.25*x[4] - 2.25*x[5] - 1.0;      dp_dx[1][4] = -1.25*x[1] + 1.25*x[3] + 1.25;      dp_dx[1][5] = -2.25*x[1] - 2.25*x[2] - 2.25*x[3] + 2.25;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.0;      dp_dx[2][2] = 0.0;      dp_dx[2][3] = -1.0;      dp_dx[2][4] = 0.0;      dp_dx[2][5] = 0.0;      
    dp_dx[3][0] = x[1] + x[2] + x[3];      dp_dx[3][1] = x[0] - 0.25*x[4] - 1.25*x[5];      dp_dx[3][2] = x[0] - 1.25*x[5];      dp_dx[3][3] = x[0] + 0.25*x[4] - 1.25*x[5];      dp_dx[3][4] = -0.25*x[1] + 0.25*x[3] + 0.25;      dp_dx[3][5] = -1.25*x[1] - 1.25*x[2] - 1.25*x[3] + 1.25;      
    dp_dx[4][0] = -x[1] -x[2] -x[3] + 1.0;      dp_dx[4][1] = -x[0] + x[5];      dp_dx[4][2] = -x[0] + x[5];      dp_dx[4][3] = -x[0] + x[5];      dp_dx[4][4] = 0.0;      dp_dx[4][5] = x[1] + x[2] + x[3] - 1.0;      
    dp_dx[5][0] = -2.0*x[1] -x[2] + 1.0;      dp_dx[5][1] = -2.0*x[0] + 1.25*x[4] + 1.25*x[5];      dp_dx[5][2] = -x[0] + 1.25*x[5];      dp_dx[5][3] = -1.25*x[4] + 1.25*x[5];      dp_dx[5][4] = 1.25*x[1] - 1.25*x[3] - 1.25;      dp_dx[5][5] = 1.25*x[1] + 1.25*x[2] + 1.25*x[3] - 1.25;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 1.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = 0.0;      
}



/**
    Update dpdx matrix of anth
*/
void dpdx_um_anth(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.00;      dp_dx[0][1] = x[3] - 1.0;      dp_dx[0][2] = -1.50;      dp_dx[0][3] = x[1] - 1.0;      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = 1.00;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      
    dp_dx[2][0] = 1.0 - x[1];      dp_dx[2][1] =-x[0] + 2.0*x[3];      dp_dx[2][2] = -2.50;      dp_dx[2][3] = 2.0*x[1] - 2.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] =-x[3];      dp_dx[3][2] = 2.50;      dp_dx[3][3] = 1.0 - x[1];      
    dp_dx[4][0] = x[1];      dp_dx[4][1] = x[0] - 2.0*x[3];      dp_dx[4][2] = 1.50;      dp_dx[4][3] = 2.0 - 2.0*x[1];      
}


/**
    Update dpdx matrix of spi
*/
void dpdx_um_spi(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 1.0;      dp_dx[0][1] = 1.0;      
    dp_dx[1][0] = -1.0;      dp_dx[1][1] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = -1.0;      
}


/**
    Update dpdx matrix of opx
*/
void dpdx_um_opx(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.0;      dp_dx[0][1] = -1.0;      dp_dx[0][2] = -1.0;      dp_dx[0][3] = -0.50;      
    dp_dx[1][0] = -x[1] -x[2] + 1.0;      dp_dx[1][1] = -x[0];      dp_dx[1][2] = -x[0];      dp_dx[1][3] = -0.50;      
    dp_dx[2][0] = x[1] + x[2];      dp_dx[2][1] = x[0];      dp_dx[2][2] = x[0];      dp_dx[2][3] = 1.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 1.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 1.0;      dp_dx[4][3] = 0.0;      
}


/**
    Update dpdx matrix of po
*/
void dpdx_um_po(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 8.0;      
    dp_dx[1][0] = -8.0;      
}

//-----------------------px for ev------------------- 
//--------------------------------------------------- 
/**
    Endmember fraction of fluid
*/
void px_um_fluid(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[0];
        p[1]           = 1.0 -x[0];
}

    
/**
    Endmember fraction of ol
*/
void px_um_ol(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 1.0 -x[0];
        p[1]           = x[0];
}

    
/**
    Endmember fraction of br
*/
void px_um_br(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 1.0 -x[0];
        p[1]           = x[0];
}

    
/**
    Endmember fraction of ch
*/
void px_um_ch(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 1.0 -x[0];
        p[1]           = x[0];
}

    
/**
    Endmember fraction of atg
*/
void px_um_atg(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[0]*x[1] + x[0]*x[2] -x[0] -x[1]*x[3] -x[1] -x[2]*x[3] -x[2] + x[3] + 1.0;
        p[1]           = x[0] - 0.5*x[1]*x[3] - 0.5*x[2]*x[3] + 0.5*x[3];
        p[2]           = -x[0]*x[1] -x[0]*x[2] + 1.5*x[1]*x[3] + 1.5*x[2]*x[3] - 1.5*x[3];
        p[3]           = x[1];
        p[4]           = x[2];
}

    
/**
    Endmember fraction of g
*/
void px_um_g(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 1.0 -x[0];
        p[1]           = x[0];
}

    
/**
    Endmember fraction of ta
*/
void px_um_ta(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[0]*x[1] + x[0]*x[2] -x[0] - 2.0*x[1] - 2.0*x[2] + 0.5*x[3]*x[4] + x[3] - 0.5*x[4] + 1.0;
        p[1]           = -x[0]*x[3] + x[0] + x[3]*x[4] -x[4];
        p[2]           = -x[0]*x[1] -x[0]*x[2] + x[0]*x[3] - 1.5*x[3]*x[4] + 1.5*x[4];
        p[3]           = 2.0*x[1] - 2.0*x[3];
        p[4]           = 2.0*x[2];
        p[5]           = x[3];
}

    
/**
    Endmember fraction of chl
*/
void px_um_chl(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + 0.25*x[1]*x[4] + 1.25*x[1]*x[5] + 1.25*x[2]*x[5] - 0.25*x[3]*x[4] + 1.25*x[3]*x[5] + 2.0*x[3] - 0.25*x[4] - 1.25*x[5];
        p[1]           = 3.0*x[0]*x[1] + 2.0*x[0]*x[2] + x[0]*x[3] - 2.0*x[0] - 1.25*x[1]*x[4] - 2.25*x[1]*x[5] -x[1] - 2.25*x[2]*x[5] -x[2] + 1.25*x[3]*x[4] - 2.25*x[3]*x[5] -x[3] + 1.25*x[4] + 2.25*x[5] + 1.0;
        p[2]           = x[1] -x[3];
        p[3]           = x[0]*x[1] + x[0]*x[2] + x[0]*x[3] - 0.25*x[1]*x[4] - 1.25*x[1]*x[5] - 1.25*x[2]*x[5] + 0.25*x[3]*x[4] - 1.25*x[3]*x[5] + 0.25*x[4] + 1.25*x[5];
        p[4]           = -x[0]*x[1] -x[0]*x[2] -x[0]*x[3] + x[0] + x[1]*x[5] + x[2]*x[5] + x[3]*x[5] -x[5];
        p[5]           = -2.0*x[0]*x[1] -x[0]*x[2] + x[0] + 1.25*x[1]*x[4] + 1.25*x[1]*x[5] + 1.25*x[2]*x[5] - 1.25*x[3]*x[4] + 1.25*x[3]*x[5] - 1.25*x[4] - 1.25*x[5];
        p[6]           = x[2];
}

    
/**
    Endmember fraction of anth
*/
void px_um_anth(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           =-x[0] + x[1]*x[3] - x[1] - 1.5*x[2] - x[3] + 1.0;
        p[1]           = x[1];
        p[2]           =-x[0]*x[1] + x[0] + 2.0*x[1]*x[3] - 2.5*x[2] - 2.0*x[3];
        p[3]           =-x[1]*x[3] + 2.5*x[2] + x[3];
        p[4]           = x[0]*x[1] - 2.0*x[1]*x[3] + 1.5*x[2] + 2.0*x[3];
}

/**
    Endmember fraction of spi
*/
void px_um_spi(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[0] + x[1] - 1.0;
        p[1]           = 1.0 -x[0];
        p[2]           = 1.0 -x[1];
}

    
/**
    Endmember fraction of opx
*/
void px_um_opx(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[0] -x[1] -x[2] - 0.5*x[3] + 1.0;
        p[1]           = -x[0]*x[1] -x[0]*x[2] + x[0] - 0.5*x[3];
        p[2]           = x[0]*x[1] + x[0]*x[2] + x[3];
        p[3]           = x[1];
        p[4]           = x[2];
}

    
/**
    Endmember fraction of po
*/
void px_um_po(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 8.0*x[0];
        p[1]           = 1.0 - 8.0*x[0];
}


//-----------------------obj for ev------------------- 
//--------------------------------------------------- 

    
/**
    Objective function of fluid
*/
double obj_um_fluid(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_fluid(SS_ref_db,x);

    
    sf[0]          = x[0];
    sf[1]          = 1.0 - x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]));
    mu[1]          = gb[1] + R*T*creal(clog(sf[1]));
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_fluid(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of ol
*/
double obj_um_ol(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_ol(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0 - x[0];
    sf[1]          = x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(cpow(sf[0], 2.0))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(cpow(sf[1], 2.0))) + mu_Gex[1];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_ol(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of br
*/
double obj_um_br(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_br(SS_ref_db,x);

    
    sf[0]          = 1.0 - x[0];
    sf[1]          = x[0];
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]));
    mu[1]          = gb[1] + R*T*creal(clog(sf[1]));
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_br(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of ch
*/
double obj_um_ch(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_ch(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0 - x[0];
    sf[1]          = x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(cpow(sf[0], 9.0))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(cpow(sf[1], 9.0))) + mu_Gex[1];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_ch(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of atg
*/
double obj_um_atg(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_atg(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[0]*x[1] + x[0]*x[2] - x[0] - x[1]*x[3] - x[1] - x[2]*x[3] - x[2] + x[3] + 1.0;
    sf[1]          = -x[0]*x[1] - x[0]*x[2] + x[0] + x[1]*x[3] + x[2]*x[3] - x[3];
    sf[2]          = x[2];
    sf[3]          = x[1];
    sf[4]          = -x[0] + 0.5*x[1]*x[3] + 0.5*x[2]*x[3] - 0.5*x[3] + 1.0;
    sf[5]          = x[0] - 0.5*x[1]*x[3] - 0.5*x[2]*x[3] + 0.5*x[3];
    sf[6]          = -0.5*x[1] - 0.5*x[2] + 1.0;
    sf[7]          = 0.5*x[1] + 0.5*x[2];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*cpow(sf[4], 2.0)*cpow(sf[6], 2.0))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[1]*cpow(sf[5], 2.0)*cpow(sf[6], 2.0))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[1]*cpow(sf[4], 2.0)*cpow(sf[6], 2.0))) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(4.0*sf[3]*cpow(sf[4], 2.0)*sf[6]*sf[7])) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(4.0*sf[2]*cpow(sf[4], 2.0)*sf[6]*sf[7])) + mu_Gex[4];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_atg(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of g
*/
double obj_um_g(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_g(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0 - x[0];
    sf[1]          = x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(cpow(sf[0], 3.0))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(cpow(sf[1], 3.0))) + mu_Gex[1];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_g(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of ta
*/
double obj_um_ta(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_ta(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[0]*x[3] - x[0] - x[3]*x[4] - x[3] + x[4] + 1.0;
    sf[1]          = -x[0]*x[3] + x[0] + x[3]*x[4] - x[4];
    sf[2]          = x[3];
    sf[3]          = x[0]*x[1] + x[0]*x[2] - x[0] - x[1] - x[2] + 0.5*x[3]*x[4] - 0.5*x[4] + 1.0;
    sf[4]          = -x[0]*x[1] - x[0]*x[2] + x[0] - 0.5*x[3]*x[4] + 0.5*x[4];
    sf[5]          = x[2];
    sf[6]          = x[1];
    sf[7]          = -x[1] - x[2] + x[3] + 1.0;
    sf[8]          = x[1] + x[2] - x[3];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*cpow(sf[3], 2.0)*cpow(sf[7], 2.0))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[1]*cpow(sf[4], 2.0)*cpow(sf[7], 2.0))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[0]*cpow(sf[4], 2.0)*cpow(sf[7], 2.0))) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(16.0*sf[0]*sf[3]*sf[6]*sf[7]*sf[8])) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(16.0*sf[0]*sf[3]*sf[5]*sf[7]*sf[8])) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(sf[2]*cpow(sf[6], 2.0)*cpow(sf[7], 2.0))) + mu_Gex[5];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_ta(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of chl
*/
double obj_um_chl(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_chl(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[0]*x[1] - x[0]*x[3] - x[0] - x[1]*x[4] - x[1] + x[3]*x[4] + x[3] + x[4] + 1.0;
    sf[1]          = -x[0]*x[1] + x[0]*x[3] + x[0] + x[1]*x[4] - x[3]*x[4] - x[4];
    sf[2]          = x[1] - x[3];
    sf[3]          = -x[0] + 0.25*x[1]*x[4] + 0.25*x[1]*x[5] + 0.25*x[2]*x[5] - 0.25*x[3]*x[4] + 0.25*x[3]*x[5] - 0.25*x[4] - 0.25*x[5] + 1.0;
    sf[4]          = x[0] - 0.25*x[1]*x[4] - 0.25*x[1]*x[5] - 0.25*x[2]*x[5] + 0.25*x[3]*x[4] - 0.25*x[3]*x[5] + 0.25*x[4] + 0.25*x[5];
    sf[5]          = x[0]*x[1] + x[0]*x[2] + x[0]*x[3] - x[0] - x[1]*x[5] - x[1] - x[2]*x[5] - x[2] - x[3]*x[5] - x[3] + x[5] + 1.0;
    sf[6]          = -x[0]*x[1] - x[0]*x[2] - x[0]*x[3] + x[0] + x[1]*x[5] + x[2]*x[5] + x[3]*x[5] - x[5];
    sf[7]          = x[2];
    sf[8]          = x[1] + x[3];
    sf[9]          = -x[1] - 0.5*x[2] + 1.0;
    sf[10]          = x[1] + 0.5*x[2];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(4.0*sf[0]*sf[10]*cpow(sf[3], 4.0)*sf[8]*sf[9])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[0]*cpow(sf[3], 4.0)*sf[5]*cpow(sf[9], 2.0))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(cpow(sf[10], 2.0)*sf[2]*cpow(sf[3], 4.0)*sf[8])) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(4.0*sf[10]*sf[1]*cpow(sf[4], 4.0)*sf[8]*sf[9])) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(sf[0]*cpow(sf[4], 4.0)*sf[6]*cpow(sf[9], 2.0))) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(sf[1]*cpow(sf[3], 4.0)*sf[5]*cpow(sf[9], 2.0))) + mu_Gex[5];
    mu[6]          = gb[6] + R*T*creal(clog(4.0*sf[0]*sf[10]*cpow(sf[3], 4.0)*sf[7]*sf[9])) + mu_Gex[6];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_chl(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of anth
*/
double obj_um_anth(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_anth(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          = -x[0] + x[1]*x[3] - 1.5*x[2] - x[3] + 1.0;
    sf[1]          = x[0] - x[1]*x[3] + 1.5*x[2] + x[3];
    sf[2]          = -x[0] + x[2] + 1.0;
    sf[3]          = x[0] - x[2];
    sf[4]          = x[1];
    sf[5]          = x[0]*x[1] - x[0] - x[1]*x[3] - x[1] + x[3] + 1.0;
    sf[6]          = -x[0]*x[1] + x[0] + x[1]*x[3] - x[3];
    sf[7]          = 0.5*x[1];
    sf[8]          = 1.0 - 0.5*x[1];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(cpow(sf[0], 2.0)*cpow(sf[2], 3.0)*cpow(sf[5], 2.0)*sf[8])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(2.0*cpow(sf[0], 2.0)*cpow(sf[2], 3.0)*cpow(sf[4], 2.0)*sqrt(sf[7])*sqrt(sf[8]))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(cpow(sf[1], 2.0)*cpow(sf[3], 3.0)*cpow(sf[6], 2.0)*sf[8])) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(cpow(sf[1], 2.0)*cpow(sf[2], 3.0)*cpow(sf[6], 2.0)*sf[8])) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(cpow(sf[1], 2.0)*cpow(sf[3], 3.0)*cpow(sf[5], 2.0)*sf[8])) + mu_Gex[4];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_anth(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}   
    
/**
    Objective function of spi
*/
double obj_um_spi(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_spi(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[1];
    sf[1]          = 1.0 - x[1];
    sf[2]          = 1.0 - x[0];
    sf[3]          = x[0];
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*sf[3])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[0]*sf[2])) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[1]*sf[3])) + mu_Gex[2];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_spi(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of opx
*/
double obj_um_opx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_opx(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[0]*x[1] + x[0]*x[2] - x[0] - x[1] - x[2] + 0.5*x[3] + 1.0;
    sf[1]          = -x[0]*x[1] - x[0]*x[2] + x[0] - 0.5*x[3];
    sf[2]          = x[2];
    sf[3]          = x[1];
    sf[4]          = -x[0] - 0.5*x[3] + 1.0;
    sf[5]          = x[0] + 0.5*x[3];
    sf[6]          = 0.5*x[1] + 0.5*x[2];
    sf[7]          = -0.5*x[1] - 0.5*x[2] + 1.0;
    
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*sf[4]*sqrt(sf[7]))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[1]*sf[5]*sqrt(sf[7]))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[0]*sf[5]*sqrt(sf[7]))) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(1.4142*sf[3]*sf[4]*cpow(sf[6], 0.25)*cpow(sf[7], 0.25))) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(1.4142*sf[2]*sf[4]*cpow(sf[6], 0.25)*cpow(sf[7], 0.25))) + mu_Gex[4];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_opx(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of po
*/
double obj_um_po(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_um_po(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0 - x[0];
    sf[1]          = x[0];
    
    mu[0]          = gb[0] + R*T*creal(clog(1.4576*cpow(sf[0], 0.875)*cpow(sf[1], 0.125))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[0])) + mu_Gex[1];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_um_po(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}



//-----------------------p2x for ev------------------- 
//--------------------------------------------------- 

/**
    Endmember to xeos for fluid
*/
void p2x_um_fluid(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = d->p[0];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for ol
*/
void p2x_um_ol(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = d->p[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}


/**
    Endmember to xeos for br
*/
void p2x_um_br(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = d->p[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for ch
*/
void p2x_um_ch(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = d->p[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for atg
*/
void p2x_um_atg(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = (3.0*d->p[1] + d->p[2])/(d->p[0] + d->p[1] + d->p[2] + 2.0);
    d->iguess[1]  = d->p[3];
    d->iguess[2]  = d->p[4];
    d->iguess[3]  = d->iguess[0] - (d->p[1]+d->p[2])/(d->p[1]+d->p[2]+d->p[0]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for g
*/
void p2x_um_g(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = d->p[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for ta
*/
void p2x_um_ta(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = (3*d->p[1]+2*d->p[2])/(1+2*(d->p[0]+d->p[1]+d->p[2])+d->p[3]+d->p[4]-d->p[5]);
    d->iguess[1]  = d->p[5]+0.5*d->p[3];
    d->iguess[2]  = 0.5*d->p[4];
    d->iguess[3]  = d->p[5];
    d->iguess[4]  = d->iguess[0]-d->p[1]/(1-d->p[5]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for chl
*/
void p2x_um_chl(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]   = d->p[6];
    d->iguess[3] = d->p[0]/2.0 + d->p[3]/2.0;
    d->iguess[1]   = d->p[2] + d->iguess[3];
    d->iguess[0]  = (5.0*(d->p[3] + d->p[4]) + d->p[5])/(5.0 + d->p[1] - d->p[2] + d->p[4] + d->p[5]);
    d->iguess[5] = d->iguess[0] - d->p[4]/(d->p[1] + d->p[4] + d->p[5]);
    d->iguess[4] = d->iguess[0] - (d->p[3] + d->p[5])/(1 - d->p[2]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}


/**
    Endmember to xeos for anth
*/
void p2x_um_anth(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[1] = d->p[1];
    d->iguess[0] = (2.0*(1.0-d->p[0]-d->p[1])+5.0*d->p[2]+2.0*d->p[3]+3.0*d->p[4])/(7.0-2.0*(d->p[1]));
    d->iguess[2] = d->iguess[0] - d->p[2] -d->p[4];
    d->iguess[3] = (d->p[3] - d->p[4] - d->iguess[2] + d->iguess[0]*d->iguess[1])/(d->iguess[1] - 1.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}



/**
    Endmember to xeos for spi
*/
void p2x_um_spi(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[1]  = d->p[0]+d->p[1];
    d->iguess[0]  = d->iguess[1] - d->p[1] + d->p[2];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for opx
*/
void p2x_um_opx(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = (2*d->p[1]+d->p[2])/(2-d->p[3]-d->p[4]);
    d->iguess[1]  = d->p[3];
    d->iguess[2]  = d->p[4];
    d->iguess[3]  = 2*(d->p[1]+d->p[2]-d->iguess[0]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for po
*/
void p2x_um_po(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = (1.0-d->p[1])/8.0;
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}


/**************************************************************************************/
/**************************************************************************************/
/*********************METAPELITE DATABASE (White et al., 2014)*************************/
/**************************************************************************************/
/**************************************************************************************/
/**
    Update dpdx matrix of liq_mp
*/
void dpdx_mp_liq(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 1.0;      dp_dx[0][1] = 0.0;      dp_dx[0][2] = 0.0;      dp_dx[0][3] = 0.0;      dp_dx[0][4] = 0.0;      dp_dx[0][5] = 0.0;      dp_dx[0][6] = 0.0;      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = x[2];      dp_dx[1][2] = x[1];      dp_dx[1][3] = 0.0;      dp_dx[1][4] = 0.0;      dp_dx[1][5] = 0.0;      dp_dx[1][6] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.0 - x[2];      dp_dx[2][2] = -x[1];      dp_dx[2][3] = 0.0;      dp_dx[2][4] = 0.0;      dp_dx[2][5] = 0.0;      dp_dx[2][6] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 1.0;      dp_dx[3][4] = 0.0;      dp_dx[3][5] = 0.0;      dp_dx[3][6] = 0.0;      
    dp_dx[4][0] = -1.0;      dp_dx[4][1] = -1.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = -1.0;      dp_dx[4][4] = -1.0;      dp_dx[4][5] = 0.0;      dp_dx[4][6] = -1.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 1.0 - x[5];      dp_dx[5][5] = -x[4];      dp_dx[5][6] = 0.0;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 0.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = x[5];      dp_dx[6][5] = x[4];      dp_dx[6][6] = 0.0;      
    dp_dx[7][0] = 0.0;      dp_dx[7][1] = 0.0;      dp_dx[7][2] = 0.0;      dp_dx[7][3] = 0.0;      dp_dx[7][4] = 0.0;      dp_dx[7][5] = 0.0;      dp_dx[7][6] = 1.0;      
}


/**
    Update dpdx matrix of pl4tr_mp
*/
void dpdx_mp_pl4tr(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.0;      dp_dx[0][1] = -1.0;      
    dp_dx[1][0] = 1.0;      dp_dx[1][1] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.0;      
}


/**
    Update dpdx matrix of bi_mp
*/
void dpdx_mp_bi(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[3] + 3.0*x[1] + x[4] + x[2] - 1.0;      dp_dx[0][1] = 3.0*x[0] - 1.0;      dp_dx[0][2] = x[0] - 1.0;      dp_dx[0][3] = x[0] - 1.0;      dp_dx[0][4] = x[0] - 1.0;      dp_dx[0][5] = -2.0/3.0;      
    dp_dx[1][0] = 1.0;      dp_dx[1][1] = 0.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      dp_dx[1][4] = 0.0;      dp_dx[1][5] = -1.0/3.0;      
    dp_dx[2][0] = -x[3] - 3.0*x[1] - x[4] - x[2];      dp_dx[2][1] = -3.0*x[0];      dp_dx[2][2] = -x[0];      dp_dx[2][3] = -x[0];      dp_dx[2][4] = -x[0];      dp_dx[2][5] = 1.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      dp_dx[3][5] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = 1.0;      dp_dx[4][5] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 1.0;      dp_dx[5][4] = 0.0;      dp_dx[5][5] = 0.0;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 1.0;      dp_dx[6][2] = 0.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = 0.0;      
}



/**
    Update dpdx matrix of g_mp
*/
void dpdx_mp_g(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[2] + x[1] - 1.0;      dp_dx[0][1] = x[0] - 1.0;      dp_dx[0][2] = x[0] - 1.0;      dp_dx[0][3] = -1.0;      
    dp_dx[1][0] = -x[2] - x[1] + 1.0;      dp_dx[1][1] = -x[0];      dp_dx[1][2] = -x[0];      dp_dx[1][3] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 0.0;      dp_dx[2][2] = 1.0;      dp_dx[2][3] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 1.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 1.0;      
}



/**
    Update dpdx matrix of ep_mp
*/
void dpdx_mp_ep(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.0;      dp_dx[0][1] = -1.0;      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = 2.0;      
    dp_dx[2][0] = 1.0;      dp_dx[2][1] = -1.0;      
}

/**
    Update dpdx matrix of ma_mp
*/
void dpdx_mp_ma(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = 1.0;      dp_dx[0][2] = -1.0;      dp_dx[0][3] = -1.0;      dp_dx[0][4] = -1.0;      
    dp_dx[1][0] = x[1] - 1.0;      dp_dx[1][1] = x[0] - 1.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      dp_dx[1][4] = 0.0;      
    dp_dx[2][0] = 1.0 - x[1];      dp_dx[2][1] = -x[0];      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 0.0;      dp_dx[2][4] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 1.0;      dp_dx[3][4] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = 1.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 1.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      
}

/**
    Update dpdx matrix of mu_mp
*/
void dpdx_mp_mu(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = 1.0;      dp_dx[0][2] = -1.0;      dp_dx[0][3] = -1.0;      dp_dx[0][4] = -1.0;      
    dp_dx[1][0] = x[1] - 1.0;      dp_dx[1][1] = x[0] - 1.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      dp_dx[1][4] = 0.0;      
    dp_dx[2][0] = 1.0 - x[1];      dp_dx[2][1] = -x[0];      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 0.0;      dp_dx[2][4] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 1.0;      dp_dx[3][4] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = 1.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 1.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      
}

/**
    Update dpdx matrix of opx_mp
*/
void dpdx_mp_opx(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[4] + x[1] - 1.0;      dp_dx[0][1] = 0.5*x[5] + x[0] - 1.0;      dp_dx[0][2] = -1.0;      dp_dx[0][3] = -1.0;      dp_dx[0][4] = 0.5*x[5] + x[0] - 1.0;      dp_dx[0][5] = 0.5*x[4] + 0.5*x[1] - 0.5;      
    dp_dx[1][0] = -x[3] - x[1] - x[2] + 1.0;      dp_dx[1][1] = 0.5*x[5] - x[0];      dp_dx[1][2] = -x[0];      dp_dx[1][3] = -x[0];      dp_dx[1][4] = 0.5*x[5];      dp_dx[1][5] = 0.5*x[4] + 0.5*x[1] - 0.5;      
    dp_dx[2][0] = -x[4] + x[3] + x[2];      dp_dx[2][1] = -x[5];      dp_dx[2][2] = x[0];      dp_dx[2][3] = x[0];      dp_dx[2][4] = -x[5] - x[0];      dp_dx[2][5] = -x[4] - x[1] + 1.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      dp_dx[3][5] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 1.0;      dp_dx[4][4] = 0.0;      dp_dx[4][5] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 1.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      dp_dx[5][5] = 0.0;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 0.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 1.0;      dp_dx[6][5] = 0.0;      
}

/**
    Update dpdx matrix of sa_mp
*/
void dpdx_mp_sa(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.0;      dp_dx[0][1] = -1.0;      dp_dx[0][2] = -1.0;      dp_dx[0][3] = -0.25;      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = 1.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      
    dp_dx[2][0] = -x[2] - x[1] + 1.0;      dp_dx[2][1] = -x[0];      dp_dx[2][2] = -x[0];      dp_dx[2][3] = -0.75;      
    dp_dx[3][0] = x[2] + x[1];      dp_dx[3][1] = x[0];      dp_dx[3][2] = x[0];      dp_dx[3][3] = 1.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 1.0;      dp_dx[4][3] = 0.0;      
}

/**
    Update dpdx matrix of cd_mp
*/
void dpdx_mp_cd(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[1] - 1.0;      dp_dx[0][1] = x[0] - 1.0;      dp_dx[0][2] = -1.0;      
    dp_dx[1][0] = 1.0 - x[1];      dp_dx[1][1] = -x[0];      dp_dx[1][2] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 0.0;      dp_dx[2][2] = 1.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 1.0;      dp_dx[3][2] = 0.0;      
}

/**
    Update dpdx matrix of st_mp
*/
void dpdx_mp_st(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[1] - 1.0;      dp_dx[0][1] = x[0] - 1.0;      dp_dx[0][2] = -1.0;      dp_dx[0][3] = -4.0/3.0;      
    dp_dx[1][0] = 1.0 - x[1];      dp_dx[1][1] = -x[0];      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.0;      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.0;      dp_dx[3][3] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 4.0/3.0;      
}

/**
    Update dpdx matrix of chl_mp
*/
void dpdx_mp_chl(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -x[2] + x[3] - x[4] - x[1];      dp_dx[0][1] = 0.25*x[5] + 1.25*x[6] - x[0];      dp_dx[0][2] = 1.25*x[6] - x[0];      dp_dx[0][3] = 0.25*x[5] + x[0] - 1.0;      dp_dx[0][4] = -0.25*x[5] + 1.25*x[6] - x[0] + 2.0;      dp_dx[0][5] = 0.25*x[3] - 0.25*x[4] + 0.25*x[1] - 0.25;      dp_dx[0][6] = 1.25*x[2] + 1.25*x[4] + 1.25*x[1] - 1.25;      
    dp_dx[1][0] = 2.0*x[2] + x[4] + 3.0*x[1] - 2.0;      dp_dx[1][1] = -1.25*x[5] - 2.25*x[6] + 3.0*x[0] - 1.0;      dp_dx[1][2] = -2.25*x[6] + 2.0*x[0] - 1.0;      dp_dx[1][3] = -1.25*x[5];      dp_dx[1][4] = 1.25*x[5] - 2.25*x[6] + x[0] - 1.0;      dp_dx[1][5] = -1.25*x[3] + 1.25*x[4] - 1.25*x[1] + 1.25;      dp_dx[1][6] = -2.25*x[2] - 2.25*x[4] - 2.25*x[1] + 2.25;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.0;      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 0.0;      dp_dx[2][4] = -1.0;      dp_dx[2][5] = 0.0;      dp_dx[2][6] = 0.0;      
    dp_dx[3][0] = x[2] - x[3] + x[4] + x[1];      dp_dx[3][1] = -0.25*x[5] - 1.25*x[6] + x[0];      dp_dx[3][2] = -1.25*x[6] + x[0];      dp_dx[3][3] = -0.25*x[5] - x[0];      dp_dx[3][4] = 0.25*x[5] - 1.25*x[6] + x[0];      dp_dx[3][5] = -0.25*x[3] + 0.25*x[4] - 0.25*x[1] + 0.25;      dp_dx[3][6] = -1.25*x[2] - 1.25*x[4] - 1.25*x[1] + 1.25;      
    dp_dx[4][0] = -x[2] - x[4] - x[1] + 1.0;      dp_dx[4][1] = x[6] - x[0];      dp_dx[4][2] = x[6] - x[0];      dp_dx[4][3] = 0.0;      dp_dx[4][4] = x[6] - x[0];      dp_dx[4][5] = 0.0;      dp_dx[4][6] = x[2] + x[4] + x[1] - 1.0;      
    dp_dx[5][0] = -x[2] - 2.0*x[1] + 1.0;      dp_dx[5][1] = 1.25*x[5] + 1.25*x[6] - 2.0*x[0];      dp_dx[5][2] = 1.25*x[6] - x[0];      dp_dx[5][3] = 1.25*x[5];      dp_dx[5][4] = -1.25*x[5] + 1.25*x[6];      dp_dx[5][5] = 1.25*x[3] - 1.25*x[4] + 1.25*x[1] - 1.25;      dp_dx[5][6] = 1.25*x[2] + 1.25*x[4] + 1.25*x[1] - 1.25;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 1.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = 0.0;      dp_dx[6][6] = 0.0;      
    dp_dx[7][0] = 0.0;      dp_dx[7][1] = 0.0;      dp_dx[7][2] = 0.0;      dp_dx[7][3] = 1.0;      dp_dx[7][4] = 0.0;      dp_dx[7][5] = 0.0;      dp_dx[7][6] = 0.0;      
}
/**
    Update dpdx matrix of ctd_mp
*/
void dpdx_mp_ctd(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[1] - 1.0;      dp_dx[0][1] = x[0] - 1.0;      dp_dx[0][2] = -1.0;      
    dp_dx[1][0] = 1.0 - x[1];      dp_dx[1][1] = -x[0];      dp_dx[1][2] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.0;      dp_dx[2][2] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.0;      
}

/**
    Update dpdx matrix of sp_mp
*/
void dpdx_mp_sp(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[2] + 1.0;      dp_dx[0][1] = 1.0;      dp_dx[0][2] = x[0] - 1.0;      
    dp_dx[1][0] = -x[2] - 1.0;      dp_dx[1][1] = 0.0;      dp_dx[1][2] = 1.0 - x[0];      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = -1.0;      dp_dx[2][2] = -1.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.0;      
}

/**
    Update dpdx matrix of ilm_mp
*/
void dpdx_mp_ilm(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = 0.0;      dp_dx[0][2] = 0.0;      dp_dx[0][3] = 1.0;      
    dp_dx[1][0] = 1.0;      dp_dx[1][1] = -1.0;      dp_dx[1][2] = -1.0;      dp_dx[1][3] = -1.0;      
    dp_dx[2][0] = -1.0;      dp_dx[2][1] = 0.0;      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 1.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 1.0;      dp_dx[4][3] = 0.0;      
}

/**
    Update dpdx matrix of mt_mp
*/
void dpdx_mp_mt(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -2.0;      dp_dx[0][1] = 3.0;      
    dp_dx[1][0] = 3.0;      dp_dx[1][1] = -3.0;      
    dp_dx[2][0] = -1.0;      dp_dx[2][1] = 0.0;      
}

/**
    Endmember fraction of liq_mp
*/
void px_mp_liq(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[0];
        p[1]           = x[1]*x[2];
        p[2]           = x[1]*(1.0 - x[2]);
        p[3]           = x[3];
        p[4]           = -x[3] - x[1] - x[6] - x[4] - x[0] + 1.0;
        p[5]           = x[4]*(1.0 - x[5]);
        p[6]           = x[4]*x[5];
        p[7]           = x[6];
}
    
/**
    Endmember fraction of pl4tr_mp
*/
void px_mp_pl4tr(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[0] - x[1] + 1.0;
        p[1]           = x[0];
        p[2]           = x[1];
}
/**
    Endmember fraction of bi_mp
*/
void px_mp_bi(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[3]*x[0] - x[3] + 3.0*x[1]*x[0] - x[1] - 2.0/3.0*x[5] + x[4]*x[0] - x[4] + x[0]*x[2] - x[0] - x[2] + 1.0;
        p[1]           = -1.0/3.0*x[5] + x[0];
        p[2]           = -x[3]*x[0] - 3.0*x[1]*x[0] + x[5] - x[4]*x[0] - x[0]*x[2];
        p[3]           = x[2];
        p[4]           = x[4];
        p[5]           = x[3];
        p[6]           = x[1];
}
    
/**
    Endmember fraction of g_mp
*/
void px_mp_g(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[3] + x[2]*x[0] - x[2] + x[0]*x[1] - x[0] - x[1] + 1.0;
        p[1]           = -x[2]*x[0] - x[0]*x[1] + x[0];
        p[2]           = x[2];
        p[3]           = x[1];
        p[4]           = x[3];
}

/**
    Endmember fraction of ep_mp
*/
void px_mp_ep(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[0] - x[1] + 1.0;
        p[1]           = 2.0*x[1];
        p[2]           = x[0] - x[1];
}
    
/**
    Endmember fraction of ma_mp
*/
void px_mp_ma(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[4] - x[2] - x[3] + x[1];
        p[1]           = x[0]*x[1] - x[0] - x[1] + 1.0;
        p[2]           = -x[0]*x[1] + x[0];
        p[3]           = x[3];
        p[4]           = x[4];
        p[5]           = x[2];
}
    
/**
    Endmember fraction of mu_mp
*/
void px_mp_mu(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[4] - x[2] - x[3] + x[1];
        p[1]           = x[0]*x[1] - x[0] - x[1] + 1.0;
        p[2]           = -x[0]*x[1] + x[0];
        p[3]           = x[3];
        p[4]           = x[4];
        p[5]           = x[2];
}
    
/**
    Endmember fraction of opx_mp
*/
void px_mp_opx(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 0.5*x[4]*x[5] + x[4]*x[0] - x[4] - x[3] + 0.5*x[1]*x[5] + x[1]*x[0] - x[1] - 0.5*x[5] - x[0] - x[2] + 1.0;
        p[1]           = 0.5*x[4]*x[5] - x[3]*x[0] + 0.5*x[1]*x[5] - x[1]*x[0] - 0.5*x[5] - x[0]*x[2] + x[0];
        p[2]           = -x[4]*x[5] - x[4]*x[0] + x[3]*x[0] - x[1]*x[5] + x[5] + x[0]*x[2];
        p[3]           = x[2];
        p[4]           = x[3];
        p[5]           = x[1];
        p[6]           = x[4];
}
    
/**
    Endmember fraction of sa_mp
*/
void px_mp_sa(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[2] - 0.25*x[3] - x[0] - x[1] + 1.0;
        p[1]           = x[1];
        p[2]           = -x[2]*x[0] - 0.75*x[3] - x[0]*x[1] + x[0];
        p[3]           = x[2]*x[0] + x[3] + x[0]*x[1];
        p[4]           = x[2];
}
    
/**
    Endmember fraction of cd_mp
*/
void px_mp_cd(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[2] + x[1]*x[0] - x[1] - x[0] + 1.0;
        p[1]           = -x[1]*x[0] + x[0];
        p[2]           = x[2];
        p[3]           = x[1];
}
    
/**
    Endmember fraction of st_mp
*/
void px_mp_st(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[2] + x[1]*x[0] - x[1] - 4.0/3.0*x[3] - x[0] + 1.0;
        p[1]           = -x[1]*x[0] + x[0];
        p[2]           = x[1];
        p[3]           = x[2];
        p[4]           = 4.0/3.0*x[3];
}

/**
    Endmember fraction of chl_mp
*/
void px_mp_chl(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 1.25*x[2]*x[6] - x[2]*x[0] + 0.25*x[3]*x[5] + x[3]*x[0] - x[3] - 0.25*x[5]*x[4] + 0.25*x[5]*x[1] - 0.25*x[5] + 1.25*x[6]*x[4] + 1.25*x[6]*x[1] - 1.25*x[6] - x[4]*x[0] + 2.0*x[4] - x[0]*x[1];
        p[1]           = -2.25*x[2]*x[6] + 2.0*x[2]*x[0] - x[2] - 1.25*x[3]*x[5] + 1.25*x[5]*x[4] - 1.25*x[5]*x[1] + 1.25*x[5] - 2.25*x[6]*x[4] - 2.25*x[6]*x[1] + 2.25*x[6] + x[4]*x[0] - x[4] + 3.0*x[0]*x[1] - 2.0*x[0] - x[1] + 1.0;
        p[2]           = -x[4] + x[1];
        p[3]           = -1.25*x[2]*x[6] + x[2]*x[0] - 0.25*x[3]*x[5] - x[3]*x[0] + 0.25*x[5]*x[4] - 0.25*x[5]*x[1] + 0.25*x[5] - 1.25*x[6]*x[4] - 1.25*x[6]*x[1] + 1.25*x[6] + x[4]*x[0] + x[0]*x[1];
        p[4]           = x[2]*x[6] - x[2]*x[0] + x[6]*x[4] + x[6]*x[1] - x[6] - x[4]*x[0] - x[0]*x[1] + x[0];
        p[5]           = 1.25*x[2]*x[6] - x[2]*x[0] + 1.25*x[3]*x[5] - 1.25*x[5]*x[4] + 1.25*x[5]*x[1] - 1.25*x[5] + 1.25*x[6]*x[4] + 1.25*x[6]*x[1] - 1.25*x[6] - 2.0*x[0]*x[1] + x[0];
        p[6]           = x[2];
        p[7]           = x[3];
}
     
/**
    Endmember fraction of ctd_mp
*/
void px_mp_ctd(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[2] + x[1]*x[0] - x[1] - x[0] + 1.0;
        p[1]           = -x[1]*x[0] + x[0];
        p[2]           = x[1];
        p[3]           = x[2];
}

/**
    Endmember fraction of sp_mp
*/
void px_mp_sp(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[1] + (x[0] - 1.0)*(x[2] + 1.0);
        p[1]           = (1.0 - x[0])*(x[2] + 1.0);
        p[2]           = -x[1] - x[2] + 1.0;
        p[3]           = x[2];
}
    
/**
    Endmember fraction of ilm_mp
*/
void px_mp_ilm(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[3];
        p[1]           = -x[1] + x[0] - x[2] - x[3];
        p[2]           = 1.0 - x[0];
        p[3]           = x[1];
        p[4]           = x[2];
}
    
/**
    Endmember fraction of mt_mp
*/
void px_mp_mt(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 3.0*x[1] - 2.0*x[0];
        p[1]           = -3.0*x[1] + 3.0*x[0];
        p[2]           = 1.0 - x[0];
}


/**
    Objective function of liq_mp
*/
double obj_mp_liq(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mp_liq(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0 - x[6];
    sf[1]          = x[0];
    sf[2]          = x[1]*x[2];
    sf[3]          = x[1]*(1.0 - x[2]);
    sf[4]          = x[3];
    sf[5]          = -x[3] - x[1] - x[6] - x[4] - x[0] + 1.0;
    sf[6]          = x[4];
    sf[7]          = x[5];
    sf[8]          = 1.0 - x[5];
    sf[9]          = x[6];
    
    
    mu[0]          = R*T*creal(clog(sf[0]*sf[1])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(sf[0]*sf[2])) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(sf[0]*sf[3])) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(sf[0]*sf[4])) + gb[3] + mu_Gex[3];
    mu[4]          = R*T*creal(clog(sf[0]*sf[5])) + gb[4] + mu_Gex[4];
    mu[5]          = R*T*creal(clog(sf[0]*sf[6]*pow(sf[8], 5.0))) + gb[5] + mu_Gex[5];
    mu[6]          = R*T*creal(clog(sf[0]*sf[6]*pow(sf[7], 5.0))) + gb[6] + mu_Gex[6];
    mu[7]          = R*T*creal(clog(pow(sf[9], 2.0))) + gb[7] + mu_Gex[7];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_liq(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}

/**
    Objective function of pl4tr_mp
*/
double obj_mp_pl4tr(unsigned n, const double *x, double *grad, void *SS_ref_db){

	SS_ref *d  = (SS_ref *) SS_ref_db;

	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;

	px_mp_pl4tr(SS_ref_db,x);

	d->sum_v = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->sum_v += d->p[i]*d->v[i];
	}
	for (int i = 0; i < d->n_em; i++){
		d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
	}

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
				it += 1;
			}
		}
	}
	
    sf[0]           = -x[0] - x[1] + 1.0;
    sf[1]           = x[0];
    sf[2]           = x[1];
    sf[3]           = 0.25*x[0] + 0.25;
    sf[4]           = 0.75 - 0.25*x[0];

	mu[0]          = R*T*creal(clog(1.7548*sf[0]*cpow(sf[3], 0.25)*cpow(sf[4], 0.75))) 	+ gb[0] + mu_Gex[0];
	mu[1]          = R*T*creal(clog(2.0*sf[1]*csqrt(sf[3])*csqrt(sf[4]))) 				+ gb[1] + mu_Gex[1];
	mu[2]          = R*T*creal(clog(1.7548*sf[2]*cpow(sf[3], 0.25)*cpow(sf[4], 0.75))) 	+ gb[2] + mu_Gex[2];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_mp_pl4tr(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
}
    
/**
    Objective function of bi_mp
*/
double obj_mp_bi(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_bi(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[3]*x[0] - x[3] + 3.0*x[1]*x[0] - x[1] - 2./3.*x[5] + x[4]*x[0] - x[4] + x[0]*x[2] - x[0] - x[2] + 1.0;
    sf[1]          = x[1];
    sf[2]          = -x[3]*x[0] - 3.0*x[1]*x[0] + 2./3.*x[5] - x[4]*x[0] - x[0]*x[2] + x[0];
    sf[3]          = x[3];
    sf[4]          = x[4];
    sf[5]          = x[2];
    sf[6]          = -x[1] + 1./3.*x[5] - x[0] + 1.0;
    sf[7]          = x[1];
    sf[8]          = -1./3.*x[5] + x[0];
    sf[9]          = -0.5*x[3] - 0.5*x[2] + 0.5;
    sf[10]          = 0.5*x[3] + 0.5*x[2] + 0.5;
    sf[11]          = 1.0 - x[4];
    sf[12]          = x[4];
    
    
    mu[0]          = R*T*creal(clog(4.0*sf[10]*pow(sf[6], 2.0)*sf[0]*pow(sf[11], 2.0)*sf[9])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(4.0*sf[10]*pow(sf[8], 2.0)*sf[2]*pow(sf[11], 2.0)*sf[9])) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(4.0*sf[10]*sf[2]*pow(sf[6], 2.0)*pow(sf[11], 2.0)*sf[9])) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(sf[5]*pow(sf[10], 2.0)*pow(sf[6], 2.0)*pow(sf[11], 2.0))) + gb[3] + mu_Gex[3];
    mu[4]          = R*T*creal(clog(4.0*sf[10]*pow(sf[6], 2.0)*pow(sf[12], 2.0)*sf[9]*sf[4])) + gb[4] + mu_Gex[4];
    mu[5]          = R*T*creal(clog(pow(sf[10], 2.0)*sf[3]*pow(sf[6], 2.0)*pow(sf[11], 2.0) + d_em[5])) + gb[5] + mu_Gex[5];
    mu[6]          = R*T*creal(clog(4.0*sf[10]*pow(sf[7], 2.0)*sf[1]*pow(sf[11], 2.0)*sf[9])) + gb[6] + mu_Gex[6];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_bi(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    

/**
    Objective function of g_mp
*/
double obj_mp_g(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_g(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          = x[2]*x[0] - x[2] + x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[1]          = -x[2]*x[0] - x[0]*x[1] + x[0];
    sf[2]          = x[2];
    sf[3]          = x[1];
    sf[4]          = 1.0 - x[3];
    sf[5]          = x[3];
    
    mu[0]          = R*T*creal(clog(pow(sf[4], 2.0)*pow(sf[0], 3.0))) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(pow(sf[4], 2.0)*pow(sf[1], 3.0))) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(pow(sf[4], 2.0)*pow(sf[2], 3.0))) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(pow(sf[4], 2.0)*pow(sf[3], 3.0))) + gb[3] + mu_Gex[3];
    mu[4]          = R*T*creal(clog(pow(sf[5], 2.0)*pow(sf[0], 3.0)+ d_em[4])) + gb[4] + mu_Gex[4];
    

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_g(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
        // printf("grad \n");
        // for (int i = 0; i < (d->n_xeos); i++){
        //     printf(" %+5f", grad[i]);
        // }
    }

    // printf("d->df %+10f d->df_raw %+10f d->factor %+10f \n",d->df,d->df_raw,d->factor);

    // printf("\n");
    // printf(" %6s %6s %6s %6s %6s %6s\n","GB","MAT_PHI","MU_GEX","SF","MU","P");      
    // for (int i = 0; i < d->n_em; i++){
    //     printf(" %+5f %+5f %+5f %+5f %+5f %+5f\n",gb[i],mat_phi[i],mu_Gex[i],sf[i],mu[i],d->p[i]);    
    // }


    return d->df;
}
    
/**
    Objective function of ep_mp
*/
double obj_mp_ep(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_ep(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[0] - x[1];
    sf[1]          = -x[0] + x[1] + 1.0;
    sf[2]          = x[0] + x[1];
    sf[3]          = -x[0] - x[1] + 1.0;
    
    
    mu[0]          = R*T*creal(clog(sf[1]*sf[3])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(sf[1]*sf[2] + d_em[1])) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(sf[0]*sf[2] + d_em[2])) + gb[2] + mu_Gex[2];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_ep(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of ma_mp
*/
double obj_mp_ma(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_ma(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          = -x[4] - x[3] + 1.0;
    sf[1]          = x[3];
    sf[2]          = x[4];
    sf[3]          = x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[4]          = -x[0]*x[1] + x[0];
    sf[5]          = x[1];
    sf[6]          = 1.0 - x[2];
    sf[7]          = x[2];
    sf[8]          = -0.5*x[4] - 0.5*x[1] + 1.0;
    sf[9]          = 0.5*x[4] + 0.5*x[1];
    
    
    mu[0]          = R*T*creal(clog(4.0*sf[5]*sf[6]*sf[9]*sf[0]*sf[8])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(sf[6]*sf[0]*sf[3]*cpow(sf[8], 2.0))) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(sf[6]*sf[4]*sf[0]*cpow(sf[8], 2.0))) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(4.0*sf[5]*sf[6]*sf[9]*sf[1]*sf[8])) + gb[3] + mu_Gex[3];
    mu[4]          = R*T*creal(clog(sf[5]*sf[6]*cpow(sf[9], 2.0)*sf[2])) + gb[4] + mu_Gex[4];
    mu[5]          = R*T*creal(clog(4.0*sf[5]*sf[9]*sf[7]*sf[0]*sf[8] + d_em[5])) + gb[5] + mu_Gex[5];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_ma(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of mu_mp
*/
double obj_mp_mu(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_mu(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          = -x[4] - x[3] + 1.0;
    sf[1]          = x[3];
    sf[2]          = x[4];
    sf[3]          = x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[4]          = -x[0]*x[1] + x[0];
    sf[5]          = x[1];
    sf[6]          = 1.0 - x[2];
    sf[7]          = x[2];
    sf[8]          = -0.5*x[4] - 0.5*x[1] + 1.0;
    sf[9]          = 0.5*x[4] + 0.5*x[1];
    
    
    mu[0]          = R*T*creal(clog(4.0*sf[5]*sf[6]*sf[9]*sf[0]*sf[8])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(sf[6]*sf[0]*sf[3]*pow(sf[8], 2.0))) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(sf[6]*sf[4]*sf[0]*pow(sf[8], 2.0))) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(4.0*sf[5]*sf[6]*sf[9]*sf[1]*sf[8])) + gb[3] + mu_Gex[3];
    mu[4]          = R*T*creal(clog(sf[5]*sf[6]*pow(sf[9], 2.0)*sf[2])) + gb[4] + mu_Gex[4];
    mu[5]          = R*T*creal(clog(4.0*sf[5]*sf[9]*sf[7]*sf[0]*sf[8] + d_em[5])) + gb[5] + mu_Gex[5];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_mu(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of opx_mp
*/
double obj_mp_opx(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_opx(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          = -0.5*x[4]*x[5] + x[3]*x[0] - x[3] - 0.5*x[1]*x[5] + x[1]*x[0] - x[1] + 0.5*x[5] + x[0]*x[2] - x[0] - x[2] + 1.0;
    sf[1]          = 0.5*x[4]*x[5] - x[3]*x[0] + 0.5*x[1]*x[5] - x[1]*x[0] - 0.5*x[5] - x[0]*x[2] + x[0];
    sf[2]          = x[1];
    sf[3]          = x[3];
    sf[4]          = x[2];
    sf[5]          = 0.5*x[4]*x[5] + x[4]*x[0] - x[4] + 0.5*x[1]*x[5] + x[1]*x[0] - x[1] - 0.5*x[5] - x[0] + 1.0;
    sf[6]          = -0.5*x[4]*x[5] - x[4]*x[0] - 0.5*x[1]*x[5] - x[1]*x[0] + 0.5*x[5] + x[0];
    sf[7]          = x[1];
    sf[8]          = x[4];
    sf[9]          = 0.5*x[3] + 0.5*x[2];
    sf[10]         = -0.5*x[3] - 0.5*x[2] + 1.0;
    
    
    mu[0]          = R*T*creal(clog(sf[0]*sf[5]*csqrt(sf[10]))) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(sf[1]*sf[6]*csqrt(sf[10]))) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(sf[6]*sf[0]*csqrt(sf[10]))) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(1.4142*sf[4]*cpow(sf[9], 0.25)*sf[5]*cpow(sf[10], 0.25))) + gb[3] + mu_Gex[3];
    mu[4]          = R*T*creal(clog(1.4142*cpow(sf[9], 0.25)*sf[3]*sf[5]*cpow(sf[10], 0.25) + d_em[4])) + gb[4] + mu_Gex[4];
    mu[5]          = R*T*creal(clog(sf[2]*sf[7]*csqrt(sf[10]))) + gb[5] + mu_Gex[5];
    mu[6]          = R*T*creal(clog(sf[8]*sf[0]*csqrt(sf[10]))) + gb[6] + mu_Gex[6];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_opx(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of sa_mp
*/
double obj_mp_sa(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_sa(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[2]*x[0] - x[2] + 0.75*x[3] + x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[1]          = -x[2]*x[0] - 0.75*x[3] - x[0]*x[1] + x[0];
    sf[2]          = x[2];
    sf[3]          = x[1];
    sf[4]          = -0.25*x[3] - x[0] + 1.0;
    sf[5]          = 0.25*x[3] + x[0];
    sf[6]          = -x[2] - x[1] + 1.0;
    sf[7]          = x[2] + x[1];
    
    
    mu[0]          = R*T*creal(clog(sf[0]*pow(sf[4], 3.0)*sf[6])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(sf[3]*sf[7]*pow(sf[4], 3.0))) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(sf[1]*pow(sf[5], 3.0)*sf[6])) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(pow(sf[5], 3.0)*sf[0]*sf[6])) + gb[3] + mu_Gex[3];
    mu[4]          = R*T*creal(clog(sf[7]*sf[2]*pow(sf[4], 3.0) + d_em[4])) + gb[4] + mu_Gex[4];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_sa(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of cd_mp
*/
double obj_mp_cd(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_mp_cd(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[0];
    sf[1]          = 1.0 - x[0];
    sf[2]          = x[1];
    sf[3]          = x[2];
    sf[4]          = 1.0 - x[2];
    
    
    mu[0]          = R*T*creal(clog(pow(sf[1], 2.0)*sf[4])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(pow(sf[0], 2.0)*sf[4])) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(sf[3]*pow(sf[1], 2.0))) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(pow(sf[2], 2.0)*sf[4])) + gb[3] + mu_Gex[3];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_cd(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of st_mp
*/
double obj_mp_st(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_st(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[1]*x[0] - x[1] - x[0] + 1.0;
    sf[1]          = -x[1]*x[0] + x[0];
    sf[2]          = x[1];
    sf[3]          = -x[2] - 1.33333333333333*x[3] + 1.0;
    sf[4]          = x[2];
    sf[5]          = x[3];
    sf[6]          = 1./3.*x[3];
    
    
    mu[0]          = R*T*creal(clog(pow(sf[3], 2.0)*pow(sf[0], 4.0))) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(pow(sf[3], 2.0)*pow(sf[1], 4.0))) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(pow(sf[3], 2.0)*pow(sf[2], 4.0))) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(pow(sf[4], 2.0)*pow(sf[0], 4.0) + d_em[3])) + gb[3] + mu_Gex[3];
    mu[4]          = R*T*creal(clog(3.0792*pow(sf[0], 4.0)*pow(sf[5], 1.5)*csqrt(sf[6]))) + gb[4] + mu_Gex[4];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_st(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of chl_mp
*/
double obj_mp_chl(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_chl(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = -x[3]*x[5] + x[3]*x[0] - x[3] + x[5]*x[4] - x[5]*x[1] + x[5] - x[4]*x[0] + x[4] + x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[1]          = x[3]*x[5] - x[3]*x[0] - x[5]*x[4] + x[5]*x[1] - x[5] + x[4]*x[0] - x[0]*x[1] + x[0];
    sf[2]          = -x[4] + x[1];
    sf[3]          = 0.25*x[2]*x[6] + 0.25*x[3]*x[5] + x[3]*x[0] - x[3] - 0.25*x[5]*x[4] + 0.25*x[5]*x[1] - 0.25*x[5] + 0.25*x[6]*x[4] + 0.25*x[6]*x[1] - 0.25*x[6] - x[0] + 1.0;
    sf[4]          = x[3];
    sf[5]          = -0.25*x[2]*x[6] - 0.25*x[3]*x[5] - x[3]*x[0] + 0.25*x[5]*x[4] - 0.25*x[5]*x[1] + 0.25*x[5] - 0.25*x[6]*x[4] - 0.25*x[6]*x[1] + 0.25*x[6] + x[0];
    sf[6]          = -x[2]*x[6] + x[2]*x[0] - x[2] - x[6]*x[4] - x[6]*x[1] + x[6] + x[4]*x[0] - x[4] + x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[7]          = x[2]*x[6] - x[2]*x[0] + x[6]*x[4] + x[6]*x[1] - x[6] - x[4]*x[0] - x[0]*x[1] + x[0];
    sf[8]          = x[2];
    sf[9]          = x[4] + x[1];
    sf[10]          = -0.5*x[2] - x[1] + 1.0;
    sf[11]          = 0.5*x[2] + x[1];
    
    
    mu[0]          = R*T*creal(clog(4.0*sf[9]*sf[11]*sf[0]*pow(sf[3], 4.0)*sf[10])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(sf[0]*pow(sf[3], 4.0)*sf[6]*pow(sf[10], 2.0))) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(sf[2]*sf[9]*pow(sf[11], 2.0)*pow(sf[3], 4.0))) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(4.0*sf[9]*sf[11]*sf[1]*pow(sf[5], 4.0)*sf[10])) + gb[3] + mu_Gex[3];
    mu[4]          = R*T*creal(clog(pow(sf[5], 4.0)*sf[7]*sf[0]*pow(sf[10], 2.0))) + gb[4] + mu_Gex[4];
    mu[5]          = R*T*creal(clog(sf[1]*pow(sf[3], 4.0)*sf[6]*pow(sf[10], 2.0))) + gb[5] + mu_Gex[5];
    mu[6]          = R*T*creal(clog(4.0*sf[11]*sf[8]*sf[0]*pow(sf[3], 4.0)*sf[10] + d_em[6])) + gb[6] + mu_Gex[6];
    mu[7]          = R*T*creal(clog(4.0*sf[9]*sf[11]*pow(sf[4], 5.0)*sf[10])) + gb[7] + mu_Gex[7];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_chl(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of ctd_mp
*/
double obj_mp_ctd(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_ctd(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0 - x[2];
    sf[1]          = x[2];
    sf[2]          = -x[1]*x[0] + x[0];
    sf[3]          = x[1]*x[0] - x[1] - x[0] + 1.0;
    sf[4]          = x[1];
    
    
    mu[0]          = R*T*creal(clog(sqrt(sf[0])*sf[3])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(sqrt(sf[0])*sf[2])) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(sqrt(sf[0])*sf[4])) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(sqrt(sf[1])*sf[3] + d_em[3])) + gb[3] + mu_Gex[3];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_ctd(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of sp_mp
*/
double obj_mp_sp(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_sp(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[1];
    sf[1]          = -x[1] - x[2] + 1.0;
    sf[2]          = x[2];
    sf[3]          = 1.0 - x[0];
    sf[4]          = x[0];
    
    
    mu[0]          = R*T*creal(clog(sf[0]*sf[4])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(sf[0]*sf[3])) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(sf[4]*sf[1] + d_em[2])) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(sf[4]*sf[2])) + gb[3] + mu_Gex[3];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_sp(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of ilm_mp
*/
double obj_mp_ilm(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_ilm(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = -0.5*x[2] + 0.5*x[3] + 0.5*x[0] - 0.5*x[1];
    sf[1]          = -0.5*x[2] - 0.5*x[3] + 0.5*x[0] - 0.5*x[1];
    sf[2]          = x[1];
    sf[3]          = x[2];
    sf[4]          = 1.0 - x[0];
    sf[5]          = 0.5*x[2] + 0.5*x[3] + 0.5*x[0] + 0.5*x[1];
    
    
    mu[0]          = R*T*creal(clog(sf[0]*sf[5])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(4.0*csqrt(sf[0])*sf[1]*csqrt(sf[5]))) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(pow(sf[4], 2.0) + d_em[2])) + gb[2] + mu_Gex[2];
    mu[3]          = R*T*creal(clog(sf[2]*sf[5])) + gb[3] + mu_Gex[3];
    mu[4]          = R*T*creal(clog(sf[3]*sf[5])) + gb[4] + mu_Gex[4];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_ilm(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
    
/**
    Objective function of mt_mp
*/
double obj_mp_mt(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_mp_mt(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 0.5 - 0.5*x[0];
    sf[1]          = -0.5*x[1] + x[0];
    sf[2]          = 0.5*x[1] - 0.5*x[0] + 0.5;
    sf[3]          = x[1];
    sf[4]          = 1.0 - x[1];
    
    
    mu[0]          = R*T*creal(clog(4.0*sf[1]*sf[3]*sf[2] + d_em[0])) + gb[0] + mu_Gex[0];
    mu[1]          = R*T*creal(clog(6.75*pow(sf[1], 4.0/3.0)*cpow(sf[3], 2.0/3.0)*cpow(sf[2], 2.0/3.0)*cpow(sf[4], 1.0/3.0) + d_em[1])) + gb[1] + mu_Gex[1];
    mu[2]          = R*T*creal(clog(4.0*sf[2]*sf[4]*sf[0])) + gb[2] + mu_Gex[2];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_mp_mt(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
/**************************************************************************************/
/**************************************************************************************/
/*********************IGNEOUS DATABASE (Holland et al., 2018)**************************/
/**************************************************************************************/
/**************************************************************************************/
/**
    Endmember to xeos for fper_S11
*/
void p2x_ig_fper(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = d->p[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

void p2x_ig_bi(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	d->iguess[0]  = (d->p[0]-2.0*d->p[1]+d->p[5]+d->p[4]+d->p[3] -1.0)/(d->p[3]+d->p[4]+d->p[5]-3.0);
	d->iguess[1]  = d->p[3];
	d->iguess[2]  = d->p[5];
	d->iguess[3]  = d->p[4];
	d->iguess[4]  = 3.0*( (d->p[0]-2.0*d->p[1]+d->p[5]+d->p[4]+d->p[3] -1.0)/(d->p[3]+d->p[4]+d->p[5]-3.0) -d->p[1]);


	if (d->z_em[4]  == 0.0){ d->iguess[3]  = eps;}
	if (d->z_em[5]  == 0.0){ d->iguess[2]  = eps;}
		
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (cordierite)
*/
void p2x_ig_cd(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;
	
	d->iguess[0]  = d->p[1];
	d->iguess[1]  = d->p[2];
	
	if (d->z_em[2]  == 0.0){ d->iguess[1]  = eps;}
		
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (clinopyroxene)
*/
void p2x_ig_cpx(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	d->iguess[0]  = (2.0*d->p[1] + d->p[8])/(d->p[1] - d->p[2] - d->p[3] - d->p[4] - 0.5*d->p[5] - d->p[6] + d->p[7] + d->p[8] - d->p[9] + 1.0);	
	d->iguess[1]  = d->p[2] + d->p[3] + d->p[4] + d->p[5];
	d->iguess[2]  = d->p[1]+d->p[7]+d->p[8];
	d->iguess[3]  = d->p[6];
	d->iguess[4]  = (d->p[7] + ((2.0*d->p[1] + d->p[8])/(d->p[1] - d->p[2] - d->p[3] - d->p[4] - 0.5*d->p[5] - d->p[6] + d->p[7] + d->p[8] - d->p[9] + 1.0) - 1.0)*(d->p[1] + d->p[7] + d->p[8]))/(-d->p[2] - d->p[3] - d->p[4] - 0.5*d->p[5] - d->p[6] - d->p[9] + 1.0);
	d->iguess[5]  = d->p[4];
	d->iguess[6]  = d->p[3];
	d->iguess[7]  = d->p[5]/2.0;	
	d->iguess[8]  = d->p[9];	

	if (d->z_em[3]  == 0.0){ d->iguess[6]  = eps;}
	if (d->z_em[4]  == 0.0){ d->iguess[5]  = eps;}
	if (d->z_em[5]  == 0.0){ d->iguess[7]  = eps;}
		
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (epidote)
*/
void p2x_ig_ep(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	d->iguess[0]  = d->p[2] + d->p[1]/2.0;
	d->iguess[1]  = d->p[1]/2.0;

	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (fluid)
*/
void p2x_ig_fl(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;
	
	d->iguess[0]  = d->p[2];
	d->iguess[1]  = d->p[1];
	d->iguess[2]  = d->p[3];
	d->iguess[3]  = d->p[4];
	d->iguess[4]  = d->p[5];
	d->iguess[5]  = d->p[6];
	d->iguess[6]  = d->p[7];
	d->iguess[7]  = d->p[8];
	d->iguess[8]  = d->p[9];
	d->iguess[9]  = d->p[10];
	
	if (d->z_em[10] == 0.0){ d->iguess[9]  = eps;}
	if (d->z_em[7]  == 0.0){ d->iguess[6]  = eps;}
	if (d->z_em[6]  == 0.0){ d->iguess[5]  = eps;}
	if (d->z_em[8]  == 0.0){ d->iguess[7]  = eps;}
		
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (garnet)
*/
void p2x_ig_g(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	d->iguess[0] = (d->p[0]+d->p[4]+d->p[5]+d->p[2]+d->p[3]-1.0)/(-1.0+d->p[2]+d->p[3]);
	d->iguess[1] = d->p[2]+d->p[3];
	d->iguess[2] = d->p[3];
	d->iguess[3] = d->p[4];
	d->iguess[4] = d->p[5]/4.0;

	if (d->z_em[3]  == 0.0){ d->iguess[2]  = eps;}
	if (d->z_em[4]  == 0.0){ d->iguess[3]  = eps;}
	if (d->z_em[5]  == 0.0){ d->iguess[4]  = eps;}
		
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (hornblende)
*/
void p2x_ig_hb(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	d->iguess[0] = (-3.5*d->p[5] - 2.0*d->p[6] - 2.5*d->p[7])/(-0.5*d->p[0] + 0.5*d->p[1] + 0.5*d->p[10] + 0.5*d->p[3] - 1.5*d->p[4] - 1.5*d->p[5] - 1.5*d->p[6] - 1.5*d->p[7] + 0.5*d->p[8] - 2.0);
	d->iguess[1] = (d->p[1]-d->p[0] + 1.0-d->p[3]-d->p[8]-d->p[4]-d->p[6]-d->p[5]-d->p[7] -2*d->p[8] - d->p[10] + 2*(d->p[3] + d->p[8]))/2.0;
	d->iguess[2] = d->p[3] + d->p[8];
	d->iguess[3] = d->p[2] + d->p[9];
	d->iguess[4] = d->p[9]/(d->p[2]+d->p[9]);
	d->iguess[5] = 1.0-d->p[3]-d->p[8]-d->p[4]-d->p[6]-d->p[5]-d->p[7];
	d->iguess[6] = d->p[8];
	d->iguess[7] = d->p[10];
	d->iguess[8] = (-3.5*d->p[5] - 2.0*d->p[6] - 2.5*d->p[7])/(-0.5*d->p[0] + 0.5*d->p[1] + 0.5*d->p[10] + 0.5*d->p[3] - 1.5*d->p[4] - 1.5*d->p[5] - 1.5*d->p[6] - 1.5*d->p[7] + 0.5*d->p[8] - 2.0) -d->p[5] -d->p[7];
	d->iguess[9] = (d->p[5] + d->p[6] - (-3.5*d->p[5] - 2.0*d->p[6] - 2.5*d->p[7])*(0.5*d->p[0] - 0.5*d->p[1] - 0.5*d->p[10] - 0.5*d->p[3] + 0.5*d->p[4] + 0.5*d->p[5] + 0.5*d->p[6] + 0.5*d->p[7] - 0.5*d->p[8] + 0.5)/(-0.5*d->p[0] + 0.5*d->p[1] + 0.5*d->p[10] + 0.5*d->p[3] - 1.5*d->p[4] - 1.5*d->p[5] - 1.5*d->p[6] - 1.5*d->p[7] + 0.5*d->p[8] - 2.0))/(-0.5*d->p[0] + 0.5*d->p[1] + 0.5*d->p[10] + 0.5*d->p[3] - 0.5*d->p[4] - 0.5*d->p[5] - 0.5*d->p[6] - 0.5*d->p[7] + 0.5*d->p[8] - 0.5);

	if (d->z_em[8]  == 0){ d->iguess[6]  = eps;}
	if (d->z_em[10] == 0){ d->iguess[7]  = eps;}
		
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (ilm)
*/
/** DEPRECATED */
void p2x_ig_ilm(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	d->iguess[0] = d->p[1]+d->p[0];
	d->iguess[1] = d->p[0];
		
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}


/** 
  endmembers to xeos (liquid)
*/
void p2x_ig_liq(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;
		
	d->iguess[0]  = (d->p[2]+d->p[10])/(1.0+3./4.*d->p[10]);	
	d->iguess[1]  = (d->p[1]+d->p[10])/(1.0+3./4.*d->p[10]);
	d->iguess[2]  = d->p[3]/(1.0+3./4.*d->p[10]);
	d->iguess[3]  = d->p[4]/(1.0+3./4.*d->p[10]);
	d->iguess[4]  = d->p[5]/(1.0+3./4.*d->p[10]);
	d->iguess[5]  = d->p[6]/(1.0+3./4.*d->p[10]);
	d->iguess[6]  = d->p[7]/(1.0+3./4.*d->p[10]);
	d->iguess[7]  = d->p[8]/(1.0+3./4.*d->p[10]);
	d->iguess[8]  = d->p[9]/(1.0+3./4.*d->p[10]);
	d->iguess[9]  = d->p[10];
	d->iguess[10] = d->p[11]/(1.0+3./4.*d->p[10]);
		
	if (d->z_em[11] == 0.0){ d->iguess[10] = eps;}
	if (d->z_em[8]  == 0.0){ d->iguess[7]  = eps;}
	if (d->z_em[7]  == 0.0){ d->iguess[6]  = eps;}
	if (d->z_em[6]  == 0.0){ d->iguess[5]  = eps;}
		
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (olivine)
*/
void p2x_ig_ol(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	d->iguess[0]  = (2.0*d->p[1]+d->p[3])/(2.0-d->p[0]);
	d->iguess[1]  = d->p[0];
	d->iguess[2]  = -d->p[0] - d->p[2] + 1.0 + (d->p[0] - 1.0)*(2.0*d->p[1] + d->p[3])/(2.0 - d->p[0]);
	
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (orthopyroxene)
*/
void p2x_ig_opx(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	d->iguess[0] = (2.0*d->p[1] + d->p[2])/(d->p[0] + d->p[1] + d->p[2] + 0.5*d->p[6] - d->p[8] + 1.0);
	d->iguess[1] = 1.0 - d->p[3] - d->p[8] - d->p[0] - d->p[1] - d->p[2];
	d->iguess[2] = d->p[3];
	d->iguess[3] = (d->p[1] + d->p[2] + (2.0*d->p[1] + d->p[2])*(d->p[3] + d->p[8] - 1.0)/(d->p[0] + d->p[1] + d->p[2] + 0.5*d->p[6] - d->p[8] + 1.0))/(-d->p[0] - d->p[1] - d->p[2] - d->p[3] - 0.5*d->p[6]);
	d->iguess[4] = d->p[7];
	d->iguess[5] = d->p[6]/2.0;
	d->iguess[6] = d->p[5];
	d->iguess[7] = d->p[8];
	
	if (d->z_em[5]  == 0.0){ d->iguess[6]  = eps;}
	if (d->z_em[4]  == 0.0){ d->iguess[4]  = eps;}
	if (d->z_em[6]  == 0.0){ d->iguess[5]  = eps;}
		
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (plagioclase)
*/
void p2x_ig_pl4T(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;
	
	d->iguess[0] = d->p[1];
	d->iguess[1] = d->p[2];
		
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/** 
  endmembers to xeos (spinel)
*/
void p2x_ig_spn(void *SS_ref_db, double eps){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	d->iguess[0]  = (1.0 - d->p[6] - d->p[7] - d->p[0] - d->p[1])/(d->p[7] + 1.0);
	d->iguess[1]  = (d->p[4] + d->p[5])/(1.0 - d->p[6] - d->p[7]);
	d->iguess[2]  = d->p[6];
	d->iguess[3]  = d->p[7];
	d->iguess[4]  = 3./2.*d->p[0] - 1./2. + 3./2.*d->p[6] + d->p[7] + ((1.0 - d->p[6] - d->p[7] - d->p[0] - d->p[1])/(d->p[7] + 1.0))/2.*(1.0+d->p[7]);
	d->iguess[5]  = ((1.0 - d->p[6] - d->p[7] - d->p[0] - d->p[1])/(d->p[7] + 1.0))*(d->p[7] + 1.0) - 3./2.*d->p[3] - 3./2.*d->p[5];
	d->iguess[6]  = -3./2.*d->p[4] + ((d->p[4] + d->p[5])/(1.0 - d->p[6] - d->p[7]))*(1./2. -1./2.*d->p[6] - 1./2.*d->p[7]);

	if (d->z_em[6]  == 0.0){ d->iguess[2]  = eps;}
	if (d->z_em[7]  == 0.0){ d->iguess[3]  = eps;}
	if (d->z_em[4]  == 0.0){ d->iguess[6]  = eps;}
	if (d->z_em[5]  == 0.0){ d->iguess[6]  = eps;}
	if (d->z_em[4]  == 0.0){ d->iguess[1]  = eps;}
	if (d->z_em[5]  == 0.0){ d->iguess[1]  = eps;}
	
	for (int i = 0; i < d->n_xeos; i++){
		if (d->iguess[i] < d->bounds[i][0]){
			d->iguess[i] = d->bounds[i][0];
		}
		if (d->iguess[i] > d->bounds[i][1]){
			d->iguess[i] = d->bounds[i][1];
		}
	}
}

/**
    Endmember to xeos for bi_mp
*/
void p2x_mp_bi(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[1]  = d->p[6];
    d->iguess[3]  = d->p[5];
    d->iguess[4]  = d->p[4];
    d->iguess[2]  = d->p[3];
    d->iguess[5]  = 3.0*(-d->iguess[3]*d->p[1] + d->iguess[3] - 3.0*d->iguess[1]*d->p[1] + d->iguess[1] + d->p[0] - d->p[1]*d->iguess[4] - d->p[1]*d->iguess[2] + d->p[1] + d->iguess[4] + d->iguess[2] - 1.0)/(d->iguess[3] + 3.0*d->iguess[1] + d->iguess[4] + d->iguess[2] - 3.0);
    d->iguess[0]  = (-d->p[2] + d->iguess[5])/(d->iguess[3] + 3.0*d->iguess[1] + d->iguess[4] + d->iguess[2]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for cd_mp
*/
void p2x_mp_cd(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[1]  = d->p[3];
    d->iguess[2]  = d->p[2];
    d->iguess[0] = d->p[1]/(1.0 - d->iguess[1]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for chl_mp
*/
void p2x_mp_chl(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]   =  d->p[6];
    d->iguess[3]   =  d->p[7];
    d->iguess[4] = (d->iguess[3] + 1.0 - d->iguess[2] -(d->p[1] - d->p[3]+d->p[5]-d->p[0] + d->p[2]+d->p[4]))/4.0;
    d->iguess[1]   =  d->p[2] + d->iguess[4];
    d->iguess[0]  = (-2.0*d->iguess[4] + d->iguess[3] + d->p[0] - 4.0*d->p[3] - 5.0*d->p[4] - d->p[5])/(d->iguess[2] + 5.0*d->iguess[3] + 2.0*d->iguess[1] - 6.0);
    d->iguess[6]  = (-2.0*pow(d->iguess[4],2.0) - 2.0*d->iguess[4]*d->iguess[2] + d->iguess[4]*d->iguess[3] + d->iguess[4]*d->p[0] - 4.0*d->iguess[4]*d->p[3] - 5.0*d->iguess[4]*d->p[4] - d->iguess[4]*d->p[5] - 2.0*d->iguess[4]*d->iguess[1] + 2.0*d->iguess[4] + d->iguess[2]*d->iguess[3] + d->iguess[2]*d->p[0] - 4.0*d->iguess[2]*d->p[3] - 4.0*d->iguess[2]*d->p[4] - d->iguess[2]*d->p[5] + 5.0*d->iguess[3]*d->p[4] + d->iguess[3]*d->iguess[1] - d->iguess[3] + d->p[0]*d->iguess[1] - d->p[0] - 4.0*d->p[3]*d->iguess[1] + 4.0*d->p[3] - 3.0*d->p[4]*d->iguess[1] - d->p[4] - d->p[5]*d->iguess[1] + d->p[5])/(d->iguess[4]*d->iguess[2] + 5.0*d->iguess[4]*d->iguess[3] + 2.0*d->iguess[4]*d->iguess[1] - 6.0*d->iguess[4] + pow(d->iguess[2],2.0) + 5.0*d->iguess[2]*d->iguess[3] + 3.0*d->iguess[2]*d->iguess[1] - 7.0*d->iguess[2] + 5.0*d->iguess[3]*d->iguess[1] - 5.0*d->iguess[3] + 2.0*pow(d->iguess[1],2.0) - 8.0*d->iguess[1] + 6.0);
    d->iguess[5]  = (10.0*pow(d->iguess[4],2.0) - 2.0*d->iguess[4]*d->iguess[2] - 25.0*d->iguess[4]*d->iguess[3] - 5.0*d->iguess[4]*d->p[0] + 20.0*d->iguess[4]*d->p[3] + 25.0*d->iguess[4]*d->p[4] + 5.0*d->iguess[4]*d->p[5] - 14.0*d->iguess[4]*d->iguess[1] + 22.0*d->iguess[4] - 4.0*pow(d->iguess[2],2.0) - 21.0*d->iguess[2]*d->iguess[3] - d->iguess[2]*d->p[0] - 4.0*d->iguess[2]*d->p[1] + 4.0*d->iguess[2]*d->p[3] - 4.0*d->iguess[2]*d->p[4] + d->iguess[2]*d->p[5] - 12.0*d->iguess[2]*d->iguess[1] + 28.0*d->iguess[2] - 20.0*d->iguess[3]*d->p[1] - 45.0*d->iguess[3]*d->p[4] - 17.0*d->iguess[3]*d->iguess[1] + 21.0*d->iguess[3] + 3.0*d->p[0]*d->iguess[1] + d->p[0] - 8.0*d->p[1]*d->iguess[1] + 24.0*d->p[1] - 12.0*d->p[3]*d->iguess[1] - 4.0*d->p[3] - 33.0*d->p[4]*d->iguess[1] + 49.0*d->p[4] - 3.0*d->p[5]*d->iguess[1] - d->p[5] - 8.0*pow(d->iguess[1],2.0) + 32.0*d->iguess[1] - 24.0)/(5.0*(-d->iguess[4]*d->iguess[2] - 5.0*d->iguess[4]*d->iguess[3] - 2.0*d->iguess[4]*d->iguess[1] + 6.0*d->iguess[4] + d->iguess[2]*d->iguess[3] + d->iguess[2]*d->iguess[1] - d->iguess[2] + 5.0*pow(d->iguess[3],2.0) + 7.0*d->iguess[3]*d->iguess[1] - 11.0*d->iguess[3] + 2.0*pow(d->iguess[1],2.0) - 8.0*d->iguess[1] + 6.0));
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for ctd_mp
*/
void p2x_mp_ctd(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]  = d->p[3];
    d->iguess[1]  = d->p[2];
    d->iguess[0] = d->p[1]/(1.0 - d->iguess[1]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for ep_mp
*/
void p2x_mp_ep(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[1]   =  d->p[1]/2.0;
    d->iguess[0]   =  d->p[2] + d->iguess[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for g_mp
*/
void p2x_mp_g(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[3]   =  d->p[4];
    d->iguess[1]   =  d->p[3];
    d->iguess[2]   =  d->p[2];
    d->iguess[0]  =  d->p[1]/(1.0 - d->iguess[2] - d->iguess[1]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}
  
/**
    Endmember to xeos for ilm_mp
*/
void p2x_mp_ilm(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]  = d->p[4];
    d->iguess[1]  = d->p[3];
    d->iguess[3]  = d->p[0];
    d->iguess[0]  = 1.0 - d->p[2];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}
 
/**
    Endmember to xeos for liq_mp
*/
void p2x_mp_liq(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]   = d->p[0];
    d->iguess[3]  = d->p[3];
    d->iguess[6] = d->p[7];
    d->iguess[1] = d->p[1] + d->p[2];
    d->iguess[2]  = d->p[1]/d->iguess[1];
    d->iguess[4]  = 1.0 - d->iguess[0] - d->iguess[1] - d->iguess[3] - d->iguess[6] - d->p[4];
    d->iguess[5]   = d->p[6]/d->iguess[4];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}
 
/**
    Endmember to xeos for ma_mp
*/
void p2x_mp_ma(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]  = d->p[5];
    d->iguess[4]  = d->p[4];
    d->iguess[3]  = d->p[3];
    d->iguess[1]  = d->p[0] + d->iguess[4] + d->iguess[3] + d->iguess[2];
    d->iguess[0] = d->p[2]/(1.0-d->iguess[1]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}
 
/**
    Endmember to xeos for mt_mp
*/
void p2x_mp_mt(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[0]  = 1.0 - d->p[2];
    d->iguess[1]   = (3.0*d->iguess[0] - d->p[1])/3.0;
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}
  
/**
    Endmember to xeos for mu_mp
*/
void p2x_mp_mu(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]  = d->p[5];
    d->iguess[4]  = d->p[4];
    d->iguess[3]  = d->p[3];
    d->iguess[1]  = d->p[0] + d->iguess[4] + d->iguess[3] + d->iguess[2];
    d->iguess[0] = d->p[2]/(1.0-d->iguess[1]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for opx_mp
*/
void p2x_mp_opx(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[4]  = d->p[6];
    d->iguess[1]  = d->p[5];
    d->iguess[3]  = d->p[4];
    d->iguess[2]  = d->p[3];
    d->iguess[0] = (-2.0*d->p[1] - d->p[2])/(d->iguess[4] + d->iguess[3] + 2.0*d->iguess[1] + d->iguess[2] - 2.0);
    d->iguess[5]  = 2.0*(pow(d->iguess[4], 2) + 2.0*d->iguess[4]*d->iguess[3] + 3.0*d->iguess[4]*d->iguess[1] + d->iguess[4]*d->p[0] + 2.0*d->iguess[4]*d->p[1] + d->iguess[4]*d->p[2] + 2.0*d->iguess[4]*d->iguess[2] - 3.0*d->iguess[4] + pow(d->iguess[3], 2) + 3.0*d->iguess[3]*d->iguess[1] + d->iguess[3]*d->p[0] + 2.0*d->iguess[3]*d->iguess[2] - 3.0*d->iguess[3] + 2.0*pow(d->iguess[1], 2) + 2.0*d->iguess[1]*d->p[0] + 2.0*d->iguess[1]*d->p[1] + d->iguess[1]*d->p[2] + 3.0*d->iguess[1]*d->iguess[2] - 4.0*d->iguess[1] + d->p[0]*d->iguess[2] - 2.0*d->p[0] - 2.0*d->p[1] - d->p[2] + pow(d->iguess[2], 2) - 3.0*d->iguess[2] + 2.0)/(pow(d->iguess[4], 2) + d->iguess[4]*d->iguess[3] + 3.0*d->iguess[4]*d->iguess[1] + d->iguess[4]*d->iguess[2] - 3.0*d->iguess[4] + d->iguess[3]*d->iguess[1] - d->iguess[3] + 2.0*pow(d->iguess[1], 2) + d->iguess[1]*d->iguess[2] - 4.0*d->iguess[1] - d->iguess[2] + 2.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}
    
/**
    Endmember to xeos for pl4tr_mp
*/
void p2x_mp_pl4tr(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[1]   = d->p[2];
    d->iguess[0]  = d->p[1];
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for sa_mp
*/
void p2x_mp_sa(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]  = d->p[4];
    d->iguess[1]  = d->p[1];
    d->iguess[0] = (4.0*d->iguess[2] + 4.0*d->p[0] + d->p[3] + 4.0*d->iguess[1] - 4.0)/(d->iguess[2] + d->iguess[1] - 4.0);
    d->iguess[3]  = 4.0/3.0*(-4.0*pow(d->iguess[2], 2) - 4.0*d->iguess[2]*d->p[0] - d->iguess[2]*d->p[2] - d->iguess[2]*d->p[3] - 8.0*d->iguess[2]*d->iguess[1] + 8.0*d->iguess[2] - 4.0*d->p[0]*d->iguess[1] + 4.0*d->p[0] - d->p[2]*d->iguess[1] + 4.0*d->p[2] - d->p[3]*d->iguess[1] + d->p[3] - 4.0*pow(d->iguess[1], 2) + 8.0*d->iguess[1] - 4.0)/(d->iguess[2] + d->iguess[1] - 4.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}


/**
    Endmember to xeos for sp_mp
*/
void p2x_mp_sp(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[2]  = d->p[3];
    d->iguess[1]  = 1.0 - d->p[2] - d->iguess[2];
    d->iguess[0] = (-d->p[1] + d->iguess[2] + 1.0)/(d->iguess[2] + 1.0);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Endmember to xeos for st_mp
*/
void p2x_mp_st(void *SS_ref_db, double eps){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    
    d->iguess[3]  = d->p[4]/(4.0/3.0);
    d->iguess[2]  = d->p[3];
    d->iguess[1]  = d->p[2];
    d->iguess[0] = d->p[1]/(1.0 - d->iguess[1]);
    
    for (int i = 0; i < d->n_xeos; i++){
        if (d->iguess[i] < d->bounds[i][0]){
            d->iguess[i] = d->bounds[i][0];
        }
        if (d->iguess[i] > d->bounds[i][1]){
            d->iguess[i] = d->bounds[i][1];
        }
    }
}

/**
    Update dpdx matrix of fper_S11
*/
void dpdx_ig_fper(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.0;      
    dp_dx[1][0] = 1.0;      
}

/** 
  update dpdpx matrix (biotite)
*/
void dpdx_ig_bi(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[2] + x[3] + x[1] - 1.0;      dp_dx[0][1] = x[0] - 1.0;      dp_dx[0][2] = x[0] - 1.0;      dp_dx[0][3] = x[0] - 1.0;      dp_dx[0][4] = -2.0/3.0;      
    dp_dx[1][0] = 1.0;      dp_dx[1][1] = 0.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      dp_dx[1][4] = -1.0/3.0;      
    dp_dx[2][0] = -x[2] - x[3] - x[1];      dp_dx[2][1] = -x[0];      dp_dx[2][2] = -x[0];      dp_dx[2][3] = -x[0];      dp_dx[2][4] = 1.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 1.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 1.0;      dp_dx[4][4] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 1.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      
}

/** 
  update dpdpx matrix (cordierite)
*/
void dpdx_ig_cd(void *SS_ref_db, const double *x){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	double **dp_dx = d->dp_dx;
   
	dp_dx[0][0]  = -1.;	   dp_dx[0][1]   = -1.;	   
	dp_dx[1][0]  =  1.;	   dp_dx[1][1]   =  0.;	   
	dp_dx[2][0]  =  0.;	   dp_dx[2][1]   =  1.;	   
}

/** 
  update dpdpx matrix (clinopyroxene)
*/
void dpdx_ig_cpx(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = -1.0;      dp_dx[0][2] = -1.0;      dp_dx[0][3] = -1.0;      dp_dx[0][4] = 0.0;      dp_dx[0][5] = 0.0;      dp_dx[0][6] = 0.0;      dp_dx[0][7] = 0.0;      dp_dx[0][8] = -1.0;      
    dp_dx[1][0] = -x[8] - x[3] + x[7] - x[1] + 1.0;      dp_dx[1][1] = -x[4] - x[0];      dp_dx[1][2] = 0.0;      dp_dx[1][3] = -x[4] - x[0];      dp_dx[1][4] = -x[8] - x[3] + x[7] - x[1] + 1.0;      dp_dx[1][5] = 0.0;      dp_dx[1][6] = 0.0;      dp_dx[1][7] = x[4] + x[0];      dp_dx[1][8] = -x[4] - x[0];      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.0;      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 0.0;      dp_dx[2][4] = 0.0;      dp_dx[2][5] = -1.0;      dp_dx[2][6] = -1.0;      dp_dx[2][7] = -2.0;      dp_dx[2][8] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 0.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      dp_dx[3][5] = 0.0;      dp_dx[3][6] = 1.0;      dp_dx[3][7] = 0.0;      dp_dx[3][8] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = 0.0;      dp_dx[4][5] = 1.0;      dp_dx[4][6] = 0.0;      dp_dx[4][7] = 0.0;      dp_dx[4][8] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      dp_dx[5][5] = 0.0;      dp_dx[5][6] = 0.0;      dp_dx[5][7] = 2.0;      dp_dx[5][8] = 0.0;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 0.0;      dp_dx[6][3] = 1.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = 0.0;      dp_dx[6][6] = 0.0;      dp_dx[6][7] = 0.0;      dp_dx[6][8] = 0.0;      
    dp_dx[7][0] = -x[2];      dp_dx[7][1] = -x[4];      dp_dx[7][2] = 1.0 - x[0];      dp_dx[7][3] = -x[4];      dp_dx[7][4] = -x[8] - x[3] + x[7] - x[1] + 1.0;      dp_dx[7][5] = 0.0;      dp_dx[7][6] = 0.0;      dp_dx[7][7] = x[4];      dp_dx[7][8] = -x[4];      
    dp_dx[8][0] = x[8] + x[3] + x[2] - x[7] + x[1] - 1.0;      dp_dx[8][1] = 2.0*x[4] + x[0];      dp_dx[8][2] = x[0];      dp_dx[8][3] = 2.0*x[4] + x[0];      dp_dx[8][4] = 2.0*x[8] + 2.0*x[3] - 2.0*x[7] + 2.0*x[1] - 2.0;      dp_dx[8][5] = 0.0;      dp_dx[8][6] = 0.0;      dp_dx[8][7] = -2.0*x[4] - x[0];      dp_dx[8][8] = 2.0*x[4] + x[0];      
    dp_dx[9][0] = 0.0;      dp_dx[9][1] = 0.0;      dp_dx[9][2] = 0.0;      dp_dx[9][3] = 0.0;      dp_dx[9][4] = 0.0;      dp_dx[9][5] = 0.0;      dp_dx[9][6] = 0.0;      dp_dx[9][7] = 0.0;      dp_dx[9][8] = 1.0;      
}

/** 
  update dpdpx matrix (epidote)
*/
void dpdx_ig_ep(void *SS_ref_db, const double *x){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	double **dp_dx = d->dp_dx;
	
	dp_dx[0][0] = -1.0;         dp_dx[0][1] = -1.0;     
	dp_dx[1][0] =  0.0;         dp_dx[1][1] =  2.0;     
	dp_dx[2][0] =  1.0;         dp_dx[2][1] = -1.0;     
}

/** 
  update dpdpx matrix (fluid)
*/
void dpdx_ig_fl(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -1.0;      dp_dx[0][1] = -1.0;      dp_dx[0][2] = -1.0;      dp_dx[0][3] = -1.0;      dp_dx[0][4] = -1.0;      dp_dx[0][5] = -1.0;      dp_dx[0][6] = -1.0;      dp_dx[0][7] = -1.0;      dp_dx[0][8] = -1.0;      dp_dx[0][9] = -1.0;      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = 1.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      dp_dx[1][4] = 0.0;      dp_dx[1][5] = 0.0;      dp_dx[1][6] = 0.0;      dp_dx[1][7] = 0.0;      dp_dx[1][8] = 0.0;      dp_dx[1][9] = 0.0;      
    dp_dx[2][0] = 1.0;      dp_dx[2][1] = 0.0;      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 0.0;      dp_dx[2][4] = 0.0;      dp_dx[2][5] = 0.0;      dp_dx[2][6] = 0.0;      dp_dx[2][7] = 0.0;      dp_dx[2][8] = 0.0;      dp_dx[2][9] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      dp_dx[3][5] = 0.0;      dp_dx[3][6] = 0.0;      dp_dx[3][7] = 0.0;      dp_dx[3][8] = 0.0;      dp_dx[3][9] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 1.0;      dp_dx[4][4] = 0.0;      dp_dx[4][5] = 0.0;      dp_dx[4][6] = 0.0;      dp_dx[4][7] = 0.0;      dp_dx[4][8] = 0.0;      dp_dx[4][9] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 1.0;      dp_dx[5][5] = 0.0;      dp_dx[5][6] = 0.0;      dp_dx[5][7] = 0.0;      dp_dx[5][8] = 0.0;      dp_dx[5][9] = 0.0;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 0.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = 1.0;      dp_dx[6][6] = 0.0;      dp_dx[6][7] = 0.0;      dp_dx[6][8] = 0.0;      dp_dx[6][9] = 0.0;      
    dp_dx[7][0] = 0.0;      dp_dx[7][1] = 0.0;      dp_dx[7][2] = 0.0;      dp_dx[7][3] = 0.0;      dp_dx[7][4] = 0.0;      dp_dx[7][5] = 0.0;      dp_dx[7][6] = 1.0;      dp_dx[7][7] = 0.0;      dp_dx[7][8] = 0.0;      dp_dx[7][9] = 0.0;      
    dp_dx[8][0] = 0.0;      dp_dx[8][1] = 0.0;      dp_dx[8][2] = 0.0;      dp_dx[8][3] = 0.0;      dp_dx[8][4] = 0.0;      dp_dx[8][5] = 0.0;      dp_dx[8][6] = 0.0;      dp_dx[8][7] = 1.0;      dp_dx[8][8] = 0.0;      dp_dx[8][9] = 0.0;      
    dp_dx[9][0] = 0.0;      dp_dx[9][1] = 0.0;      dp_dx[9][2] = 0.0;      dp_dx[9][3] = 0.0;      dp_dx[9][4] = 0.0;      dp_dx[9][5] = 0.0;      dp_dx[9][6] = 0.0;      dp_dx[9][7] = 0.0;      dp_dx[9][8] = 1.0;      dp_dx[9][9] = 0.0;      
    dp_dx[10][0] = 0.0;      dp_dx[10][1] = 0.0;      dp_dx[10][2] = 0.0;      dp_dx[10][3] = 0.0;      dp_dx[10][4] = 0.0;      dp_dx[10][5] = 0.0;      dp_dx[10][6] = 0.0;      dp_dx[10][7] = 0.0;      dp_dx[10][8] = 0.0;      dp_dx[10][9] = 1.0;      
}

/** 
  update dpdpx matrix (garnet)
*/
void dpdx_ig_g(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[1] - 1.0;      dp_dx[0][1] = x[0] - 1.0;      dp_dx[0][2] = 0.0;      dp_dx[0][3] = -1.0;      dp_dx[0][4] = -4.0;      
    dp_dx[1][0] = 1.0 - x[1];      dp_dx[1][1] = -x[0];      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      dp_dx[1][4] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 1.0;      dp_dx[2][2] = -1.0;      dp_dx[2][3] = 0.0;      dp_dx[2][4] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 1.0;      dp_dx[4][4] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 4.0;      
}

/** 
  update dpdpx matrix (hornblende)
*/
void dpdx_ig_hb(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = -1.0;      dp_dx[0][2] = 1.0;      dp_dx[0][3] = -0.50;      dp_dx[0][4] = 0.0;      dp_dx[0][5] = 1.0;      dp_dx[0][6] = -1.0;      dp_dx[0][7] = -1.0;      dp_dx[0][8] = 0.0;      dp_dx[0][9] = 0.0;      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = 1.0;      dp_dx[1][2] = -1.0;      dp_dx[1][3] = -0.50;      dp_dx[1][4] = 0.0;      dp_dx[1][5] = 0.0;      dp_dx[1][6] = 1.0;      dp_dx[1][7] = 0.0;      dp_dx[1][8] = 0.0;      dp_dx[1][9] = 0.0;      
    dp_dx[2][0] = 0.0;      dp_dx[2][1] = 0.0;      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 1.0 - x[4];      dp_dx[2][4] = -x[3];      dp_dx[2][5] = 0.0;      dp_dx[2][6] = 0.0;      dp_dx[2][7] = 0.0;      dp_dx[2][8] = 0.0;      dp_dx[2][9] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      dp_dx[3][5] = 0.0;      dp_dx[3][6] = -1.0;      dp_dx[3][7] = 0.0;      dp_dx[3][8] = 0.0;      dp_dx[3][9] = 0.0;      
    dp_dx[4][0] = x[5] + x[2] - 1.0;      dp_dx[4][1] = x[9];      dp_dx[4][2] = x[0] - 1.0;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = 0.0;      dp_dx[4][5] = x[0] - 1.0;      dp_dx[4][6] = x[9];      dp_dx[4][7] = x[9];      dp_dx[4][8] = -1.5;      dp_dx[4][9] = x[6] + x[7] + x[1] - 1.0;      
    dp_dx[5][0] = x[5] - x[6] - x[7] - x[1] + x[2] + 1.0;      dp_dx[5][1] = 2.0*x[9] - x[0];      dp_dx[5][2] = x[0];      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      dp_dx[5][5] = x[0];      dp_dx[5][6] = 2.0*x[9] - x[0];      dp_dx[5][7] = 2.0*x[9] - x[0];      dp_dx[5][8] = -2.5;      dp_dx[5][9] = 2.0*x[6] + 2.0*x[7] + 2.0*x[1] - 2.0;      
    dp_dx[6][0] = -x[5] - x[2];      dp_dx[6][1] = -x[9];      dp_dx[6][2] = -x[0];      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = -x[0];      dp_dx[6][6] = -x[9];      dp_dx[6][7] = -x[9];      dp_dx[6][8] = 2.5;      dp_dx[6][9] = -x[6] - x[7] - x[1] + 1.0;      
    dp_dx[7][0] = -x[5] + x[6] + x[7] + x[1] - x[2];      dp_dx[7][1] = -2.0*x[9] + x[0];      dp_dx[7][2] = -x[0];      dp_dx[7][3] = 0.0;      dp_dx[7][4] = 0.0;      dp_dx[7][5] = -x[0];      dp_dx[7][6] = -2.0*x[9] + x[0];      dp_dx[7][7] = -2.0*x[9] + x[0];      dp_dx[7][8] = 1.5;      dp_dx[7][9] = -2.0*x[6] - 2.0*x[7] - 2.0*x[1] + 2.0;      
    dp_dx[8][0] = 0.0;      dp_dx[8][1] = 0.0;      dp_dx[8][2] = 0.0;      dp_dx[8][3] = 0.0;      dp_dx[8][4] = 0.0;      dp_dx[8][5] = 0.0;      dp_dx[8][6] = 1.0;      dp_dx[8][7] = 0.0;      dp_dx[8][8] = 0.0;      dp_dx[8][9] = 0.0;      
    dp_dx[9][0] = 0.0;      dp_dx[9][1] = 0.0;      dp_dx[9][2] = 0.0;      dp_dx[9][3] = x[4];      dp_dx[9][4] = x[3];      dp_dx[9][5] = 0.0;      dp_dx[9][6] = 0.0;      dp_dx[9][7] = 0.0;      dp_dx[9][8] = 0.0;      dp_dx[9][9] = 0.0;      
    dp_dx[10][0] = 0.0;      dp_dx[10][1] = 0.0;      dp_dx[10][2] = 0.0;      dp_dx[10][3] = 0.0;      dp_dx[10][4] = 0.0;      dp_dx[10][5] = 0.0;      dp_dx[10][6] = 0.0;      dp_dx[10][7] = 1.0;      dp_dx[10][8] = 0.0;      dp_dx[10][9] = 0.0;      
}

/** 
  update dpdpx matrix (ilm)
*/
void dpdx_ig_ilm(void *SS_ref_db, const double *x){
	SS_ref *d  = (SS_ref *) SS_ref_db;
	double **dp_dx = d->dp_dx;
	
	dp_dx[0][0] = 0.0;          dp_dx[0][1] = 1.0;      
	dp_dx[1][0] = 1.0;          dp_dx[1][1] = -1.0;      
	dp_dx[2][0] = -1.0;         dp_dx[2][1] = 0.0;   
}

/** 
  update dpdpx matrix (liquid)
*/
void dpdx_ig_liq(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = -3.0*x[9]/4.0 - 1.0;      dp_dx[0][1] = -3.0*x[9]/4.0 - 1.0;      dp_dx[0][2] = -3.0*x[9]/4.0 - 1.0;      dp_dx[0][3] = -3.0*x[9]/4.0 - 1.0;      dp_dx[0][4] = -3.0*x[9]/4.0 - 1.0;      dp_dx[0][5] = -3.0*x[9]/4.0 - 1.0;      dp_dx[0][6] = -3.0*x[9]/4.0 - 1.0;      dp_dx[0][7] = -3.0*x[9]/4.0 - 1.0;      dp_dx[0][8] = -3.0*x[9]/4.0 - 1.0;      dp_dx[0][9] = -3.0*x[6]/4.0 - 3.0*x[3]/4.0 - 3.0*x[2]/4.0 - 3.0*x[10]/4.0 - 3.0*x[5]/4.0 - 3.0*x[4]/4.0 - 3.0*x[8]/4.0 - 3.0*x[1]/4.0 - 3.0*x[7]/4.0 - 3.0*x[0]/4.0 + 1.0;      dp_dx[0][10] = -3.0*x[9]/4.0 - 1.0;      
    dp_dx[1][0] = 0.0;      dp_dx[1][1] = 3.0*x[9]/4.0 + 1.0;      dp_dx[1][2] = 0.0;      dp_dx[1][3] = 0.0;      dp_dx[1][4] = 0.0;      dp_dx[1][5] = 0.0;      dp_dx[1][6] = 0.0;      dp_dx[1][7] = 0.0;      dp_dx[1][8] = 0.0;      dp_dx[1][9] = 3.0*x[1]/4.0 - 1.0;      dp_dx[1][10] = 0.0;      
    dp_dx[2][0] = 3.0*x[9]/4.0 + 1.0;      dp_dx[2][1] = 0.0;      dp_dx[2][2] = 0.0;      dp_dx[2][3] = 0.0;      dp_dx[2][4] = 0.0;      dp_dx[2][5] = 0.0;      dp_dx[2][6] = 0.0;      dp_dx[2][7] = 0.0;      dp_dx[2][8] = 0.0;      dp_dx[2][9] = 3.0*x[0]/4.0 - 1.0;      dp_dx[2][10] = 0.0;      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 3.0*x[9]/4.0 + 1.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      dp_dx[3][5] = 0.0;      dp_dx[3][6] = 0.0;      dp_dx[3][7] = 0.0;      dp_dx[3][8] = 0.0;      dp_dx[3][9] = 3.0*x[2]/4.0;      dp_dx[3][10] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 0.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 3.0*x[9]/4.0 + 1.0;      dp_dx[4][4] = 0.0;      dp_dx[4][5] = 0.0;      dp_dx[4][6] = 0.0;      dp_dx[4][7] = 0.0;      dp_dx[4][8] = 0.0;      dp_dx[4][9] = 3.0*x[3]/4.0;      dp_dx[4][10] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 3.0*x[9]/4.0 + 1.0;      dp_dx[5][5] = 0.0;      dp_dx[5][6] = 0.0;      dp_dx[5][7] = 0.0;      dp_dx[5][8] = 0.0;      dp_dx[5][9] = 3.0*x[4]/4.0;      dp_dx[5][10] = 0.0;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 0.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = 3.0*x[9]/4.0 + 1.0;      dp_dx[6][6] = 0.0;      dp_dx[6][7] = 0.0;      dp_dx[6][8] = 0.0;      dp_dx[6][9] = 3.0*x[5]/4.0;      dp_dx[6][10] = 0.0;      
    dp_dx[7][0] = 0.0;      dp_dx[7][1] = 0.0;      dp_dx[7][2] = 0.0;      dp_dx[7][3] = 0.0;      dp_dx[7][4] = 0.0;      dp_dx[7][5] = 0.0;      dp_dx[7][6] = 3.0*x[9]/4.0 + 1.0;      dp_dx[7][7] = 0.0;      dp_dx[7][8] = 0.0;      dp_dx[7][9] = 3.0*x[6]/4.0;      dp_dx[7][10] = 0.0;      
    dp_dx[8][0] = 0.0;      dp_dx[8][1] = 0.0;      dp_dx[8][2] = 0.0;      dp_dx[8][3] = 0.0;      dp_dx[8][4] = 0.0;      dp_dx[8][5] = 0.0;      dp_dx[8][6] = 0.0;      dp_dx[8][7] = 3.0*x[9]/4.0 + 1.0;      dp_dx[8][8] = 0.0;      dp_dx[8][9] = 3.0*x[7]/4.0;      dp_dx[8][10] = 0.0;      
    dp_dx[9][0] = 0.0;      dp_dx[9][1] = 0.0;      dp_dx[9][2] = 0.0;      dp_dx[9][3] = 0.0;      dp_dx[9][4] = 0.0;      dp_dx[9][5] = 0.0;      dp_dx[9][6] = 0.0;      dp_dx[9][7] = 0.0;      dp_dx[9][8] = 3.0*x[9]/4.0 + 1.0;      dp_dx[9][9] = 3.0*x[8]/4.0;      dp_dx[9][10] = 0.0;      
    dp_dx[10][0] = 0.0;      dp_dx[10][1] = 0.0;      dp_dx[10][2] = 0.0;      dp_dx[10][3] = 0.0;      dp_dx[10][4] = 0.0;      dp_dx[10][5] = 0.0;      dp_dx[10][6] = 0.0;      dp_dx[10][7] = 0.0;      dp_dx[10][8] = 0.0;      dp_dx[10][9] = 1.0;      dp_dx[10][10] = 0.0;      
    dp_dx[11][0] = 0.0;      dp_dx[11][1] = 0.0;      dp_dx[11][2] = 0.0;      dp_dx[11][3] = 0.0;      dp_dx[11][4] = 0.0;      dp_dx[11][5] = 0.0;      dp_dx[11][6] = 0.0;      dp_dx[11][7] = 0.0;      dp_dx[11][8] = 0.0;      dp_dx[11][9] = 3.0*x[10]/4.0;      dp_dx[11][10] = 3.0*x[9]/4.0 + 1.0;      
}

/** 
  update dpdpx matrix (olivine)
*/
void dpdx_ig_ol(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = 0.0;      dp_dx[0][1] = 1.0;      dp_dx[0][2] = 0.0;      
    dp_dx[1][0] = 1.0;      dp_dx[1][1] = 0.0;      dp_dx[1][2] = -1.0;      
    dp_dx[2][0] = x[1] - 1.0;      dp_dx[2][1] = x[0] - 1.0;      dp_dx[2][2] = -1.0;      
    dp_dx[3][0] = -x[1];      dp_dx[3][1] = -x[0];      dp_dx[3][2] = 2.0;      
}

/** 
  update dpdpx matrix (orthopyroxene)
*/
void dpdx_ig_opx(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double **dp_dx = d->dp_dx;

    dp_dx[0][0] = x[2] + x[7] - 1.0;      dp_dx[0][1] = -x[3] - 1.0;      dp_dx[0][2] = x[0] - 1.0;      dp_dx[0][3] = -x[7] + x[5] - x[1] + 1.0;      dp_dx[0][4] = 0.0;      dp_dx[0][5] = x[3];      dp_dx[0][6] = 0.0;      dp_dx[0][7] = -x[3] + x[0] - 1.0;      
    dp_dx[1][0] = -x[7] + x[5] - x[1] + 1.0;      dp_dx[1][1] = -x[3] - x[0];      dp_dx[1][2] = 0.0;      dp_dx[1][3] = -x[7] + x[5] - x[1] + 1.0;      dp_dx[1][4] = 0.0;      dp_dx[1][5] = x[3] + x[0];      dp_dx[1][6] = 0.0;      dp_dx[1][7] = -x[3] - x[0];      
    dp_dx[2][0] = -x[2] - x[5] + x[1];      dp_dx[2][1] = 2.0*x[3] + x[0];      dp_dx[2][2] = -x[0];      dp_dx[2][3] = 2.0*x[7] - 2.0*x[5] + 2.0*x[1] - 2.0;      dp_dx[2][4] = 0.0;      dp_dx[2][5] = -2.0*x[3] - x[0];      dp_dx[2][6] = 0.0;      dp_dx[2][7] = 2.0*x[3];      
    dp_dx[3][0] = 0.0;      dp_dx[3][1] = 0.0;      dp_dx[3][2] = 1.0;      dp_dx[3][3] = 0.0;      dp_dx[3][4] = 0.0;      dp_dx[3][5] = 0.0;      dp_dx[3][6] = 0.0;      dp_dx[3][7] = 0.0;      
    dp_dx[4][0] = 0.0;      dp_dx[4][1] = 1.0;      dp_dx[4][2] = 0.0;      dp_dx[4][3] = 0.0;      dp_dx[4][4] = -1.0;      dp_dx[4][5] = -2.0;      dp_dx[4][6] = -1.0;      dp_dx[4][7] = 0.0;      
    dp_dx[5][0] = 0.0;      dp_dx[5][1] = 0.0;      dp_dx[5][2] = 0.0;      dp_dx[5][3] = 0.0;      dp_dx[5][4] = 0.0;      dp_dx[5][5] = 0.0;      dp_dx[5][6] = 1.0;      dp_dx[5][7] = 0.0;      
    dp_dx[6][0] = 0.0;      dp_dx[6][1] = 0.0;      dp_dx[6][2] = 0.0;      dp_dx[6][3] = 0.0;      dp_dx[6][4] = 0.0;      dp_dx[6][5] = 2.0;      dp_dx[6][6] = 0.0;      dp_dx[6][7] = 0.0;      
    dp_dx[7][0] = 0.0;      dp_dx[7][1] = 0.0;      dp_dx[7][2] = 0.0;      dp_dx[7][3] = 0.0;      dp_dx[7][4] = 1.0;      dp_dx[7][5] = 0.0;      dp_dx[7][6] = 0.0;      dp_dx[7][7] = 0.0;      
    dp_dx[8][0] = 0.0;      dp_dx[8][1] = 0.0;      dp_dx[8][2] = 0.0;      dp_dx[8][3] = 0.0;      dp_dx[8][4] = 0.0;      dp_dx[8][5] = 0.0;      dp_dx[8][6] = 0.0;      dp_dx[8][7] = 1.0;      
}

/** 
  update dpdpx matrix (plagioclase)
*/
void dpdx_ig_pl4T(void *SS_ref_db, const double *x){
	SS_ref *d  = (SS_ref *) SS_ref_db;

	double **dp_dx = d->dp_dx;

	dp_dx[0][0] = -1.0;         dp_dx[0][1] = -1.0;      
	dp_dx[1][0] =  1.0;         dp_dx[1][1] =  0.0;      
	dp_dx[2][0] =  0.0;         dp_dx[2][1] =  1.0;  
}

/** 
  update dpdpx matrix (spinel)
*/
void dpdx_ig_spn(void *SS_ref_db, const double *x){
	SS_ref *d  = (SS_ref *) SS_ref_db;
	
	double **dp_dx = d->dp_dx;

	dp_dx[0][0] = -x[3]/3. - 1./3.;        dp_dx[0][1] = 0.0;         					dp_dx[0][2] = -1.0;        dp_dx[0][3] = -x[0]/3. - 2./3.;        dp_dx[0][4] = 2./3.;       dp_dx[0][5] = 0.0;         dp_dx[0][6] = 0.0;     
	dp_dx[1][0] = -2.*x[3]/3. - 2./3.;     dp_dx[1][1] = 0.0;         					dp_dx[1][2] = 0.0;         dp_dx[1][3] = -2.*x[0]/3. - 1./3.;     dp_dx[1][4] = -2./3.;      dp_dx[1][5] = 0.0;         dp_dx[1][6] = 0.0;     
	dp_dx[2][0] = x[3]/3. + 1./3.;         dp_dx[2][1] = x[2]/3. + x[3]/3. - 1./3.;         	dp_dx[2][2] = x[1]/3.;        dp_dx[2][3] = x[0]/3. + x[1]/3.;         	dp_dx[2][4] = 0.0;         dp_dx[2][5] = 2./3.;       dp_dx[2][6] = 2./3.;     
	dp_dx[3][0] = 2.*x[3]/3. + 2./3.;      dp_dx[3][1] = 2.*x[2]/3. + 2.*x[3]/3. - 2./3.;    dp_dx[3][2] = 2.*x[1]/3.;     dp_dx[3][3] = 2.*x[0]/3. + 2.*x[1]/3.;    dp_dx[3][4] = 0.0;         dp_dx[3][5] = -2./3.;      dp_dx[3][6] = -2./3.;     
	dp_dx[4][0] = 0.0;         			dp_dx[4][1] = -x[2]/3. - x[3]/3. + 1./3.;         dp_dx[4][2] = -x[1]/3.;       dp_dx[4][3] = -x[1]/3.;         		dp_dx[4][4] = 0.0;         dp_dx[4][5] = 0.0;         dp_dx[4][6] = -2./3.;     
	dp_dx[5][0] = 0.0;         			dp_dx[5][1] = -2.*x[2]/3. - 2.*x[3]/3. + 2./3.;   dp_dx[5][2] = -2.*x[1]/3.;    dp_dx[5][3] = -2.*x[1]/3.;         		dp_dx[5][4] = 0.0;         dp_dx[5][5] = 0.0;         dp_dx[5][6] = 2./3.;     
	dp_dx[6][0] = 0.0;         			dp_dx[6][1] = 0.0;         					dp_dx[6][2] = 1.0;         dp_dx[6][3] = 0.0;         			dp_dx[6][4] = 0.0;         dp_dx[6][5] = 0.0;         dp_dx[6][6] = 0.0;     
	dp_dx[7][0] = 0.0;         			dp_dx[7][1] = 0.0;        					dp_dx[7][2] = 0.0;         dp_dx[7][3] = 1.0;         			dp_dx[7][4] = 0.0;         dp_dx[7][5] = 0.0;         dp_dx[7][6] = 0.0;     
}


/** 
  update dpdpx matrix (biotite)
*/
void px_ig_bi(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -2.*x[4]/3.0 + x[2]*x[0] - x[2] + x[3]*x[0] - x[3] + x[0]*x[1] - x[0] - x[1] + 1.0;
        p[1]           = -x[4]/3.0 + x[0];
        p[2]           = x[4] - x[2]*x[0] - x[3]*x[0] - x[0]*x[1];
        p[3]           = x[1];
        p[4]           = x[3];
        p[5]           = x[2];
}
/**
    Endmember fraction of fper_S11
*/
void px_ig_fper(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = 1.0 - 1.0*x[0];
        p[1]           = x[0];
}

/** 
  update px matrix (cordierite)
*/
void px_ig_cd(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -x[1] - x[0] + 1.0;
        p[1]           = x[0];
        p[2]           = x[1];
}

/** 
  update px matrix (clinopyroxene)
*/
void px_ig_cpx(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = - x[8] - x[3] - x[2] - x[1] + 1.0;
        p[1]           = - x[4]*x[8] - x[4]*x[3] + x[4]*x[7] - x[4]*x[1] + x[4] - x[8]*x[0] - x[3]*x[0] + x[7]*x[0] - x[0]*x[1] + x[0];
        p[2]           = - x[6] - x[5] - 2.0*x[7] + x[1];
        p[3]           = x[6];
        p[4]           = x[5];
        p[5]           = 2.0*x[7];
        p[6]           = x[3];
        p[7]           = - x[4]*x[8] - x[4]*x[3] + x[4]*x[7] - x[4]*x[1] + x[4] - x[2]*x[0] + x[2];
        p[8]           = 2.0*x[4]*x[8] + 2.0*x[4]*x[3] - 2.0*x[4]*x[7] + 2.0*x[4]*x[1] - 2.0*x[4] + x[8]*x[0] + x[3]*x[0] + x[2]*x[0] - x[7]*x[0] + x[0]*x[1] - x[0];
        p[9]           = x[8];
}

/** 
  update px matrix (epidote)
*/
void px_ig_ep(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = - x[1] - x[0] + 1.0;
        p[1]           = 2.0*x[1];
        p[2]           = - x[1] + x[0];
}

/** 
  update px matrix (fluid)
*/
void px_ig_fl(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = - x[6] - x[3] - x[2] - x[9] - x[5] - x[4] - x[8] - x[1] - x[7] - x[0] + 1.0;
        p[1]           = x[1];
        p[2]           = x[0];
        p[3]           = x[2];
        p[4]           = x[3];
        p[5]           = x[4];
        p[6]           = x[5];
        p[7]           = x[6];
        p[8]           = x[7];
        p[9]           = x[8];
        p[10]           = x[9];
}

/** 
  update px matrix (garnet)
*/
void px_ig_g(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[1]*x[0] - x[1] - x[3] - 4.0*x[4] - x[0] + 1.0;
        p[1]           = - x[1]*x[0] + x[0];
        p[2]           = x[1] - x[2];
        p[3]           = x[2];
        p[4]           = x[3];
        p[5]           = 4.0*x[4];
}

/** 
  update px matrix (hornblende)
*/
void px_ig_hb(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = -0.5*x[3] + x[5] - x[6] - x[7] - x[1] + x[2];
        p[1]           = -0.5*x[3] + x[6] + x[1] - x[2];
        p[2]           = - x[3]*x[4] + x[3];
        p[3]           = - x[6] + x[2];
        p[4]           = -1.5*x[8] + x[9]*x[6] + x[9]*x[7] + x[9]*x[1] - x[9] + x[5]*x[0] - x[5] + x[0]*x[2] - x[0] - x[2] + 1.0;
        p[5]           = -2.5*x[8] + 2.0*x[9]*x[6] + 2.0*x[9]*x[7] + 2.0*x[9]*x[1] - 2.0*x[9] + x[5]*x[0] - x[6]*x[0] - x[7]*x[0] - x[0]*x[1] + x[0]*x[2] + x[0];
        p[6]           = 2.5*x[8] - x[9]*x[6] - x[9]*x[7] - x[9]*x[1] + x[9] - x[5]*x[0] - x[0]*x[2];
        p[7]           = 1.5*x[8] - 2.0*x[9]*x[6] - 2.0*x[9]*x[7] - 2.0*x[9]*x[1] + 2.0*x[9] - x[5]*x[0] + x[6]*x[0] + x[7]*x[0] + x[0]*x[1] - x[0]*x[2];
        p[8]           = x[6];
        p[9]           = x[3]*x[4];
        p[10]           = x[7];
}

/** 
  update px matrix (ilm)
*/
void px_ig_ilm(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[1];
        p[1]           = x[0] - x[1];
        p[2]           = 1.0 - x[0];
}

/** 
  update px matrix (liquid)
*/
void px_ig_liq(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = - x[6] - x[3] - x[2] - x[10] - x[5] - x[4] - x[8] - x[1] - x[7] - x[0] + 0.25*x[9]*(-3.0*x[6] - 3.0*x[3] - 3.0*x[2] - 3.0*x[10] - 3.0*x[5] - 3.0*x[4] - 3.0*x[8] - 3.0*x[1] - 3.0*x[7] - 3.0*x[0] + 4.0) + 1.0;
        p[1]           = 3.0*x[1]*x[9]/4.0 + x[1] - x[9];
        p[2]           = 3.0*x[0]*x[9]/4.0 + x[0] - x[9];
        p[3]           = 3.0*x[2]*x[9]/4.0 + x[2];
        p[4]           = 3.0*x[3]*x[9]/4.0 + x[3];
        p[5]           = 3.0*x[4]*x[9]/4.0 + x[4];
        p[6]           = 3.0*x[5]*x[9]/4.0 + x[5];
        p[7]           = 3.0*x[6]*x[9]/4.0 + x[6];
        p[8]           = 3.0*x[7]*x[9]/4.0 + x[7];
        p[9]           = 3.0*x[8]*x[9]/4.0 + x[8];
        p[10]          = x[9];
        p[11]          = 3.0*x[10]*x[9]/4.0 + x[10];
}

/** 
  update px matrix (olivine)
*/
void px_ig_ol(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = x[1];
        p[1]           = - x[2] + x[0];
        p[2]           = - x[2] + x[1]*x[0] - x[1] - x[0] + 1.0;
        p[3]           = 2.0*x[2] - x[1]*x[0];
}

/** 
  update px matrix (orthopyroxene)
*/
void px_ig_opx(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = - x[3]*x[7] + x[3]*x[5] - x[3]*x[1] + x[3] + x[2]*x[0] - x[2] + x[7]*x[0] - x[7] - x[0] - x[1] + 1.0;
        p[1]           = - x[3]*x[7] + x[3]*x[5] - x[3]*x[1] + x[3] - x[7]*x[0] + x[5]*x[0] - x[0]*x[1] + x[0];
        p[2]           = 2.0*x[3]*x[7] - 2.0*x[3]*x[5] + 2.0*x[3]*x[1] - 2.0*x[3] - x[2]*x[0] - x[5]*x[0] + x[0]*x[1];
        p[3]           = x[2];
        p[4]           = - x[6] - x[4] - 2.0*x[5] + x[1];
        p[5]           = x[6];
        p[6]           = 2.0*x[5];
        p[7]           = x[4];
        p[8]           = x[7];
}

/** 
  update px matrix (plagioclase 4T)
*/
void px_ig_pl4T(void *SS_ref_db, const double *x){
    SS_ref *d  = (SS_ref *) SS_ref_db;
    double *p = d->p;
        p[0]           = - x[0] - x[1] + 1.0;
        p[1]           = x[0];
        p[2]           = x[1];
}

/** 
  update px matrix (spinel)
*/
void px_ig_spn(void *SS_ref_db, const double *x){
	SS_ref *d  = (SS_ref *) SS_ref_db;
	double *p = d->p;

	p[0]           =  2.*x[4]/3. - x[2] - x[3]*x[0]/3. - 2.*x[3]/3. - x[0]/3. + 1./3.;
	p[1]           = -2.*x[4]/3. - 2.*x[3]*x[0]/3. - x[3]/3. - 2.*x[0]/3. + 2./3.;
	p[2]           =  2.*x[5]/3. + 2.*x[6]/3. + x[2]*x[1]/3. + x[3]*x[0]/3. + x[3]*x[1]/3. + x[0]/3. - x[1]/3.;
	p[3]           = -2.*x[5]/3. - 2.*x[6]/3. + 2.*x[2]*x[1]/3. + 2.*x[3]*x[0]/3. + 2.*x[3]*x[1]/3. + 2.*x[0]/3. - 2.*x[1]/3.;
	p[4]           = -2.*x[6]/3. - x[2]*x[1]/3. - x[3]*x[1]/3. + x[1]/3.;
	p[5]           =  2.*x[6]/3. - 2.*x[2]*x[1]/3. - 2.*x[3]*x[1]/3. + 2.*x[1]/3.;
	p[6]           =  x[2];
	p[7]           =  x[3];
}

   
/**
    Objective function of fper_S11
*/
double obj_ig_fper(unsigned n, const double *x, double *grad, void *SS_ref_db){
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    px_ig_fper(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0*x[0];
    sf[1]          = 1.0 - x[0];
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[1])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[0])) + mu_Gex[1];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_ig_fper(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}


/** 
  objective function of biotite
*/
double obj_ig_bi(unsigned  n, const double *x, double *grad, void *SS_ref_db) {
	SS_ref *d  = (SS_ref *) SS_ref_db;

	
	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;
    double *d_em   = d->d_em;
	px_ig_bi(SS_ref_db,x);

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
				it += 1;
			}
		}
	}
	
    sf[0]           = x[2]*x[0] - x[2] - 2.0/3.0*x[4] + x[3]*x[0] - x[3] + x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[1]           = -x[2]*x[0] + 2.0/3.0*x[4] - x[3]*x[0] - x[0]*x[1] + x[0];
    sf[2]           = x[2];
    sf[3]           = x[3];
    sf[4]           = x[1];
    sf[5]           = 1.0/3.0*x[4] - x[0] + 1.0;
    sf[6]           = -1.0/3.0*x[4] + x[0];
    sf[7]           = -0.5*x[2] - 0.5*x[1] + 0.5;
    sf[8]           = 0.5*x[2] + 0.5*x[1] + 0.5;
    sf[9]           = 1.0 - x[3];
    sf[10]          = x[3];

	mu[0]          = R*T*creal(clog( 4.0*sf[0]*pow(sf[5], 2.0)*sf[7]*sf[8]*pow(sf[9], 2.0))) + gb[0] + mu_Gex[0];
	mu[1]          = R*T*creal(clog( 4.0*sf[1]*pow(sf[6], 2.0)*sf[7]*sf[8]*pow(sf[9], 2.0))) + gb[1] + mu_Gex[1];
	mu[2]          = R*T*creal(clog( 4.0*sf[1]*pow(sf[5], 2.0)*sf[7]*sf[8]*pow(sf[9], 2.0))) + gb[2] + mu_Gex[2];
	mu[3]          = R*T*creal(clog( sf[4]*pow(sf[5], 2.0)*pow(sf[8], 2.0)*pow(sf[9], 2.0))) + gb[3] + mu_Gex[3];
	mu[4]          = R*T*creal(clog( 4.0*sf[3]*pow(sf[5], 2.0)*sf[7]*sf[8])* pow(sf[10], 2.0)) + gb[4] + mu_Gex[4];
	mu[5]          = R*T*creal(clog( sf[2]*pow(sf[5], 2.0)*pow(sf[8], 2.0)*pow(sf[9], 2.0) + d_em[5])) + gb[5] + mu_Gex[5];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_bi(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
};
/** 
  objective function of cordierite
*/
double obj_ig_cd(unsigned  n, const double *x, double *grad, void *SS_ref_db) {

	SS_ref *d  = (SS_ref *) SS_ref_db;


	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;

	px_ig_cd(SS_ref_db,x);

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
				it += 1;
			}
		}
	}
	
    sf[0]           = x[0];
    sf[1]           = 1.0 - x[0];
    sf[2]           = x[1];
    sf[3]           = 1.0 - x[1];
    
	mu[0]            = R*T*creal(clog( pow(sf[1], 2.0)*sf[3])) + gb[0]  + mu_Gex[0];
	mu[1]            = R*T*creal(clog( pow(sf[0], 2.0)*sf[3])) + gb[1]  + mu_Gex[1];
	mu[2]            = R*T*creal(clog( pow(sf[1], 2.0)*sf[2])) + gb[2]  + mu_Gex[2];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_cd(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
};

/** 
  objective function of clinopyroxene
*/
double obj_ig_cpx(unsigned n, const double *x, double *grad, void *SS_ref_db) {

	SS_ref *d  = (SS_ref *) SS_ref_db;


	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *p      = d->p;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;
    double *d_em   = d->d_em;

	px_ig_cpx(SS_ref_db,x);

	d->sum_v = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->sum_v += d->p[i]*d->v[i];
	}
	for (int i = 0; i < d->n_em; i++){
		d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
	}

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
				it += 1;
			}
		}
	}
	
    sf[0]           = x[8]*x[4] + x[8]*x[0] - x[8] + x[3]*x[4] + x[3]*x[0] - x[3] - x[4]*x[7] + x[4]*x[1] - x[4] - x[7]*x[0] + x[7] + x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[1]           = -x[8]*x[4] - x[8]*x[0] - x[3]*x[4] - x[3]*x[0] + x[4]*x[7] - x[4]*x[1] + x[4] + x[7]*x[0] - x[0]*x[1] + x[0];
    sf[2]           = -x[6] - x[5] + x[8] + x[3] - 2.0*x[7] + x[1];
    sf[3]           = x[5];
    sf[4]           = x[6];
    sf[5]           = x[7];
    sf[6]           = -x[8]*x[4] - x[3]*x[4] - x[2]*x[0] + x[2] + x[4]*x[7] - x[4]*x[1] + x[4];
    sf[7]           = x[8]*x[4] + x[3]*x[4] + x[2]*x[0] - x[4]*x[7] + x[4]*x[1] - x[4];
    sf[8]           = -x[8] - x[3] - x[2] + 1.0;
    sf[9]           = x[3];
    sf[10]           = x[8];
    sf[11]           = 1.0 - 0.5*x[1];
    sf[12]           = 0.5*x[1];

	mu[0]          = R*T*creal(clog(sf[0]*csqrt(sf[11])*sf[8])) + gb[0] + mu_Gex[0];
	mu[1]          = R*T*creal(clog(csqrt(sf[11])*sf[1]*sf[7])) + gb[1] + mu_Gex[1];
	mu[2]          = R*T*creal(clog(1.4142*cpow(sf[11], 0.25)*cpow(sf[12], 0.25)*sf[2]*sf[8])) + gb[2] + mu_Gex[2];
	mu[3]          = R*T*creal(clog(1.4142*cpow(sf[11], 0.25)*cpow(sf[12], 0.25)*sf[4]*sf[8] + d_em[3])) + gb[3] + mu_Gex[3];
	mu[4]          = R*T*creal(clog(1.4142*cpow(sf[11], 0.25)*cpow(sf[12], 0.25)*sf[3]*sf[8] + d_em[4])) + gb[4] + mu_Gex[4];
	mu[5]          = R*T*creal(clog(2.8284*csqrt(sf[0])*cpow(sf[11], 0.25)*cpow(sf[12], 0.25)*csqrt(sf[5])*sf[8])) + gb[5] + mu_Gex[5];
	mu[6]          = R*T*creal(clog(csqrt(sf[11])*sf[2]*sf[9])) + gb[6] + mu_Gex[6];
	mu[7]          = R*T*creal(clog(sf[0]*csqrt(sf[11])*sf[6])) + gb[7] + mu_Gex[7];
	mu[8]          = R*T*creal(clog(sf[0]*csqrt(sf[11])*sf[7])) + gb[8] + mu_Gex[8];
	mu[9]          = R*T*creal(clog(sf[10]*csqrt(sf[11])*sf[2])) + gb[9] + mu_Gex[9];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_cpx(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
};


/** 
  objective function of epidote
*/
double obj_ig_ep(unsigned  n, const double *x, double *grad, void *SS_ref_db) {
    SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_ig_ep(SS_ref_db,x);

    for (int i = 0; i < n_em; i++){
        mu_Gex[i] = 0.0;
        int it    = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
                it += 1;
            }
        }
    }
    
    sf[0]          = x[0] - x[1];
    sf[1]          = -x[0] + x[1] + 1.0;
    sf[2]          = x[0] + x[1];
    sf[3]          = -x[0] - x[1] + 1.0;
    
    mu[0]          = gb[0] + R*T*creal(clog(sf[1]*sf[3])) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(sf[1]*sf[2] + d_em[1])) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(sf[0]*sf[2] + d_em[2])) + mu_Gex[2];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_ig_ep(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
};

/** 
  objective function of fluid
*/
double obj_ig_fl(unsigned  n, const double *x, double *grad, void *SS_ref_db) {

	SS_ref *d  = (SS_ref *) SS_ref_db;


	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;
    double *d_em   = d->d_em;

	px_ig_fl(SS_ref_db,x);

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
				it += 1;
			}
		}
	}
	
    sf[0]           = -x[6] - x[3] - x[2] - x[9] - x[5] - x[4] - x[8] - x[1] - x[7] - x[0] + 1.0;
    sf[1]           = x[1];
    sf[2]           = x[0];
    sf[3]           = x[2];
    sf[4]           = x[3];
    sf[5]           = x[4];
    sf[6]           = x[5];
    sf[7]           = x[6];
    sf[8]           = x[7];
    sf[9]           = x[8];
    sf[10]           = x[9];
    sf[11]           = 1.0 - x[9];

	mu[0]          = R*T*creal(clog(sf[0]*sf[11]))  + gb[0]  + mu_Gex[0];
	mu[1]          = R*T*creal(clog(sf[11]*sf[1]))  + gb[1]  + mu_Gex[1];
	mu[2]          = R*T*creal(clog(sf[11]*sf[2]))  + gb[2]  + mu_Gex[2];
	mu[3]          = R*T*creal(clog(sf[11]*sf[3]))  + gb[3]  + mu_Gex[3];
	mu[4]          = R*T*creal(clog(sf[11]*sf[4]))  + gb[4]  + mu_Gex[4];
	mu[5]          = R*T*creal(clog(sf[11]*sf[5]))  + gb[5]  + mu_Gex[5];
	mu[6]          = R*T*creal(clog(sf[11]*sf[6] + d_em[6]))  + gb[6]  + mu_Gex[6];
	mu[7]          = R*T*creal(clog(sf[11]*sf[7] + d_em[7]))  + gb[7]  + mu_Gex[7];
	mu[8]          = R*T*creal(clog(sf[11]*sf[8]))  + gb[8]  + mu_Gex[8];
	mu[9]          = R*T*creal(clog(sf[11]*sf[9]))  + gb[9]  + mu_Gex[9];
	mu[10]         = R*T*creal(clog( pow(sf[10], 2.0))) + gb[10] + mu_Gex[10];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_fl(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
};


/** 
  objective function of garnet
*/
double obj_ig_g(unsigned   n, const double *x, double *grad, void *SS_ref_db) {

	SS_ref *d  = (SS_ref *) SS_ref_db;


	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;
    double *d_em   = d->d_em;

	px_ig_g(SS_ref_db,x);

	d->sum_v = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->sum_v += d->p[i]*d->v[i];
	}
	for (int i = 0; i < d->n_em; i++){
		d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
	}

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
				it += 1;
			}
		}
	}
	
    sf[0]           = x[1]*x[0] - x[1] - x[0] + 1.0;
    sf[1]           = -x[1]*x[0] + x[0];
    sf[2]           = x[1];
    sf[3]           = -x[3] - x[2] - 2.0*x[4] + 1.0;
    sf[4]           = x[3];
    sf[5]           = x[2];
    sf[6]           = x[4];

	mu[0]          = R*T*creal(clog( pow(sf[0], 3.0)* pow(sf[3], 2.0))) + gb[0] + mu_Gex[0];
	mu[1]          = R*T*creal(clog( pow(sf[1], 3.0)* pow(sf[3], 2.0))) + gb[1] + mu_Gex[1];
	mu[2]          = R*T*creal(clog( pow(sf[2], 3.0)* pow(sf[3], 2.0))) + gb[2] + mu_Gex[2];
	mu[3]          = R*T*creal(clog( pow(sf[2], 3.0)* pow(sf[5], 2.0) + d_em[3])) + gb[3] + mu_Gex[3];
	mu[4]          = R*T*creal(clog( pow(sf[0], 3.0)* pow(sf[4], 2.0) + d_em[4])) + gb[4] + mu_Gex[4];
	mu[5]          = R*T*creal(clog(8.0* pow(sf[0], 3.0)*sf[3]*sf[6])) + gb[5] + mu_Gex[5];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_g(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
};


/** 
  objective function of hornblende
*/
double obj_ig_hb(unsigned  n, const double *x, double *grad, void *SS_ref_db) {
     SS_ref *d         = (SS_ref *) SS_ref_db;

    int n_em          = d->n_em;
    double P          = d->P;
    double T          = d->T;
    double R          = d->R;

    double *gb        = d->gb_lvl;
    double *mat_phi   = d->mat_phi;
    double *mu_Gex    = d->mu_Gex;
    double *sf        = d->sf;
    double *mu        = d->mu;
    double *d_em      = d->d_em;
    px_ig_hb(SS_ref_db,x);

    d->sum_v = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_v += d->p[i]*d->v[i];
    }
    for (int i = 0; i < n_em; i++){
        d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
    }
    
    for (int i = 0; i < d->n_em; i++){
        mu_Gex[i] = 0.0;
        int it = 0;
        for (int j = 0; j < d->n_xeos; j++){
            for (int k = j+1; k < d->n_em; k++){
                mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
                it += 1;
            }
        }
    }
    
    sf[0]          = 1.0 - x[3];
    sf[1]          = -x[3]*x[4] + x[3];
    sf[2]          = x[3]*x[4];
    sf[3]          = -x[0] + x[8] + 1.0;
    sf[4]          = x[0] - x[8];
    sf[5]          = x[0]*x[1] + x[0]*x[6] + x[0]*x[7] - x[0] - x[1]*x[9] - x[1] - x[6]*x[9] - x[6] - x[7]*x[9] - x[7] + x[9] + 1.0;
    sf[6]          = -x[0]*x[1] - x[0]*x[6] - x[0]*x[7] + x[0] + x[1]*x[9] + x[6]*x[9] + x[7]*x[9] - x[9];
    sf[7]          = x[1];
    sf[8]          = x[6];
    sf[9]          = x[7];
    sf[10]          = x[5];
    sf[11]          = x[0]*x[2] + x[0]*x[5] - x[0] + x[1]*x[9] - x[2] - x[5] + x[6]*x[9] + x[7]*x[9] - 1.5*x[8] - x[9] + 1.0;
    sf[12]          = -x[0]*x[2] - x[0]*x[5] + x[0] - x[1]*x[9] - x[6]*x[9] - x[7]*x[9] + 1.5*x[8] + x[9];
    sf[13]          = x[2];
    sf[14]          = -0.5*x[1] + 0.5*x[2] - 0.25*x[3] - 0.5*x[6] - 0.5*x[7] + 1.0;
    sf[15]          = 0.5*x[1] - 0.5*x[2] + 0.25*x[3] + 0.5*x[6] + 0.5*x[7];
    sf[16]          = 1.0 - x[7];

    mu[0]          = gb[0] + R*T*creal(clog(sf[0]*cpow(sf[10], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[5], 2.0))) + mu_Gex[0];
    mu[1]          = gb[1] + R*T*creal(clog(2.0*sf[0]*cpow(sf[10], 2.0)*sqrt(sf[14])*sqrt(sf[15])*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[7], 2.0))) + mu_Gex[1];
    mu[2]          = gb[2] + R*T*creal(clog(8.0*cpow(sf[10], 2.0)*sqrt(sf[14])*sqrt(sf[15])*cpow(sf[16], 2.0)*sf[1]*cpow(sf[3], 3.0)*sf[5]*sf[7])) + mu_Gex[2];
    mu[3]          = gb[3] + R*T*creal(clog(sf[0]*cpow(sf[13], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[7], 2.0))) + mu_Gex[3];
    mu[4]          = gb[4] + R*T*creal(clog(sf[0]*cpow(sf[11], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[5], 2.0))) + mu_Gex[4];
    mu[5]          = gb[5] + R*T*creal(clog(sf[0]*cpow(sf[12], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[4], 3.0)*cpow(sf[6], 2.0))) + mu_Gex[5];
    mu[6]          = gb[6] + R*T*creal(clog(sf[0]*cpow(sf[12], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[6], 2.0))) + mu_Gex[6];
    mu[7]          = gb[7] + R*T*creal(clog(sf[0]*cpow(sf[12], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[4], 3.0)*cpow(sf[5], 2.0))) + mu_Gex[7];
    mu[8]          = gb[8] + R*T*creal(clog(sf[0]*cpow(sf[13], 2.0)*sf[14]*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[8], 2.0) + d_em[8])) + mu_Gex[8];
    mu[9]          = gb[9] + R*T*creal(clog(8.0*cpow(sf[10], 2.0)*sqrt(sf[14])*sqrt(sf[15])*cpow(sf[16], 2.0)*sf[2]*cpow(sf[3], 3.0)*sf[5]*sf[7])) + mu_Gex[9];
    mu[10]         = gb[10] + R*T*creal(clog(2.0*sf[0]*cpow(sf[10], 2.0)*sqrt(sf[14])*sqrt(sf[15])*cpow(sf[16], 2.0)*cpow(sf[3], 3.0)*cpow(sf[9], 2.0))) + mu_Gex[10];
    
    d->sum_apep = 0.0;
    for (int i = 0; i < n_em; i++){
        d->sum_apep += d->ape[i]*d->p[i];
    }
    d->factor = d->fbc/d->sum_apep;

    d->df_raw = 0.0;
    for (int i = 0; i < n_em; i++){
        d->df_raw += mu[i]*d->p[i];
    }
    d->df = d->df_raw * d->factor;

    if (grad){
        double *dfx    = d->dfx;
        double **dp_dx = d->dp_dx;
        dpdx_ig_hb(SS_ref_db,x);
        for (int i = 0; i < (d->n_xeos); i++){
            dfx[i] = 0.0;
            for (int j = 0; j < n_em; j++){
                dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
            }
            grad[i] = creal(dfx[i]);
        }
    }

    return d->df;
}
   

  
/** 
  objective function of ilmenite
*/
double obj_ig_ilm(unsigned n, const double *x, double *grad, void *SS_ref_db) {

	SS_ref *d  = (SS_ref *) SS_ref_db;


	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;
    double *d_em   = d->d_em;

	px_ig_ilm(SS_ref_db,x);

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
				it += 1;
			}
		}
	}
	
    sf[0]           = 0.5*x[1] + 0.5*x[0];
    sf[1]           = -0.5*x[1] + 0.5*x[0];
    sf[2]           = 1.0 - x[0];
    sf[3]           = -0.5*x[1] + 0.5*x[0];
    sf[4]           = 0.5*x[1] + 0.5*x[0];
    sf[5]           = 1.0 - x[0];

	mu[0]         = R*T*creal(clog(csqrt(sf[0])*csqrt(sf[4]))) + gb[0] + mu_Gex[0];
	mu[1]         = R*T*creal(clog(2.0*cpow(sf[0], 0.25)*cpow(sf[1], 0.25)*cpow(sf[3], 0.25)*cpow(sf[4], 0.25))) + gb[1] + mu_Gex[1];
	mu[2]         = R*T*creal(clog(csqrt(sf[2])*csqrt(sf[5]) + d_em[2])) + gb[2] + mu_Gex[2];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_ilm(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
};
   
/** 
  objective function of liquid
*/
double obj_ig_liq(unsigned n, const double *x, double *grad, void *SS_ref_db) {

	SS_ref *d  = (SS_ref *) SS_ref_db;

    double muGex;
	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;
    double *d_em   = d->d_em;
	px_ig_liq(SS_ref_db,x);
	
	d->sum_v = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->sum_v += d->p[i]*d->v[i];
	}
	for (int i = 0; i < d->n_em; i++){
		d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
	}

	for (int i = 0; i < d->n_em; i++){
		muGex = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				muGex -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
				it += 1;
			}
		}
        mu_Gex[i] = muGex;
	}

	sf[0]           = -x[6] - x[3] - x[2] - x[10] - x[5] - x[4] - x[8] - x[1] - x[7] - x[0] + 0.25*x[9]*(-3.0*x[6] - 3.0*x[3] - 3.0*x[2] - 3.0*x[10] - 3.0*x[5] - 3.0*x[4] - 3.0*x[8] - 3.0*x[1] - 3.0*x[7] - 3.0*x[0] + 4.0) + 1.0;
	sf[1]           = 0.75*x[1]*x[9] + x[1] - x[9];
	sf[2]           = 0.75*x[0]*x[9] + x[0] - x[9];
	sf[3]           = 0.75*x[4]*x[9] + x[4];
	sf[4]           = 0.75*x[5]*x[9] + x[5];
	sf[5]           = 0.75*x[6]*x[9] + x[6];
	sf[6]           = 0.75*x[7]*x[9] + x[7];
	sf[7]           = 0.75*x[8]*x[9] + x[8];
	sf[8]           = x[9];
	sf[9]           = x[3] + x[2] + 0.75*x[9]*(x[3] + x[2]);
	sf[10]          = -0.75*x[10]*x[9] - x[10] + 1.0;
	sf[11]          = 4.0*x[2];
	sf[12]          = 4.0*x[3];
	sf[13]          = x[0];
	sf[14]          = x[1];
	sf[15]          = 4.0*x[3] + 4.0*x[2] + x[1] + x[0];
	sf[16]          = x[10];
	sf[17]          = 1.0 - x[10];

	mu[0]         = R*T*creal(clog( sf[0]*1.0/sf[10]*pow(sf[17], 2.0))) 					+ gb[0] + mu_Gex[0];
	mu[1]         = R*T*creal(clog( 1.0/sf[10]*sf[14]*1.0/sf[15]*pow(sf[17], 2.0)*sf[1])) 	+ gb[1] + mu_Gex[1];
	mu[2]         = R*T*creal(clog( 1.0/sf[10]*sf[13]*1.0/sf[15]*pow(sf[17], 2.0)*sf[2])) 	+ gb[2] + mu_Gex[2];
	mu[3]         = R*T*creal(clog( 1.0/sf[10]*pow(sf[11], 4.0)* (1./pow(sf[15], 4.0))*pow(sf[17], 2.0)*sf[9])) + gb[3] + mu_Gex[3];
	mu[4]         = R*T*creal(clog( 1.0/sf[10]*pow(sf[12], 4.0)* (1./pow(sf[15], 4.0))*pow(sf[17], 2.0)*sf[9])) + gb[4] + mu_Gex[4];
	mu[5]         = R*T*creal(clog( 1.0/sf[10]*pow(sf[17], 2.0)*sf[3])) 					+ gb[5] + mu_Gex[5];
	mu[6]         = R*T*creal(clog( 1.0/sf[10]*pow(sf[17], 2.0)*sf[4] + d_em[6])) 			+ gb[6] + mu_Gex[6];
	mu[7]         = R*T*creal(clog( 1.0/sf[10]*pow(sf[17], 2.0)*sf[5] + d_em[7])) 			+ gb[7] + mu_Gex[7];
	mu[8]         = R*T*creal(clog( 1.0/sf[10]*pow(sf[17], 2.0)*sf[6])) 					+ gb[8] + mu_Gex[8];
	mu[9]         = R*T*creal(clog( 1.0/sf[10]*pow(sf[17], 2.0)*sf[7])) 					+ gb[9] + mu_Gex[9];
	mu[10]        = R*T*creal(clog( 1.0/sf[10]*pow(sf[17], 2.0)*sf[8])) 					+ gb[10] + mu_Gex[10];
	mu[11]        = R*T*creal(clog( pow(sf[16], 2.0))) 										+ gb[11] + mu_Gex[11];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_liq(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
};


/** 
  objective function of olivine
*/
double obj_ig_ol(unsigned  n, const double *x, double *grad, void *SS_ref_db) {

	SS_ref *d  = (SS_ref *) SS_ref_db;


	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;

	px_ig_ol(SS_ref_db,x);

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
				it += 1;
			}
		}
	}
    sf[0]          =  x[2] - x[0] + 1.0;
    sf[1]          = -x[2] + x[0];
    sf[2]          =  x[1]*x[0] - x[1] - x[2] - x[0] + 1.0;
    sf[3]          = -x[1]*x[0] + x[2] + x[0];
    sf[4]          =  x[1];
    
	mu[0]          = R*T*creal(clog(sf[0]*sf[4])) + gb[0] + mu_Gex[0];
	mu[1]          = R*T*creal(clog(sf[1]*sf[3])) + gb[1] + mu_Gex[1];
	mu[2]          = R*T*creal(clog(sf[0]*sf[2])) + gb[2] + mu_Gex[2];
	mu[3]          = R*T*creal(clog(sf[0]*sf[3])) + gb[3] + mu_Gex[3];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_ol(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
};

/** 
  objective function of orthopyroxene
*/
double obj_ig_opx(unsigned n, const double *x, double *grad, void *SS_ref_db) {

	SS_ref *d  = (SS_ref *) SS_ref_db;


	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;
    double *d_em   = d->d_em;

	px_ig_opx(SS_ref_db,x);

	d->sum_v = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->sum_v += d->p[i]*d->v[i];
	}
	for (int i = 0; i < d->n_em; i++){
		d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
	}

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
				it += 1;
			}
		}
	}
	
    sf[0]           = x[7]*x[3] + x[7]*x[0] - x[7] - x[3]*x[5] + x[3]*x[1] - x[3] - x[5]*x[0] + x[5] + x[0]*x[1] - x[0] - x[1] + 1.0;
    sf[1]           = -x[7]*x[3] - x[7]*x[0] + x[3]*x[5] - x[3]*x[1] + x[3] + x[5]*x[0] - x[0]*x[1] + x[0];
    sf[2]           = -x[6] - x[4] + x[7] - 2.0*x[5] + x[1];
    sf[3]           = x[4];
    sf[4]           = x[6];
    sf[5]           = x[5];
    sf[6]           = x[2]*x[0] - x[2] - x[7]*x[3] + x[7]*x[0] - x[7] + x[3]*x[5] - x[3]*x[1] + x[3] - x[0] + 1.0;
    sf[7]           = -x[2]*x[0] + x[7]*x[3] - x[7]*x[0] - x[3]*x[5] + x[3]*x[1] - x[3] + x[0];
    sf[8]           = x[2];
    sf[9]           = x[7];
    sf[10]           = 1.0 - 0.5*x[1];
    sf[11]           = 0.5*x[1];

	mu[0]          = R*T*creal(clog(sf[0]*csqrt(sf[10])*sf[6])) + gb[0] + mu_Gex[0];
	mu[1]          = R*T*creal(clog(csqrt(sf[10])*sf[1]*sf[7])) + gb[1] + mu_Gex[1];
	mu[2]          = R*T*creal(clog(sf[0]*csqrt(sf[10])*sf[7])) + gb[2] + mu_Gex[2];
	mu[3]          = R*T*creal(clog(sf[0]*csqrt(sf[10])*sf[8])) + gb[3] + mu_Gex[3];
	mu[4]          = R*T*creal(clog(1.4142*cpow(sf[10], 0.25)*cpow(sf[11], 0.25)*sf[2]*sf[6])) + gb[4] + mu_Gex[4];
	mu[5]          = R*T*creal(clog(1.4142*cpow(sf[10], 0.25)*cpow(sf[11], 0.25)*sf[4]*sf[6] + d_em[5])) + gb[5] + mu_Gex[5];
	mu[6]          = R*T*creal(clog(2.8284*csqrt(sf[0])*cpow(sf[10], 0.25)*cpow(sf[11], 0.25)*csqrt(sf[5])*sf[6])) + gb[6] + mu_Gex[6];
	mu[7]          = R*T*creal(clog(1.4142*cpow(sf[10], 0.25)*cpow(sf[11], 0.25)*sf[3]*sf[6] + d_em[7])) + gb[7] + mu_Gex[7];
	mu[8]          = R*T*creal(clog(csqrt(sf[10])*sf[2]*sf[9])) + gb[8] + mu_Gex[8];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_opx(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
};


/** 
  objective function of plagioclase 4T
*/
double obj_ig_pl4T(unsigned  n, const double *x, double *grad, void *SS_ref_db) {

	SS_ref *d  = (SS_ref *) SS_ref_db;


	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;

	px_ig_pl4T(SS_ref_db,x);

	d->sum_v = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->sum_v += d->p[i]*d->v[i];
	}
	for (int i = 0; i < d->n_em; i++){
		d->mat_phi[i] = (d->p[i]*d->v[i])/d->sum_v;
	}

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->mat_phi[j])*(d->eye[i][k] - d->mat_phi[k])*(d->W[it]*2.0*d->v[i]/(d->v[j]+d->v[k]));
				it += 1;
			}
		}
	}
	
    sf[0]           = -x[0] - x[1] + 1.0;
    sf[1]           = x[0];
    sf[2]           = x[1];
    sf[3]           = 0.25*x[0] + 0.25;
    sf[4]           = 0.75 - 0.25*x[0];

	mu[0]          = R*T*creal(clog(1.7548*sf[0]*cpow(sf[3], 0.25)*cpow(sf[4], 0.75))) 	+ gb[0] + mu_Gex[0];
	mu[1]          = R*T*creal(clog(2.0*sf[1]*csqrt(sf[3])*csqrt(sf[4]))) 				+ gb[1] + mu_Gex[1];
	mu[2]          = R*T*creal(clog(1.7548*sf[2]*cpow(sf[3], 0.25)*cpow(sf[4], 0.75))) 	+ gb[2] + mu_Gex[2];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_pl4T(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
};

/** 
  objective function of spinel
*/
double obj_ig_spn(unsigned n, const double *x, double *grad, void *SS_ref_db) {
 
	SS_ref *d  = (SS_ref *) SS_ref_db;

	int n_em   = d->n_em;
	double P   = d->P;
	double T   = d->T;
	double R   = d->R;

	double *gb     = d->gb_lvl;
	double *mat_phi= d->mat_phi;
	double *mu_Gex = d->mu_Gex;
	double *sf     = d->sf;
	double *mu     = d->mu;
    double *d_em   = d->d_em;

	px_ig_spn(SS_ref_db,x);

	for (int i = 0; i < d->n_em; i++){
		mu_Gex[i] = 0.0;
		int it = 0;
		for (int j = 0; j < d->n_xeos; j++){
			for (int k = j+1; k < d->n_em; k++){
				mu_Gex[i] -= (d->eye[i][j] - d->p[j])*(d->eye[i][k] - d->p[k])*(d->W[it]);
				it += 1;
			}
		}
	}
	
    sf[0]           = 2.0*x[4]/3.0 -x[3]*x[0]/3.0 +x[3]/3.0 -x[0]/3.0 + 1.0/3.0;
    sf[1]           = 2.0*x[5]/3.0 +x[3]*x[0]/3.0 +x[0]/3.0;
    sf[2]           = -2.0*x[4]/3.0 - 2.0*x[5]/3.0 - 2.0*x[6]/3.0 + 2.0*x[2]*x[1]/3.0 + 2.0*x[3]*x[1]/3.0 - x[3]/3.0 - 2.0*x[1]/3.0 + 2.0/3.0;
    sf[3]           = 2.0*x[6]/3.0 - 2.0*x[2]*x[1]/3.0 - 2.0*x[3]*x[1]/3.0 + 2.0*x[1]/3.0;
    sf[4]           = -x[4]/3.0 -x[3]*x[0]/3.0 +x[3]/3.0 -x[0]/3.0 + 1.0/3.0;
    sf[5]           = -x[5]/3.0 +x[3]*x[0]/3.0 +x[0]/3.0;
    sf[6]           = x[4]/3.0 +x[5]/3.0 +x[6]/3.0 + 2.0*x[2]*x[1]/3.0 -x[2] + 2.0*x[3]*x[1]/3.0 - 5.0*x[3]/6.0 - 2.0*x[1]/3.0 + 2.0/3.0;
    sf[7]           = -x[6]/3.0 - 2.0*x[2]*x[1]/3.0 - 2.0*x[3]*x[1]/3.0 + 2.0*x[1]/3.0;
    sf[8]           = x[2];
    sf[9]           = 0.5*x[3];

	mu[0]          = R*T*creal(clog(sf[0]*sf[6])) + gb[0] + mu_Gex[0];
	mu[1]          = R*T*creal(clog(2.0*sf[2]*csqrt(sf[4])*csqrt(sf[6]))) + gb[1] + mu_Gex[1];
	mu[2]          = R*T*creal(clog(sf[1]*sf[6])) + gb[2] + mu_Gex[2];
	mu[3]          = R*T*creal(clog(2.0*sf[2]*csqrt(sf[5])*csqrt(sf[6]))) + gb[3] + mu_Gex[3];
	mu[4]          = R*T*creal(clog(sf[1]*sf[7] + d_em[4])) + gb[4] + mu_Gex[4];
	mu[5]          = R*T*creal(clog(2.0*sf[3]*csqrt(sf[5])*csqrt(sf[7]) + d_em[5])) + gb[5] + mu_Gex[5];
	mu[6]          = R*T*creal(clog(sf[0]*sf[8] + d_em[6])) + gb[6] + mu_Gex[6];
	mu[7]          = R*T*creal(clog(2.0*sf[0]*csqrt(sf[4])*csqrt(sf[9]))) + gb[7] + mu_Gex[7];

	d->sum_apep = 0.0;
	for (int i = 0; i < n_em; i++){
	   d->sum_apep += d->ape[i]*d->p[i];
	}
	d->factor = d->fbc/d->sum_apep;

	d->df_raw = 0.0;
	for (int i = 0; i < d->n_em; i++){
		d->df_raw += mu[i]*d->p[i];
	}
	d->df = d->df_raw * d->factor;
	
	if (grad){
	double *dfx    = d->dfx;
	double **dp_dx = d->dp_dx;
		dpdx_ig_spn(SS_ref_db,x);
		for (int i = 0; i < (d->n_xeos); i++){
		   dfx[i] = 0.0;
		   for (int j = 0; j < n_em; j++){
			   dfx[i] += (mu[j] - (d->ape[j]/d->sum_apep)*d->df_raw)*d->factor*dp_dx[j][i];
		   }
		   grad[i] = creal(dfx[i]);
		}
	}

	return d->df;
}

SS_ref P2X(					global_variable 	 gv,
							SS_ref 				 SS_ref_db, 
							bulk_info 			 z_b,
							char    			*name				){

	double eps = gv.bnd_val;

	/* Associate the right solid-solution data */
     if(gv.EM_database == 0){
		if 	(strcmp( name, "liq") == 0 ){
			p2x_mp_liq(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "pl4tr") == 0){
			p2x_mp_pl4tr(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "bi") == 0){
			p2x_mp_bi(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "g") == 0){
			p2x_mp_g(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "ep") == 0){
			p2x_mp_ep(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "ma") == 0){
			p2x_mp_ma(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "mu") == 0){
			p2x_mp_mu(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "opx") == 0){
			p2x_mp_opx(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "sa") == 0){
			p2x_mp_sa(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "cd") == 0){
			p2x_mp_cd(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "st") == 0){
			p2x_mp_st(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "chl") == 0){
			p2x_mp_chl(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "ctd") == 0){
			p2x_mp_ctd(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "sp") == 0){
			p2x_mp_sp(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "ilm") == 0){
			p2x_mp_ilm(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "mt") == 0){
			p2x_mp_mt(&SS_ref_db, eps);	
		}
		else{
			printf("\nsolid solution '%s' is not in the database\n",name);		
		}	
     }
     else if (gv.EM_database == 1){
      if (strcmp( name, "liq") == 0 ){
         p2x_mb_liq(&SS_ref_db, eps);
      }
      else if (strcmp( name, "hb") == 0 ){
         p2x_mb_hb(&SS_ref_db, eps);
      }
      else if (strcmp( name, "aug") == 0 ){
         p2x_mb_aug(&SS_ref_db, eps);
      }
      else if (strcmp( name, "dio") == 0 ){
         p2x_mb_dio(&SS_ref_db, eps);
      }
      else if (strcmp( name, "opx") == 0 ){
         p2x_mb_opx(&SS_ref_db, eps);
      }
      else if (strcmp( name, "g") == 0 ){
         p2x_mb_g(&SS_ref_db, eps);
      }
      else if (strcmp( name, "ol") == 0 ){
         p2x_mb_ol(&SS_ref_db, eps);
      }
      else if (strcmp( name, "pl4tr") == 0 ){
         p2x_mb_pl4tr(&SS_ref_db, eps);
      }
      else if (strcmp( name, "abc") == 0 ){
         p2x_mb_abc(&SS_ref_db, eps);
      }
      else if (strcmp( name, "k4tr") == 0 ){
         p2x_mb_k4tr(&SS_ref_db, eps);
      }
      else if (strcmp( name, "sp") == 0 ){
         p2x_mb_sp(&SS_ref_db, eps);
      }
      else if (strcmp( name, "ilm") == 0 ){
         p2x_mb_ilm(&SS_ref_db, eps);
      }
      else if (strcmp( name, "ilmm") == 0 ){
         p2x_mb_ilmm(&SS_ref_db, eps);
      }
      else if (strcmp( name, "ep") == 0 ){
         p2x_mb_ep(&SS_ref_db, eps);
      }
      else if (strcmp( name, "bi") == 0 ){
         p2x_mb_bi(&SS_ref_db, eps);
      }
      else if (strcmp( name, "mu") == 0 ){
         p2x_mb_mu(&SS_ref_db, eps);
      }
      else if (strcmp( name, "chl") == 0 ){
         p2x_mb_chl(&SS_ref_db, eps);
      }
      else{
         printf("\nsolid solution '%s' is not in the database\n",name);
      }
   }
	else if(gv.EM_database == 2){
		if 	(strcmp( name, "bi") == 0 ){
			p2x_ig_bi(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "fper")  == 0){
			p2x_ig_fper(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "cd")  == 0){
			p2x_ig_cd(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "cpx") == 0){
			p2x_ig_cpx(&SS_ref_db, eps);
		}	
		else if (strcmp( name, "ep")  == 0){
			p2x_ig_ep(&SS_ref_db, eps);
		}
		else if (strcmp( name, "fl")  == 0){
			p2x_ig_fl(&SS_ref_db, eps);
		}		
		else if (strcmp( name, "g")   == 0){
			p2x_ig_g(&SS_ref_db, eps);
		}
		else if (strcmp( name, "hb")  == 0){
			p2x_ig_hb(&SS_ref_db, eps);
		}	
		else if (strcmp( name, "ilm") == 0){
			p2x_ig_ilm(&SS_ref_db, eps);
		}
		else if (strcmp( name, "liq") == 0){
			p2x_ig_liq(&SS_ref_db, eps);
		}	
		else if (strcmp( name, "ol")  == 0){
			p2x_ig_ol(&SS_ref_db, eps);
		}
		else if (strcmp( name, "opx") == 0){
			p2x_ig_opx(&SS_ref_db, eps);
		}
		else if (strcmp( name, "pl4T")  == 0){
			p2x_ig_pl4T(&SS_ref_db, eps);
		}	
		else if (strcmp( name, "spn") == 0){
			p2x_ig_spn(&SS_ref_db, eps);	
		}
		else{
			printf("\nsolid solution '%s' is not in the database\n",name);		
		}	
	}
	else if(gv.EM_database == 4){
		if 	(strcmp( name, "fluid") == 0 ){
			p2x_um_fluid(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "ol")  == 0){
			p2x_um_ol(&SS_ref_db, eps);	
		}
		else if (strcmp( name, "br") == 0){
			p2x_um_br(&SS_ref_db, eps);
		}	
		else if (strcmp( name, "ch")  == 0){
			p2x_um_ch(&SS_ref_db, eps);
		}
		else if (strcmp( name, "atg")  == 0){
			p2x_um_atg(&SS_ref_db, eps);
		}		
		else if (strcmp( name, "g")   == 0){
			p2x_um_g(&SS_ref_db, eps);
		}
		else if (strcmp( name, "ta")  == 0){
			p2x_um_ta(&SS_ref_db, eps);
		}	
		else if (strcmp( name, "chl") == 0){
			p2x_um_chl(&SS_ref_db, eps);
		}
		else if (strcmp( name, "anth") == 0){
			p2x_um_anth(&SS_ref_db, eps);
		}
		else if (strcmp( name, "spi")  == 0){
			p2x_um_spi(&SS_ref_db, eps);	
		}	
		else if (strcmp( name, "opx")  == 0){
			p2x_um_opx(&SS_ref_db, eps);
		}
		else if (strcmp( name, "po") == 0){
			p2x_um_po(&SS_ref_db, eps);
		}
		else{
			printf("\nsolid solution '%s' is not in the database\n",name);		
		}	
	}

	return SS_ref_db;
};

SS_ref PC_function(		global_variable 	 gv,
						SS_ref 				 SS_ref_db, 
						bulk_info 	 		 z_b,
						char    			*name				){

	double G0 = 0.0;

	/* Associate the right solid-solution data */
    if(gv.EM_database == 0){
		if 	(strcmp( name, "liq") == 0 ){
			G0 = obj_mp_liq(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "pl4tr") == 0){	
			G0 = obj_mp_pl4tr(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "bi") == 0){	
			G0 = obj_mp_bi(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "g") == 0){	
			G0 = obj_mp_g(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "ep") == 0){	
			G0 = obj_mp_ep(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "ma") == 0){	
			G0 = obj_mp_ma(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "mu") == 0){	
			G0 = obj_mp_mu(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "opx") == 0){	
			G0 = obj_mp_opx(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "sa") == 0){	
			G0 = obj_mp_sa(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "cd") == 0){	
			G0 = obj_mp_cd(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "st") == 0){	
			G0 = obj_mp_st(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "chl") == 0){	
			G0 = obj_mp_chl(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "ctd") == 0){	
			G0 = obj_mp_ctd(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "sp") == 0){	
			G0 = obj_mp_sp(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "ilm") == 0){	
			G0 = obj_mp_ilm(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "mt") == 0){	
			G0 = obj_mp_mt(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else{
			printf("\nsolid solution '%s' is not in the database\n",name);		
		}	
	}
    else if (gv.EM_database == 1){
      if (strcmp( name,"liq") == 0 ){
         G0 = obj_mb_liq(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "hb") == 0 ){
         G0 = obj_mb_hb(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "aug") == 0 ){
         G0 = obj_mb_aug(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "dio") == 0 ){
         G0 = obj_mb_dio(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "opx") == 0 ){
         G0 = obj_mb_opx(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "g") == 0 ){
         G0 = obj_mb_g(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "ol") == 0 ){
         G0 = obj_mb_ol(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "pl4tr") == 0 ){
         G0 = obj_mb_pl4tr(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "abc") == 0 ){
         G0 = obj_mb_abc(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "k4tr") == 0 ){
         G0 = obj_mb_k4tr(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "sp") == 0 ){
         G0 = obj_mb_sp(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "ilm") == 0 ){
         G0 = obj_mb_ilm(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "ilmm") == 0 ){
         G0 = obj_mb_ilmm(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "ep") == 0 ){
         G0 = obj_mb_ep(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "bi") == 0 ){
         G0 = obj_mb_bi(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "mu") == 0 ){
         G0 = obj_mb_mu(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else if (strcmp( name, "chl") == 0 ){
         G0 = obj_mb_chl(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
      }
      else{
         printf("\nsolid solution '%s' is not in the database\n",name);
      }
   }
	else if(gv.EM_database == 2){
		if 	(strcmp( name, "bi") == 0 ){
			G0 = obj_ig_bi(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "fper")  == 0){
			G0 = obj_ig_fper(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "cd")  == 0){
			G0 = obj_ig_cd(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "cpx") == 0){	
			G0 = obj_ig_cpx(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
				}	
		else if (strcmp( name, "ep")  == 0){
			G0 = obj_ig_ep(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "fl")  == 0){
			G0 = obj_ig_fl(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}		
		else if (strcmp( name, "g")   == 0){
			G0 = obj_ig_g(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "hb")  == 0){
			G0 = obj_ig_hb(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}	
		else if (strcmp( name, "ilm") == 0){
			G0 = obj_ig_ilm(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}	
		else if (strcmp( name, "liq") == 0){
			G0 = obj_ig_liq(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "ol")  == 0){
			G0 = obj_ig_ol(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "opx") == 0){
			G0 = obj_ig_opx(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "pl4T")  == 0){
			G0 = obj_ig_pl4T(SS_ref_db.n_xeos, SS_ref_db.iguess,SS_ref_db.dfx, &SS_ref_db);
		}	
		else if (strcmp( name, "spn") == 0){	
			G0 = obj_ig_spn(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else{
			printf("\nsolid solution '%s' is not in the database\n",name);		
		}	
	}
	else if(gv.EM_database == 4){
		if 	(strcmp( name, "fluid") == 0 ){
			G0 = obj_um_fluid(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "ol")  == 0){
			G0 = obj_um_ol(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "br") == 0){	
			G0 = obj_um_br(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
				}	
		else if (strcmp( name, "ch")  == 0){
			G0 = obj_um_ch(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "atg")  == 0){
			G0 = obj_um_atg(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}		
		else if (strcmp( name, "g")   == 0){
			G0 = obj_um_g(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "ta")  == 0){
			G0 = obj_um_ta(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}	
		else if (strcmp( name, "chl") == 0){
			G0 = obj_um_chl(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}	
		else if (strcmp( name, "anth") == 0){
			G0 = obj_um_anth(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "spi")  == 0){
			G0 = obj_um_spi(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}	
		else if (strcmp( name, "opx")  == 0){
			G0 = obj_um_opx(SS_ref_db.n_xeos, SS_ref_db.iguess, 	SS_ref_db.dfx, &SS_ref_db);
		}
		else if (strcmp( name, "po") == 0){
			G0 = obj_um_po(SS_ref_db.n_xeos, SS_ref_db.iguess, SS_ref_db.dfx, &SS_ref_db);
		}

		else{
			printf("\nsolid solution '%s' is not in the database\n",name);		
		}	
	}

	/** get driving force for simplex pseudocompounds */
	SS_ref_db.df = G0;
	
	/** initialize composition */
	for (int j = 0; j < gv.len_ox; j++){
	   SS_ref_db.ss_comp[j] = 0.0;
	}
	
	/* set mu = 0 for absent oxides */
	for (int j = 0; j < SS_ref_db.n_em; j++){
	   SS_ref_db.mu[j] *= SS_ref_db.z_em[j];
	} 
	
	/* find solution phase composition*/
	for (int i = 0; i < SS_ref_db.n_em; i++){
	   for (int j = 0; j < gv.len_ox; j++){
		   SS_ref_db.ss_comp[j] += SS_ref_db.Comp[i][j]*SS_ref_db.p[i]*SS_ref_db.z_em[i];
	   } 
	}
	
	/* check if site fractions are satisfied */
	SS_ref_db.sf_ok = 1;
	for (int i = 0; i < SS_ref_db.n_sf; i++){
		if (SS_ref_db.sf[i] < gv.eps_sf_pc || isnan(SS_ref_db.sf[i]) == 1|| isinf(SS_ref_db.sf[i]) == 1){
			SS_ref_db.sf_ok = 0;	
			break;
		}
	}

	return SS_ref_db;
};

int get_phase_id(		global_variable 	 gv,
						char    			*name				){
	int id = 0;
	for (int i = 0; i < gv.len_ss; i++){
		if 	(strcmp( name, gv.SS_list[i]) == 0 ){
			id = i;
			break;
		}
	}

	return id;
};

