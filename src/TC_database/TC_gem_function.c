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
/**
  Function to calculate chemical potential of endmembers/pure phases for thermocalc database
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <complex.h> 

#include "../MAGEMin.h"
#include "../toolkit.h"
#include "TC_gem_function.h"


#define eps 1e-8
#define max_ox 15

typedef struct Helmholtz_WP {
    double R;
    double no[9];
    double gammao[5];
    double c[55];
    double d[55];
    double t[55];
    double n[57];
    double alpha[3];
    double beta[3];
    double gamma[3];
    double epsilon[3];

    double a[2];
    double b[2];
    double A[2];
    double B[2];
    double C[2];
    double F[2];
    double E[2];

    // output
    double helmholtz;   
    double helmholtzD;   
    double helmholtzDD;  

} HelmholtzWP;

HelmholtzWP helm_WP = {
    461.51805,
    {0.0,-8.32044648201, 6.6832105268, 3.00632, 0.012436, 0.97315, 1.27950, 0.96956, 0.24873},
    {1.28728967, 3.53734222, 7.74073708, 9.24437796, 27.5075105},

    {0.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 6.0, 6.0, 6.0, 6.0, 0.0, 0.0, 0},
    {0.0,1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 4.0, 4.0, 5.0, 7.0, 9.0, 10.0, 11.0, 13.0, 15.0, 1.0, 2.0, 2.0, 2.0, 3.0, 4.0, 4.0, 4.0, 5.0, 6.0, 6.0, 7.0, 9.0, 9.0, 9.0, 9.0, 9.0, 10.0, 10.0, 12.0, 3.0, 4.0, 4.0, 5.0, 14.0, 3.0, 6.0, 6.0, 6.0, 3.0, 3.0, 3.0},
    {0.0,-0.5, 0.875, 1.0, 0.5, 0.75, 0.375, 1.0, 4.0, 6.0, 12.0, 1.0, 5.0, 4.0, 2.0, 13.0, 9.0, 3.0, 4.0, 11.0, 4.0, 13.0, 1.0, 7.0, 1.0, 9.0, 10.0, 10.0, 3.0, 7.0, 10.0, 10.0, 6.0, 10.0, 10.0, 1.0, 2.0, 3.0, 4.0, 8.0, 6.0, 9.0, 8.0, 16.0, 22.0, 23.0, 23.0, 10.0, 50.0, 44.0, 46.0, 50.0, 0.0, 1.0, 4.0},
    {0.0,0.12533547935523e-01,0.78957634722828e+01,-0.87803203303561e+01, 0.31802509345418,-0.26145533859358,-0.78199751687981e-02, 0.88089493102134e-02,-0.66856572307965, 0.20433810950965,-0.66212605039687e-04,-0.19232721156002,-0.25709043003438, 0.16074868486251,-0.40092828925807e-01, 0.39343422603254e-06,-0.75941377088144e-05, 0.56250979351888e-03,-0.15608652257135e-04, 0.11537996422951e-08, 0.36582165144204e-06,-0.13251180074668e-11,-0.62639586912454e-09,-0.10793600908932, 0.17611491008752e-01, 0.22132295167546,-0.40247669763528, 0.58083399985759, 0.49969146990806e-02,-0.31358700712549e-01,-0.74315929710341, 0.47807329915480, 0.20527940895948e-01,-0.13636435110343, 0.14180634400617e-01, 0.83326504880713e-02,-0.29052336009585e-01, 0.38615085574206e-01,-0.20393486513704e-01,-0.16554050063734e-02, 0.19955571979541e-02, 0.15870308324157e-03,-0.16388568342530e-04, 0.43613615723811e-01, 0.34994005463765e-01,-0.76788197844621e-01, 0.22446277332006e-01,-0.62689710414685e-04,-0.55711118565645e-09,-0.19905718354408, 0.31777497330738,-0.11841182425981,-0.31306260323435e+02,0.31546140237781e+02,-0.25213154341695e+04,-0.14874640856724,0.31806110878444},
    {20.0,20.0,20.0} ,
    {150.0, 150.0, 150.0},
    {1.21,1.21,1.25},
    {1.0,1.0,1.0},
    {3.5,3.5},
    {0.85,0.95},
    {0.32,0.32},
    {0.2,0.2},
    {28.0,32.0},
    {700.0,800.0},
    {0.3,0.3},

    0.0,
    0.0,
    0.0
};

typedef struct Helmholtz_HGK {
    //  constant
    double refT;
    double refrho;
    double refV;
    double refP;
    double refVis;
    double refCond;
    double refSurfTension;
    double refF;
    double refS;
    double refSoundSpeed;

    double A0[18];
    double A1[5];
    double A20;
    double yc[4];
    double z0;
    double ki[36];
    double li[36];
    double A3[36];
    double mi[4];
    double ni[4];
    double alpha[4];
    double beta[4];

    double ri[4];
    double ti[4];
    double A4[4];

    // output
    double helmholtz;   
    double helmholtzD;   
    double helmholtzDD;  

} HelmholtzHGK;

HelmholtzHGK helm_HGK = {
    647.27,
    317.763,
    1.0/317.763,
    22.115e+06,
    55.071e-06,
    0.49450,
    235.8e-03,
    69595.89,
    107.5222,
    263.810,
    {-0.130840393653E+2,-0.857020420940E+2,0.765192919131E-2,-0.620600116069E+0,-0.106924329402E+2,-0.280671377296E+1,0.119843634845E+3,-0.823907389256E+2,0.555864146443E+2,-0.310698122980E+2,0.136200239305E+2,-0.457116129409E+1,0.115382128188E+1,-0.214242224683E+0,0.282800597384E-1,-0.250384152737E-2,0.132952679669E-3,-0.319277411208E-5},
    {0.15383053E+1,-0.81048367E+0,-0.68305748E+1,0.00000000E+0,0.86756271E+0},
    0.42923415E+1,
    {0.59402227E-1,-0.28128238E-1,0.56826674E-3,-0.27987451E-3},
    0.317763E+0,
    {1.0,1.0,1.0,1.0,2.0,2.0,2.0,2.0,3.0,3.0,3.0,3.0,4.0,4.0,4.0,4.0,5.0,5.0,5.0,5.0,6.0,6.0,6.0,6.0,7.0,7.0,7.0,7.0,9.0,9.0,9.0,9.0,3.0,3.0,1.0,5.0},
    {1.0,2.0,4.0,6.0,1.0,2.0,4.0,6.0,1.0,2.0,4.0,6.0,1.0,2.0,4.0,6.0,1.0,2.0,4.0,6.0,1.0,2.0,4.0,6.0,1.0,2.0,4.0,6.0,1.0,2.0,4.0,6.0,0.0,3.0,3.0,3.0},
    {-0.76221190138079E+1,0.32661493707555E+2,0.11305763156821E+2,-0.10015404767712E+1,0.12830064355028E+3,-0.28371416789846E+3,0.24256279839182E+3,-0.99357645626725E+2,-0.12275453013171E+4,0.23077622506234E+4,-0.16352219929859E+4,0.58436648297764E+3,0.42365441415641E+4,-0.78027526961828E+4,0.38855645739589E+4,-0.91225112529381E+3,-0.90143895703666E+4,0.15196214817734E+5,-0.39616651358508E+4,-0.72027511617558E+3,0.11147126705990E+5,-0.17412065252210E+5,0.99918281207782E+3,0.33504807153854E+4,-0.64752644922631E+4,0.98323730907847E+4,0.83877854108422E+3,-0.27919349903103E+4,0.11112410081192E+4,-0.17287587261807E+4,-0.36233262795423E+3,0.61139429010144E+3,0.32968064728562E+2,0.10411239605066E+3,-0.38225874712590E+2,-0.20307478607599E+3},
    {2.0,2.0,2.0,4.0},
    {0.0,2.0,0.0,0.0},
    {34.0,40.0,30.0,1050.0},
    {20000.0,20000.0,40000.0,25.0},
    {0.10038928E+1, 0.10038928E+1, 0.10038928E+1, 0.48778492E+1},
    {0.98876821E+0, 0.98876821E+0, 0.99124013E+0, 0.41713659E+0},
    {-0.32329494E-2, -0.24139355E-1, 0.79027651E-3, -0.13362857E+1},

    0.0,
    0.0,
    0.0
};

void HelmholtzHGK_calc( HelmholtzHGK    *HGK,
                        double           TK,
                        double           D        ){
    HelmholtzHGK *in  = (HelmholtzHGK *) HGK;

    //Dimensionless temperature and density
    double t = TK/in->refT;
    double d = D/in->refrho;

    /**--------------------------------------------------------------------------
     HGK0 - Auxilliary 0
    --------------------------------------------------------------------------*/

    double aux0helmholtz    = (in->A0[0] + in->A0[1] * t) * log(t);
    double aux = 0.0;
    for (int i = 2; i < 18; i++){
        aux0helmholtz    +=   in->A0[i] * pow(t,((double)i-4.0));
    }

    /**--------------------------------------------------------------------------
     HGK1 - Auxilliary 1
    --------------------------------------------------------------------------*/
    double aux1helmholtz    = 0.0;

    for (int i = 0; i < 5; i++){
        aux1helmholtz    += + d * in->A1[i] * pow(t,(1.0-(double)i));
    }

    double aux1helmholtzD   = aux1helmholtz/d;

    /**--------------------------------------------------------------------------
     HGK2 - Auxilliary 2
    --------------------------------------------------------------------------*/
    double t3   = pow(t,-3.0);
    double t5   = pow(t,-5.0);

    double y     = d * (in->yc[0] + in->yc[1]*log(t) + in->yc[2]*t3 + in->yc[3]*t5);
    double y_r   = y/d;
    double y_rr  = 0.0;

    double x    = 1.0/(1.0 - y);
    double x2   = x * x;
    double x_r  = y_r * x2;
    double x_rr = y_rr * x2 + 2.0 * y_r * x_r * x;

    double u     = log(d * x);
    double u_r   = x_r/x + 1.0/d;
    double u_rr  = x_rr/x - x_r*x_r/(x*x) - 1.0/(d*d);

    double c1 = -130.0/3.0;
    double c2 =  169.0/6.0;
    double c3 = -14.0;

    double aux2helmholtz    = in->A20 * t * (u + c1*x  + c2*x*x + c3*y);
    double aux2helmholtzD   = in->A20 * t * (u_r + c1*x_r + 2.0*c2*x*x_r + c3*y_r);
    double aux2helmholtzDD  = in->A20 * t * (u_rr + c1*x_rr + 2.0*c2*(x*x_rr + x_r*x_r) + c3*y_rr);

    /**--------------------------------------------------------------------------
     HGK3 - Auxilliary 3
    --------------------------------------------------------------------------*/
    double z     =  1.0 - exp(-in->z0 * d);
    double z_r   =  in->z0 * (1.0 - z);
    double z_rr  = -in->z0 * z_r;

    double aux3helmholtz    = 0.0;
    double aux3helmholtzD   = 0.0;
    double aux3helmholtzDD  = 0.0;

    double lambda     = 0.0; 
    double lambda_r   = 0.0;
    double lambda_rr  = 0.0;

    for (int i = 0; i < 36; ++i){

        lambda     =  in->A3[i] * pow(t,(-in->li[i])) * pow(z,(in->ki[i]));
        lambda_r   =  in->ki[i]*z_r*lambda/z;
        lambda_rr  =  lambda_r*(z_rr/z_r + lambda_r/lambda - z_r/z);

        aux3helmholtz    +=  lambda;
        aux3helmholtzD   +=  lambda_r;
        aux3helmholtzDD  +=  lambda_rr;
    }

    /*--------------------------------------------------------------------------
     HGK4 - Auxilliary 4
    --------------------------------------------------------------------------*/
    double aux4helmholtz    = 0.;
    double aux4helmholtzD   = 0.;
    double aux4helmholtzDD  = 0.;

    double delta     = 0.0;
    double tau       = 0.0; 
    double delta_r   = 0.0; 

    double delta_m   = 0.0; 
    double delta_n   = 0.0; 
    
    double psi       = 0.0; 
    double psi_r     = 0.0; 
    double psi_rr    = 0.0; 
    
    double theta     = 0.0;  
    double theta_r   = 0.0;  
    double theta_rr  = 0.0;  

    for(int i = 0; i < 4; ++i){
        delta   = (d - in->ri[i])/in->ri[i];
        tau     = (t - in->ti[i])/in->ti[i];
        delta_r = 1.0/in->ri[i];

        delta_m = pow(delta,in->mi[i]);
        delta_n = pow(delta,in->ni[i]);
        
        psi    =  (in->ni[i] - in->alpha[i]*in->mi[i]*delta_m)*(delta_r/delta);
        psi_r  = -(in->ni[i] + in->alpha[i]*in->mi[i]*(in->mi[i] - 1.0)*delta_m)* pow((delta_r/delta),2.0);
        psi_rr = (2.0*in->ni[i] - in->alpha[i]*in->mi[i]*(in->mi[i] - 1.0)*(in->mi[i] - 2.0)*delta_m)*pow((delta_r/delta),3.0);
        
        theta     =  in->A4[i]*delta_n*exp(-in->alpha[i]*delta_m - in->beta[i]*tau*tau);
        theta_r   =  psi*theta;
        theta_rr  =  psi_r*theta + psi*theta_r;

        aux4helmholtz    +=  theta;
        aux4helmholtzD   +=  theta_r;
        aux4helmholtzDD  +=  theta_rr;
    }

    /*--------------------------------------------------------------------------
     FINAL WATER HELMHOLTZ STATE
    --------------------------------------------------------------------------*/
    // Assemble the contributions from each auxiliary Helmholtz state
    in->helmholtz    = aux0helmholtz    + aux1helmholtz    + aux2helmholtz    + aux3helmholtz    + aux4helmholtz;
    in->helmholtzD   =                    aux1helmholtzD   + aux2helmholtzD   + aux3helmholtzD   + aux4helmholtzD;
    in->helmholtzDD  =                                       aux2helmholtzDD  + aux3helmholtzDD  + aux4helmholtzDD;

    // Convert the Helmholtz free energy of water and its derivatives to dimensioned form
    in->helmholtz    *=  in->refF;
    in->helmholtzD   *=  in->refF/in->refrho;
    in->helmholtzDD  *=  in->refF/in->refrho/in->refrho;
}

void HelmholtzWP_calc( HelmholtzWP      *WP,
                        double           TK,
                        double           D,
                        double           Tcr,
                        double           Dcr        ){
    HelmholtzWP *in  = (HelmholtzWP *) WP;
    int j;
    double tau   = Tcr/TK;
    double delta = D/Dcr;

    double phio     =  log(delta) + in->no[1] + in->no[2]*tau + in->no[3]*log(tau);
    double phio_d   =  1.0/delta;
    double phio_dd  = -1.0/pow(delta,2.0);

    double ee = 0.0;

    for(int i = 4; i < 9; ++i){
        j = i - 4;
        ee = exp(in->gammao[j] * tau);
        phio     = phio + in->no[i] * log(1.0 - 1.0/ee);
    }

    double phir     = 0.0;
    double phir_d   = 0.0;
    double phir_dd  = 0.0;

    double xA,xA_d,xA_dd,dci;
    for(int i = 1; i < 8; ++i){

        xA     = in->n[i]*pow(delta,in->d[i])*pow(tau,in->t[i]);
        xA_d   = in->d[i]/delta * xA;
        xA_dd  = (in->d[i] - 1.0)/delta * xA_d;

        phir     +=  xA;
        phir_d   +=  xA_d;
        phir_dd  +=  xA_dd;
    }

    double     xB,xB_d,xB_dd;
    for(int i = 8; i < 52; ++i){

        dci = pow(delta,in->c[i]);

        xB     =  in->n[i]*pow(delta,in->d[i])*pow(tau,in->t[i])*exp(-dci);
        xB_d   = (in->d[i] - in->c[i]*dci)/delta * xB;
        xB_dd  = (in->d[i] - in->c[i]*dci - 1.0)/delta * xB_d - dci*pow((in->c[i]/delta),2.0) * xB;

        phir     = phir + xB;
        phir_d   = phir_d + xB_d;
        phir_dd  = phir_dd + xB_dd;
    }

    double aux1d,aux2d,xC,xC_d,xC_dd;

    for(int i = 52; i < 55; ++i){
        j = i - 52;

        aux1d = (in->d[i]/delta - 2.0*in->alpha[j]*(delta - in->epsilon[j]));
        aux2d = (in->d[i]/pow(delta,2.0) + 2.0*in->alpha[j]);

        xC     = in->n[i] * pow(delta,in->d[i]) * pow(tau,in->t[i]) * exp(-in->alpha[j]*pow((delta - in->epsilon[j]),2.0) - in->beta[j]*pow((tau - in->gamma[j]),2.0));
        xC_d   = aux1d * xC;
        xC_dd  = aux1d * xC_d - aux2d * xC;

        phir     = phir + xC;
        phir_d   = phir_d + xC_d;
        phir_dd  = phir_dd + xC_dd;
    }

    double dd,tt,theta,theta_d,theta_dd,psi,psi_d,psi_dd;
    double Delta,Delta_d,Delta_dd;
    double DeltaPow,DeltaPow_d,DeltaPow_dd;
    double xD,xD_d,xD_dd;
        
    for(int i = 55; i < 57; ++i){
        j = i - 55;

        dd = pow((delta - 1.0),2.0);
        tt = pow((tau - 1.0),2.0);
        
        theta     = (1.0 - tau) + in->A[j]*pow(dd,(0.5/in->E[j]));
        theta_d   = (theta + tau - 1.0)/(delta - 1.0)/in->E[j];
        theta_dd  = (1.0/in->E[j] - 1.0) * theta_d/(delta - 1.0);
        
        psi     = exp(-in->C[j]*dd - in->F[j]*tt);
        psi_d   = -2.0*in->C[j]*(delta - 1.0) * psi;
        psi_dd  = -2.0*in->C[j]*(psi + (delta - 1.0) * psi_d);
        
        Delta     = theta*theta + in->B[j]*pow(dd,in->a[j]);
        Delta_d   = 2.0*(theta*theta_d + in->a[j]*(Delta - theta*theta)/(delta - 1));
        Delta_dd  = 2.0*(theta_d*theta_d + theta*theta_dd + in->a[j] * ((Delta_d - 2.0*theta*theta_d)/(delta - 1.0) - (Delta - theta*theta)/pow((delta - 1.0),2.0)));
        
        DeltaPow     =  pow(Delta,in->b[j]);
        DeltaPow_d   =  in->b[j]*Delta_d/Delta * DeltaPow;
        DeltaPow_dd  = (in->b[j]*Delta_dd/Delta + in->b[j]*(in->b[j] - 1.0)*pow((Delta_d/Delta),2.0)) * DeltaPow;
        
        xD     = in->n[i]*DeltaPow*delta*psi;
        xD_d   = in->n[i]*(DeltaPow*(psi + delta*psi_d) + DeltaPow_d*delta*psi);
        xD_dd  = in->n[i]*(DeltaPow*(2.0*psi_d + delta*psi_dd) + 2.0*DeltaPow_d*(psi + delta*psi_d) + DeltaPow_dd*delta*psi);
    
        phir     +=  xD;
        phir_d   +=  xD_d;
        phir_dd  +=  xD_dd;
    }

    double phi     = phio     + phir    ;
    double phi_d   = phio_d   + phir_d  ;
    double phi_dd  = phio_dd  + phir_dd ;

    double dD   =  1.0/Dcr;

    double phiD   = phi_d*dD;
    double phiDD  = phi_dd*dD*dD;

    /*--------------------------------------------------------------------------
     FINAL WATER HELMHOLTZ STATE
    --------------------------------------------------------------------------*/
    in->helmholtz    = in->R*TK*phi;
    in->helmholtzD   = in->R*TK*phiD;
    in->helmholtzDD  = in->R*TK*phiDD;
}


// typedef struct PropSub_datas {
//     double R;
//     double no[8];

// } PropSub_data;

// PropSub_data PS_data = {
    
//     0.0,
//     0.0
// };


// function [density, densityT, densityP, densityTT, densityTP, densityPP] = Rho(Pbar,TK,Rhocalc)
void rho_wat_calc(      solvent_prop    *wat,
                        double           Pbar,
                        double           TK,
                        char            *opt        ){

    solvent_prop *d     = (solvent_prop *) wat;

    HelmholtzWP   WP    = helm_WP;
    HelmholtzHGK  HGK   = helm_HGK;

    double Tcr = 647.096;                                       //waterCriticalTemperature;
    double Dcr = 322.0;                                         //waterCriticalDensity;
    
    if (strcmp( opt, "HGK") == 0 || strcmp( opt, "WP") == 0){
        double D        = 0.0;
        double b1       =  1.99274064;
        double b2       =  1.09965342;
        double b3       = -0.510839303;
        double b4       = -1.75493479;
        double b5       = -45.5170352;
        double b6       = -6.74694450e+05;
        double t        = 1.0 - TK/Tcr;
        double t13      = pow(t,1.0/3.0);
        double t23      = t13 * t13;
        double t53      = t13 * t23 * t23;
        double t163     = t13 * t53 * t53 * t53;
        double t433     = t163 * t163 * t53 * t * t;
        double t1103    = t433 * t433 * t163 * t53 * t;
        
        if (TK > Tcr){
            D = 0.99*Dcr;
        }   
        else{
            D = Dcr * (1.0 + b1*t13 + b2*t23 + b3*t53 + b4*t163 + b5*t433 + b6*t1103);
        }
        
        double Pcr = 22.064e6;                                      //waterCriticalPressure; - In Pa - so 220.64 bars

        double f,df;
        //Apply Newton's method to the pressure-density equation
        for (int i = 0; i < 100; i++){
            if (strcmp( opt, "HGK") == 0 ){                         //Haar, Gallagher, and Kell (1984)
                HelmholtzHGK_calc( &HGK, TK, D );
                f  = (D*D*HGK.helmholtzD - (Pbar*100000.0))/Pcr;
                df = (2.0*D*HGK.helmholtzD + D*D*HGK.helmholtzDD)/Pcr;
                if (D > f/df){
                    D = D - f/df;   
                }
                else{
                   D = (Pbar*100000.0)/(D*HGK.helmholtzD);     
                }
            }
            else{                                                   //Wagner and Pruss (1995)
                HelmholtzWP_calc( &WP, TK, D, Tcr, Dcr);

                f  = (D*D*WP.helmholtzD - (Pbar*100000.0))/Pcr;
                df = (2.0*D*WP.helmholtzD + D*D*WP.helmholtzDD)/Pcr;     
                if (D > f/df){
                    D = D - f/df;   
                }
                else{
                   D = (Pbar*100000.0)/(D*WP.helmholtzD);     
                }
            } 

            if(fabs(f) < 1e-8){
                break;
            }
        }
        
        //Set the density and its partial derivatives of the thermodynamic state of water
        d->density   =  D;
    }
}


void propSolvent_JN91_calc(     solvent_prop    *wat,
                                double           TK     ){

    solvent_prop *in  = (solvent_prop *) wat;


    double Tr = 298.15;
    double Dr = 1000.0;
    double t  = TK/Tr;
    double r  = in->density/Dr;
    double a[] = {0.0000000000e+00,0.1470333593e+02,0.2128462733e+03,-0.1154445173e+03,0.1955210915e+02,-0.8330347980e+02,0.3213240048e+02,-0.6694098645e+01,-0.3786202045e+02,0.6887359646e+02,-0.2729401652e+02};

    /*--------------------------------------------------------------------------
    CALCULATION
    --------------------------------------------------------------------------*/
    double k1 = 1.0;  // k0
    double k2 = a[1]/t;   // k1
    double k3 = a[2]/t + a[3] + a[4]*t;   // k2
    double k4 = a[5]/t + a[6]*t + a[7]*t*t;   // k3
    double k5 = a[8]/t/t + a[9]/t + a[10];   // k4

    double k[] = {k1,k2,k3,k4,k5};

    double epsilon = 0.0;
    double ri;
    for (int i = 0; i < 5; i++){
        ri = pow(r,((double)i));
        epsilon += k[i]*ri;
    }

    //--------------------------------------------------------------------------
    // FINAL BORN FUNCTIONS
    //--------------------------------------------------------------------------
    in->epsilon = epsilon;
    in->Z = -1.0/epsilon;;
}

void propSolvent_SV14_calc(     solvent_prop    *wat,
                                double           Pbar,
                                double           TK       ){

    solvent_prop *in  = (solvent_prop *) wat;
    // Constants
    double a[] = {0.0000000000, -1.576377e-03, 6.810288e-02, 7.548755e-01};
    double b[] = {0.0000000000, -8.016651e-05, -6.871618e-02, 4.747973};

    double TC = TK - 273.15;
    double density = in->density/1000.0;

    double epsilon = exp(b[1]*TC + b[2]*pow(TC,0.5) + b[3])*pow(density,(a[1]*TC + a[2]*pow(TC,0.5) + a[3]));

    in->epsilon = epsilon;
    in->Z = -1.0/epsilon;

}


void propSolvent_FE97_calc(     solvent_prop    *wat,
                                double           Pbar,
                                double           TK       ){

    solvent_prop *in  = (solvent_prop *) wat;
    //--------------------------------------------------------------------------
    // Constants
    double II[] = {1.0, 1.0, 1.0, 2.0, 3.0, 3.0, 4.0, 5.0, 6.0, 7.0, 10.0};
    double J[] = {0.25, 1.0, 2.5, 1.5, 1.5, 2.5, 2.0, 2.0, 5.0, 0.5, 10.0};
    double n[] = {0.978224486826, -0.957771379375, 0.237511794148, 0.714692244396, -0.298217036956, -0.108863472196, 0.949327488264e-1, -0.980469816509e-2, 0.165167634970e-4, 0.937359795772e-4, -0.12317921872e-9, 0.196096504426e-2};

    double Tcr = 647.096; //waterCriticalTemperature;
    double Dcr = 322.0; //waterCriticalDensity;
    double Pcr = 22.064; //waterCriticalPressure; - In MPa - Do I need to convert? 22.064e6 Pa = 220.64 bars

    //--------------------------------------------------------------------------
    // CALCULATION
    //--------------------------------------------------------------------------
    // Constants
    double k        = 1.380658e-23;
    double Na       = 6.0221367e23;
    double alfa     = 1.636e-40;
    double epsilon0 = 8.854187817e-12;
    double mu       = 6.138e-30;
    double M        = 0.018015268;

    double g = 1.0 + n[11]*(in->density/Dcr)/(pow((Tcr/228.0/(Tcr/TK)-1.0),1.2));
    double epsilon, A, B;
    for (int i = 0; i < 11; i++){
        g = g + n[i]*pow((in->density/Dcr),II[i])*pow((Tcr/TK),J[i]);
    }
    A = Na * mu*mu * in->density * g / M / epsilon0 / k / TK;
    B = Na * alfa * in->density / 3.0 / M / epsilon0;
    epsilon = (1.0+A+5.0*B + pow((9.0+2.0*A+18.0*B+A*A+10.0*A*B+9.0*B*B),0.5)) / 4.0 / (1.0-B);

    //--------------------------------------------------------------------------
    // FINAL EPSILON AND BORN FUNCTIONS
    //--------------------------------------------------------------------------
    in->epsilon = epsilon;
    in->Z = -1.0/epsilon;

}

PP_ref TC_G_EM_function(	int 		 EM_dataset, 
							int 		 len_ox,
							int         *id,
							double 		*bulk_rock, 
							double 		*apo, 
							double 		 P, 
							double 		 T, 
							char 		*name, 
							char		*state			
){
	/* Get thermodynamic data */
	EM_db EM_return;
	int i, p_id = find_EM_id(name);
	EM_return   = Access_EM_DB(p_id, EM_dataset);
	
	/* Get composition (in molar amount) */
	double composition[len_ox];
	for (i = 0; i < len_ox; i ++){
		composition[i] = EM_return.Comp[id[i]];
	}

	double t0, 	p0, 	R;
	double pth, theta, 	vv;
	double enthalpy, 	entropy, volume;
	double cpa, cpb, cpc, cpd;
	double alpha0, kappa0, kappa0p, kappa0pp, dkappa0dT;
	double cpterms, n;
	double vterm 	= 0.0;
	double ta 		= 0.0;
	double tb 		= 0.0;
	double tc 		= 0.0;
	double kbar2bar = 1e3;
	double RTlnf 	= 0.0;
	t0 				= 298.15;
	p0 				= 0.001;
	R  				= 0.0083144; 
	
	enthalpy 		= EM_return.input_1[0];
	entropy  		= EM_return.input_1[1];
	volume   		= EM_return.input_1[2];
	
	cpa      		= EM_return.input_2[0];
	cpb      		= EM_return.input_2[1];
	cpc      		= EM_return.input_2[2];
	cpd      		= EM_return.input_2[3];
	
	alpha0   		= EM_return.input_3[0];
	kappa0   		= EM_return.input_3[1];
	kappa0p  		= EM_return.input_3[2];
	kappa0pp 		= EM_return.input_3[3];	

	cpterms  		= cpa* (T - t0) +   cpb* (pow(T,2.0) - pow(t0,2.0))/2.0 - 
                                        cpc* (1.0/T - 1.0/t0) + 
                                   2.0* cpd* (pow(T,0.5) - pow(t0,0.5))     - 
							   T* (2.0* cpa* (log(pow(T,0.5)) - log(pow(t0,0.5))) 
                               + cpb* (T - t0) - 
							   cpc/2.0* (pow(T,-2.) - pow(t0,-2.0)) - 2.0* cpd* (pow(T,-0.5) - pow(t0,-0.5)));
							   
	n        		= EM_return.Comp[max_ox];
	
	char liq_tail[] = "L";
	if ( EndsWithTail(name, liq_tail) == 1 ) {
		dkappa0dT        = EM_return.input_3[4];
		pth              = 0.0;
		vv               = volume * exp(alpha0 * (T-t0));
		kappa0           = kappa0 + (dkappa0dT * (T-t0));
	}
	else {
		theta = (double)(round(10636/(entropy*1e3/n + 6.44)));
		pth   = theta* alpha0* kappa0 / (exp(theta/t0) * pow(theta/t0,2.0) / pow(exp(theta/t0) - 1.,2.0)) * (1./(exp(theta/(T)) - 1.) - 1./(exp(theta/t0) - 1.));
		vv    = volume;
	}
	/* EOS After Pitzer and Sterner, 1994 - API, The Journal of Chemical Physics */
	if (strcmp( name, "H2O") == 0 || strcmp( name, "CO2") == 0 ){

		double p_bar = 1000.*P; //in bar
		double c1, c2, c3, c4, c5, c6, c7, c8, c9, c10; 

		if (strcmp( name, "H2O") == 0){
			c1  =  0.24657688e6 / T + 0.51359951e2;
			c2  =  0.58638965e0 / T - 0.28646939e-2 + 0.31375577e-4 * T;
			c3  = -0.62783840e1 / T + 0.14791599e-1 + 0.35779579e-3 * T +  0.15432925e-7 * pow(T,2.0);
			c4  = -0.42719875e0 - 0.16325155e-4 * T;
			c5  =  0.56654978e4 / T - 0.16580167e2 + 0.76560762e-1 * T;
			c6  =  0.10917883e0;
			c7  =  0.38878656e13 / pow(T,4.0) - 0.13494878e9 / pow(T,2.0) + 0.30916564e6 / T + 0.75591105e1;
			c8  = -0.65537898e5 / T + 0.18810675e3;
			c9  = -0.14182435e14 / pow(T,4.0) + 0.18165390e9 / pow(T,2.0) - 0.19769068e6 / T - 0.23530318e2;
			c10 =  0.92093375e5 / T + 0.12246777e3;
		}
		else {	// can only be CO2
			c1  =  0.18261340e7 / T + 	0.79224365e2;
			c2  =  						0.66560660e-4 	+ 0.57152798e-5 * T + 0.30222363e-9 * pow(T,2.0);
			c3 	= 						0.59957845e-2 	+ 0.71669631e-4 * T + 0.62416103e-8 * pow(T,2.0);
			c4  = -0.13270279e1 / T +  -0.15210731e0  	+ 0.53654244e-3 * T - 0.71115142e-7 * pow(T,2.0);
			c5  =  0.12456776e0 / T +   0.49045367e1    + 0.98220560e-2 * T + 0.55962121e-5 * pow(T,2.0);
			c6  = 				     	0.75522299e0;
			c7  = -0.39344644e+12 / pow(T,4.0) + 0.90918237e8 / pow(T,2.0) + 0.42776716e6 / T - 0.22347856e2;
			c8 	=  0.40282608e3 / T +   0.11971627e3;
			c9  =  0.22995650e8 / pow(T,2.0) - 0.78971817e5 / T - 0.63376456e2;
			c10 =  0.95029765e5 / T + 0.18038071e2;
		}

		/* solve for volume at P, T */
		int    err,  k;
		double vsub, yr;
		double R1     = 83.144;
		double data[] = {R1,T,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,p_bar};
		
		double x1     = 3.0;
		double x2     = R1*T/P;
		
		double e      = 1e-14;
		int maxiter   = 500;
		int mode      = 0;               														/** Mode is used to send the right *data (see root_finding.c) */
		
		vsub          =  BrentRoots(x1,x2,data,e,mode,maxiter, &yr, &k, &err);
		
		double r      =   1.0/vsub;
		double Ares   =   R1*T*( c1*r + (1.0/(c2 + c3*r + c4*pow(r, 2.0) + c5*pow(r, 3.0) + c6*pow(r, 4.0)) - 1.0/c2) - c7/c8*(exp(-c8*r) - 1.0) - c9/c10*(exp(-c10*r) - 1.0) );
		vterm         =   (Ares + p_bar*vsub + R1*T*(log( R1*T / vsub ) - 1.0)) * 1e-4;	

        // printf("vterm: %f\n", vterm);
        // printf("Ares: %f\n", Ares);
	}
	/**
	 	here we use the CORK EOS to calculate G_O2 (see Holland & Powell, 1991) 
		Critical Temperature and pressures are taken from Holland & Powell 1991 for H2
	*/
	else if(strcmp( name, "H2") == 0){
		double Tc, Pc;
		// if(strcmp( name, "O2") == 0){
		// 	Tc     =  154.75;
		// 	Pc     =  0.05;	
		// }
		// else{ //(strcmp( name, "H2") == 0){	
			Tc     =  41.2;
			Pc     =  0.0211;	
		// }

		double a0 	  =  5.45963e-5;
		double a1     = -8.63920e-6;
		double b0     =  9.18301e-4;
		double c0     = -3.30558e-5;
		double c1 	  =  2.30524e-6;
		double d0     =  6.93054e-7;
		double d1     = -8.38293e-8;

		double a	  = a0*pow(Tc,5.0/2.0)/Pc + a1*pow(Tc,3.0/2.0)/Pc*T;
		double b      = b0*Tc/Pc;
		double c      = c0*Tc/pow(Pc,3.0/2.0) + c1/pow(Pc,3.0/2.0)*T;
		double d      = d0*Tc/pow(Pc,2.0) + d1/pow(Pc,2.0)*T;

		vterm 		  =(R*T/P + b - (a*R*sqrt(T))/((R*T+b*P)*(R*T + 2.0*b*P)) + c*sqrt(P) + d*P)/100.0; 
        // RTlnf        = R*T*log(1000.0 * P) + b*P + a/(b*sqrt(T)) * (log(R*T + b*P) - log(R*T + 2.0*P*b)) + 2.0/3.0*c*P*sqrt(P) + d/2.0*P*P;
	}
    else if(strcmp( name, "O2") == 0){
        vterm = 0.0;
    }
	else {
		ta     = (1. + kappa0p)/(1. + kappa0p + kappa0 * kappa0pp);
		tb     = (kappa0p + pow(kappa0p,2.0) - (kappa0 * kappa0pp))/(kappa0 * (1. + kappa0p));
		tc     = (1. + kappa0p + kappa0 * kappa0pp)/(kappa0p + pow(kappa0p,2.0) - kappa0 * kappa0pp);
		vterm  = vv*((P-p0)*(1.-ta)+ta*(-pow(1.+tb*(P-pth),(1.0-tc))+pow(1.0 + tb * (p0 - pth),(1.0 - tc)))/(tb* (tc - 1.)))/((1. - ta) + ta* pow(1. + tb * p0,(-tc)));
	}

    
	double gbase = (enthalpy - T*entropy + cpterms + vterm + RTlnf);	
	double landaut, smax, vmax, sfdh, sfdhv, sfw, sfwv, sfn, sffac; 
	double god, sod, q, v;
	double tc0, q20, q2;
	
	sfn 	= 0.0;
	vmax 	= 0.0;
	smax    = 0.0;
	landaut = 0.0;
	god     = 0.0;
	
	if ( EndsWithTail(name, liq_tail) != 1 ) {
		if (EM_return.input_3[4] == 0.0){
			landaut = 0.; smax = -0.001; vmax = 0.; sfdh = 0.; sfdhv = 0.;
			sfw     = 0.; sfwv = 0.;     sfn  = 0.; sffac = 0.;
		}
		else if (EM_return.input_3[4] == 1.0){
			landaut = EM_return.input_3[5]; smax = EM_return.input_3[6]; vmax = EM_return.input_3[7];
			sfdh = 0.; sfdhv = 0.; sfw = 0.; sfwv = 0.; sfn = 0.; sffac = 0.;
		}
		else if (EM_return.input_3[4] == 2.0){		
			landaut = 0.; smax = -0.001; vmax = 0.;
			sfdh = EM_return.input_3[5]; sfdhv = EM_return.input_3[6]; sfw = EM_return.input_3[7];
			sfwv = EM_return.input_3[8]; sfn   = EM_return.input_3[9]; sffac = EM_return.input_3[10];
		}
		if (sfn > 0.){
			if (strcmp( state, "ordered") == 0 ){
				god = 0.;
			}
			else if (strcmp( state, "disordered") == 0 ){
				if (sffac < 0.){
					sod = sffac * R * (log(1./(sfn+1.)) + sfn*log(sfn/(sfn+1.)))*(1./sffac-sfn)/(sfn+1.);
				}
				else {
					sod = sffac * R * (log(1./(sfn+1.)) + sfn*log(sfn/(sfn+1.)));
				}
				god = sfdh + P*sfdhv + T*sod;
			}
			else if (strcmp( state, "equilibrium") == 0 ){
				if (sffac < 0.){
					/* solve for volume at P, T */
					int    err,  k;
					double vsub, yr;
					double data[] = {sfdh,P,sfdhv,sfw,T,sfwv,sfn,R,sffac};
					double x1     = eps;
					double x2     = 1.0-eps;
					double e      = 1e-12;
					int maxiter   = 500;
					int mode      = 1;               /** Mode is used to send the right *data (see root_finding.c) */
					
					q             =  BrentRoots(x1,x2,data,e,mode,maxiter,&yr,&k,&err);
					sod           = (((1. + sfn*q)*log((1. + sfn*q)/(sfn+1.)) + sfn*(1.-q)*log(sfn*(1.-q)/(sfn+1.)) - sffac*(sfn*(1.-q)*log((1.-q)/(sfn+1.)) + sfn*(sfn+q)*log((sfn+q)/(sfn+1.)) ))/(sfn+1.));
				}
				else {
					/* Test function to define min/max */
					v = eps;
					double v1 = ( sfdh + P*sfdhv + (sfw + P*sfwv)*(2.*v - 1.) + sffac*sfn/(sfn + 1.)*R*T * log(sfn*pow(1. - v,2.0)/((1. + sfn*v)*(sfn + v))) );
					v = 1.0-eps;
					double v2 = ( sfdh + P*sfdhv + (sfw + P*sfwv)*(2.*v - 1.) + sffac*sfn/(sfn + 1.)*R*T * log(sfn*pow(1. - v,2.0)/((1. + sfn*v)*(sfn + v))) );
					
					/* solve for volume at P, T */
					double x1, x2;
					int    err,  k;
					double vsub, yr;
					double data[] = {sfdh,P,sfdhv,sfw,sfwv,sffac,sfn,R,T};
					
					if (check_sign(v1, v2) == 1) {	x1     = eps;	x2     = 1.0 - eps;	}
					else {							x1     = 0.;	x2     = 1.0 - eps;	}
					
					double e     = 1e-12;
					int maxiter  = 500;
					int mode     = 2;               /** Mode is used to send the right *data (see root_finding.c) */
					q            =  BrentRoots(x1,x2,data,e,mode,maxiter,&yr,&k,&err);
					sod          =  (sffac*((1.+sfn*q)*log((1. + sfn*q)/(sfn + 1.)) + sfn*(1. - q)*log((1. - q)/(sfn + 1.)) + sfn*(1. - q)*log(sfn*(1. - q)/(sfn + 1.)) + sfn*(sfn + q)*log((sfn + q)/(sfn + 1.))) / (sfn + 1.));
				}
				god              = sfdh + P*sfdhv + q*(sfw - sfdh + P*(sfwv - sfdhv)) - pow(q,2.0)*(sfw + P*sfwv) + R*T*sod;
			}
			else {
				printf("wrong state (HAS TO BE: ordered, disordered or equilibrium)");
			}
		}
		else if (smax > 0.0){
			tc0 = landaut;
			q20 = sqrt(1.0 - t0 / tc0);
			if (strcmp( state, "ordered") == 0 ){
				god = smax*tc0*(-(2./3.) + q20*(1.0 - pow(q20,2.)/3.)) - T*smax*(q20 - 1.0) + P*vmax*(q20 - 1.0);
			}
			else if (strcmp( state, "disordered") == 0 ){
				god = smax*tc0*q20*(1.0 - pow(q20,2.)/3.) - T*smax*q20 + P*vmax*q20;
			}
			else if (strcmp( state, "equilibrium") == 0 ){
				if (vmax == 0){
					tc  = tc0;
				}
				else{
					tc  = tc0 + P * vmax / smax;
				}
				if(T >  tc){
					q2  = 0.0;
				}
				else{
					q2  = pow((tc - T) / tc0, 0.5);
				}
				god = smax*(tc0*(q20*(1.0 - (1./3.)*pow(q20, 2.0)) + (1./3.)*pow(q2, 3.0)) - q2*tc) - T*smax*(q20 - q2) + P*vmax*q20;
			}
			else{
				printf("wrong state (HAS TO BE: ordered, disordered or equilibrium)");
			}
		}
		else{
			god = 0.0;
		}
		gbase = gbase + god;
		
 	}
 	
	/* fill structure to send back to main */
	PP_ref PP_ref_db;
	
	/* Calculate normalizing factor using bulk-rock composition */
	double factor  = 0.0;
	
	/* Calculate the number of atoms in the bulk-rock composition */
	double fbc     = 0.0;
	for (i = 0; i < len_ox; i++){
		fbc += bulk_rock[i]*apo[i];
	}
	
	/* Calculate the number of atom in the solution */
	double ape = 0.0;
	for (i = 0; i < len_ox; i++){
		ape += composition[i]*apo[i];
	}
	
	/* Calculate normalizing factor */
	factor = fbc/ape;

	strcpy(PP_ref_db.Name, name);
	for (i = 0; i < len_ox; i++){
		PP_ref_db.Comp[i] = composition[i];
	}
	PP_ref_db.gbase   =  gbase;
	PP_ref_db.factor  =  factor;
	PP_ref_db.phase_shearModulus  =  (EM_return.input_4[0]*kbar2bar + (P - p0)*(EM_return.input_4[1])*kbar2bar + (T - t0)*(EM_return.input_4[2]))/kbar2bar;

	return (PP_ref_db);
}

PP_ref G_FS_function(		int 		 len_ox,
							solvent_prop *wat,
							int         *id,
							double 		*bulk_rock, 
                            double      *ElH,
							double 		*apo, 
							double 		 P, 
							double 		 T, 
							char 		*name, 
							char		*state			
){
	solvent_prop *in = (solvent_prop *) wat;


	/* Get thermodynamic data */
	FS_db FS_return;
	int i, p_id = find_FS_id(name);
	FS_return   = Access_FS_DB(p_id);
	
	/* Get composition (in molar amount) */
	double composition[len_ox];
	for (i = 0; i < len_ox; i ++){
		composition[i] = FS_return.Comp[id[i]];
	}
	double ZPrTr 	= -0.1278034682e-1;
	double YPrTr 	= -0.5798650444e-4;
	double Pref 	= 1.0;
	double Tref 	= 298.15;
	double eta 		= 0.166027e6;
	double theta 	= 228.0;
	double psi 		= 2600.0;

	double Volume 	= FS_return.input_1[0];
	double Hr       = FS_return.input_1[2]/4.184;
	double Sr 		= FS_return.input_1[1]/4.184;
	double Gf 		= FS_return.input_1[3]/4.184;
	double a1 		= FS_return.input_2[0];
	double a2 		= FS_return.input_2[1];
	double a3 		= FS_return.input_2[2];
	double a4 		= FS_return.input_2[3];
	double c1 		= FS_return.input_2[4];
	double c2 		= FS_return.input_2[5];
	double wr 		= FS_return.input_2[6];
	double charge	= FS_return.input_3[0];


	/* compute properties of the solute */
	double TK 		= T;
	double Pbar 	= P*1000.0;
	double Z 		= in->Z;
	double TC 		= T - 273.15;
	double ag1 		= -2.037662;
	double ag2 		=  5.747000e-3;
	double ag3 		= -6.557892e-6;
	double bg1 		=  6.107361;
	double bg2 		= -1.074377e-2;
	double bg3 		=  1.268348e-5;
	double ag 		= ag1 + ag2*TC + ag3*TC*TC;
	double bg 		= bg1 + bg2*TC + bg3*TC*TC;
	double r 		= in->density/1000.0;
	double g   		=  ag * pow((1.0 - r),bg);	
    in->g           = g;

	/* Born coefficient */
	double w, z, re, reref, X1, X2, G;
	if (charge == 0.0){
		w   = wr;
	} 
	else{
		z       = charge;
		reref   = z*z/(wr/eta + z/3.082);
		re      = reref + fabs(z) * g;
		X1      = -eta * (fabs(z*z*z)/(re*re) - z/pow((3.082 + g),2.0));
		X2      = 2.0*eta * (z*z*z*z/(re*re*re) - z/pow((3.082 + g),3.0));
		w       = eta * (z*z/re - z/(3.082 + g));
	}

	G  = 4.184 * (Gf - Sr * (TK - Tref) - c1 * (TK * log(TK / Tref) - TK + Tref) + a1 * (Pbar - Pref) + a2 * log((psi + Pbar) / (psi + Pref)) - c2 * ((1.0 / (TK - theta) - 1.0 / (Tref - theta)) * ((theta - TK) / theta) - TK / pow(theta,2.0) * log((Tref * (TK - theta)) / (TK * (Tref - theta)))) + (1.0 / (TK - theta)) * (a3 * (Pbar - Pref) + a4 * log((psi + Pbar) / (psi + Pref))) + (w * (-Z - 1.0) - wr * (-ZPrTr - 1.0) + wr * YPrTr * (TK - Tref)));
	G /= 1000.0; // turn to kJ


	/* fill structure to send back to main */
	PP_ref PP_ref_db;
	
	/* Calculate normalizing factor using bulk-rock composition */
	double factor  = 0.0;
	
	/* Calculate the number of atoms in the bulk-rock composition */
	double fbc     = 0.0;
	for (i = 0; i < len_ox; i++){
		fbc += bulk_rock[i]*apo[i];
	}
	
	/* Calculate the number of atom in the solution */
	double ape = 0.0;
	for (i = 0; i < len_ox; i++){
		ape += composition[i]*apo[i];
	}
	
	/* Calculate normalizing factor */
	factor = fbc/ape;

	strcpy(PP_ref_db.Name, name);
	for (i = 0; i < len_ox; i++){
		PP_ref_db.Comp[i] = composition[i];
	}

    /* convert Gibbs energy from SUPCRT to HSC convention */
    // double cor        =  SUPCRT_to_HSC(ElH, composition, len_ox);

	PP_ref_db.charge  =  charge;
	PP_ref_db.gbase   =  G;// + cor;
	PP_ref_db.factor  =  factor;
	PP_ref_db.phase_shearModulus  = 0.0;
    
	printf(" %4s %+10f %+10f | factor: %+10f\n",name,G,PP_ref_db.gbase,PP_ref_db.factor);
	for (i = 0; i < len_ox; i++){
		printf("%+10f",PP_ref_db.Comp[i]*PP_ref_db.factor); 
	}
	printf("\n\n");

	return (PP_ref_db);
}
