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
/* Source: Pourteau et al. (2014), Contrib Mineral Petrol; Berman (1988)
   EOS formalism via the THERIAK JUN92.bs lineage. */
#include <math.h>
#include <string.h>

#include "../MAGEMin.h"
#include "BR_PP_endmembers.h"
#include "BR_gem_function.h"

#define BR_Tr    298.15
#define BR_Pr    1.0
#define BR_kbar2bar 1000.0

void BR_berman_HS( double T, double H0, double S0, const double cp[7],
                    double *H_out, double *S_out               ){
    double k0 = cp[0], k1 = cp[1], k2 = cp[2], k3 = cp[3];
    double k4 = cp[4], k5 = cp[5], k6 = cp[6];

    double H = H0 + k0*(T-BR_Tr) + 2.0*k1*(sqrt(T)-sqrt(BR_Tr))
             - k2*(1.0/T - 1.0/BR_Tr) - 0.5*k3*(1.0/(T*T) - 1.0/(BR_Tr*BR_Tr))
             + k4*log(T/BR_Tr) + 0.5*k5*(T*T-BR_Tr*BR_Tr) + (1.0/3.0)*k6*(T*T*T-BR_Tr*BR_Tr*BR_Tr);
    double S = S0 + k0*log(T/BR_Tr)
             - 2.0*k1*(1.0/sqrt(T) - 1.0/sqrt(BR_Tr))
             - 0.5*k2*(1.0/(T*T) - 1.0/(BR_Tr*BR_Tr))
             - (1.0/3.0)*k3*(1.0/(T*T*T) - 1.0/(BR_Tr*BR_Tr*BR_Tr))
             - k4*(1.0/T - 1.0/BR_Tr) + k5*(T-BR_Tr) + 0.5*k6*(T*T-BR_Tr*BR_Tr);

    *H_out = H;
    *S_out = S;
}

double BR_landa_G(double T, double P, const double landa[4]){
    double tq1b = landa[0], teq = landa[1], l1 = landa[2], l2 = landa[3];
    if (tq1b == 0.0){
        return 0.0;
    }
    double ctr = teq*(P - BR_Pr);
    double tqp = tq1b + ctr;
    double T2  = fmin(T, tqp);
    double T1  = (fabs(l2) > 1e-20) ? -l1/l2 + ctr : T2;
    if (T < T1){
        return 0.0;
    }
    double A = -ctr*T1*T1;
    double B = T1*(2.0*ctr + T1);
    double C = -(ctr + 2.0*T1);
    double dH = A*(T2-T1) + B*(T2*T2-T1*T1)/2.0 + C*(T2*T2*T2-T1*T1*T1)/3.0
              + (T2*T2*T2*T2-T1*T1*T1*T1)/4.0;
    double dS = A*log(T2/T1) + B*(T2-T1) + C*(T2*T2-T1*T1)/2.0
              + (T2*T2*T2-T1*T1*T1)/3.0;
    return (dH - T*dS)*l2*l2;
}

double BR_disorder_G(double T, double P, const double dis[10]){
    double d1 = dis[0], d2 = dis[1], d3 = dis[2], d4 = dis[3], d5 = dis[4];
    double d6 = dis[5], d8 = dis[6], td0 = dis[7], tdmax = dis[8], vadj = dis[9];
    if (fabs(td0) <= 1e-20 || fabs(tdmax) <= 1e-20 || T <= td0){
        return 0.0;
    }
    double td  = fmin(T, tdmax);
    double dH  = d1*(td-td0) + d2*(td*td-td0*td0)/2.0
               - d3*(1.0/td-1.0/td0) + 2.0*d4*(sqrt(td)-sqrt(td0))
               + d5*(td*td*td-td0*td0*td0)/3.0 + d6*log(td/td0)
               - d8*(1.0/(td*td)-1.0/(td0*td0))/2.0;
    double dS  = d1*log(td/td0) + d2*(td-td0)
               - d3*(1.0/(td*td)-1.0/(td0*td0))/2.0
               - 2.0*d4*(1.0/sqrt(td)-1.0/sqrt(td0))
               + d5*(td*td-td0*td0)/2.0 - d6*(1.0/td-1.0/td0)
               - d8*(1.0/(td*td*td)-1.0/(td0*td0*td0))/3.0;
    double vd  = (fabs(vadj) > 10.0) ? dH/(10.0*vadj) : 0.0;
    return dH - T*dS + vd*(P - BR_Pr);
}

double BR_berman_EOS_dG(double T, double P, double V0, const double eos[4]){
    double v1 = eos[0], v2 = eos[1], v3 = eos[2], v4 = eos[3];
    return V0*( (v1/2.0-v2)*(P*P-BR_Pr*BR_Pr) + v2*(P*P*P-BR_Pr*BR_Pr*BR_Pr)/3.0
              + (1.0-v1+v2+v3*(T-BR_Tr)+v4*(T-BR_Tr)*(T-BR_Tr))*(P-BR_Pr) );
}

double BR_albite_Salje_G(double T, double P){
    if (T >= 1290.0){ return 0.0; }

    double A0   = 5.479;
    double B    = 6854.0;
    double BOD  = -9301.0;
    double COD  = 43600.0;
    double DD0  = -2.171;
    double DD1  = -3.043;
    double DD2  = -0.001569;
    double DD3  = 0.000002109;
    double AOD0 = 41.620;
    double TC   = 1251.0;
    double TOD  = 824.1;

    double D    = DD0 + DD1*T + DD2*T*T + DD3*T*T*T;
    double T1   = T - TC;
    double T2   = T - TOD;
    double D3   = D*D*D;
    double D5   = D3*D*D;

    double E1   = A0*AOD0*T1*T2/D - D;
    double E2   = AOD0*B*T2/D + BOD*pow(A0,3)*pow(T1,3)/D3;
    double E3   = 3.0*BOD*A0*A0*B*T1*T1/D3 + COD*pow(A0,5)*pow(T1,5)/D5;
    double E4   = 3.0*BOD*A0*B*B*T1/D3 + 5.0*pow(A0,4)*COD*B*pow(T1,4)/D5;
    double E5   = BOD*pow(B,3)/D3 + 10.0*COD*pow(A0,3)*B*B*pow(T1,3)/D5;
    double E6   = 10.0*COD*A0*A0*pow(B,3)*T1*T1/D5;
    double E7   = 5.0*COD*A0*pow(B,4)*T1/D5;
    double E8   = COD*pow(B,5)/D5;

    double q_int  = 0.1;
    double f_last = -1.0;
    double q      = (T < 1100.0) ? 1.1 : 0.6;
    double q_last = q;
    int    ite    = 0;

    q -= q_int;
    if (q < 0.00001){ q = 0.0001; }
    while (1){
        double q2 = q*q;
        double f  = q*(E1 + q2*(E2 + q2*(E3 + q2*(E4 + q2*(E5 + q2*(E6 + q2*(E7 + q2*E8)))))));
        ite += 1;
        if (fabs(q - q_last) < 0.0001 || ite > 30){ break; }
        if (f*f_last < 0.001){
            q_int *= 0.5;
            q     += q_int;
        }
        else{
            q_last = q;
            f_last = f;
            q     -= q_int;
            if (q < 0.00001){ q = 0.0001; }
        }
    }

    double q_c = q;
    for (int k = 0; k < 50; k++){
        double q2 = q_c*q_c;
        double f  = q_c*(E1 + q2*(E2 + q2*(E3 + q2*(E4 + q2*(E5 + q2*(E6 + q2*(E7 + q2*E8)))))));
        double df = E1 + q2*(3.0*E2 + q2*(5.0*E3 + q2*(7.0*E4 + q2*(9.0*E5 + q2*(11.0*E6 + q2*(13.0*E7 + q2*15.0*E8))))));
        if (df == 0.0){ break; }
        double dq = f/df;
        q_c      -= dq;
        if (!(q_c > 0.0 && q_c < 1.5)){ break; }
        if (fabs(dq) < 1e-15){ break; }
    }
    if (q_c > 0.0 && q_c < 1.5 && fabs(q_c - q) < 1e-2){ q = q_c; }

    double A    = A0*T1;
    double Qod  = (-A*q - B*q*q*q)/D;

    double S    = -( 0.5*A0*q*q + 0.5*AOD0*Qod*Qod + (DD1 + 2.0*DD2*T + 3.0*DD3*T*T)*q*Qod );
    double H    = -0.5*A0*TC*q*q + 0.25*B*pow(q,4) - 0.5*AOD0*TOD*Qod*Qod + 0.25*BOD*pow(Qod,4)
                + COD/6.0*pow(Qod,6) + (DD0 - DD2*T*T - 2.0*DD3*T*T*T)*q*Qod;
    double V    = H/335282.925;

    return H - T*S + (P - BR_Pr)*V;
}

static void BR_cubic(double b, double c, double d, double *x1, double *x2, double *x2i, double *x3){
    double pi = 3.14159263538979;
    *x2  = 0.0;
    *x2i = 0.0;
    *x3  = 0.0;
    if (fabs(c) < 1e-12 && fabs(d) < 1e-12){
        *x1 = -b;
        return;
    }
    double q  = ((2.0*b*b*b)/27.0 - (b*c)/3.0 + d)/2.0;
    double p  = (3.0*c - b*b)/9.0;
    double r  = sqrt(fabs(p));
    if (r*q < 0.0){ r = -r; }
    double ff = q/(r*r*r);
    if (p > 0.0){
        double phi3 = log(ff + sqrt(ff*ff + 1.0))/3.0;
        *x1  = -r*(exp(phi3) - exp(-phi3)) - b/3.0;
        *x2i = 1.0;
    }
    else if (q*q + p*p*p > 0.0){
        double phi3 = log(ff + sqrt(ff*ff - 1.0))/3.0;
        *x1  = -r*(exp(phi3) + exp(-phi3)) - b/3.0;
        *x2i = 1.0;
    }
    else{
        double phi3 = atan(sqrt(1.0 - ff*ff)/ff)/3.0;
        *x1  = -2.0*r*cos(phi3) - b/3.0;
        *x2  = 2.0*r*cos(pi/3.0 - phi3) - b/3.0;
        *x2i = 0.0;
        *x3  = 2.0*r*cos(pi/3.0 + phi3) - b/3.0;
    }
}

static double BR_psat_h2o(double T){
    static const double A[8] = {-7.8889166, 2.5514255, -6.716169, 33.239495, -105.38479, 174.35319, -148.39348, 48.631602};
    if (T <= 314.0){
        return exp(6.3573118 - 8858.843/T + 607.56335/pow(T, 0.6));
    }
    double v   = T/647.25;
    double w   = fabs(1.0 - v);
    double wsq = sqrt(w);
    double ff  = 0.0;
    for (int i = 0; i < 8; i++){
        ff += A[i]*w;
        w  *= wsq;
    }
    return 220.93*exp(ff/v);
}

static double BR_haar_pressure(double rh, double T, double B, double BB, const double *taui, double *dpr_out){
    static const double GI[40] = {
        -.53062968529023e4, .22744901424408e5, .78779333020687e4, -.69830527374994e3,
         .17863832875422e6,-.39514731563338e6, .33803884280753e6, -.13855050202703e6,
        -.25637436613260e7, .48212575981415e7,-.34183016969660e7,  .12223156417448e7,
         .11797433655832e8,-.21734810110373e8, .10829952168620e8, -.25441998064049e7,
        -.31377774947767e8, .52911910757704e8,-.13802577177877e8, -.25109914369001e7,
         .46561826115608e8,-.72752773275387e8, .41774246148294e7,  .14016358244614e8,
        -.31555231392127e8, .47929666384584e8, .40912664781209e7, -.13626369388386e8,
         .69625220862664e7,-.10834900096447e8,-.22722827401688e7,  .38365486000660e7,
         .68833257944332e5, .21757245522644e6,-.26627944829770e5, -.70730418082074e6,
        -.225e1, -1.68e1, .055e1, -93.0e1 };
    static const int    KI[40]   = {1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5, 6,6,6,6, 7,7,7,7, 9,9,9,9, 3,3, 1, 5, 2,2,2, 4};
    static const int    LI[40]   = {1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 0,3,3,3, 0,2,0,0};
    static const double RHOI[4]  = {0.319, 0.310, 0.310, 1.550};
    static const double TTTI[4]  = {640.0, 640.0, 641.6, 270.0};
    static const double ALPI[4]  = {34.0, 40.0, 30.0, 1050.0};
    static const double BETI[4]  = {2.0e4, 2.0e4, 4.0e4, 25.0};

    double RT    = 4.6152*T;
    double rh2   = rh*rh;
    double y     = rh*B/4.0;
    double er    = exp(-rh);
    double y3    = (1.0 - y)*(1.0 - y)*(1.0 - y);
    double aly   = 11.0*y;
    double bety  = 44.33333333333333*y*y;
    double f1    = (1.0 + aly + bety)/y3;
    double f2    = 4.0*y*(BB/B - 3.5);
    double ermi[10];
    ermi[0] = 1.0;
    ermi[1] = 1.0 - er;
    for (int i = 2; i < 10; i++){ ermi[i] = ermi[i-1]*ermi[1]; }

    double pr  = 0.0;
    double dpr = 0.0;
    for (int i = 0; i < 36; i++){
        double s = GI[i]/taui[LI[i]]*ermi[KI[i]-1];
        pr  += s;
        dpr += (2.0 + rh*(KI[i]*er - 1.0)/ermi[1])*s;
    }
    for (int i = 36; i < 40; i++){
        int    j     = i - 36;
        double del   = rh/RHOI[j] - 1.0;
        double rhoi2 = RHOI[j]*RHOI[j];
        double tau   = T/TTTI[j] - 1.0;
        double qhex  = -ALPI[j]*pow(del, LI[i]) - BETI[j]*tau*tau;
        double q10   = (qhex > -150.0) ? GI[i]*pow(del, LI[i])*exp(-ALPI[j]*pow(del, LI[i]) - BETI[j]*tau*tau) : 0.0;
        double qm    = LI[i]/del - KI[i]*ALPI[j]*pow(del, KI[i]-1);
        double s     = q10*qm*rh2/RHOI[j];
        pr  += s;
        dpr += s*(2.0/rh + qm/RHOI[j]) - rh2/rhoi2*q10*(LI[i]/del/del + KI[i]*(KI[i]-1)*ALPI[j]*pow(del, KI[i]-2));
    }
    pr   = rh*(rh*er*pr + RT*(f1 + f2));
    dpr  = rh*er*dpr + RT*((1.0 + 2.0*aly + 3.0*bety)/y3 + 3.0*y*f1/(1.0 - y) + 2.0*f2);
    *dpr_out = dpr;
    return pr;
}

double BR_water_Haar_G(double T, double P){
    static const double GI[40] = {
        -.53062968529023e4, .22744901424408e5, .78779333020687e4, -.69830527374994e3,
         .17863832875422e6,-.39514731563338e6, .33803884280753e6, -.13855050202703e6,
        -.25637436613260e7, .48212575981415e7,-.34183016969660e7,  .12223156417448e7,
         .11797433655832e8,-.21734810110373e8, .10829952168620e8, -.25441998064049e7,
        -.31377774947767e8, .52911910757704e8,-.13802577177877e8, -.25109914369001e7,
         .46561826115608e8,-.72752773275387e8, .41774246148294e7,  .14016358244614e8,
        -.31555231392127e8, .47929666384584e8, .40912664781209e7, -.13626369388386e8,
         .69625220862664e7,-.10834900096447e8,-.22722827401688e7,  .38365486000660e7,
         .68833257944332e5, .21757245522644e6,-.26627944829770e5, -.70730418082074e6,
        -.225e1, -1.68e1, .055e1, -93.0e1 };
    static const int    KI[40]   = {1,1,1,1, 2,2,2,2, 3,3,3,3, 4,4,4,4, 5,5,5,5, 6,6,6,6, 7,7,7,7, 9,9,9,9, 3,3, 1, 5, 2,2,2, 4};
    static const int    LI[40]   = {1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 1,2,4,6, 0,3,3,3, 0,2,0,0};
    static const double CI[18]   = {.19730271018e2, .209662681977e2, -.483429455355e0, .605743189245e1, 22.56023885, -9.87532442,
                                    -.43135538513e1, .458155781e0, -.47754901883e-1, .41238460633e-2, -.27929052852e-3, .14481695261e-4,
                                    -.56473658748e-6, .16200446e-7, -.3303822796e-9, .451916067368e-11, -.370734122708e-13, .137546068238e-15};
    static const double RHOI[4]  = {0.319, 0.310, 0.310, 1.550};
    static const double TTTI[4]  = {640.0, 640.0, 641.6, 270.0};
    static const double ALPI[4]  = {34.0, 40.0, 30.0, 1050.0};
    static const double BETI[4]  = {2.0e4, 2.0e4, 4.0e4, 25.0};

    double RT    = 4.6152*T;
    double GREF  = -54955.2356146119409;
    double T0    = 647.073;
    double taui[7];
    taui[0] = 1.0;
    taui[1] = T/T0;
    for (int i = 2; i < 7; i++){ taui[i] = taui[i-1]*taui[1]; }

    double B  = -0.3540782*log(taui[1]) + 0.7478629 + 0.007159876/taui[3] - 0.003528426/taui[5];
    double BB = 1.1278334 - 0.5944001/taui[1] - 5.010996/taui[2] + 0.63684256/taui[4];

    double PS = 220.55;
    if (T <= 647.25){ PS = BR_psat_h2o(T); }

    double ark = 1.279186e8 - 2.241415e4*T;
    double brk = 1.428062e1 + 6.092237e-4*T;
    double oft = ark/(P*sqrt(T));
    double buk = -10.0*8.31441*T/P;
    double cuk = oft - brk*brk + brk*buk;
    double duk = -brk*oft;
    double x1, x2, x2i, x3, vol;
    BR_cubic(buk, cuk, duk, &x1, &x2, &x2i, &x3);
    if (fabs(x2i) > 1e-12){
        vol = x1;
    }
    else if (P < PS){
        vol = fmax(x1, fmax(x2, x3));
    }
    else{
        vol = fmin(x1, fmin(x2, x3));
    }
    double rhn = (vol <= 0.0) ? 1.9 : 1.0/vol*18.0152;

    double rh  = rhn;
    for (int loo = 0; loo < 100; loo++){
        rh = rhn;
        if (rh <= 1e-8){ rh = 1e-8; }
        if (rh > 1.9){ rh = 1.9; }
        double dpr;
        double pr = BR_haar_pressure(rh, T, B, BB, taui, &dpr);
        if (dpr <= 0.0){
            rhn = (P <= PS) ? rhn*0.95 : rhn*1.05;
        }
        else{
            if (dpr < 0.01){ dpr = 0.01; }
            double s = (P - pr)/dpr;
            if (fabs(s) > 0.1){ s = 0.1*s/fabs(s); }
            rhn = rh + s;
        }
        double dp = fabs(1.0 - pr/P);
        double dr = fabs(1.0 - rhn/rh);
        if (dp < 1e-5 && dr < 1e-5){ break; }
    }
    rh = rhn;

    double y  = rh*B/4.0;
    double x  = 1.0 - y;
    double er = exp(-rh);
    double ermi[10];
    ermi[0] = 1.0;
    ermi[1] = 1.0 - er;
    for (int i = 2; i < 10; i++){ ermi[i] = ermi[i-1]*ermi[1]; }

    double aa = RT*(-log(x) - 43.33333333333333/x + 28.16666666666667/x/x + 4.0*y*(BB/B - 3.5) + 15.16666666666667 + log(rh*RT/1.01325));
    for (int i = 0; i < 36; i++){
        aa += GI[i]/KI[i]/taui[LI[i]]*ermi[KI[i]];
    }
    for (int i = 36; i < 40; i++){
        int    j    = i - 36;
        double del  = rh/RHOI[j] - 1.0;
        double tau  = T/TTTI[j] - 1.0;
        double qhex = -ALPI[j]*pow(del, KI[i]) - BETI[j]*tau*tau;
        if (qhex > -150.0){
            aa += GI[i]*pow(del, LI[i])*exp(qhex);
        }
    }

    double tr  = T/1.0e2;
    double w   = pow(tr, -3);
    double aid = 1.0 + (CI[0]/tr + CI[1])*log(tr);
    for (int i = 2; i < 18; i++){
        aid += CI[i]*w;
        w   *= tr;
    }
    aa -= RT*aid;

    return (aa + P/rh)*1.80152 - GREF - 306685.5925;
}

static PP_ref BR_pack_PP_ref(char *name, int len_ox, const double *Comp,
                              double *bulk_rock, double *apo, double gbase_J){
    PP_ref PP_ref_db;
    strcpy(PP_ref_db.Name, name);
    for (int i = 0; i < len_ox; i++){
        PP_ref_db.Comp[i] = Comp[i];
    }
    double fbc = 0.0;
    for (int i = 0; i < len_ox; i++){
        fbc += bulk_rock[i]*apo[i];
    }
    double ape = 0.0;
    for (int i = 0; i < len_ox; i++){
        ape += Comp[i]*apo[i];
    }
    PP_ref_db.gbase              = gbase_J/BR_kbar2bar;
    PP_ref_db.factor             = fbc/ape;
    PP_ref_db.phase_shearModulus = 0.0;
    PP_ref_db.phase_bulkModulus  = 0.0;
    PP_ref_db.phase_expansivity  = 0.0;
    PP_ref_db.phase_cp           = 0.0;
    return PP_ref_db;
}

PP_ref BR_G_EM_function(   int          EM_dataset,
                            int          len_ox,
                            int         *id,
                            double      *bulk_rock,
                            double      *apo,
                            double       Pkbar,
                            double       T,
                            char        *name,
                            char        *state          ){
    (void) EM_dataset;
    (void) id;
    (void) state;

    double P = Pkbar * BR_kbar2bar;

    int pp_id = BR_find_PP_id(name);
    if (pp_id >= 0){
        PP_db_br PP_return = Access_BR_PP_DB(pp_id);
        double H_T, S_T;
        BR_berman_HS(T, PP_return.H, PP_return.S, PP_return.cp_berman, &H_T, &S_T);
        double gbase_J = (H_T - T*S_T) + BR_berman_EOS_dG(T, P, PP_return.V, PP_return.eos_berman)
                       + BR_landa_G(T, P, PP_return.landa) + BR_disorder_G(T, P, PP_return.dis);
        if (strcmp(name, "ab") == 0){
            gbase_J += BR_albite_Salje_G(T, P);
        }
        if (strcmp(name, "h2o") == 0){
            gbase_J = BR_water_Haar_G(T, P);
        }
        return BR_pack_PP_ref(name, len_ox, PP_return.Comp, bulk_rock, apo, gbase_J);
    }

    PP_ref PP_ref_db;
    memset(&PP_ref_db, 0, sizeof(PP_ref_db));
    strcpy(PP_ref_db.Name, name);
    return PP_ref_db;
}
