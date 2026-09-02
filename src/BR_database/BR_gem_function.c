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

void BR_berman_HS( double T, double H0, double S0, const double cp[11],
                    double *H_out, double *S_out               ){
    double k0 = cp[0], k1 = cp[1], k2 = cp[2], k3 = cp[3];
    double k4 = cp[4], k5 = cp[5], k6 = cp[6];
    double Tt = cp[7], deltaH = cp[8], l1 = cp[9], l2 = cp[10];

    double H = H0 + k0*(T-BR_Tr) + 2.0*k1*(sqrt(T)-sqrt(BR_Tr))
             - k2*(1.0/T - 1.0/BR_Tr) - 0.5*k3*(1.0/(T*T) - 1.0/(BR_Tr*BR_Tr))
             + k4*log(T/BR_Tr) + 0.5*k5*(T*T-BR_Tr*BR_Tr) + (1.0/3.0)*k6*(T*T*T-BR_Tr*BR_Tr*BR_Tr);
    double S = S0 + k0*log(T/BR_Tr)
             - 2.0*k1*(1.0/sqrt(T) - 1.0/sqrt(BR_Tr))
             - 0.5*k2*(1.0/(T*T) - 1.0/(BR_Tr*BR_Tr))
             - (1.0/3.0)*k3*(1.0/(T*T*T) - 1.0/(BR_Tr*BR_Tr*BR_Tr))
             - k4*(1.0/T - 1.0/BR_Tr) + k5*(T-BR_Tr) + 0.5*k6*(T*T-BR_Tr*BR_Tr);

    if (Tt != 0.0){
        if (T > Tt){
            H += deltaH + 0.5*l1*l1*(Tt*Tt-BR_Tr*BR_Tr)
               + (2.0/3.0)*l1*l2*(Tt*Tt*Tt-BR_Tr*BR_Tr*BR_Tr)
               + 0.25*l2*l2*(Tt*Tt*Tt*Tt-BR_Tr*BR_Tr*BR_Tr*BR_Tr);
            S += deltaH/Tt + l1*l1*(Tt-BR_Tr) + l1*l2*(Tt*Tt-BR_Tr*BR_Tr)
               + (1.0/3.0)*l2*l2*(Tt*Tt*Tt-BR_Tr*BR_Tr*BR_Tr);
        }
        else{
            H += 0.5*l1*l1*(T*T-BR_Tr*BR_Tr)
               + (2.0/3.0)*l1*l2*(T*T*T-BR_Tr*BR_Tr*BR_Tr)
               + 0.25*l2*l2*(T*T*T*T-BR_Tr*BR_Tr*BR_Tr*BR_Tr);
            S += l1*l1*(T-BR_Tr) + l1*l2*(T*T-BR_Tr*BR_Tr)
               + (1.0/3.0)*l2*l2*(T*T*T-BR_Tr*BR_Tr*BR_Tr);
        }
    }

    *H_out = H;
    *S_out = S;
}

double BR_berman_EOS_dG(double T, double P, double V0, const double eos[4]){
    double v1 = eos[0], v2 = eos[1], v3 = eos[2], v4 = eos[3];
    return V0*( (v1/2.0-v2)*(P*P-BR_Pr*BR_Pr) + v2*(P*P*P-BR_Pr*BR_Pr*BR_Pr)/3.0
              + (1.0-v1+v2+v3*(T-BR_Tr)+v4*(T-BR_Tr)*(T-BR_Tr))*(P-BR_Pr) );
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
        double gbase_J = (H_T - T*S_T) + BR_berman_EOS_dG(T, P, PP_return.V, PP_return.eos_berman);
        return BR_pack_PP_ref(name, len_ox, PP_return.Comp, bulk_rock, apo, gbase_J);
    }

    PP_ref PP_ref_db;
    memset(&PP_ref_db, 0, sizeof(PP_ref_db));
    strcpy(PP_ref_db.Name, name);
    return PP_ref_db;
}
