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
#ifndef __BR_NLOPT_OPT_FUNCTION_H_
#define __BR_NLOPT_OPT_FUNCTION_H_

#include "../MAGEMin.h"

typedef SS_ref (*NLopt_type) (      global_variable      gv,
                                     SS_ref               SS_ref_db      );

SS_ref NLopt_opt_br_ctd_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_car_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_chl_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_mica_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_talc_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_ilm_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_bt_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_ol_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_ep_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_opx_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_amph_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_spl_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_stau_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_crd_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_grt_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_omph_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_amphx_function(global_variable gv, SS_ref SS_ref_db);
SS_ref NLopt_opt_br_fsp_function(global_variable gv, SS_ref SS_ref_db);

void BR_NLopt_opt_init(NLopt_type *NLopt_opt, global_variable gv);

#endif
