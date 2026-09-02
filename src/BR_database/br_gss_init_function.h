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
#ifndef __BR_GSS_INIT_FUNCTION_H_
#define __BR_GSS_INIT_FUNCTION_H_

#include "../initialize.h"

SS_ref G_SS_br_ctd_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_car_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_chl_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_mica_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_talc_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_ilm_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_bt_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_ol_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_ep_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_opx_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_amph_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_spl_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_stau_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_crd_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_grt_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_omph_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_amphx_init_function(SS_ref SS_ref_db, global_variable gv);
SS_ref G_SS_br_fsp_init_function(SS_ref SS_ref_db, global_variable gv);

void BR_SS_init(SS_init_type *SS_init, global_variable gv);

#endif
