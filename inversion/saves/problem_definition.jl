#=
    Problem definition
    Here is define the range of parameters to be adjusted during calibration

    For W's parameters the structure is of the form:

    phase{
        W1 -> active val P_dep Tep scalar_range P_dep_range T_dep_range
        W2 -> active val P_dep Tep scalar_range P_dep_range T_dep_range
        ...
        Wn -> active val P_dep Tep scalar_range P_dep_range T_dep_range
    }

    1. here range is [lb,ub] 
    2. only the phase parameters listed above will overwrite the default values. 
        if you want to fix one or more parameters just set the bounds to be equal to desired vallue
    3. make sure the name of the phase fits the one from the database you currently want to use!
=#
mutable struct def_Ws{ _T  } 
    phase       ::  String
    ss_num      ::  Int64                       #this is the id of the ss_num
    W_num       ::  Int64                       #id of the W parameter
    val         ::  _T 
    P_dep       ::  _T 
    T_dep       ::  _T 

    val_range   ::  Vector{Float64}
    P_dep_range ::  Vector{Float64} 
    T_dep_range ::  Vector{Float64} 
end


#= IGNEOUS DATABASE INFOS
    Here the informations needed to fill the definition problem are displayed for convenience

    spl: 1; n_W's, 28; n_em's, 8
    bi: 2; n_W's, 15; n_em's, 6
    cd: 3; n_W's, 3; n_em's, 3
    cpx: 4; n_W's, 45; n_em's, 10
    ep: 5; n_W's, 3; n_em's, 3
    g: 6; n_W's, 15; n_em's, 6
    hb: 7; n_W's, 55; n_em's, 11
    ilm: 8; n_W's, 10; n_em's, 5
    liq: 9; n_W's, 66; n_em's, 12
    ol: 10; n_W's, 6; n_em's, 4
    opx: 11; n_W's, 36; n_em's, 9
    pl4T: 12; n_W's, 3; n_em's, 3
    fl: 13; n_W's, 55; n_em's, 11
=#

                   # ph  id   n act val Pdep Tdep   valRange    Prange      Trange
def_W  = [
    def_Ws{Float64}("ol", 10, 1, 0.0, 0.0, 0.0, [0.0,40.0], [-1.0,1.0], [-1.0,1.0]),
    def_Ws{Float64}("ol", 10, 2, 0.0, 0.0, 0.0, [0.0,40.0], [-1.0,1.0], [-1.0,1.0]),
    def_Ws{Float64}("ol", 10, 3, 0.0, 0.0, 0.0, [0.0,40.0], [-1.0,1.0], [-1.0,1.0]),
    def_Ws{Float64}("ol", 10, 4, 0.0, 0.0, 0.0, [0.0,40.0], [-1.0,1.0], [-1.0,1.0]),
    def_Ws{Float64}("ol", 10, 5, 0.0, 0.0, 0.0, [0.0,40.0], [-1.0,1.0], [-1.0,1.0]),
    def_Ws{Float64}("ol", 10, 6, 0.0, 0.0, 0.0, [0.0,40.0], [-1.0,5.0], [-1.0,1.0]),
]