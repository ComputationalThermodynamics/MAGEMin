function [xEOS] = xEOS_From_NN(SS,XY_vec)
% This returns the xEOS of a solution model as cmputed by a NN for a given
% list of T,P points (and later chemistry)
%

% Load mat-file
load('NN_test1');

% Compute NN
switch SS
    case 'spn'
        xEOS    = net_spn([XY_vec]')';
        lb      =   [0 0 0 0 -1 -1 -1];
        ub      =   [1 1 1 1  1  1  1];
    case 'cpx'
        xEOS    =   net_cpx([XY_vec]')';
        lb      =   [0 0 0 0 -1 0 0 0 0];
        ub      =   [1 2 1 1  1 1 1 1 1];

    case 'opx'  
        xEOS    =   net_opx([XY_vec]')';
        lb      =   [0 0 0 -1 0 0 0 0];
        ub      =   [1 1 1 2  1 1 1 1];
    case 'pli'
        xEOS    =   net_pli([XY_vec]')';
        lb      =   [0 0];
        ub      =   [1 1];
        
    case 'g'
        xEOS    =   net_g([XY_vec]')';
        lb      =   [0 0 0 0 0];
        ub      =   [1 1 1 1 1];
        
    case 'ol'
        xEOS    =   net_ol([XY_vec]')';
        lb      =   [0 0 0];
        ub      =   [1 1 1];
    case 'liq'
        xEOS    =   net_liq([XY_vec]')';
        lb      =   [0 0 0 0 0 0 0 0 0 0 0];
        ub      =   [1 1 1 1 1 1 1 1 1 1 1];
        
    otherwise
        warning(['xEOS NN not implemented for phase ',SS])
        xEOS = [];
end


if ~isempty(xEOS)
    % ensure that bounds are satisfied
    for i=1:length(lb)
        ind_lb = find(xEOS(:,i)<lb(i)); xEOS(ind_lb,i)=lb(i) + 1e-10;
        ind_ub = find(xEOS(:,i)>ub(i)); xEOS(ind_ub,i)=ub(i) - 1e-10;
    end
end
