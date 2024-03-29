function [Gamma, xEOS] = Interpolate_Gamma_xEOS_SinglePoint(Data, PseudoSectionData, Method, num_Neighbours)
% This interpolates Gamma and x_EOS on a single point given by Data
%
% Mostly used for debugging



XY_vec  = PseudoSectionData.XY_vec;
TP_norm = max(XY_vec)-min(XY_vec);


% Find current point:
[newPoints] = nearestneighbour([[Data.T Data.P]./TP_norm]', [XY_vec./TP_norm]','n', 1);

Computation = PseudoSectionData.Computation;

if isempty(Method)
    Method = 'Average surrounding points';
end
if isempty(num_Neighbours)
    num_Neighbours =9;
end

Computation.Gamma_Method	=   Method;
Computation.num_Neighbours  =   num_Neighbours;
Computation.EOS_Method      =   Computation.EOS_Method; 


% Update Gamma on points
elements    =   PseudoSectionData.elements;
PhaseData   =   Update_Gamma_onPoints(PseudoSectionData.PhaseData, XY_vec, elements, newPoints,[],Computation);
Gamma       =   PhaseData{newPoints}.Gamma;       % Gamma interpolated from surrounding points


% Update xEOS on points
PhaseData   = Update_xEOS_onPoints(PseudoSectionData.PhaseData,XY_vec, newPoints,[],Computation);
xEOS        = PhaseData{newPoints}.StartingValues_xEOS;       % xEOS combined from surrounding points

