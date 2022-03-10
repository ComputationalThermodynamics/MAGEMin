function [EM_Fraction, SS_Names] = xEOS_From_NN(P,T)
% This returns the EM fractions of a solution model as computed by a pre-trained
% NN for a given list of T,P points (@ some stage perhaps chemistry as well)
%

minFrac = 5e-3; % below this, it is assumed to be zero

% Load mat-file
load('NN_EMFractions_KLB1_v1');

PT_vec          =   [P(:),T(:)]';
SolutionNames   = 	net_EMFraction.SolutionNames;

% Predict the mass fractions for current [P,T] point
Fractions       =       net_EMFraction.EMFractions(PT_vec)';
Fractions(Fractions<minFrac) = 0;

NN_SS_Names     =       fieldnames(net_EMFraction);
id = contains(NN_SS_Names,'EMFractions');   NN_SS_Names(id)=[];
id = contains(NN_SS_Names,'SolutionNames'); NN_SS_Names(id)=[];


% Compute all fractions
SS_Names    =   [];
for i = 1:length(NN_SS_Names)
    SS              =   NN_SS_Names{i};
    
    EM_Frac{:,i}    =   num2cell(net_EMFraction.(SS)(PT_vec)',2);
end

ind_SS       =  contains(SolutionNames,NN_SS_Names);
Fractions_SS =  Fractions(:,ind_SS);
for iPoint=1:length(P)
    ind         =   find(Fractions_SS(iPoint,:)>0);
    EM_Frac1D   =   [];
    
    SS_names_1D = [];
    for j=1:length(ind)
        EM_Frac1D{j}    =   EM_Frac{ind(j)}{iPoint};
        SS_names_1D{j}  =   NN_SS_Names{ind(j)};
    end
    EM_Fraction{iPoint} =   EM_Frac1D;
    SS_Names{iPoint}    =   SS_names_1D;
end