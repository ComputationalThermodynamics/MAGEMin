% function ComputeReactionLines(PseudoSectionData)
% This computes reaction lines from a AMR mesh with phase assemblages
% These lines are rather coarse (and should really be refined with a
% bisection & running MAGEMin again)


load ../test_KLB_lowPres

data = PseudoSectionData

% determine 
NumAssemblage = data.NumAssemblage;
T             = data.TP_vec(:,1);
P             = data.TP_vec(:,2);

% Every AMR element in the Phase diagram grid is ordered as:
%
%  4     3
%  o --- o
%  |     |
%  o --- o 
%  1     2
%
%
% And at every point we know the number of the stable assemblage 


elem        = data.elements;

reaction1_2 = abs(diff(NumAssemblage(elem(:,1:2  )),[],2));
reaction2_3 = abs(diff(NumAssemblage(elem(:,2:3  )),[],2));
reaction3_4 = abs(diff(NumAssemblage(elem(:,3:4  )),[],2));
reaction4_1 = abs(diff(NumAssemblage(elem(:,[4 1])),[],2));


% Determine the number of reactions:
reactions       = unique([reaction1_2; reaction2_3;reaction3_4;  reaction4_1])
reactions       = reactions(2:end);     % not counting zero
numReact        = length(reactions);
ReactionNumber  = 0:numReact+1;         % Give a number to each reaction on diagram      


% Determine the points of the reaction:
[~,reacN_12]=   ismember(reaction1_2, reactions);
ind         =   find(reaction1_2>0);
reac1_2_pts =   [mean(T(elem(ind,1:2  )),2), mean(P(elem(ind,1:2  )),2)];
reac1_2     =   reacN_12(ind);

[~,reacN_23]=   ismember(reaction2_3, reactions);
ind         =   find(reaction2_3>0);
reac2_3_pts =   [mean(T(elem(ind,2:3  )),2), mean(P(elem(ind,2:3  )),2)];
reac2_3     =   reacN_23(ind);

[~,reacN_34]=   ismember(reaction3_4, reactions);
ind         = find(reaction3_4>0);
reac3_4_pts = [mean(T(elem(ind,3:4  )),2), mean(P(elem(ind,3:4  )),2)];
reac3_4     =   reacN_34(ind);

[~,reacN_41]=   ismember(reaction4_1, reactions);
ind         = find(reaction4_1>0);
reac4_1_pts = [mean(T(elem(ind,[4 1])),2), mean(P(elem(ind,[4 1])),2)];
reac4_1     =   reacN_41(ind);


% Reconstruct line pieces 
ReactionNumbers = [reacN_12 reacN_23 reacN_34 reacN_41];
Reactions       = [reaction1_2 reaction2_3 reaction3_4 reaction4_1];
Reactions_T     = [mean(T(elem(:,1:2)),2) mean(T(elem(:,2:3)),2) mean(T(elem(:,3:4)),2) mean(T(elem(:,[4 1])),2)];
Reactions_P     = [mean(P(elem(:,1:2)),2) mean(P(elem(:,2:3)),2) mean(P(elem(:,3:4)),2) mean(P(elem(:,[4 1])),2)];


% This contains two points that each form line segments
REACT           =   10;
ind             =   find(sum(ReactionNumbers==REACT,2)==2);   % reaction 1

id              =   (ReactionNumbers(ind,:)==REACT);
Tpts            =   Reactions_T(ind,:);
Ppts            =   Reactions_P(ind,:);
Tvec            =   zeros(length(Tpts),2);
Pvec            =   zeros(length(Tpts),2);
for i=1:size(Tvec,1)
    Tvec(i,:) = Tpts(i,id(i,:));
    Pvec(i,:) = Ppts(i,id(i,:));
end


% NOTE: there appears to be a big, perhaps because we only check the
% differences between two fields and not the sum at the same time, so we
% may not actuall have unique phase fields here