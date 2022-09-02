% MATLAB script to compute silicate melt viscosity.
%
%  Citation: Giordano D, Russell JK, & Dingwell DB (2008)
%  Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science
%  Letters, v. 271, 123-134.
%
% ________________________________________________________________
% INPUT: Chemical compositions of silicate melts (wt% oxides) as:
% SiO2 TiO2 Al2O3 FeO(T) MnO MgO CaO Na2O K2O P2O5 H2O F2O-1
% One line for each melt composition
% _________________________________________________________________

% ________________________________________________________________
% OUTPUT: eta values at T(C) values (1 line per melt)
% ________________________________________________________________

% VFT Multicomponent-Model Coedfficients
% -----------------------------------------------------------------

% modified and optimised by Tobias Keller, 10. June, 2022

function   [eta] = grdmodel08(wt,TC)

AT  =  -4.55;
bb  =  [159.56  -173.34 72.13 75.69 -38.98 -84.08 141.54 -2.43 -0.91 17.62];
cc  =  [2.75 15.72 8.32 10.2 -12.29 -99.54 0.3 ];

% molar weights
mw  = [ 60.0843, 79.8658, 101.961276, 71.8444, 70.937449,40.3044,56.0774, 61.97894, 94.1960, 141.9446,18.01528, 18.9984];

% convert wt% to mol%, note all normalized on anhydrous components first
wtn  = [wt(:,1:10).*(100-wt(:,11))./(sum(wt(:,1:10),2)+wt(:,12)) wt(:,11) 0.5.*wt(:,12).*(100-wt(:,11))./(sum(wt(:,1:10),2)+wt(:,12))];
mp   = wtn./mw;
mol  = 100.*(mp./sum(mp,2));

% Load composition-basis matrix for multiplication against model-coefficients
% Result is two matrices bcf[nx by 10] and ccf[nx by 7]
siti  =  mol(:,1) + mol(:,2);
tial  =  mol(:,2)+mol(:,3);
fmm   =  mol(:,4) + mol(:,5) + mol(:,6);
nak   =  mol(:,8) + mol(:,9);
b1    =  siti;
b2    =  mol(:,3);
b3    =  mol(:,4) + mol(:,5) + mol(:,10);
b4    =  mol(:,6);
b5    =  mol(:,7);
b6    =  mol(:,8) + mol(:,11) + mol(:,12);
b7    =  mol(:,11) + mol(:,12) + log(1+mol(:,11));
b12   =  siti.*fmm;
b13   =  (siti + mol(:,3) + mol(:,10)).*( nak + mol(:,11) );
b14   =  mol(:,3).*nak;

c1    =  mol(:,1);
c2    =  tial;
c3    =  fmm;
c4    =  mol(:,7);
c5    =  nak;
c6    =  log(1+mol(:,11) + mol(:,12));
c11   =  mol(:,3) + fmm + mol(:,7) - mol(:,10);
c11   =  c11.*(nak + mol(:,11) + mol(:,12));
bcf   =  [b1 b2 b3 b4 b5 b6 b7 b12 b13 b14];
ccf   =  [c1 c2 c3 c4 c5 c6 c11];

BT    =  sum(bb.*bcf,2);
CT    =  sum(cc.*ccf,2);

TK    =  TC + 273.15;
eta   =  10.^(min(18,max(-6,AT + BT./(TK(:)-CT))));

