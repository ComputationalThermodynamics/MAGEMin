function [] = ExportData_LaMEM(PseudoSectionData, FileName, FileDir, GridResolution)
% Exports the pseudosection data to a format that can be used by LaMEM
% 
% In particular, we save: 
%
%   Density (average)
%   Density (liquid)
%   Density (solid)
%   MeltFraction
%
% To a file & we interpolate data from the AMR grid to a regular grid
%
%  Usage:
%       [] = ExportData_LaMEM(PseudoSectionData, FileName, FileDir, GridResolution, info)
%
%       PseudoSectionData: structure, created with the GUI
%       FileName:           Filename of phase diagram
%       FileDir:            Directory of phase diagram
%       GridResolution:     # of grid points in Temperature and Pressure
%                           direction
%
%

nt      =   GridResolution(1); % # of gridpoints Temp
np      =   GridResolution(2); % # of grid points P

%========================================================= 
% Create regular T,P grid
%========================================================= 
Tmin    =   min(PseudoSectionData.TP_vec(:,1));
Pmin    =   min(PseudoSectionData.TP_vec(:,2));
Tmax    =   max(PseudoSectionData.TP_vec(:,1));
Pmax    =   max(PseudoSectionData.TP_vec(:,2));

dT      =   (Tmax-Tmin)/(nt-1);
dP      =   (Pmax-Pmin)/(np-1);
[T,P]   =   meshgrid(Tmin:dT:Tmax, Pmin:dP:Pmax);

%========================================================= 
% Interpolate to regular grid 
%========================================================= 
Rho_sol                 =   PseudoSectionData.Rho_sol;
ind                     =   find(Rho_sol>0);
Rho_sol(Rho_sol==0)     =   min(Rho_sol(ind));
Rho_sol(isnan(Rho_sol)) =   min(Rho_sol(ind));

Rho_liq                 =   PseudoSectionData.Rho_liq;
ind                     =   find(Rho_liq>0);
Rho_liq(Rho_liq==0)     =   min(Rho_liq(ind));
Rho_liq(isnan(Rho_liq)) =   min(Rho_liq(ind));

% shift bounds slightly to avoid NaN's
Pi      = P;   Pi(end,:) = Pi(end,:) - 1e-6; Pi(1,:) = Pi(1,:) + 1e-6;
Ti      = T;   Ti(:,end) = Ti(:,end) - 1e-6; Ti(:,1) = Ti(:,1) + 1e-6;

% Interpolate:
rho_sol = Interpolate_AMR_grid(PseudoSectionData.elements, PseudoSectionData.TP_vec, Rho_sol,               Ti, Pi);
rho_liq = Interpolate_AMR_grid(PseudoSectionData.elements, PseudoSectionData.TP_vec, Rho_liq,               Ti, Pi);
melt    = Interpolate_AMR_grid(PseudoSectionData.elements, PseudoSectionData.TP_vec, PseudoSectionData.liq, Ti, Pi);
Vp      = Interpolate_AMR_grid(PseudoSectionData.elements, PseudoSectionData.TP_vec, PseudoSectionData.Vp,  Ti, Pi);
Vs      = Interpolate_AMR_grid(PseudoSectionData.elements, PseudoSectionData.TP_vec, PseudoSectionData.Vs,  Ti, Pi);

VpVs      = Vp./Vs;
ind       = find(VpVs>3);
VpVs(ind) = NaN;


%========================================================= 
% Transfer to format that LaMEM expects
%========================================================= 
T_K     =   T(:)+273.15;    % in K
P_bar   =   P(:)*1e3;       % in bar
Pmin    =   Pmin*1e3;
dP      =   dP*1e3;
Tmin    =   Tmin+273.15;
dT      =   dT*1;



%========================================================= 
% Write output 
%========================================================= 
fid = fopen([FileDir,filesep,FileName],'w+');

% Write header info
fprintf(fid,'8 \n');
fprintf(fid,'  \n');
fprintf(fid,'Phase diagram version: 0.11 \n');
fprintf(fid,'DATA:  \n');
fprintf(fid,'Phase diagram always needs this 5 columns: \n');
fprintf(fid,'1                     2              3         4       5         6        7      8      \n');
fprintf(fid,'rho_fluid [kg/m3]   melt []    rho [kg/m3]   T [K]   P [bar]  Vp[km/s] Vs[km/s] Vp/Vs[] \n');
fprintf(fid,'  \n');
fprintf(fid,'  \n');
fprintf(fid,'LINES:    \n'); 
fprintf(fid,'1:     Number of columns (to be able to increase this later)   \n'); 
fprintf(fid,'1-50:  Comments     \n');
fprintf(fid,'50:    Lowest T [K] \n');
fprintf(fid,'51:    T increment     \n');
fprintf(fid,'52:    # of T values     \n');
fprintf(fid,'53-55: Same for pressure \n');


% Write some specific info
for i=1:5
    fprintf(fid,'  \n');    
end
fprintf(fid,'Phase diagram computed with MAGEMin  \n');   
fprintf(fid,'Phase diagram name         :  %s \n',PseudoSectionData.Name);
fprintf(fid,'This file generated on     :  %s \n',date);
fprintf(fid,'Chemistry:  \n');

Chem = table2cell(PseudoSectionData.Chemistry.MolProp);
for i=1:length(Chem)
    fprintf(fid,'%8s %10.8f \n',Chem{i,1},Chem{i,2});
end

for i=1:13
    fprintf(fid,'  \n');    
end

fprintf(fid,'%f \n', Tmin);
fprintf(fid,'%f \n', dT);
fprintf(fid,'%i \n', nt);
fprintf(fid,'%f \n', Pmin);
fprintf(fid,'%f \n', dP);
fprintf(fid,'%i \n', np);


rho_liq =   rho_liq(:);
melt    =   melt(:);
rho_sol =   rho_sol(:);

for i=1:length(T_K)
    fprintf(fid,'%1.2f %1.5f %1.2f %1.4f %1.3f %1.2f %1.2f %1.2f\n', rho_liq(i), melt(i), rho_sol(i), T_K(i), P_bar(i), Vp(i), Vs(i), VpVs(i));
end
fclose(fid);



%========================================================= 
% Plot (using same as usual plotting routine)
%========================================================= 
rho         =   reshape(rho_sol,np,nt);
rho_fluid   =   reshape(rho_liq,np,nt);
melt        =   reshape(melt   ,np,nt);
T           =   reshape(T_K    ,np,nt);       % degree C
P           =   reshape(P_bar  ,np,nt)./1e3;  % kbar  - LaMEM input must be bar!



figure(1),clf
subplot(221)
pcolor(T,P,rho)
shading interp
colorbar
ylabel('P [kbar]')
xlabel('T [Celcius]')
title('density rock w/out melt!')

subplot(222)
pcolor(T,P,melt)
shading interp
colorbar
ylabel('P [kbar]')
xlabel('T [Celcius]')
title('Melt content')

subplot(223)
pcolor(T,P,rho_fluid)
ylabel('P [kbar]')
xlabel('T')
title('density of fluid')
shading interp
colorbar

subplot(224)
rho_average = (1-melt).*rho + melt.*rho_fluid;
pcolor(T,P,rho_average)
ylabel('P [kbar]')
xlabel('T')
title('Combined density of fluid + solid (should NOT be passed to LaMEM!)')
shading interp
colorbar




 


    

figure(2),clf
subplot(221)
pcolor(T,P,Vp)
shading interp
colorbar
ylabel('P [kbar]')
xlabel('T [Celcius]')
title('Vp [km/s]')

subplot(222)
pcolor(T,P,Vs)
shading interp
colorbar
ylabel('P [kbar]')
xlabel('T [Celcius]')
title('Vs [km/s]')

subplot(223)
pcolor(T,P,VpVs)
shading interp
colorbar
ylabel('P [kbar]')
xlabel('T [Celcius]')
title('Vp/Vs []')



