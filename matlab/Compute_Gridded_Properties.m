function [T2D, P2D, Rho, Gibbs, Cp] = Compute_Gridded_Properties(varargin)
% This interpolates the properties of the AMR pseudosection grid to a 
% regular grid. 
%

PseudoSectionData = varargin{1};
if nargin==1
    nT = 100;
    nP = 100;
elseif nargin==3
    nT = varargin{2};
    nP = varargin{3};
else
    error('Wrong number of input arguments')
end


% 2D grid
Tmin        = min(PseudoSectionData.TP_vec(:,1));
Tmax        = max(PseudoSectionData.TP_vec(:,1));
dT          = (Tmax-Tmin)/(nT-1);

Pmin        = min(PseudoSectionData.TP_vec(:,2));
Pmax        = max(PseudoSectionData.TP_vec(:,2));
dP          = (Pmax-Pmin)/(nP-1);

[T2D,P2D]   = meshgrid(Tmin:dT:Tmax,Pmin:dP:Pmax);

% Interpolate data on grid
Rho         =   Interpolate_AMR_grid(PseudoSectionData.elements, PseudoSectionData.TP_vec, PseudoSectionData.Rho_tot, T2D(:), P2D(:));
Rho         =   reshape(Rho,size(T2D));

Gibbs       =   Interpolate_AMR_grid(PseudoSectionData.elements, PseudoSectionData.TP_vec, PseudoSectionData.Gibbs,     T2D(:), P2D(:));
Gibbs       =   reshape(Gibbs,size(T2D));

ind         =   find(isnan(Gibbs));
if ~isempty(ind)
%     Gibbs(ind)  =   griddata(PseudoSectionData.TP_vec(:,1), PseudoSectionData.TP_vec(:,2), PseudoSectionData.Gibbs,     T2D(ind), P2D(ind));
end



T0          =   273.15;
Cp          =   -(T2D(:,2:end-1)+T0).*diff(Gibbs,2,2)/dT^2;  
Cp          =   [Cp(:,1), Cp, Cp(:,end)];


