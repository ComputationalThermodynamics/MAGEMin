function [surface] = Grid_Gibbs_Surface(Local, ind_x, ind_y)
% This creates a 3D plot of the driving force surface from a range of
% samples stored in Local
%

AddSmooth   =    true;

xEOS        =   Local.xEOS_end(:,[ind_x ind_y]);
dF          =   Local.df;



% Round coordinates a bit
% xEOS     	=   round(xEOS*5,1)/5;  % round points to be on top
xEOS     	=   round(xEOS*2,1)/2;  % round points to be on top


x           = unique(xEOS(:,1));
y           = unique(xEOS(:,2));

[x2d,y2d]   = meshgrid(x,y);

for i=1:size(x2d,1)
    for j=1:size(x2d,2)
        ind = find(xEOS(:,1)==x2d(i,j) & xEOS(:,2)==y2d(i,j));
        
        if ~isempty(ind)
            z2d(i,j) = min(dF(ind));
        else
            z2d(i,j) = NaN;
        end
    end
end

surface.x = x2d;
surface.y = y2d;
surface.z = z2d;


% Use griddata to extrapolate points that are covered by data
ind                 =   find(~isnan(z2d));
surface.z_Tn        =   griddata(x2d(ind),y2d(ind),z2d(ind),x2d,y2d, 'nearest');       % has values everywhere


if AddSmooth
    
    % Use gridfit to fit a smooth surface
    numP        =   50;
    surface.xS   =   linspace(min(xEOS(:,1)),max(xEOS(:,1)), numP);
    surface.yS   =   linspace(min(xEOS(:,2)),max(xEOS(:,2)), numP);
    surface.zS   =   gridfit(x2d(:),y2d(:),surface.z_Tn(:),surface.xS,surface.yS);
    
    % MASK areas not covered by data
    [x,y]               =   meshgrid(surface.xS,surface.yS);
    ind                 =   find(~isnan(z2d));
    zT                  =   griddata(x2d(ind),y2d(ind),z2d(ind),x,y, 'linear');
    ind                 =   isnan(zT);
    surface.zS(ind)     =   NaN;
end

%==========================================================================

% Do the same but now for the proportions of the endmembers rather than for
% x-eos (this is how we actually discretized the space)
Prop_end    =   Local.Prop_end(:,[ind_x ind_y]);
Prop_end    =   round(Prop_end*2,1)/2;  % round points to be on top



xP          = unique(Prop_end(:,1));
yP          = unique(Prop_end(:,2));

[xP2d,yP2d]   = meshgrid(xP,yP);

for i=1:size(xP2d,1)
    for j=1:size(xP2d,2)
        ind = find(Prop_end(:,1)==xP2d(i,j) & Prop_end(:,2)==yP2d(i,j));
        
        if ~isempty(ind)
            zP2d(i,j) = min(dF(ind));
        else
            zP2d(i,j) = NaN;
        end
    end
end

surface.xP = xP2d;
surface.yP = yP2d;
surface.zP = zP2d;



% Use griddata to extrapolate points that are covered by data
ind                 =   find(~isnan(zP2d));
surface.zP_Tn       =   griddata(xP2d(ind),yP2d(ind),zP2d(ind),xP2d,yP2d, 'nearest');       % has values everywhere
surface.zP_Tl       =   griddata(xP2d(ind),yP2d(ind),zP2d(ind),xP2d,yP2d, 'linear');        % has NaN outside covered range    
surface.zP          =   surface.zP_Tl;

% Use gridfit to fit a smooth surface
numP            =   50;
surface.xP_S    =   linspace(min(Prop_end(:,1)),max(Prop_end(:,1)), numP);
surface.yP_S    =   linspace(min(Prop_end(:,2)),max(Prop_end(:,2)), numP);

if AddSmooth
    % Use nearest neigbour to create smooth surface
    ind                 =   find(~isnan(surface.zP_Tn));
    surface.zP_S        =   gridfit(xP2d(ind),yP2d(ind),surface.zP_Tn(ind),surface.xP_S,surface.yP_S);
    
    
    
    
    % determine points outside covered rangfe
    [x,y]               =   meshgrid(surface.xP_S,surface.yP_S);
    ind                 =   find(~isnan(zP2d));
    zP_out              =   griddata(xP2d(ind),yP2d(ind),zP2d(ind),x,y, 'linear');
    ind                 =   isnan(zP_out);
    surface.zP_S(ind)   =   NaN;
end

