    % Function to create uniform grid to plot isocontours (e.g., isentropic lines)
    function [PseudoSectionData] = Update_interpolated_grids(PseudoSectionData)

    XY_vec = PseudoSectionData.XY_vec;
  
    np      = 50;
    Tstep   = (max(XY_vec(:,1))-min(XY_vec(:,1)))/(np-1);
    Pstep   = (max(XY_vec(:,2))-min(XY_vec(:,2)))/(np-1);
    Tx      = min(XY_vec(:,1)):Tstep:max(XY_vec(:,1));
    Px      = min(XY_vec(:,2)):Pstep:max(XY_vec(:,2));

    [X,Y] = meshgrid(Tx,Px);

    PseudoSectionData.ugrid_Tx               = Tx;
    PseudoSectionData.ugrid_Px               = Px;
    PseudoSectionData.ugrid_entropy          = griddata(XY_vec(:,1),XY_vec(:,2),PseudoSectionData.entropy(:,1),X,Y);
    PseudoSectionData.ugrid_Density          = griddata(XY_vec(:,1),XY_vec(:,2),PseudoSectionData.Rho_tot(:,1),X,Y);
    PseudoSectionData.ugrid_Density_liquid   = griddata(XY_vec(:,1),XY_vec(:,2),PseudoSectionData.Rho_liq(:,1),X,Y);
    PseudoSectionData.ugrid_Density_solid    = griddata(XY_vec(:,1),XY_vec(:,2),PseudoSectionData.Rho_sol(:,1),X,Y);
    PseudoSectionData.ugrid_Vp               = griddata(XY_vec(:,1),XY_vec(:,2),PseudoSectionData.Vp(:,1),X,Y);
    PseudoSectionData.ugrid_Vs               = griddata(XY_vec(:,1),XY_vec(:,2),PseudoSectionData.Vs(:,1),X,Y);
    PseudoSectionData.ugrid_liquid_fraction  = griddata(XY_vec(:,1),XY_vec(:,2),PseudoSectionData.liq(:,1),X,Y);
 
end


