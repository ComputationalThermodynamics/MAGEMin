    % Function to create uniform grid to plot isocontours (e.g., isentropic lines)
    function [PseudoSectionData] = Update_interpolated_grids(PseudoSectionData)

    TP_vec = PseudoSectionData.TP_vec;
  
    np      = 100;
    Tstep   = (max(TP_vec(:,1))-min(TP_vec(:,1)))/(np-1);
    Pstep   = (max(TP_vec(:,2))-min(TP_vec(:,2)))/(np-1);
    Tx      = min(TP_vec(:,1)):Tstep:max(TP_vec(:,1));
    Px      = min(TP_vec(:,2)):Pstep:max(TP_vec(:,2));

    [X,Y] = meshgrid(Tx,Px);

    PseudoSectionData.ugrid_Tx               = Tx;
    PseudoSectionData.ugrid_Px               = Px;
    PseudoSectionData.ugrid_entropy          = griddata(TP_vec(:,1),TP_vec(:,2),PseudoSectionData.entropy(:,1),X,Y);
    PseudoSectionData.ugrid_Density          = griddata(TP_vec(:,1),TP_vec(:,2),PseudoSectionData.Rho_tot(:,1),X,Y);
    PseudoSectionData.ugrid_Density_liquid   = griddata(TP_vec(:,1),TP_vec(:,2),PseudoSectionData.Rho_liq(:,1),X,Y);
    PseudoSectionData.ugrid_Density_solid    = griddata(TP_vec(:,1),TP_vec(:,2),PseudoSectionData.Rho_sol(:,1),X,Y);
    PseudoSectionData.ugrid_Vp               = griddata(TP_vec(:,1),TP_vec(:,2),PseudoSectionData.Vp(:,1),X,Y);
    PseudoSectionData.ugrid_Vs               = griddata(TP_vec(:,1),TP_vec(:,2),PseudoSectionData.Vs(:,1),X,Y);
    PseudoSectionData.ugrid_liquid_fraction  = griddata(TP_vec(:,1),TP_vec(:,2),PseudoSectionData.liq(:,1),X,Y);
 
end


