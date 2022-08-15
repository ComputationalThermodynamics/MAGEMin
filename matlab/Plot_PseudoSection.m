function Plot_PseudoSection(PseudoSectionData, Data)

fontsiz=16;

patch('Faces',  PseudoSectionData.elements, ...
     'Vertices', PseudoSectionData.TP_vec,...
     'FaceVertexCData',Data,...
     'FaceColor','interp', ...
     'Linestyle','-');

 axis('tight');
 colorbar;
 
 xlabel('Temperature [C]','Fontsize',fontsiz);
 ylabel('Pressure [kbar]','Fontsize',fontsiz);
 