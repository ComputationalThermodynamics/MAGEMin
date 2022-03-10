function [coordinates, elements4, irregular] = GenerateRegularMesh_AMR(T_vec1D,P_vec1D)
% Generates a regular mesh for use with AMR methods (initial mesh)



[T_2D,P_2D]     =   meshgrid(T_vec1D,P_vec1D);
id              =   ones(size(P_2D));
id(find(id))    =   find(id);

id1             =   id(1:end-1,1:end-1);
id2             =   id(1:end-1,2:end  );
id3             =   id(2:end  ,2:end  );
id4             =   id(2:end  ,1:end-1);

elements4       =   [id1(:), id2(:), id3(:), id4(:)];
coordinates     =   [T_2D(:), P_2D(:)];
irregular       =   zeros(0,3);

