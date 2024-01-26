function [edge2nodes,element3edges,element4edges,varargout] ...
             = provideGeometricData(elements3,elements4,varargin)
% provideGeometricData: returns edge-based geometric data for mesh
%
%Usage:
% for triangular meshes:
%   [edges2nodes,element3edges,~,dirichlet2edges,neumann2edges] ...
%    = provideGeometricData(elements3,zeros(0,4),dirichlet,neumann)
% for quadrialteral meshes:
%   [edges2nodes,~,elemen4edges,dirichlet2edges,neumann2edges] ...
%    = provideGeometricData(zeros(0,3),elements4,dirichlet,neumann)
% for mixed meshes:
%   [edges2nodes,element3edges,element4edges,dirichlet2edges,neumann2edges] ...
%    = provideGeometricData(elements3,elements4,dirichlet,neumann)
%
%Comments:
%
%    provideGeometricData expects as input a mesh described
%    by the fields elements3, elements4, dirichlet (optional), and neumann
%    (optional). The function chooses a numbering of the edges and 
%    returns the relations to nodes, edges, and the boundary conditions.
%
%    edges2nodes(k) returns the indices of the two nodes of the k-th edge.
%    element3edges(j,k) provides the edge number of the edge between the two 
%    nodes elements3(j,k) and elements3(j,k+1). Analougously, for elements4.
%    dirichlet2edges(k) provides the number of the k-th Dirichlet edge given
%    by dirichlet(k,:). The same applies for neumann2edges.
%
%Remark:
%
%    This program is a supplement to the paper 
%    >> Adaptive Mesh Refinement in 2D - An Efficient Implementation in Matlab <<
%    by S. Funken, and A. Schmidt. The reader should 
%    consult that paper for more information.   
%
%Authors:
% 
%    S. Funken, A. Schmidt  20-08-18

%*** Obtain geometric information on edges
%*** Collect all edges
edges = [reshape(elements3(:,[1:3,2:3,1]),[],2); ...
         reshape(elements4(:,[1:4,2:4,1]),[],2)] ;  
ptr = [3*size(elements3,1),4*size(elements4,1),zeros(1,nargin-2)]; 
for j = 1:nargin-2
  ptr(j+2) = size(varargin{j},1);
  edges = [edges;varargin{j}];
end   
ptr = cumsum(ptr);
%*** Create numbering of edges
[edge2nodes,~,ie] = unique(sort(edges,2),'rows');
element3edges = reshape(ie(       1:ptr(1)),[],3);
element4edges = reshape(ie(ptr(1)+1:ptr(2)),[],4);
%*** Provide boundary2edges
for j = 1:nargin-2
  varargout{j} = ie(ptr(j+1)+1:ptr(j+2));
end


