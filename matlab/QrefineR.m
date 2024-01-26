function [coordinates,newElements,newIrregular,prPts, varargout,r_coor,r_elem] ...
    = QrefineR(coordinates,elements,irregular,varargin)

%QrefineR: local refinement of quadrilateral mesh by red refinement, 
%          where marked elements are refined by bisecting all edges 
%          of the element
%
%Usage:
%
% [coordinates,elements4,irregular,dirichlet,neumann] ...
%    = QrefineR(coordinates,elements4,irregular,dirichlet,neumann,marked)
% or
%
% [coordinates,elements4,irregular] ...
%    = QrefineR(coordinates,elements4,irregular,marked)
%
%Comments:
%
%    QrefineR expects as input a mesh described by the 
%    fields coordinates, elements4, irregular, dirichlet (optional) and neumann 
%    (optional). The vector marked contains the indices of elements which
%    are refined by refining all edges of the element.
%    Further elements will be refined by a red refinement to obtain. 
%    1-Irregularity of the mesh is ensured by the 1-Irregular Rule.
% 
%    The function returns the refined mesh in terms of the same data as
%    for the input.
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

coordinates_initial = coordinates;
nE = size(elements,1);
markedElements = varargin{end};

%*** Obtain geometric information on edges
[edge2nodes,irregular2edges,element2edges,...
    boundary2edges{1:nargin-4}] ...
    = provideGeometricData(irregular,elements,varargin{1:end-1});
%*** Mark edges for refinement and existing hanging nodes
edge2newNode = zeros(1,size(edge2nodes,1));
edge2newNode(element2edges(markedElements,:)) = 1;
edge2newNode(irregular2edges(:,1)) = 1;
kdx = 1;
while ~isempty(kdx) || ~isempty(swap)
    markedEdge = edge2newNode(element2edges);
    %*** Change flags for elements
    kdx = find(sum(abs(markedEdge),2)<4 & ...
        (sum(abs(markedEdge),2)>2 | min(markedEdge,[],2)<0));
    [idx,jdx] = find(~markedEdge(kdx,:));
    edge2newNode(element2edges(kdx(idx)+(jdx-1)*nE)) = 1;
    %*** Change flags for irregular marker elements
    markedEdge = edge2newNode(irregular2edges);
    flag = irregular2edges(any(markedEdge(:,2:end),2),1);
    swap = find(edge2newNode(flag)~=-1);
    edge2newNode(flag(swap)) = -1;
end
%*** Generate new nodes on edges
edge2newNode(irregular2edges(:,1)) = -1;
idx = edge2newNode>0;
edge2newNode(idx) = size(coordinates,1) + (1:nnz(idx));
coordinates(edge2newNode(idx),:)=(coordinates(edge2nodes(idx,1),:)...
                            +coordinates(edge2nodes(idx,2),:))/2;
                        
                        
prPts = [  edge2nodes(idx,1),     edge2nodes(idx,2), edge2nodes(idx,2)*0, edge2nodes(idx,2)*0];

%*** Refine boundary conditions
varargout = cell(nargout-4,1);
for j = 1:nargout-4
    boundary = varargin{j};
    if ~isempty(boundary)
        newNodes = edge2newNode(boundary2edges{j})';
        markedEdges = find(newNodes);
        if ~isempty(markedEdges)
            boundary = [boundary(~newNodes,:); ...
                boundary(markedEdges,1),newNodes(markedEdges); ...
                newNodes(markedEdges),boundary(markedEdges,2)];
        end
    end
    varargout{j} = boundary;
end
%*** Provide new nodes for refinement of elements
edge2newNode(irregular2edges(:,1)) = irregular(:,3);
newNodes = reshape(edge2newNode(element2edges),[],4);
%*** Determine type of refinement for each element
reftyp = (newNodes~=0)*2.^(0:3)';
none   = reftyp < 15;
red    = reftyp == 15;
%*** Generate new interior nodes if red elements are refined
idx = find(red);
midNodes = zeros(nE,1);
midNodes(idx) = size(coordinates,1)+(1:length(idx));
coordinates = [coordinates;...
    ( coordinates(elements(idx,1),:) ... 
    + coordinates(elements(idx,2),:) ...
    + coordinates(elements(idx,3),:) ...
    + coordinates(elements(idx,4),:) )/4];


prPts = [prPts; [  elements(idx,1),     elements(idx,2), elements(idx,3), elements(idx,4)]];


%*** Generate element numbering for refined mesh
idx             = zeros(nE,1);
idx(none)       = 1;
idx(red)        = 4;
idx             = [1;1+cumsum(idx)];
%*** Generate new elements
newElements = zeros(idx(end)-1,4);
newElements(idx(none),:) = elements(none,:);
newElements([idx(red),1+idx(red),2+idx(red),3+idx(red)],:)=...
   [elements(red,1),newNodes(red,1),midNodes(red),newNodes(red,4);...
    elements(red,2),newNodes(red,2),midNodes(red),newNodes(red,1);...
    elements(red,3),newNodes(red,3),midNodes(red),newNodes(red,2);...
    elements(red,4),newNodes(red,4),midNodes(red),newNodes(red,3)];

%*** Generate new irregularity data
kdx = find(reftyp >0 & reftyp < 15);
[idx,jdx,val] = find(newNodes(kdx,:));
edx = element2edges(kdx(idx)+(jdx-1)*nE);
newIrregular = [edge2nodes(edx,:),val(:)];
newNodes = reshape(edge2newNode(irregular2edges(:,2:3)),[],2);
kdx = find(sum(newNodes,2)~=0);
[idx,jdx,val] = find(newNodes(kdx,:));
edx = irregular2edges(kdx(idx)+(jdx-1+1)*size(irregular2edges,1));
newIrregular = [newIrregular;[edge2nodes(edx(:),:),val(:)]];

%% Determine duplicate coordinates and remove them (causes plotting issues later)
indDuplicates         =   1:size(coordinates,1);
[coordUnique,iA,iC]   =   unique(coordinates,'rows','stable');
if length(iA)<length(indDuplicates)
    indDuplicates(iA)=[];        % this will now have all duplicate points
    
    for i=1:length(indDuplicates)
       ind = find( (coordinates_initial(:,1)==coordinates(indDuplicates(i),1)) & ...
                                 (coordinates_initial(:,2)==coordinates(indDuplicates(i),2)) );
        Duplicate_num(i) = ind(1);
    end


%     Duplicate_num = iC(indDuplicates);  % indices of the duplicates
    
    % to not mix up the list of points, we only replace it in the
    % "elements"
    for i=1:length(indDuplicates)
        ind                 = find(newElements==indDuplicates(i));
        if ~isempty(ind)
            newElements(ind)    = Duplicate_num(i);
        end
    end

    % what we should do is 


end


 







