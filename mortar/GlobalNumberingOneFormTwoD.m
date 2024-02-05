function gn = GlobalNumberingOneFormTwoD(n,p)
%
%   GlobalNumberingOneFormTwoD computes the global numbering of 1-form
%   degrees of freedom on a primal mesh.
%
%   gn = GlobalNumberingOneFormTwoD(n,p)
%
%   input:
%       n   :: the number of elements; if n is a number it has n x n
%              elements, if it is a vector, then it has n(1) x n(2) elements
%       p   :: the order of the 0-form approximation; if p is a number then 
%              px=py=p, if it is a vector, then px!=py
%
%   output:
%       gn  :: global numbering of the lines of the mesh (gn is a N x (2*p*(p+1) 
%              matrix where N is the total number of elements)
%
%   Local (per element) numbering:
%
%       1 - surfaces are numbered from left to right and bottom to top
%       2 - horizontal edges are numbered from left to right and bottom to top    
%       3 - vertical edges are numbered from bottom to top and left to right
%       4 - nodes are numbered from left to right and bottom to top
%
%   Global numbering:
%
%       note: within each element local numbering rules apply
%       elements:
%       1.1 - elements are numbered from left to right and bottom to top 
%       surfaces
%       2.1 - next the element surfaces are numbered
%       edges
%       3.1 - internal horizontal edges of the elements
%       3.2 - internal horizontal element boundary edges
%       3.3 - internal vertical edges of the elements
%       3.4 - internal vertical element boundary edges
%       3.5 - bottom/top domain boundaries
%       3.6 - left/right domain boundaries
%       nodes
%       4.1 - internal element nodes
%       4.2 - internal element corner nodes
%       4.3 - internal horizontal element boundary nodes
%       4.4 - internal vertical element boundary nodes
%       4.5 - bottom/top domain boundary nodes (corners excluded)
%       4.6 - left/right domain boundary nodes (corners excluded)
%       4.7 - bottom/top domain corner nodes
%       4.8 - bottom/top internal domain boundary nodes
%       4.9 - left/right internal domain boundary nodes
%
%   Copyright 2010 Artur Palha
%   $Revision: 1.2 $  $Date: 03/03/2010 $
%   Revised by: Peter Kuystermans
%   $Revision: 2.0 $  $Date: 13/07/2011 $    
%   Revised by: Peter Kuystermans
%   $Revision: 3.0 $  $Date: 17/10/2011 $ 
    
%-------------------------------------------------------------------------%
% input checks                                                            %
%-------------------------------------------------------------------------%       
 
    if length(n)>2
        fprintf(':: n :: has too many values \n');
        return
    elseif length(n)==1
        n = [n n];
    end   
    if prod(n)==0
        fprintf(':: n :: has at least one zero \n');
        return
    end
    
    if length(p)>2
        fprintf(':: p :: has too many values \n');
        return
    elseif length(p)==1
        p = [p p];
    end    
    if prod(p)==0
        fprintf(':: p :: has at least one zero \n');
        return
    end      
    
%-------------------------------------------------------------------------%
% parameters and storage                                                  %
%-------------------------------------------------------------------------%   
    
    % compute the total number of elements
    nElements = n(1)*n(2);
    
    % compute the total number of internal (horizontal and vertical) edges
    nInternalEdgesX = p(1)*(p(2)-1)*nElements;
    nInternalEdgesY = p(2)*(p(1)-1)*nElements;

    % compute the number of (horizontal and vertical) edges per element
    nElementEdgesX = p(1)*(p(2)+1);
    nElementEdgesY = p(2)*(p(1)+1);
    
    % allocate memory space for global numbering matrix
    gn = zeros(nElements,nElementEdgesX+nElementEdgesY,'uint32');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%-------------------------------------------------------------------------%
% compute the numbering for the horizontal edges (xi basis functions)     %
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------------------------------------%
% number the horizontal internal edges of the elements                    %
%-------------------------------------------------------------------------%

    % select the internal edges of the elements and number them
    internalEdgesColumns = p(1)+1:p(1)*p(2);
    gn(:,internalEdgesColumns) = reshape(1:nInternalEdgesX,[p(1)*(p(2)-1) nElements])';
    
    % keep track of last edge number
    edgesgn = nInternalEdgesX; 
    
%-------------------------------------------------------------------------%
% number the internal horizontal element boundary edges                   %
%-------------------------------------------------------------------------%
    
    % bottom and top edge internal element numbers
    bottomEdgeEdgesInternalNumber = 1:p(1);
    topEdgeEdgesInternalNumber = bottomEdgeEdgesInternalNumber+p(1)*p(2);
    
    if n(2)~=1
        edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+p(1)*n(1)*(n(2)-1));
        for k=1:(n(2)-1)
            gn((1:n(1))+n(1)*(k-1),topEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*p(1)*n(1)+1):(k*n(1)*p(1))),[p(1) n(1)])';
            gn((1:n(1))+n(1)*(k-1+1),bottomEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*p(1)*n(1)+1):(k*n(1)*p(1))),[p(1) n(1)])';
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%-------------------------------------------------------------------------%
% compute the numbering for the vertical edges (eta basis functions)      %
%-------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % store the last stored edge number
    lastNumber = edgesgn(end);
    
%-------------------------------------------------------------------------%
% number the vertical internal edges of the elements                      %
%-------------------------------------------------------------------------%

    % select the internal edges of the elements and number them
    internalEdgesColumns = p(2)+1:p(1)*p(2);
    gn(:,internalEdgesColumns+nElementEdgesX) = reshape(1:uint32(nInternalEdgesY),[p(2)*(p(1)-1) nElements])'+lastNumber;
    
    % keep track of last edge number
    edgesgn = uint32(nInternalEdgesY+lastNumber); 
    
%-------------------------------------------------------------------------%
% number the internal vertical element boundary edges                     %
%-------------------------------------------------------------------------%

    % left and right edge internal element numbers
    leftEdgeEdgesInternalNumber = 1:p(2);
    rightEdgeEdgesInternalNumber = leftEdgeEdgesInternalNumber+p(1)*p(2);
        
    if n(1)~=1
        edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+p(2)*n(2)*(n(1)-1));
        for k=1:(n(1)-1)
            gn(k:n(1):(n(1)*n(2)),rightEdgeEdgesInternalNumber+nElementEdgesX) = reshape(edgesgn(((k-1)*p(2)*n(2)+1):(k*n(2)*p(2))),[p(2) n(2)])';
            gn((k:n(1):(n(1)*n(2)))+1,leftEdgeEdgesInternalNumber+nElementEdgesX) = reshape(edgesgn(((k-1)*p(2)*n(2)+1):(k*n(2)*p(2))),[p(2) n(2)])';
        end
    end
     
%-------------------------------------------------------------------------%
% domain boundary edges                                                   %
%-------------------------------------------------------------------------%

    % lower boundary
    edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)*p(1));
    gn(1:n(1),bottomEdgeEdgesInternalNumber) = reshape(edgesgn,[p(1) n(1)])';     
    % upper boundary
    edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)*p(1));
    gn(((n(2)-1)*n(1)+1):(n(1)*n(2)),topEdgeEdgesInternalNumber) = reshape(edgesgn,[p(1) n(1)])';    
    % left boundary
    edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)*p(2));
    gn(1:n(1):((n(2)-1)*n(1))+1,leftEdgeEdgesInternalNumber+nElementEdgesX) = reshape(edgesgn,[p(2) n(2)])'; 
    % right boundary
    edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)*p(2));
    gn((1:n(1):((n(2)-1)*n(1))+1)+n(1)-1,rightEdgeEdgesInternalNumber+nElementEdgesX) = reshape(edgesgn,[p(2) n(2)])'; 
    
end