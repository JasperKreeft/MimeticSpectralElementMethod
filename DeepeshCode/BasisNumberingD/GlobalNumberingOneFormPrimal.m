function gn = GlobalNumberingOneFormPrimal(n,p)
%GlobalNumberingOneFormPrimal computes the global numbering of 1-form
%degrees of freedom on a  primal mesh.
%
%   gn = GlobalNumberingOneFormPrimal(n, p)
%
%   Where:
%       n   :: the number of elements. If n is a number it has n x n
%              elements, if it is a vector, then it has n(1) x n(2) elements.
%       p   :: the order of the 0-form approximation
%
%   Returns the global numbering of the edges of the mesh. gn is a matrix
%   N x (p+1)^2, where N is the total number of elements.
%
%   Internal numbering:
%
%       Internally edges are numbered from left to right and from bottom to
%       top. So the first node of element 1, that is gn(1,1) is the bottom
%       left node of the first element. The node gn(1,p+1) is the top left
%       node of the first element. The node gn(1,p*(p+1)+1) is the right
%       bottom node of the first element. The node gn(1, (p+1)*(p+1)) is
%       the top right node of the first element.
%
%   Global numbering:
%
%       1- internal edges of the elements, following the internal numbering;
%       2- internal corner edges, from bottom to top and left to right;
%       3- edges at vertical edges, from bottom to top and left to right;
%       4- edges at horizontal edges, from left to right and bottom to top;
%       5- corner edges at the boundary: bottom boundary from left to
%          right, top boundary, from left to right, left boundary from bottom
%          to top, right boundary from bottom to top;
%       6- edges at edges of the boundaries: bottom boundary from left to
%          right, top boundary, from left to right, left boundary from bottom
%          to top, right boundary from bottom to top;
%
%   See also: GLOBALNUMBERINGONEFORMDUAL, PLOTGLOBALNUMBERINGONEFORMPRIMAL

%   Copyright 2010 Artur Palha
%   $Revision: 1.2 $  $Date: 2010/03/03 $
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this code is based on the GlobalNumberingZeroFormPrimal, hence the
    % many comments of the code. The parts not needed are just comment out
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 
    if length(n)>2
        % check if n has more than two values
        disp(sprintf(':: n :: has too many values'));
        return
    elseif length(n)==1
            n = [n n];
    end
    
    if prod(n)==0
        % check if n has a zero value
        disp(sprintf(':: n :: has at least one zero'));
        return
    end
    % compute the number of elements
    nElements = n(1)*n(2);

    % compute the number of interior elements
    if prod(n)==1
        nInternalElements = 0;
    else
        nInternalElements = (n(1)-2)*(n(2)-2);
    end
    
%--------------------------------------------------------------------------
% compute the numbering for the xi basis functions part
%--------------------------------------------------------------------------
    
    % compute the number of internal edges
    nInternalEdges = p*(p-1)*nElements;

    % compute the number of internal edges per element
    nElementInternalEdges = p*(p-1);

    % compute the total number of edges
    nEdges = (n(1)*p+1)*(n(2)*(p-1));

    % compute the number of edges in each element
    nElementEdges = p*(p+1);

    % compute the number of edges in the boundary
    nBoundaryEdges = 2*n(1)*p;

    % allocate memory space for global numbering matrix
    % Notice that edges in elements are numbered from bottom to top and from
    % left to right.
    % Elements follow the same kind of numbering: bottom to
    % top and left to right.
    gn = zeros(nElements, 2*nElementEdges, 'uint32');

    % number internal edges

    % first select the internal edges of the elements
    internalEdgesColumns = zeros(nElementInternalEdges,1,'uint32');
    for i=1:p
        for j=2:p
            internalEdgesColumns((i-1)*(p-1)+j-1) = (i-1)*(p+1) + j;
        end
    end

    gn(:,internalEdgesColumns) = reshape(1:nInternalEdges,[p*(p-1) nElements])';
    
    edgesgn = nInternalEdges; % just to keep track in which number we are.
    
    % there are no corner edges...
    
%     % number corner boundary edges of the interior elements
% 
%     % compute the global numbering of the edges
%     edgesgn = uint32(nInternalEdges+1):uint32(nInternalEdges+(n(1)-1)*(n(2)-1));
%     
%     % compute the internal elements
%     internalElements = zeros(nInternalElements,1,'uint32');
%     for i=2:(n(1)-1)
%         for j=2:(n(2)-1)
%             internalElements((i-2)*(n(2)-2)+j-1) = (i-1)*(n(2)) + j;
%         end 
%     end
% 
%     % compute the internal numbering of the corner edges
%     cornerEdgesInternalNumber = [1 p+1 p*(p+1)+1 (p+1)*(p+1)];
% 
%     % number them
%     if nInternalElements >= 1
%         % first the ones of the inner elements
%         for k=1:(n(1)-2)
%             gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerEdgesInternalNumber(1)) = edgesgn(((k-1)*(n(2)-1)+1):(k*(n(2)-1)-1));
%             gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerEdgesInternalNumber(2)) = edgesgn(((k-1)*(n(2)-1)+1+1):(k*(n(2)-1)-1+1));
%             gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerEdgesInternalNumber(3)) = edgesgn(((k)*(n(2)-1)+1):((k+1)*(n(2)-1)-1));
%             gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerEdgesInternalNumber(4)) = edgesgn(((k)*(n(2)-1)+1+1):((k+1)*(n(2)-1)-1+1));
%         end
% 
%         if n(1)>1
%             % the ones of the elements on the left boundary
%             k=1;
%             gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))-n(2), cornerEdgesInternalNumber(3)) = edgesgn(((k-1)*(n(2)-1)+1):(k*(n(2)-1)-1));
%             gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))-n(2), cornerEdgesInternalNumber(4)) = edgesgn(((k-1)*(n(2)-1)+1+1):(k*(n(2)-1)-1+1));
% 
%             % the ones of the elements on the right boundary
%             k=n(1)-2;
%             gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))+n(2), cornerEdgesInternalNumber(1)) = edgesgn(((k)*(n(2)-1)+1):((k+1)*(n(2)-1)-1));
%             gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))+n(2), cornerEdgesInternalNumber(2)) = edgesgn(((k)*(n(2)-1)+1+1):((k+1)*(n(2)-1)-1+1));
%         end
% 
%         if n(2)>1
%             % the ones of the elements on the bottom boundary
%             elements = 1:n(2)-2:((n(1)-2)*(n(2)-2));
%             gn(internalElements(elements)-1, cornerEdgesInternalNumber(2)) = edgesgn((1:n(2)-1:((n(1)-2)*(n(2)-1))));
%             gn(internalElements(elements)-1, cornerEdgesInternalNumber(4)) = edgesgn((1:n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);
% 
%             % the ones of the elements on the top boundary
%             elements = elements+n(2)-3;
%             gn(internalElements(elements)+1, cornerEdgesInternalNumber(1)) = edgesgn((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)));
%             gn(internalElements(elements)+1, cornerEdgesInternalNumber(3)) = edgesgn(((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);
%         end
%         % the four corner elements
%         gn(1,cornerEdgesInternalNumber(4)) = edgesgn(1); % left bottom
%         if n(2)>1
%             gn(n(2),cornerEdgesInternalNumber(3)) = edgesgn(n(2)-1); % left top
%         end
%         if n(1)>1
%             gn((n(1)-1)*n(2)+1,cornerEdgesInternalNumber(2)) = edgesgn((n(1)-2)*(n(2)-1)+1); % right bottom
%         end
%         gn(n(1)*n(2),cornerEdgesInternalNumber(1)) = edgesgn(end); % right top
%     end
    
    % number the edges of the internal edges
    
    % leftEdgeEdgesInternalNumber = 2:p;
    % rightEdgeEdgesInternalNumber = (p*(p+1)+2):((p+1)*(p+1)-1);
    
    bottomEdgeEdgesInternalNumber = 1:p+1:(p*(p+1));
    topEdgeEdgesInternalNumber = bottomEdgeEdgesInternalNumber + p;
    
    % vertical edges are not in dxi one forms
    
%     % start with the vertical edges
%     if n(1)~=1
%         if nInternalElements >= 1
%             edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+(p-1)*n(2)*(n(1)-1));
%         else
%             edgesgn = uint32(nInternalEdges+1):uint32(nInternalEdges+(p-1)*n(2)*(n(1)-1));
%         end
% 
%         for k=1:(n(1)-1)
%             gn(((k-1)*n(2)+1):k*n(2),rightEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*(p-1)*n(2)+1):(k*n(2)*(p-1))),[p-1 n(2)])';
%             gn((k*n(2)+1):(k+1)*n(2),leftEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*(p-1)*n(2)+1):(k*n(2)*(p-1))),[p-1 n(2)])';
%         end
%     end
    
    % now the horizontal ones
    if n(2)~=1
        if nInternalElements >= 1
            edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+p*n(1)*(n(2)-1));
        else
            edgesgn = uint32(nInternalEdges+1):uint32(nInternalEdges+p*n(1)*(n(2)-1));
        end

        for k=1:(n(2)-1)
            gn(k:n(2):(n(1)*n(2)),topEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*p*n(1)+1):(k*n(1)*p)),[p n(1)])';
            gn((k:n(2):(n(1)*n(2)))+1,bottomEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*p*n(1)+1):(k*n(1)*p)),[p n(1)])';
        end
    end
    
    % finally number the edges at the actual boundary
    
    % there are no corner edges
%     lowerLeft = 1;
%     upperLeft = p+1;
%     lowerRight = p*(p+1)+1;
%     upperRight = (p+1)*(p+1);
% 
%     % start with the edges of the elements corners except the corner edges of
%     % the domain
%     if n(1)>1
%         edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)-1);
%         gn(1:n(2):((n(1)-1)*n(2)),lowerRight) = edgesgn'; % bottom boundary
%         gn((n(2)+1):n(2):(n(1)*n(2)),lowerLeft) = edgesgn'; % bottom boundary
%     
%         edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)-1);
%         gn((1:n(2):((n(1)-1)*n(2)))+n(2)-1,upperRight) = edgesgn'; % top boundary
%         gn(((n(2)+1):n(2):(n(1)*n(2)))+n(2)-1,upperLeft) = edgesgn'; % top boundary
%     end
%     
%     if n(2)>1
%         edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)-1);
%         gn(1:(n(2)-1),upperLeft) = edgesgn'; % left boundary
%         gn(2:n(2),lowerLeft) = edgesgn'; % left boundary
% 
%         edgesgn = (edgesgn(end)+1):(edgesgn(end)+n(2)-1);
%         gn(((n(1)-1)*n(2)+1):(n(1)*n(2)-1),upperRight) = edgesgn'; % right boundary
%         gn(((n(1)-1)*n(2)+2):(n(1)*n(2)),lowerRight) = edgesgn'; % right boundary
%     end
%     
%     if nElements>1
%         % now the edges at the corners of the boundary
%         edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+4);
%     else
%         edgesgn = uint32(nInternalEdges+1):uint32(nInternalEdges+4);
%     end
%     
%     gn(1,lowerLeft) = edgesgn(1);
%     gn(n(2),upperLeft) = edgesgn(2);
%     gn((n(1)-1)*n(2)+1,lowerRight) = edgesgn(3);
%     gn(n(1)*n(2),upperRight) = edgesgn(4);

%     edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)*(p-1));% left boundary
%     gn(1:n(2),leftEdgeEdgesInternalNumber) = reshape(edgesgn,[(p-1) n(2)])'; % left boundary
% 
%     edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)*(p-1));% right boundary
%     gn(((n(1)-1)*n(2)+1):(n(1)*n(2)),rightEdgeEdgesInternalNumber) = reshape(edgesgn,[(p-1) n(2)])'; % right boundary

    
%--------------------------------------------------------------------------
% compute the numbering for the eta basis functions part
%--------------------------------------------------------------------------
    
    % store the last number given to a dXi edge
    lastNumber = edgesgn(end);
    
    % compute the number of internal edges
    nInternalEdges = p*(p-1)*nElements;

    % compute the number of internal edges per element
    nElementInternalEdges = p*(p-1);

    % compute the total number of edges
    nEdges = (n(1)*p+1)*(n(2)*(p-1));

    % compute the number of edges in each element
    nElementEdges = p*(p+1);

    % compute the number of edges in the boundary
    nBoundaryEdges = 2*n(1)*p;

    % it was already allocated
    
%     % allocate memory space for global numbering matrix
%     % Notice that edges in elements are numbered from bottom to top and from
%     % left to right.
%     % Elements follow the same kind of numbering: bottom to
%     % top and left to right.
%     gn = zeros(nElements, nElementEdges, 'uint32');

    % number internal edges

    % first select the internal edges of the elements
    internalEdgesColumns = zeros(nElementInternalEdges,1,'uint32');
    for i=2:p
        for j=1:p
            internalEdgesColumns((i-2)*p+j) = (i-2)*p + j + p;
        end
    end
    
    gn(:,internalEdgesColumns+nElementEdges) = reshape(1:uint32(nInternalEdges),[p*(p-1) nElements])' + lastNumber;
    
    edgesgn = uint32(nInternalEdges+lastNumber); % just to keep track in which number we are.
    
    % there are no corner edges...
    
%     % number corner boundary edges of the interior elements
% 
%     % compute the global numbering of the edges
%     edgesgn = uint32(nInternalEdges+1):uint32(nInternalEdges+(n(1)-1)*(n(2)-1));
%     
%     % compute the internal elements
%     internalElements = zeros(nInternalElements,1,'uint32');
%     for i=2:(n(1)-1)
%         for j=2:(n(2)-1)
%             internalElements((i-2)*(n(2)-2)+j-1) = (i-1)*(n(2)) + j;
%         end 
%     end
% 
%     % compute the internal numbering of the corner edges
%     cornerEdgesInternalNumber = [1 p+1 p*(p+1)+1 (p+1)*(p+1)];
% 
%     % number them
%     if nInternalElements >= 1
%         % first the ones of the inner elements
%         for k=1:(n(1)-2)
%             gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerEdgesInternalNumber(1)) = edgesgn(((k-1)*(n(2)-1)+1):(k*(n(2)-1)-1));
%             gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerEdgesInternalNumber(2)) = edgesgn(((k-1)*(n(2)-1)+1+1):(k*(n(2)-1)-1+1));
%             gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerEdgesInternalNumber(3)) = edgesgn(((k)*(n(2)-1)+1):((k+1)*(n(2)-1)-1));
%             gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerEdgesInternalNumber(4)) = edgesgn(((k)*(n(2)-1)+1+1):((k+1)*(n(2)-1)-1+1));
%         end
% 
%         if n(1)>1
%             % the ones of the elements on the left boundary
%             k=1;
%             gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))-n(2), cornerEdgesInternalNumber(3)) = edgesgn(((k-1)*(n(2)-1)+1):(k*(n(2)-1)-1));
%             gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))-n(2), cornerEdgesInternalNumber(4)) = edgesgn(((k-1)*(n(2)-1)+1+1):(k*(n(2)-1)-1+1));
% 
%             % the ones of the elements on the right boundary
%             k=n(1)-2;
%             gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))+n(2), cornerEdgesInternalNumber(1)) = edgesgn(((k)*(n(2)-1)+1):((k+1)*(n(2)-1)-1));
%             gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))+n(2), cornerEdgesInternalNumber(2)) = edgesgn(((k)*(n(2)-1)+1+1):((k+1)*(n(2)-1)-1+1));
%         end
% 
%         if n(2)>1
%             % the ones of the elements on the bottom boundary
%             elements = 1:n(2)-2:((n(1)-2)*(n(2)-2));
%             gn(internalElements(elements)-1, cornerEdgesInternalNumber(2)) = edgesgn((1:n(2)-1:((n(1)-2)*(n(2)-1))));
%             gn(internalElements(elements)-1, cornerEdgesInternalNumber(4)) = edgesgn((1:n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);
% 
%             % the ones of the elements on the top boundary
%             elements = elements+n(2)-3;
%             gn(internalElements(elements)+1, cornerEdgesInternalNumber(1)) = edgesgn((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)));
%             gn(internalElements(elements)+1, cornerEdgesInternalNumber(3)) = edgesgn(((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);
%         end
%         % the four corner elements
%         gn(1,cornerEdgesInternalNumber(4)) = edgesgn(1); % left bottom
%         if n(2)>1
%             gn(n(2),cornerEdgesInternalNumber(3)) = edgesgn(n(2)-1); % left top
%         end
%         if n(1)>1
%             gn((n(1)-1)*n(2)+1,cornerEdgesInternalNumber(2)) = edgesgn((n(1)-2)*(n(2)-1)+1); % right bottom
%         end
%         gn(n(1)*n(2),cornerEdgesInternalNumber(1)) = edgesgn(end); % right top
%     end
    
    % number the edges of the internal edges
    
    leftEdgeEdgesInternalNumber = 1:p;
    rightEdgeEdgesInternalNumber = leftEdgeEdgesInternalNumber + p*p;
    
%     bottomEdgeEdgesInternalNumber = 1:p+1:(p*(p+1));
%     topEdgeEdgesInternalNumber = bottomEdgeEdgesInternalNumber + p;
    
    % vertical edges are not in dxi one forms
    
    % start with the vertical edges
    if n(1)~=1
        if nInternalElements >= 1
            edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+p*n(2)*(n(1)-1));
        else
            edgesgn = uint32(nInternalEdges+1+lastNumber):uint32(nInternalEdges+p*n(2)*(n(1)-1)+lastNumber);
        end

        for k=1:(n(1)-1)
            gn(((k-1)*n(2)+1):k*n(2),rightEdgeEdgesInternalNumber+nElementEdges) = reshape(edgesgn(((k-1)*p*n(2)+1):(k*n(2)*p)),[p n(2)])';
            gn((k*n(2)+1):(k+1)*n(2),leftEdgeEdgesInternalNumber+nElementEdges) = reshape(edgesgn(((k-1)*p*n(2)+1):(k*n(2)*p)),[p n(2)])';
        end
    end
     
    % horizontal edges are not in deta one forms
    
%     % now the horizontal ones
%     if n(2)~=1
%         if nInternalElements >= 1
%             edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+p*n(1)*(n(2)-1));
%         else
%             edgesgn = uint32(nInternalEdges+1):uint32(nInternalEdges+p*n(1)*(n(2)-1));
%         end
% 
%         for k=1:(n(2)-1)
%             gn(k:n(2):(n(1)*n(2)),topEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*p*n(1)+1):(k*n(1)*p)),[p n(1)])';
%             gn((k:n(2):(n(1)*n(2)))+1,bottomEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*p*n(1)+1):(k*n(1)*p)),[p n(1)])';
%         end
%     end
    
    % finally number the edges at the actual boundary
    
    % there are no corner edges
%     lowerLeft = 1;
%     upperLeft = p+1;
%     lowerRight = p*(p+1)+1;
%     upperRight = (p+1)*(p+1);
% 
%     % start with the edges of the elements corners except the corner edges of
%     % the domain
%     if n(1)>1
%         edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)-1);
%         gn(1:n(2):((n(1)-1)*n(2)),lowerRight) = edgesgn'; % bottom boundary
%         gn((n(2)+1):n(2):(n(1)*n(2)),lowerLeft) = edgesgn'; % bottom boundary
%     
%         edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)-1);
%         gn((1:n(2):((n(1)-1)*n(2)))+n(2)-1,upperRight) = edgesgn'; % top boundary
%         gn(((n(2)+1):n(2):(n(1)*n(2)))+n(2)-1,upperLeft) = edgesgn'; % top boundary
%     end
%     
%     if n(2)>1
%         edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)-1);
%         gn(1:(n(2)-1),upperLeft) = edgesgn'; % left boundary
%         gn(2:n(2),lowerLeft) = edgesgn'; % left boundary
% 
%         edgesgn = (edgesgn(end)+1):(edgesgn(end)+n(2)-1);
%         gn(((n(1)-1)*n(2)+1):(n(1)*n(2)-1),upperRight) = edgesgn'; % right boundary
%         gn(((n(1)-1)*n(2)+2):(n(1)*n(2)),lowerRight) = edgesgn'; % right boundary
%     end
%     
%     if nElements>1
%         % now the edges at the corners of the boundary
%         edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+4);
%     else
%         edgesgn = uint32(nInternalEdges+1):uint32(nInternalEdges+4);
%     end
%     
%     gn(1,lowerLeft) = edgesgn(1);
%     gn(n(2),upperLeft) = edgesgn(2);
%     gn((n(1)-1)*n(2)+1,lowerRight) = edgesgn(3);
%     gn(n(1)*n(2),upperRight) = edgesgn(4);

%--------------------------------------------------------------------------
% the boundary edges are left to the end in order to take care of Neumann
% boundary conditions
%--------------------------------------------------------------------------

    % now the edges at the boundary edges
    edgesgn =  uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)*p);% lower boundary
    gn(1:n(2):((n(1)-1)*n(2))+1,bottomEdgeEdgesInternalNumber) = reshape(edgesgn,[p n(1)])'; % lower boundary

    edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)*p);% upper boundary
    gn((1:n(2):((n(1)-1)*n(2))+1)+n(2)-1,topEdgeEdgesInternalNumber) = reshape(edgesgn,[p n(1)])'; % upper boundary
    edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)*p);% left boundary
    gn(1:n(2),leftEdgeEdgesInternalNumber+nElementEdges) = reshape(edgesgn,[p n(2)])'; % left boundary

    edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)*p);% right boundary
    gn(((n(1)-1)*n(2)+1):(n(1)*n(2)),rightEdgeEdgesInternalNumber+nElementEdges) = reshape(edgesgn,[p n(2)])'; % right boundary



end