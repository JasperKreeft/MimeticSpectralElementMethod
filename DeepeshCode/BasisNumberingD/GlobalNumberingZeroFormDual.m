function gn = GlobalNumberingZeroFormDual(n,p)
%GlobalNumberingZeroFormDual computes the global numbering of 0-form
%degrees of freedom on a  dual mesh.
%
%   gn = GlobalNumberingOneFormDual(n, p)
%
%   Where:
%       n   :: the number of elements. If n is a number it has n x n
%              elements, if it is a vector, then it has n(1) x n(2) elements.
%       p   :: the order of the 0-form approximation
%
%   Returns the global numbering of the nodes of the mesh. gn is a matrix
%   N x (p^2 + 4p), where N is the total number of elements.
%
%   Internal numbering:
%
%       Internally nodes are numbered from left to right and from bottom to
%       top. So the first node of element 1, that is gn(1,1) is the bottom
%       left node of the first element. The node gn(1,p+1) is the top left
%       node of the first element. The node gn(1,p*(p+1)+1) is the right
%       bottom node of the first element. The node gn(1, (p+1)*(p+1)) is
%       the top right node of the first element.
%
%   Global numbering:
%
%       1- internal nodes of the elements, following the internal numbering;
%       2- nodes at vertical edges, from bottom to top and left to right;
%       3- nodes at horizontal edges, from left to right and bottom to top;
%       4- nodes at edges of the boundaries: bottom boundary from left to
%          right, top boundary, from left to right, left boundary from bottom
%          to top, right boundary from bottom to top;
%
%   See also: GLOBALNUMBERINGZEROFORMPRIMAL, PLOTGLOBALNUMBERINZEROFORMDUAL

%   Copyright 2010 Artur Palha
%   $Revision: 1.2 $  $Date: 2010/03/08 $
    
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

    % compute the number of internal nodes
    nInternalNodes = p*p*nElements;

    % compute the number of internal nodes per element
    nElementInternalNodes = p*p;

    % compute the total number of nodes
    nNodes = (n(1)*p+1)*(n(2)*p+1) - (n(2)+1)*(n(1)+1);

    % compute the number of nodes in each element
    nElementNodes = (p+2)*(p+2);

    % compute the number of nodes in the boundary
    nBoundaryNodes = 2*n(1)*p+2*n(2)*p;

    % allocate memory space for global numbering matrix
    % Notice that nodes in elements are numbered from left to right and from
    % bottom to top.
    % Elements follow the same kind of numbering: left to right and bottom to
    % top.
    gn = zeros(nElements, nElementNodes, 'uint32');

    % number internal nodes

    % first select the internal nodes of the elements

    internalNodesColumns = uint32(1):uint32(nElementInternalNodes);

    gn(:,internalNodesColumns) = reshape(uint32(1):uint32(nInternalNodes),[p*p nElements])';
    
    % number the nodes of the internal edges
    
    leftEdgeNodesInternalNumber = (p*p+1):2:(p*p+2*p);
    rightEdgeNodesInternalNumber = (p*p+2):2:(p*p+2*p);
    bottomEdgeNodesInternalNumber = (p*p+2*p+1):2:(p*p+4*p);
    topEdgeNodesInternalNumber = (p*p+2*p+2):2:(p*p+4*p);
    
    nodesgn = uint32(nInternalNodes);
    
    % start with the vertical edges
    if n(1)~=1
        if nInternalElements >= 1
            nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p)*n(2)*(n(1)-1));
        else
            nodesgn = uint32(nInternalNodes+1):uint32(nInternalNodes+(p)*n(2)*(n(1)-1));
        end

        for k=1:(n(1)-1)
            gn(((k-1)*n(2)+1):k*n(2),rightEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p)*n(2)+1):(k*n(2)*(p))),[p n(2)])';
            gn((k*n(2)+1):(k+1)*n(2),leftEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p)*n(2)+1):(k*n(2)*(p))),[p n(2)])';
        end
    end
    
    % now the horizontal ones
    if n(2)~=1
        if n(1)~=1
            nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p)*n(1)*(n(2)-1));
        else
            nodesgn = uint32(nInternalNodes+1):uint32(nInternalNodes+(p)*n(1)*(n(2)-1));
        end

        for k=1:(n(2)-1)
            gn(k:n(2):(n(1)*n(2)),topEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p)*n(1)+1):(k*n(1)*(p))),[p n(1)])';
            gn((k:n(2):(n(1)*n(2)))+1,bottomEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p)*n(1)+1):(k*n(1)*(p))),[p n(1)])';
        end
    end

    % now the nodes at the boundary edges
    nodesgn =  uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)*(p));% lower boundary
    gn(1:n(2):((n(1)-1)*n(2))+1,bottomEdgeNodesInternalNumber) = reshape(nodesgn,[(p) n(1)])'; % lower boundary

    nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)*(p));% upper boundary
    gn((1:n(2):((n(1)-1)*n(2))+1)+n(2)-1,topEdgeNodesInternalNumber) = reshape(nodesgn,[(p) n(1)])'; % upper boundary

    nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)*(p));% left boundary
    gn(1:n(2),leftEdgeNodesInternalNumber) = reshape(nodesgn,[(p) n(2)])'; % left boundary

    nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)*(p));% right boundary
    gn(((n(1)-1)*n(2)+1):(n(1)*n(2)),rightEdgeNodesInternalNumber) = reshape(nodesgn,[(p) n(2)])'; % right boundary
    
    
    % number corner boundary nodes of the interior elements

    % compute the global numbering of the nodes
    if (n(1)-1)*(n(2)-1)~=0
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(n(1)-1)*(n(2)-1));
    end
    
    % compute the internal elements
    internalElements = zeros(nInternalElements,1,'uint32');
    for i=2:(n(1)-1)
        for j=2:(n(2)-1)
            internalElements((i-2)*(n(2)-2)+j-1) = (i-1)*(n(2)) + j;
        end 
    end

    % compute the internal numbering of the corner nodes
    cornerNodesInternalNumber = ((p+2)*(p+2)-3):((p+2)*(p+2));
    
    lowerLeft = cornerNodesInternalNumber(1);
    upperLeft = cornerNodesInternalNumber(2);
    lowerRight = cornerNodesInternalNumber(3); 
    upperRight = cornerNodesInternalNumber(4);

    % number them
    if nInternalElements >= 1
        % first the ones of the inner elements
        for k=1:(n(1)-2)
            gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerNodesInternalNumber(1)) = nodesgn(((k-1)*(n(2)-1)+1):(k*(n(2)-1)-1));
            gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerNodesInternalNumber(2)) = nodesgn(((k-1)*(n(2)-1)+1+1):(k*(n(2)-1)-1+1));
            gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerNodesInternalNumber(3)) = nodesgn(((k)*(n(2)-1)+1):((k+1)*(n(2)-1)-1));
            gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerNodesInternalNumber(4)) = nodesgn(((k)*(n(2)-1)+1+1):((k+1)*(n(2)-1)-1+1));
        end

        if n(1)>1
            % the ones of the elements on the left boundary
            k=1;
            gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))-n(2), cornerNodesInternalNumber(3)) = nodesgn(((k-1)*(n(2)-1)+1):(k*(n(2)-1)-1));
            gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))-n(2), cornerNodesInternalNumber(4)) = nodesgn(((k-1)*(n(2)-1)+1+1):(k*(n(2)-1)-1+1));

            % the ones of the elements on the right boundary
            k=n(1)-2;
            gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))+n(2), cornerNodesInternalNumber(1)) = nodesgn(((k)*(n(2)-1)+1):((k+1)*(n(2)-1)-1));
            gn(internalElements((((k-1)*(n(2)-2)+1):(k*(n(2)-2))))+n(2), cornerNodesInternalNumber(2)) = nodesgn(((k)*(n(2)-1)+1+1):((k+1)*(n(2)-1)-1+1));
        end

        if n(2)>1
            % the ones of the elements on the bottom boundary
            elements = 1:n(2)-2:((n(1)-2)*(n(2)-2));
            gn(internalElements(elements)-1, cornerNodesInternalNumber(2)) = nodesgn((1:n(2)-1:((n(1)-2)*(n(2)-1))));
            gn(internalElements(elements)-1, cornerNodesInternalNumber(4)) = nodesgn((1:n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);

            % the ones of the elements on the top boundary
            elements = elements+n(2)-3;
            gn(internalElements(elements)+1, cornerNodesInternalNumber(1)) = nodesgn((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)));
            gn(internalElements(elements)+1, cornerNodesInternalNumber(3)) = nodesgn(((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);
        end
    end
    
    if nElements>1
        % the four corner elements
        gn(1,cornerNodesInternalNumber(4)) = nodesgn(1); % left bottom
        if n(2)>1
            gn(n(2),cornerNodesInternalNumber(3)) = nodesgn(n(2)-1); % left top
        end
        if n(1)>1
            gn((n(1)-1)*n(2)+1,cornerNodesInternalNumber(2)) = nodesgn((n(1)-2)*(n(2)-1)+1); % right bottom
        end
        gn(n(1)*n(2),cornerNodesInternalNumber(1)) = nodesgn(end); % right top
    end
    
    % corner nodes of the boundary elements, the nodes in the interior of
    % the domain
    if nInternalElements == 0
        if n(1) == 2 && n(2) > 1
            gn(2:(n(2)-1),end-1) = nodesgn(1:(end-1));
            gn(2:(n(2)-1),end) = nodesgn(2:end);
            gn((n(2)+2):(2*n(2)-1),end-3) = nodesgn(1:(end-1));
            gn((n(2)+2):(2*n(2)-1),end-2) = nodesgn(2:end);
        end
        if n(2) == 2 && n(1) > 1
            gn(3:2:(2*n(1)-3),end-2) = nodesgn(1:(end-1));
            gn(3:2:(2*n(1)-3),end) = nodesgn(2:end);
            gn((3:2:(2*n(1)-3))+1,end-3) = nodesgn(1:(end-1));
            gn((3:2:(2*n(1)-3))+1,end-1) = nodesgn(2:end);
        end
        
    end
    
    % start with the nodes of the boundary elements corners except the corner nodes of
    % the domain
    if n(1)>1
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)-1);
        gn(1:n(2):((n(1)-1)*n(2)),lowerRight) = nodesgn'; % bottom boundary
        gn((n(2)+1):n(2):(n(1)*n(2)),lowerLeft) = nodesgn'; % bottom boundary
    
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)-1);
        gn((1:n(2):((n(1)-1)*n(2)))+n(2)-1,upperRight) = nodesgn'; % top boundary
        gn(((n(2)+1):n(2):(n(1)*n(2)))+n(2)-1,upperLeft) = nodesgn'; % top boundary
    end
    
    if n(2)>1
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)-1);
        gn(1:(n(2)-1),upperLeft) = nodesgn'; % left boundary
        gn(2:n(2),lowerLeft) = nodesgn'; % left boundary

        nodesgn = (nodesgn(end)+1):(nodesgn(end)+n(2)-1);
        gn(((n(1)-1)*n(2)+1):(n(1)*n(2)-1),upperRight) = nodesgn'; % right boundary
        gn(((n(1)-1)*n(2)+2):(n(1)*n(2)),lowerRight) = nodesgn'; % right boundary
    end
    
    % now the nodes at the corners of the boundary
    if nElements>1
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+4);
    else
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+4);
    end
    
    gn(1,cornerNodesInternalNumber(1)) = nodesgn(1);
    gn(n(2),cornerNodesInternalNumber(2)) = nodesgn(2);
    gn((n(1)-1)*n(2)+1,cornerNodesInternalNumber(3)) = nodesgn(3);
    gn(n(1)*n(2),cornerNodesInternalNumber(4)) = nodesgn(4);

end