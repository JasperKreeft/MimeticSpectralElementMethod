function gn = GlobalNumberingZeroFormPrimal(n,p)
%GlobalNumberingZeroFormPrimal computes the global numbering of 0-form
%degrees of freedom on a  primal mesh.
%
%   gn = GlobalNumberingOneFormPrimal(n, p)
%
%   Where:
%       n   :: the number of elements. If n is a number it has n x n
%              elements, if it is a vector, then it has n(1) x n(2) elements.
%       p   :: the order of the 0-form approximation
%
%   Returns the global numbering of the nodes of the mesh. gn is a matrix
%   N x (p+1)^2, where N is the total number of elements.
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
%       2- internal corner nodes, from bottom to top and left to right;
%       3- nodes at vertical edges, from bottom to top and left to right;
%       4- nodes at horizontal edges, from left to right and bottom to top;
%       5- corner nodes at the boundary: bottom boundary from left to
%          right, top boundary, from left to right, left boundary from bottom
%          to top, right boundary from bottom to top;
%       6- nodes at edges of the boundaries: bottom boundary from left to
%          right, top boundary, from left to right, left boundary from bottom
%          to top, right boundary from bottom to top;
%
%   See also: GLOBALNUMBERINGZEROFORMDUAL, PLOTGLOBALNUMBERINZEROFORMPRIMAL

%   Copyright 2010 Artur Palha
%   $Revision: 1.2 $  $Date: 2010/03/02 $
    
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
    nInternalNodes = (p-1)*(p-1)*nElements;

    % compute the number of internal nodes per element
    nElementInternalNodes = (p-1)*(p-1);

    % compute the total number of nodes
    nNodes = (n(1)*p+1)*(n(2)*p+1);

    % compute the number of nodes in each element
    nElementNodes = (p+1)*(p+1);

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
    internalNodesColumns = zeros(nElementInternalNodes,1,'uint32');
    for i=2:p
        for j=2:p
            internalNodesColumns((i-2)*(p-1)+j-1) = (i-1)*(p+1) + j;
        end
    end
    
    if nInternalNodes==0
        nodesgn=0;
    else
        % generate the global number
        nodesgn = 1:nInternalNodes;
        % number them
        gn(:,internalNodesColumns) = reshape(nodesgn,[(p-1)*(p-1) nElements])';
    end

    % number corner boundary nodes of the interior elements

    % compute the global numbering of the nodes
    if (n(1)-1)*(n(2)-1)>0
        nodesgn = uint32(nInternalNodes+1):uint32(nInternalNodes+(n(1)-1)*(n(2)-1));
    end
        
    % compute the internal elements
    internalElements = zeros(nInternalElements,1,'uint32');
    for i=2:(n(1)-1)
        for j=2:(n(2)-1)
            internalElements((i-2)*(n(2)-2)+j-1) = (i-1)*(n(2)) + j;
        end 
    end

    % compute the internal numbering of the corner nodes
    cornerNodesInternalNumber = [1 p+1 p*(p+1)+1 (p+1)*(p+1)];

    % number them
    if nInternalElements >= 1
        % first the ones of the inner elements
        for k=1:(n(1)-2)
            gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerNodesInternalNumber(1)) = nodesgn(((k-1)*(n(2)-1)+1):(k*(n(2)-1)-1));
            gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerNodesInternalNumber(2)) = nodesgn(((k-1)*(n(2)-1)+1+1):(k*(n(2)-1)-1+1));
            gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerNodesInternalNumber(3)) = nodesgn(((k)*(n(2)-1)+1):((k+1)*(n(2)-1)-1));
            gn(internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2))), cornerNodesInternalNumber(4)) = nodesgn(((k)*(n(2)-1)+1+1):((k+1)*(n(2)-1)-1+1));
        end     
    end
    
    if (n(1)>1) && (n(2)>2)
        % the ones of the elements on the left boundary
        gn(2:(n(2)-1), cornerNodesInternalNumber(3)) = nodesgn(1:(n(2)-2));
        gn(2:(n(2)-1), cornerNodesInternalNumber(4)) = nodesgn(2:(n(2)-1));

        % the ones of the elements on the right boundary
        gn(((n(1)-1)*n(2)+2):(n(1)*n(2)-1), cornerNodesInternalNumber(1)) = nodesgn(((n(2)-1)*(n(1)-2)+1):((n(2)-1)*(n(1)-1)-1));
        gn(((n(1)-1)*n(2)+2):(n(1)*n(2)-1), cornerNodesInternalNumber(2)) = nodesgn(((n(2)-1)*(n(1)-2)+2):((n(2)-1)*(n(1)-1)));
    end

    if (n(2)>1) && (n(1)>2)
        % the ones of the elements on the bottom boundary
        gn((n(2)+1):n(2):((n(1)-1)*n(2)), cornerNodesInternalNumber(2)) = nodesgn((1:n(2)-1:((n(1)-2)*(n(2)-1))));
        gn((n(2)+1):n(2):((n(1)-1)*n(2)), cornerNodesInternalNumber(4)) = nodesgn((1:n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);

        % the ones of the elements on the top boundary
        gn((2*n(2)):n(2):((n(1)-1)*n(2)), cornerNodesInternalNumber(1)) = nodesgn((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)));
        gn((2*n(2)):n(2):((n(1)-1)*n(2)), cornerNodesInternalNumber(3)) = nodesgn(((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);
    end 
    
    if (n(1)-1)*(n(2)-1)>0
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
    
    % number the nodes of the internal edges
    
    leftEdgeNodesInternalNumber = 2:p;
    rightEdgeNodesInternalNumber = (p*(p+1)+2):((p+1)*(p+1)-1);
    bottomEdgeNodesInternalNumber = (p+2):p+1:(p*(p+1));
    topEdgeNodesInternalNumber = bottomEdgeNodesInternalNumber + p;
    
    if p>1
        % start with the vertical edges
        if n(1)~=1
           % if nInternalElements >= 1
           %     nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p-1)*n(2)*(n(1)-1));
           % else
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p-1)*n(2)*(n(1)-1));
            %end

            for k=1:(n(1)-1)
                gn(((k-1)*n(2)+1):k*n(2),rightEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p-1)*n(2)+1):(k*n(2)*(p-1))),[p-1 n(2)])';
                gn((k*n(2)+1):(k+1)*n(2),leftEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p-1)*n(2)+1):(k*n(2)*(p-1))),[p-1 n(2)])';
            end
        end

        % now the horizontal ones
        if n(2)~=1
            %if n(1)~=1
            %    nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p-1)*n(1)*(n(2)-1));
            %else
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p-1)*n(1)*(n(2)-1));
            %end

            for k=1:(n(2)-1)
                gn(k:n(2):(n(1)*n(2)),topEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p-1)*n(1)+1):(k*n(1)*(p-1))),[p-1 n(1)])';
                gn((k:n(2):(n(1)*n(2)))+1,bottomEdgeNodesInternalNumber) = reshape(nodesgn(((k-1)*(p-1)*n(1)+1):(k*n(1)*(p-1))),[p-1 n(1)])';
            end
        end
    end
    % finally number the nodes at the actual boundary
    lowerLeft = 1;
    upperLeft = p+1;
    lowerRight = p*(p+1)+1;
    upperRight = (p+1)*(p+1);

    % start with the nodes of the elements corners except the corner nodes of
    % the domain
%     if (p==1) && ((n(1)>1) || (n(2)>1))
%         nodesgn
%     end
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
        nodesgn = uint32(nInternalNodes+1):uint32(nInternalNodes+4);
    end
    
    gn(1,lowerLeft) = nodesgn(1);
    gn(n(2),upperLeft) = nodesgn(2);
    gn((n(1)-1)*n(2)+1,lowerRight) = nodesgn(3);
    gn(n(1)*n(2),upperRight) = nodesgn(4);
    
    if p>1
        % now the nodes at the boundary edges
        nodesgn =  uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)*(p-1));% lower boundary
        gn(1:n(2):((n(1)-1)*n(2))+1,bottomEdgeNodesInternalNumber) = reshape(nodesgn,[(p-1) n(1)])'; % lower boundary

        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)*(p-1));% upper boundary
        gn((1:n(2):((n(1)-1)*n(2))+1)+n(2)-1,topEdgeNodesInternalNumber) = reshape(nodesgn,[(p-1) n(1)])'; % upper boundary

        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)*(p-1));% left boundary
        gn(1:n(2),leftEdgeNodesInternalNumber) = reshape(nodesgn,[(p-1) n(2)])'; % left boundary

        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)*(p-1));% right boundary
        gn(((n(1)-1)*n(2)+1):(n(1)*n(2)),rightEdgeNodesInternalNumber) = reshape(nodesgn,[(p-1) n(2)])'; % right boundary
    end
end