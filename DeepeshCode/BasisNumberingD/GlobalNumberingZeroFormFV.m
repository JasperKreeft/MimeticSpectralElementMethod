function gn = GlobalNumberingZeroFormFV(n,p,bcFlag)
%GlobalNumberingZeroFormPrimalPeriodic computes the global numbering of 0-form
%degrees of freedom on a  primal mesh. Note that the numbering enforces
%periodic boundary conditions on all boundaries (top/bottom and left/right).
%
%   gn = GlobalNumberingZeroFormPrimalPeriodic(n, p)
%
%   Where:
%       n   :: the number of elements. If n is a number it has n x n
%              elements, if it is a vector, then it has n(1) x n(2) elements.
%       p   :: the order of the 0-form approximation
%       bcFlag :: specifies where periodic boundary conditions are to be
%                 specified: [false true] means periodic left/right but not
%                 top/bottom; [true false] means periodic bottom/top but
%                 not left right; [true true] means periodic bottom/top and
%                 left/right
%
%   Returns the global numbering of the nodes of the mesh. gn is a matrix
%   (p+1)^2 x N, where N is the total number of elements.
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
%             GLOBALNUMBERINGZEROFORMPRIMAL

%   Copyright 2012 Deepesh Toshniwal
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

    %% XI FINITE VOLUMES
    
    pX = p;
    pY = p+1;
    
    % compute the number of internal nodes
    nInternalNodes = (pX-1)*(pY-1)*nElements;

    % compute the number of internal nodes per element
    nElementInternalNodes = (pX-1)*(pY-1);

    % compute the total number of nodes
    nNodes = (n(1)*pX+1)*(n(2)*pY+1);

    % compute the number of nodes in each element
    nElementNodes = (pX+1)*(pY+1);

    % compute the number of nodes in the boundary
    nBoundaryNodes = 2*n(1)*pX+2*n(2)*pY;

    % allocate memory space for global numbering matrix
    % Notice that nodes in elements are numbered from left to right and from
    % bottom to top.
    % Elements follow the same kind of numbering: left to right and bottom to
    % top.
    gn = zeros(nElementNodes,nElements);

    % number internal nodes

    % first select the internal nodes of the elements
    internalNodesColumns = zeros(nElementInternalNodes,1);
    for i=2:pX
        for j=2:pY
            internalNodesColumns((i-2)*(pY-1)+j-1) = (i-1)*(pY+1) + j;
        end
    end
    
    if nInternalNodes==0
        nodesgn=0;
    else
        % generate the global number
        nodesgn = 1:nInternalNodes;
        % number them
        gn(internalNodesColumns,:) = reshape(nodesgn,[(pX-1)*(pY-1) nElements]);
    end

    % number corner boundary nodes of the interior elements

    % compute the global numbering of the nodes
    if (n(1)-1)*(n(2)-1)>0
        nodesgn = uint32(nInternalNodes+1):uint32(nInternalNodes+(n(1)-1)*(n(2)-1));
    end
        
    % compute the internal elements
    internalElements = zeros(nInternalElements,1);
    for i=2:(n(1)-1)
        for j=2:(n(2)-1)
            internalElements((i-2)*(n(2)-2)+j-1) = (i-1)*(n(2)) + j;
        end 
    end

    % compute the internal numbering of the corner nodes
    cornerNodesInternalNumber = [1 pY+1 pX*(pY+1)+1 (pX+1)*(pY+1)];

    % number them
    if nInternalElements >= 1
        % first the ones of the inner elements
        for k=1:(n(1)-2)
            gn(cornerNodesInternalNumber(1),internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2)))) = nodesgn(((k-1)*(n(2)-1)+1):(k*(n(2)-1)-1));
            gn(cornerNodesInternalNumber(2),internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2)))) = nodesgn(((k-1)*(n(2)-1)+1+1):(k*(n(2)-1)-1+1));
            gn(cornerNodesInternalNumber(3),internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2)))) = nodesgn(((k)*(n(2)-1)+1):((k+1)*(n(2)-1)-1));
            gn(cornerNodesInternalNumber(4),internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2)))) = nodesgn(((k)*(n(2)-1)+1+1):((k+1)*(n(2)-1)-1+1));
        end     
    end
    
    if (n(1)>1) && (n(2)>2)
        % the ones of the elements on the left boundary
        gn(cornerNodesInternalNumber(3),2:(n(2)-1)) = nodesgn(1:(n(2)-2));
        gn(cornerNodesInternalNumber(4),2:(n(2)-1)) = nodesgn(2:(n(2)-1));

        % the ones of the elements on the right boundary
        gn(cornerNodesInternalNumber(1),((n(1)-1)*n(2)+2):(n(1)*n(2)-1)) = nodesgn(((n(2)-1)*(n(1)-2)+1):((n(2)-1)*(n(1)-1)-1));
        gn(cornerNodesInternalNumber(2),((n(1)-1)*n(2)+2):(n(1)*n(2)-1)) = nodesgn(((n(2)-1)*(n(1)-2)+2):((n(2)-1)*(n(1)-1)));
    end

    if (n(2)>1) && (n(1)>2)
        % the ones of the elements on the bottom boundary
        gn(cornerNodesInternalNumber(2),(n(2)+1):n(2):((n(1)-1)*n(2))) = nodesgn((1:n(2)-1:((n(1)-2)*(n(2)-1))));
        gn(cornerNodesInternalNumber(4),(n(2)+1):n(2):((n(1)-1)*n(2))) = nodesgn((1:n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);

        % the ones of the elements on the top boundary
        gn(cornerNodesInternalNumber(1),(2*n(2)):n(2):((n(1)-1)*n(2))) = nodesgn((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)));
        gn(cornerNodesInternalNumber(3),(2*n(2)):n(2):((n(1)-1)*n(2))) = nodesgn(((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);
    end 
    
    if (n(1)-1)*(n(2)-1)>0
       % the four corner elements
        gn(cornerNodesInternalNumber(4),1) = nodesgn(1); % left bottom
        if n(2)>1
            gn(cornerNodesInternalNumber(3),n(2)) = nodesgn(n(2)-1); % left top
        end
        if n(1)>1
            gn(cornerNodesInternalNumber(2),(n(1)-1)*n(2)+1) = nodesgn((n(1)-2)*(n(2)-1)+1); % right bottom
        end
        gn(cornerNodesInternalNumber(1),n(1)*n(2)) = nodesgn(end); % right top   
    end
    
    % number the nodes of the internal edges
    
    leftEdgeNodesInternalNumber = 2:pY;
    rightEdgeNodesInternalNumber = (pX*(pY+1)+2):((pX+1)*(pY+1)-1);
    bottomEdgeNodesInternalNumber = (pY+2):(pY+1):(pX*(pY+1));
    topEdgeNodesInternalNumber = bottomEdgeNodesInternalNumber + pY;
    
    if pY>1
        % start with the vertical edges
        if n(1)~=1
           % if nInternalElements >= 1
           %     nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p-1)*n(2)*(n(1)-1));
           % else
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(pY-1)*n(2)*(n(1)-1));
            %end

            for k=1:(n(1)-1)
                gn(rightEdgeNodesInternalNumber,((k-1)*n(2)+1):k*n(2)) = reshape(nodesgn(((k-1)*(pY-1)*n(2)+1):(k*n(2)*(pY-1))),[pY-1 n(2)]);
                gn(leftEdgeNodesInternalNumber,(k*n(2)+1):(k+1)*n(2)) = reshape(nodesgn(((k-1)*(pY-1)*n(2)+1):(k*n(2)*(pY-1))),[pY-1 n(2)]);
            end
        end
    end

    if pX>1
        % now the horizontal ones
        if n(2)~=1
            %if n(1)~=1
            %    nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p-1)*n(1)*(n(2)-1));
            %else
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(pX-1)*n(1)*(n(2)-1));
            %end

            for k=1:(n(2)-1)
                gn(topEdgeNodesInternalNumber,k:n(2):(n(1)*n(2))) = reshape(nodesgn(((k-1)*(pX-1)*n(1)+1):(k*n(1)*(pX-1))),[pX-1 n(1)]);
                gn(bottomEdgeNodesInternalNumber,(k:n(2):(n(1)*n(2)))+1) = reshape(nodesgn(((k-1)*(pX-1)*n(1)+1):(k*n(1)*(pX-1))),[pX-1 n(1)]);
            end
        end
    end
    % finally number the nodes at the actual boundary
    lowerLeft = 1;
    upperLeft = pY+1;
    lowerRight = pX*(pY+1)+1;
    upperRight = (pX+1)*(pY+1);

    % start with the nodes of the elements corners except the corner nodes of
    % the domain
%     if (p==1) && ((n(1)>1) || (n(2)>1))
%         nodesgn
%     end
    if n(1)>1
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)-1);
        gn(lowerRight,1:n(2):((n(1)-1)*n(2))) = nodesgn; % bottom boundary
        gn(lowerLeft,(n(2)+1):n(2):(n(1)*n(2))) = nodesgn; % bottom boundary
    
        if ~bcFlag(1)
            nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)-1);
        end
        
        gn(upperRight,(1:n(2):((n(1)-1)*n(2)))+n(2)-1) = nodesgn; % top boundary
        gn(upperLeft,((n(2)+1):n(2):(n(1)*n(2)))+n(2)-1) = nodesgn; % top boundary
    end
    
    if n(2)>1
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)-1);
        gn(upperLeft,1:(n(2)-1)) = nodesgn; % left boundary
        gn(lowerLeft,2:n(2)) = nodesgn; % left boundary

        if ~bcFlag(2)
            nodesgn = (nodesgn(end)+1):(nodesgn(end)+n(2)-1);
        end
        
        gn(upperRight,((n(1)-1)*n(2)+1):(n(1)*n(2)-1)) = nodesgn; % right boundary
        gn(lowerRight,((n(1)-1)*n(2)+2):(n(1)*n(2))) = nodesgn; % right boundary
    end
    
    % now the nodes at the corners of the boundary
    
    if nElements>1
        if ~bcFlag(1)
            if ~bcFlag(2)
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+4);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,n(2)) = nodesgn(2);
                gn(lowerRight,(n(1)-1)*n(2)+1) = nodesgn(3);
                gn(upperRight,n(1)*n(2)) = nodesgn(4);
            else
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+2);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,n(2)) = nodesgn(2);
                gn(lowerRight,(n(1)-1)*n(2)+1) = nodesgn(1);
                gn(upperRight,n(1)*n(2)) = nodesgn(2);
            end
        else
            if ~bcFlag(2)
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+2);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,n(2)) = nodesgn(1);
                gn(lowerRight,(n(1)-1)*n(2)+1) = nodesgn(2);
                gn(upperRight,n(1)*n(2)) = nodesgn(2);
            else
                nodesgn = uint32(nodesgn(end)+1);
                gn(lowerLeft,1) = nodesgn;
                gn(upperLeft,n(2)) = nodesgn;
                gn(lowerRight,(n(1)-1)*n(2)+1) = nodesgn;
                gn(upperRight,n(1)*n(2)) = nodesgn;
            end
        end
    else
        if ~bcFlag(1)
            if ~bcFlag(2)
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+4);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,1) = nodesgn(2);
                gn(lowerRight,1) = nodesgn(3);
                gn(upperRight,1) = nodesgn(4);
            else
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+2);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,1) = nodesgn(2);
                gn(lowerRight,1) = nodesgn(1);
                gn(upperRight,1) = nodesgn(2);
            end
        else
            if ~bcFlag(2)
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+2);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,1) = nodesgn(1);
                gn(lowerRight,1) = nodesgn(2);
                gn(upperRight,1) = nodesgn(2);
            else
                nodesgn = uint32(nodesgn(end)+1);
                gn(lowerLeft,1) = nodesgn;
                gn(upperLeft,1) = nodesgn;
                gn(lowerRight,1) = nodesgn;
                gn(upperRight,1) = nodesgn;
            end
        end
    end
    
    if pX>1
        % now the nodes at the boundary edges
        nodesgn =  uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)*(pX-1));% lower boundary
        gn(bottomEdgeNodesInternalNumber,1:n(2):((n(1)-1)*n(2))+1) = reshape(nodesgn,[(pX-1) n(1)]); % lower boundary
        
        if ~bcFlag(1)
            nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)*(pX-1));% upper boundary
        end
        gn(topEdgeNodesInternalNumber,(1:n(2):((n(1)-1)*n(2))+1)+n(2)-1) = reshape(nodesgn,[(pX-1) n(1)]); % upper boundary

        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)*(pY-1));% left boundary
        gn(leftEdgeNodesInternalNumber,1:n(2)) = reshape(nodesgn,[(pY-1) n(2)]); % left boundary
        
        if ~bcFlag(2)
            nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)*(pY-1));% right boundary
        end
        gn(rightEdgeNodesInternalNumber,((n(1)-1)*n(2)+1):(n(1)*n(2))) = reshape(nodesgn,[(pY-1) n(2)]); % right boundary
    end
    
    gnXi = gn';
    
    %% ETA FINITE VOLUMES
    
    pX = p+1;
    pY = p;
    
    % compute the number of internal nodes
    nInternalNodes = (pX-1)*(pY-1)*nElements;

    % compute the number of internal nodes per element
    nElementInternalNodes = (pX-1)*(pY-1);

    % compute the total number of nodes
    nNodes = (n(1)*pX+1)*(n(2)*pY+1);

    % compute the number of nodes in each element
    nElementNodes = (pX+1)*(pY+1);

    % compute the number of nodes in the boundary
    nBoundaryNodes = 2*n(1)*pX+2*n(2)*pY;

    % allocate memory space for global numbering matrix
    % Notice that nodes in elements are numbered from left to right and from
    % bottom to top.
    % Elements follow the same kind of numbering: left to right and bottom to
    % top.
    gn = zeros(nElementNodes,nElements);

    % number internal nodes

    % first select the internal nodes of the elements
    internalNodesColumns = zeros(nElementInternalNodes,1);
    for i=2:pX
        for j=2:pY
            internalNodesColumns((i-2)*(pY-1)+j-1) = (i-1)*(pY+1) + j;
        end
    end
    
    if nInternalNodes==0
        nodesgn=0;
    else
        % generate the global number
        nodesgn = 1:nInternalNodes;
        % number them
        gn(internalNodesColumns,:) = reshape(nodesgn,[(pX-1)*(pY-1) nElements]);
    end

    % number corner boundary nodes of the interior elements

    % compute the global numbering of the nodes
    if (n(1)-1)*(n(2)-1)>0
        nodesgn = uint32(nInternalNodes+1):uint32(nInternalNodes+(n(1)-1)*(n(2)-1));
    end
        
    % compute the internal elements
    internalElements = zeros(nInternalElements,1);
    for i=2:(n(1)-1)
        for j=2:(n(2)-1)
            internalElements((i-2)*(n(2)-2)+j-1) = (i-1)*(n(2)) + j;
        end 
    end

    % compute the internal numbering of the corner nodes
    cornerNodesInternalNumber = [1 pY+1 pX*(pY+1)+1 (pX+1)*(pY+1)];

    % number them
    if nInternalElements >= 1
        % first the ones of the inner elements
        for k=1:(n(1)-2)
            gn(cornerNodesInternalNumber(1),internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2)))) = nodesgn(((k-1)*(n(2)-1)+1):(k*(n(2)-1)-1));
            gn(cornerNodesInternalNumber(2),internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2)))) = nodesgn(((k-1)*(n(2)-1)+1+1):(k*(n(2)-1)-1+1));
            gn(cornerNodesInternalNumber(3),internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2)))) = nodesgn(((k)*(n(2)-1)+1):((k+1)*(n(2)-1)-1));
            gn(cornerNodesInternalNumber(4),internalElements(((k-1)*(n(2)-2)+1):(k*(n(2)-2)))) = nodesgn(((k)*(n(2)-1)+1+1):((k+1)*(n(2)-1)-1+1));
        end     
    end
    
    if (n(1)>1) && (n(2)>2)
        % the ones of the elements on the left boundary
        gn(cornerNodesInternalNumber(3),2:(n(2)-1)) = nodesgn(1:(n(2)-2));
        gn(cornerNodesInternalNumber(4),2:(n(2)-1)) = nodesgn(2:(n(2)-1));

        % the ones of the elements on the right boundary
        gn(cornerNodesInternalNumber(1),((n(1)-1)*n(2)+2):(n(1)*n(2)-1)) = nodesgn(((n(2)-1)*(n(1)-2)+1):((n(2)-1)*(n(1)-1)-1));
        gn(cornerNodesInternalNumber(2),((n(1)-1)*n(2)+2):(n(1)*n(2)-1)) = nodesgn(((n(2)-1)*(n(1)-2)+2):((n(2)-1)*(n(1)-1)));
    end

    if (n(2)>1) && (n(1)>2)
        % the ones of the elements on the bottom boundary
        gn(cornerNodesInternalNumber(2),(n(2)+1):n(2):((n(1)-1)*n(2))) = nodesgn((1:n(2)-1:((n(1)-2)*(n(2)-1))));
        gn(cornerNodesInternalNumber(4),(n(2)+1):n(2):((n(1)-1)*n(2))) = nodesgn((1:n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);

        % the ones of the elements on the top boundary
        gn(cornerNodesInternalNumber(1),(2*n(2)):n(2):((n(1)-1)*n(2))) = nodesgn((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)));
        gn(cornerNodesInternalNumber(3),(2*n(2)):n(2):((n(1)-1)*n(2))) = nodesgn(((1+n(2)-2):n(2)-1:((n(1)-2)*(n(2)-1)))+n(2)-1);
    end 
    
    if (n(1)-1)*(n(2)-1)>0
       % the four corner elements
        gn(cornerNodesInternalNumber(4),1) = nodesgn(1); % left bottom
        if n(2)>1
            gn(cornerNodesInternalNumber(3),n(2)) = nodesgn(n(2)-1); % left top
        end
        if n(1)>1
            gn(cornerNodesInternalNumber(2),(n(1)-1)*n(2)+1) = nodesgn((n(1)-2)*(n(2)-1)+1); % right bottom
        end
        gn(cornerNodesInternalNumber(1),n(1)*n(2)) = nodesgn(end); % right top   
    end
    
    % number the nodes of the internal edges
    
    leftEdgeNodesInternalNumber = 2:pY;
    rightEdgeNodesInternalNumber = (pX*(pY+1)+2):((pX+1)*(pY+1)-1);
    bottomEdgeNodesInternalNumber = (pY+2):(pY+1):(pX*(pY+1));
    topEdgeNodesInternalNumber = bottomEdgeNodesInternalNumber + pY;
    
    if pY>1
        % start with the vertical edges
        if n(1)~=1
           % if nInternalElements >= 1
           %     nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p-1)*n(2)*(n(1)-1));
           % else
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(pY-1)*n(2)*(n(1)-1));
            %end

            for k=1:(n(1)-1)
                gn(rightEdgeNodesInternalNumber,((k-1)*n(2)+1):k*n(2)) = reshape(nodesgn(((k-1)*(pY-1)*n(2)+1):(k*n(2)*(pY-1))),[pY-1 n(2)]);
                gn(leftEdgeNodesInternalNumber,(k*n(2)+1):(k+1)*n(2)) = reshape(nodesgn(((k-1)*(pY-1)*n(2)+1):(k*n(2)*(pY-1))),[pY-1 n(2)]);
            end
        end
    end

    if pX>1
        % now the horizontal ones
        if n(2)~=1
            %if n(1)~=1
            %    nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(p-1)*n(1)*(n(2)-1));
            %else
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+(pX-1)*n(1)*(n(2)-1));
            %end

            for k=1:(n(2)-1)
                gn(topEdgeNodesInternalNumber,k:n(2):(n(1)*n(2))) = reshape(nodesgn(((k-1)*(pX-1)*n(1)+1):(k*n(1)*(pX-1))),[pX-1 n(1)]);
                gn(bottomEdgeNodesInternalNumber,(k:n(2):(n(1)*n(2)))+1) = reshape(nodesgn(((k-1)*(pX-1)*n(1)+1):(k*n(1)*(pX-1))),[pX-1 n(1)]);
            end
        end
    end
    % finally number the nodes at the actual boundary
    lowerLeft = 1;
    upperLeft = pY+1;
    lowerRight = pX*(pY+1)+1;
    upperRight = (pX+1)*(pY+1);

    % start with the nodes of the elements corners except the corner nodes of
    % the domain
%     if (p==1) && ((n(1)>1) || (n(2)>1))
%         nodesgn
%     end
    if n(1)>1
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)-1);
        gn(lowerRight,1:n(2):((n(1)-1)*n(2))) = nodesgn; % bottom boundary
        gn(lowerLeft,(n(2)+1):n(2):(n(1)*n(2))) = nodesgn; % bottom boundary
    
        if ~bcFlag(1)
            nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)-1);
        end
        
        gn(upperRight,(1:n(2):((n(1)-1)*n(2)))+n(2)-1) = nodesgn; % top boundary
        gn(upperLeft,((n(2)+1):n(2):(n(1)*n(2)))+n(2)-1) = nodesgn; % top boundary
    end
    
    if n(2)>1
        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)-1);
        gn(upperLeft,1:(n(2)-1)) = nodesgn; % left boundary
        gn(lowerLeft,2:n(2)) = nodesgn; % left boundary

        if ~bcFlag(2)
            nodesgn = (nodesgn(end)+1):(nodesgn(end)+n(2)-1);
        end
        
        gn(upperRight,((n(1)-1)*n(2)+1):(n(1)*n(2)-1)) = nodesgn; % right boundary
        gn(lowerRight,((n(1)-1)*n(2)+2):(n(1)*n(2))) = nodesgn; % right boundary
    end
    
    % now the nodes at the corners of the boundary
    
    if nElements>1
        if ~bcFlag(1)
            if ~bcFlag(2)
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+4);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,n(2)) = nodesgn(2);
                gn(lowerRight,(n(1)-1)*n(2)+1) = nodesgn(3);
                gn(upperRight,n(1)*n(2)) = nodesgn(4);
            else
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+2);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,n(2)) = nodesgn(2);
                gn(lowerRight,(n(1)-1)*n(2)+1) = nodesgn(1);
                gn(upperRight,n(1)*n(2)) = nodesgn(2);
            end
        else
            if ~bcFlag(2)
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+2);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,n(2)) = nodesgn(1);
                gn(lowerRight,(n(1)-1)*n(2)+1) = nodesgn(2);
                gn(upperRight,n(1)*n(2)) = nodesgn(2);
            else
                nodesgn = uint32(nodesgn(end)+1);
                gn(lowerLeft,1) = nodesgn;
                gn(upperLeft,n(2)) = nodesgn;
                gn(lowerRight,(n(1)-1)*n(2)+1) = nodesgn;
                gn(upperRight,n(1)*n(2)) = nodesgn;
            end
        end
    else
        if ~bcFlag(1)
            if ~bcFlag(2)
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+4);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,1) = nodesgn(2);
                gn(lowerRight,1) = nodesgn(3);
                gn(upperRight,1) = nodesgn(4);
            else
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+2);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,1) = nodesgn(2);
                gn(lowerRight,1) = nodesgn(1);
                gn(upperRight,1) = nodesgn(2);
            end
        else
            if ~bcFlag(2)
                nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+2);
                gn(lowerLeft,1) = nodesgn(1);
                gn(upperLeft,1) = nodesgn(1);
                gn(lowerRight,1) = nodesgn(2);
                gn(upperRight,1) = nodesgn(2);
            else
                nodesgn = uint32(nodesgn(end)+1);
                gn(lowerLeft,1) = nodesgn;
                gn(upperLeft,1) = nodesgn;
                gn(lowerRight,1) = nodesgn;
                gn(upperRight,1) = nodesgn;
            end
        end
    end
    
    if pX>1
        % now the nodes at the boundary edges
        nodesgn =  uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)*(pX-1));% lower boundary
        gn(bottomEdgeNodesInternalNumber,1:n(2):((n(1)-1)*n(2))+1) = reshape(nodesgn,[(pX-1) n(1)]); % lower boundary
        
        if ~bcFlag(1)
            nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(1)*(pX-1));% upper boundary
        end
        gn(topEdgeNodesInternalNumber,(1:n(2):((n(1)-1)*n(2))+1)+n(2)-1) = reshape(nodesgn,[(pX-1) n(1)]); % upper boundary

        nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)*(pY-1));% left boundary
        gn(leftEdgeNodesInternalNumber,1:n(2)) = reshape(nodesgn,[(pY-1) n(2)]); % left boundary
        
        if ~bcFlag(2)
            nodesgn = uint32(nodesgn(end)+1):uint32(nodesgn(end)+n(2)*(pY-1));% right boundary
        end
        gn(rightEdgeNodesInternalNumber,((n(1)-1)*n(2)+1):(n(1)*n(2))) = reshape(nodesgn,[(pY-1) n(2)]); % right boundary
    end
    
    gnEta = gn';
    
    clear gn
    
    gn.Xi = gnXi;
    gn.Eta = gnEta;
    
end