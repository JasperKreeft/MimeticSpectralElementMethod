function Gn = GlobalNumberingVectorValuedOneFormPrimal(n,p,varargin)
% GlobalNumberingOneFormPrimal computes the global numbering of 1-form
% degrees of freedom on a  primal mesh.
%
%   Gn.Xi: Global numbering for 1-forms corresponding to finite-volumes
%   that envelope Xi edges on primal mesh.
%
%   Gn.Eta: Global numbering for 1-forms corresponding to finite-volumes
%   that envelope Eta edges on primal mesh.
%
%   periodic = [BottomTopFlag; LeftRightFlag]
%   
%   Copyright 2012 Deepesh Toshniwal
%   $Revision: 1.2 $  $Date: 2012/04/17 $
        
    if (size(varargin,2))
        periodic = varargin{1};
    else
        periodic = [false false];
    end
 
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
    
    %% Compute numbering for edges in Finite-Volume Type 1
    
    pY = p+1;
    pX = p;
    
    %%% XI %%%
    
    % compute the number of internal xi edges
    nInternalXiEdges = (pX)*(pY-1)*nElements;

    % compute the number of internal edges per element
    nElementInternalXiEdges = pX*(pY-1);

    % compute the total number of edges
    nXiEdges = (n(1)*pX)*(n(2)*pY+1);

    % compute the number of edges in each element
    nElementXiEdges = pX*(pY+1);

    % compute the number of edges in the boundary
    nBoundaryXiEdges = 2*n(1)*pX;
    
    %%% ETA %%%
    
    % compute the number of internal eta edges
    nInternalEtaEdges = (pY)*(pX-1)*nElements;

    % compute the number of internal edges per element
    nElementInternalEtaEdges = (pY)*(pX-1);

    % compute the total number of edges
    nEtaEdges = (n(2)*pY)*(n(1)*pX+1);

    % compute the number of edges in each element
    nElementEtaEdges = pY*(pX+1);

    % compute the number of edges in the boundary
    nBoundaryEtaEdges = 2*n(2)*pY;
    
    %%% Allocate memory space

    % allocate memory space for global numbering matrix
    % Notice that edges in elements are numbered from bottom to top and from
    % left to right.
    % Elements follow the same kind of numbering: bottom to
    % top and left to right.
    gn = zeros(nElements, nElementEtaEdges + nElementXiEdges, 'uint32');

    % number internal edges

    % first select the internal edges of the elements
    internalXiEdgesColumns = repmat((2:pY)',1,pX) + repmat(0:(pY+1):(pX-1)*(pY+1),pY-1,1);
    internalXiEdgesColumns = internalXiEdgesColumns(:);
    
    gn(:,internalXiEdgesColumns) = reshape(1:nInternalXiEdges,[nElementInternalXiEdges nElements])';
    
    edgesgn = nInternalXiEdges; % just to keep track in which number we are.
    
    % number the edges of the internal edges
    bottomEdgeEdgesInternalNumber = 1:(pY+1):(pX*(pY+1));
    topEdgeEdgesInternalNumber = bottomEdgeEdgesInternalNumber + pY;
    
    % now the horizontal ones
    if n(2)~=1
        if nInternalElements >= 1
            edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+pX*n(1)*(n(2)-1));
        else
            edgesgn = uint32(nInternalXiEdges+1):uint32(nInternalXiEdges+pX*n(1)*(n(2)-1));
        end

        for k=1:(n(2)-1)
            gn(k:n(2):(n(1)*n(2)),topEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*pX*n(1)+1):(k*n(1)*pX)),[pX n(1)])';
            gn((k:n(2):(n(1)*n(2)))+1,bottomEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*pX*n(1)+1):(k*n(1)*pX)),[pX n(1)])';
        end
    end
    
    % store the last number given to a dXi edge
    lastNumber = edgesgn(end);

    % first select the internal edges of the elements
    internalEtaEdgesColumns = pX*(pY+1) + pY + (1:(pX-1)*pY);
    gn(:,internalEtaEdgesColumns) = reshape(1:uint32(nInternalEtaEdges),[nElementInternalEtaEdges nElements])' + lastNumber;
    
    edgesgn = uint32(nInternalEtaEdges+lastNumber); % just to keep track in which number we are.
    
    % number the edges of the internal edges
    leftEdgeEdgesInternalNumber = 1:pY;
    rightEdgeEdgesInternalNumber = leftEdgeEdgesInternalNumber + pX*pY;
    
    % start with the vertical edges
    if n(1)~=1
        if nInternalElements >= 1
            edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+pY*n(2)*(n(1)-1));
        else
            edgesgn = uint32(nInternalEtaEdges+1+lastNumber):uint32(nInternalEtaEdges+pY*n(2)*(n(1)-1)+lastNumber);
        end

        for k=1:(n(1)-1)
            gn(((k-1)*n(2)+1):k*n(2),rightEdgeEdgesInternalNumber+nElementXiEdges) = reshape(edgesgn(((k-1)*pY*n(2)+1):(k*n(2)*pY)),[pY n(2)])';
            gn((k*n(2)+1):(k+1)*n(2),leftEdgeEdgesInternalNumber+nElementXiEdges) = reshape(edgesgn(((k-1)*pY*n(2)+1):(k*n(2)*pY)),[pY n(2)])';
        end
    end

    %--------------------------------------------------------------------------
    % the boundary edges are left to the end in order to take care of Neumann
    % boundary conditions
    %--------------------------------------------------------------------------

    % now the edges at the boundary edges
    edgesgn =  uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)*pX);% lower boundary
    gn(1:n(2):((n(1)-1)*n(2))+1,bottomEdgeEdgesInternalNumber) = reshape(edgesgn,[pX n(1)])'; % lower boundary

    if ~(periodic(1))
        edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)*pX);% upper boundary
    end
    gn((1:n(2):((n(1)-1)*n(2))+1)+n(2)-1,topEdgeEdgesInternalNumber) = reshape(edgesgn,[pX n(1)])'; % upper boundary
    
    edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)*pY);% left boundary
    gn(1:n(2),leftEdgeEdgesInternalNumber+nElementXiEdges) = reshape(edgesgn,[pY n(2)])'; % left boundary

    if ~(periodic(2))
        edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)*pY);% right boundary
    end
    gn(((n(1)-1)*n(2)+1):(n(1)*n(2)),rightEdgeEdgesInternalNumber+nElementXiEdges) = reshape(edgesgn,[pY n(2)])'; % right boundary
    
    Gn.Xi = gn;
    
    %% Compute numbering for edges in Finite-Volume Type 2
    
    pY = p;
    pX = p+1;
    
    %%% XI %%%
    
    % compute the number of internal xi edges
    nInternalXiEdges = (pX)*(pY-1)*nElements;

    % compute the number of internal edges per element
    nElementInternalXiEdges = pX*(pY-1);

    % compute the total number of edges
    nXiEdges = (n(1)*pX)*(n(2)*pY+1);

    % compute the number of edges in each element
    nElementXiEdges = pX*(pY+1);

    % compute the number of edges in the boundary
    nBoundaryXiEdges = 2*n(1)*pX;
    
    %%% ETA %%%
    
    % compute the number of internal eta edges
    nInternalEtaEdges = (pY)*(pX-1)*nElements;

    % compute the number of internal edges per element
    nElementInternalEtaEdges = (pY)*(pX-1);

    % compute the total number of edges
    nEtaEdges = (n(2)*pY)*(n(1)*pX+1);

    % compute the number of edges in each element
    nElementEtaEdges = pY*(pX+1);

    % compute the number of edges in the boundary
    nBoundaryEtaEdges = 2*n(2)*pY;
    
    %%% Allocate memory space

    % allocate memory space for global numbering matrix
    % Notice that edges in elements are numbered from bottom to top and from
    % left to right.
    % Elements follow the same kind of numbering: bottom to
    % top and left to right.
    gn = zeros(nElements, nElementEtaEdges + nElementXiEdges, 'uint32');

    % number internal edges

    % first select the internal edges of the elements
    internalXiEdgesColumns = repmat((2:pY)',1,pX) + repmat(0:(pY+1):(pX-1)*(pY+1),pY-1,1);
    internalXiEdgesColumns = internalXiEdgesColumns(:);
    
    gn(:,internalXiEdgesColumns) = reshape(1:nInternalXiEdges,[nElementInternalXiEdges nElements])';
    
    edgesgn = nInternalXiEdges; % just to keep track in which number we are.
    
    % number the edges of the internal edges
    bottomEdgeEdgesInternalNumber = 1:(pY+1):(pX*(pY+1));
    topEdgeEdgesInternalNumber = bottomEdgeEdgesInternalNumber + pY;
    
    % now the horizontal ones
    if n(2)~=1
        if nInternalElements >= 1
            edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+pX*n(1)*(n(2)-1));
        else
            edgesgn = uint32(nInternalXiEdges+1):uint32(nInternalXiEdges+pX*n(1)*(n(2)-1));
        end

        for k=1:(n(2)-1)
            gn(k:n(2):(n(1)*n(2)),topEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*pX*n(1)+1):(k*n(1)*pX)),[pX n(1)])';
            gn((k:n(2):(n(1)*n(2)))+1,bottomEdgeEdgesInternalNumber) = reshape(edgesgn(((k-1)*pX*n(1)+1):(k*n(1)*pX)),[pX n(1)])';
        end
    end
    
    % store the last number given to a dXi edge
    lastNumber = edgesgn(end);

    % first select the internal edges of the elements
    internalEtaEdgesColumns = pX*(pY+1) + pY + (1:(pX-1)*pY);
    gn(:,internalEtaEdgesColumns) = reshape(1:uint32(nInternalEtaEdges),[nElementInternalEtaEdges nElements])' + lastNumber;
    
    edgesgn = uint32(nInternalEtaEdges+lastNumber); % just to keep track in which number we are.
    
    % number the edges of the internal edges
    leftEdgeEdgesInternalNumber = 1:pY;
    rightEdgeEdgesInternalNumber = leftEdgeEdgesInternalNumber + pX*pY;
    
    % start with the vertical edges
    if n(1)~=1
        if nInternalElements >= 1
            edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+pY*n(2)*(n(1)-1));
        else
            edgesgn = uint32(nInternalEtaEdges+1+lastNumber):uint32(nInternalEtaEdges+pY*n(2)*(n(1)-1)+lastNumber);
        end

        for k=1:(n(1)-1)
            gn(((k-1)*n(2)+1):k*n(2),rightEdgeEdgesInternalNumber+nElementXiEdges) = reshape(edgesgn(((k-1)*pY*n(2)+1):(k*n(2)*pY)),[pY n(2)])';
            gn((k*n(2)+1):(k+1)*n(2),leftEdgeEdgesInternalNumber+nElementXiEdges) = reshape(edgesgn(((k-1)*pY*n(2)+1):(k*n(2)*pY)),[pY n(2)])';
        end
    end

    %--------------------------------------------------------------------------
    % the boundary edges are left to the end in order to take care of Neumann
    % boundary conditions
    %--------------------------------------------------------------------------

    % now the edges at the boundary edges
    edgesgn =  uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)*pX);% lower boundary
    gn(1:n(2):((n(1)-1)*n(2))+1,bottomEdgeEdgesInternalNumber) = reshape(edgesgn,[pX n(1)])'; % lower boundary

    if ~(periodic(1))
        edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(1)*pX);% upper boundary
    end
    gn((1:n(2):((n(1)-1)*n(2))+1)+n(2)-1,topEdgeEdgesInternalNumber) = reshape(edgesgn,[pX n(1)])'; % upper boundary
    
    edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)*pY);% left boundary
    gn(1:n(2),leftEdgeEdgesInternalNumber+nElementXiEdges) = reshape(edgesgn,[pY n(2)])'; % left boundary

    if ~(periodic(2))
        edgesgn = uint32(edgesgn(end)+1):uint32(edgesgn(end)+n(2)*pY);% right boundary
    end
    gn(((n(1)-1)*n(2)+1):(n(1)*n(2)),rightEdgeEdgesInternalNumber+nElementXiEdges) = reshape(edgesgn,[pY n(2)])'; % right boundary
    
    Gn.Eta = gn;
    
end