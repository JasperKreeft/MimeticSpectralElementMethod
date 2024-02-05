function gn = GlobalNumberingZeroFormPrimal1D(n,p,varargin)
% GlobalNumberingZeroFormPrimal computes the global numbering of 0-form
% degrees of freedom on a  primal mesh.
%
%   gn = GlobalNumberingOneFormPrimal(n, p, options)
%
%   Where:
%       n   :: the number of elements. n is a number.
%       p   :: the order of the 0-form approximation.
%       
%   Optional inputs:
%       circular    ::  use circular boundary conditions or not.    
%
%   Returns the global numbering of the nodes of the mesh. gn is a matrix
%   N x (p+1), where N is the total number of elements.
%
%   Internal numbering: - Need to Review this
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

%   Copyright 2010 Artur Palha & Pedro Pinto Rebelo
%   $Revision: 1 $  $Date: 2011/04/13 $

    % check the number of optional arguments
    optargin = size(varargin,2);
    
    % if there is at least one optional argument use the first one as the
    % value for the circular variable
    if optargin > 0
        circular = varargin{1};
    else
        circular = false;
    end

    if length(n)>1
        % check if n has more than two values
        disp(sprintf(':: n :: has too many values'));
        return
    end
    
    if prod(n)==0
        % check if n has a zero value
        disp(sprintf(':: n :: is zero'));
        return
    end
    % compute the number of elements
    nElements = n;
    
    % compute the number of internal edges
    nInternalEdges = (p-1)*nElements;

    % compute the number of internal edges per element
    nElementInternalEdges = (p-1);

    % compute the total number of edges
    nEdges = (n*p+1);

    % compute the number of edges in each element
    nElementEdges = (p+1);

    gn = zeros(nElements, nElementEdges,'uint32');

    % first select the internal edges of the elements
    internalEdgesColumns = zeros(nElementInternalEdges,1,'uint32');
    
    for i=1:p-1
       internalEdgesColumns(i) = i+1;
    end

    gn(:,internalEdgesColumns) = reshape(1:nInternalEdges,[(p-1) nElements])';

    % Boundary edges of interior elements
    edgesgn = uint32(nInternalEdges+1):uint32(nEdges)-2;
        for k=1:(n-1)
            gn(k,end) = edgesgn(k);
            gn(k+1,1) = edgesgn(k);
        end
  

    % now the edges at the boundary edges
    edgesgn = nEdges - 1; % left boundary
    gn(1,1) = edgesgn; % left boundary
    
    if circular
        gn(end,end) = gn(1,1); % if boundary conditions are circular, then the
                               % values on the right are the same as the
                               % ones on the left
    else
        edgesgn = nEdges; % right boundary
        gn(end,end) = edgesgn; % right boundary
    end
    
end