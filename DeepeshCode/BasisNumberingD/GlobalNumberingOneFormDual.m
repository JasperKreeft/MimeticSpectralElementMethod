function gn = GlobalNumberingOneFormDual(n,p)
%GlobalNumberingOneFormDual computes the global numbering of 1-form
%degrees of freedom on a  dual mesh.
%
%   gn = GlobalNumberingOneFormDual(n, p)
%
%   Where:
%       n   :: the number of elements. If n is a number it has n x n
%              elements, if it is a vector, then it has n(1) x n(2) elements.
%       p   :: the order of the 1-form approximation
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
%   See also: GLOBALNUMBERINGZEROFORMPRIMAL, GLOBALNUMBERINGONEFORMPRIMAL, PLOTGLOBALNUMBERINGONEFORMDUAL

%   Copyright 2010 Artur Palha
%   $Revision: 1.2 $  $Date: 2010/03/08 $
    

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
    
    % compute the total number of edges
    nEdges = n(1)*n(2)*p*(p+1);

    % compute the number of edges in each element
    nElementEdges = p*(p+1);
    
    % allocate memory space for global numbering matrix
    % Notice that edges in elements are numbered from bottom to top and from
    % left to right.
    % Elements follow the same kind of numbering: bottom to
    % top and left to right.
    gn = zeros(nElements, 2*nElementEdges, 'uint32');
    
%--------------------------------------------------------------------------
% compute the numbering for the eta basis functions part
%--------------------------------------------------------------------------
    
    % generate the numbering for the eta basis functions part
    edgesgn = 1:nEdges;
    gn(:,1:nElementEdges) = reshape(edgesgn, [nElementEdges nElements])';
    
    
%--------------------------------------------------------------------------
% compute the numbering for the xi basis functions part
%--------------------------------------------------------------------------

    % generate the numbering for the xi basis functions part
    edgesgn = (edgesgn(end)+1):(2*nEdges);
    gn(:,(nElementEdges+1):end) = reshape(edgesgn, [nElementEdges nElements])';
    

end