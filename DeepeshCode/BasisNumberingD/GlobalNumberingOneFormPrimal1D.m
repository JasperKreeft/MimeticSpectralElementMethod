function gn = GlobalNumberingOneFormPrimal1D(n,p)
%GlobalNumberingOneFormPrimal computes the global numbering of 2-form
%degrees of freedom on a  primal mesh.
%
%   gn = GlobalNumberingTwoFormPrimal(n, p)
%
%   Where:
%       n   :: the number of elements. n is the number of elements .
%       p   :: the order of the 0-form approximation
%
%   Returns the global numbering of the nodes of the mesh. gn is a matrix
%   N x p, where N is the total number of elements.
%
%   Internal numbering:
%
%       Internally nodes are numbered from left to right and from bottom to
%       top. So the first node of element 1, that is gn(1,1) is the bottom
%       left node of the first element. 
%
%   See also: GLOBALNUMBERINGZEROFORMDUAL, PLOTGLOBALNUMBERINZEROFORMPRIMAL

%   Copyright 2010 Artur Palha & Pedro Pinto Rebelo
%   $Revision: 1.2 $  $Date: 2010/04/13 $
    
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

    % compute the total number of degrees of freedom
    dof = p*nElements;
    
    % compute the global numbering
    gn = reshape(uint32(1):uint32(dof), p, nElements)';

end