function gn = GlobalNumberingTwoFormPrimal(n,p)
%GlobalNumberingTwoFormPrimal computes the global numbering of 2-form
%degrees of freedom on a  primal mesh.
%
%   gn = GlobalNumberingTwoFormPrimal(n, p)
%
%   Where:
%       n   :: the number of elements. If n is a number it has n x n
%              elements, if it is a vector, then it has n(1) x n(2) elements.
%       p   :: the order of the 0-form approximation
%
%   Returns the global numbering of the nodes of the mesh. gn is a matrix
%   N x p^2, where N is the total number of elements.
%
%   Internal numbering:
%
%       Internally nodes are numbered from left to right and from bottom to
%       top. So the first node of element 1, that is gn(1,1) is the bottom
%       left node of the first element. 
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

    % compute the total number of degrees of freedom
    dof = p*p*nElements;
    
    % compute the global numbering
    gn = reshape(uint32(1):uint32(dof), p*p, nElements)';

end