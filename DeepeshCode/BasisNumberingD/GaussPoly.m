function result = GaussPoly(x,p)
%GaussPoly Returns the p+1 Gauss Lagrange interpolant basis functions,
%evaluated at x.
%   result = GaussPoly(x,p) gives the p+1 Gauss basis interpolants 
%   evaluated at x.
%   It returns a (n+1)d matrix with the values of the Gauss basis
%   interpolants evaluated at x, where n is the dimension of x. If x is
%   a vector it returns a 2d matrix whose rows are the values of the
%   evaluated polynomial in x. If x is a 2d matrix then it returns a 3d
%   matrix whose planes (first and second index) are the evaluations and
%   the third index corresponds to the polynomial.
%
%   If x=[] then it computes the Gauss polynomial basis at the nodes,
%   that is, the result is a sparse identity matrix.
%
%   Copyright 2009 Artur Palha
%   $Revision: 1.0 $  $Date: 2009/11/25 13:38:00 $
    
    sizeOfx = size(x);
    if sizeOfx(1) == 1,
        sizeOfx = sizeOfx(2);
    elseif sizeOfx(2) == 1,
        sizeOfx = sizeOfx(1);
        x = x';
    end
    
    if length(x) == 0
        result = speye(p+1);
        return
    end
    
    % allocate memory space for the result
    result = zeros([p+1 sizeOfx]);
    
    % compute Gauss roots
    roots = GaussQuad(p);
    
    % auxiliary matrix to select roots for each polynomial
    rootSelection = true(p+1,p+1);
    rootSelection(linspace(1,(p+1)^2,p+1)) = false;
    
    % Compute each polynomial n using the formula:
    % \frac{\prod_{i=1\\i\neq n}^{p+1}(x-r_{i})}{\prod_{i=1\\i\neq
    % j}^{p+1}(r_{n}-r_{i})}
    %
    % For the top product one uses the built in function poly, based upon
    % the roots. For the bottom part, one just computed the product.
    
    if p == 0
        result(:) = 1.0;
    elseif p ~= 1,
        for n=1:p+1,
            repmatRoots = repmat(roots(rootSelection(n,:)),[1 sizeOfx]);
            result(n,:) = prod((repmat(x,[p,1])-repmatRoots)./ (roots(n)-repmatRoots));
        end
    else
        for n=1:p+1,
            repmatRoots = repmat(roots(rootSelection(n,:)),[1 sizeOfx]);
            result(n,:) = (repmat(x,[p,1])-repmatRoots)./ (roots(n)-repmatRoots);
        end
    end
end