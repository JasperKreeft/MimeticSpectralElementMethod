function result = DerivativePolyNodes(p, polyType, options)
%DerivativePolyNodes Returns the derivative of the polynomial interpolant basis
%functions at the nodal points.
%
%   result = DerivativePoly(x, p, polyType) computes the derivative of the
%   polynomial interpolating basis functions of type polyType, at the nodal
%   points.
%
%   result(i,j) :: is the value of dh_{i}(x_{k})
%   
%   The derivative of the Lagrange interpolating basis functions
%   (l_{n}^{p}(x)) are given by:
%   
%   \frac{dl_{n}(x)}{dx} = \sum_{i=1}_{i\neq n}^{p+1}\prod_{j\neq n}_{j\neq i}
%                          \frac{1}{x_{n}-x_{i}}\frac{x-x_{j}}{x_{n}-x_{j}}
%
%
%   For computation at the nodes a more efficient and accurate formula can
%   be used, see [1]:
%       
%             | \frac{c_{k}}{c_{j}}\frac{1}{x_{k}-x_{j}},          k \neq j
%             |
%   d_{kj} = <
%             | \sum_{l=1,l\neq k}^{p+1}\frac{1}{x_{k}-x_{l}},     k = j
%             |
%   
%   with
%
%       c_{k} = \prod_{l=1,l\neq k}^{p+1} (x_{k}-x_{l})   
%
%   [1] Costa, B., Don, W. S.: On the computation of high order
%       pseudospectral derivatives, Applied Numerical Mathematics, vol.33
%       (1-4), pp. 151-159

%   Copyright 2009 Artur Palha
%   $Revision: 1.0 $  $Date: 2009/11/25 13:38:00 $
    
    if nargin == 2
        options = [true true];
    end
    
    % check if polyType is a valid one
    if ~TestPolyType(polyType)
        disp(sprintf(':: %s :: is not a valid type of polynomial', polyType));
        return
    end
    
    % preallocate memory space for result
    result = zeros(p+1, p+1);
    
    % compute the nodes of the type of polynomial
    roots = eval(sprintf('%sQuad(%d)', polyType, p));
    
    % compute all (ri - rj) = xi_xj(i,j)
    xi_xj = repmat(roots, [1 p+1]) - repmat(roots', [p+1 1]);
    
    % transform all the diagonals to 1
    xi_xj(linspace(1,(p+1)^2,(p+1))) = 1.0;
    
    % compute (ci's)
    c = prod(xi_xj, 2);
    
    % compute ck/cj = ck_cj(k,j) matrix
    ci_cj = repmat(c, [1 p+1]) ./ repmat(c', [p+1 1]);
    
    % compute the off diagonal results
    result = ci_cj ./ xi_xj;
    
    % compute the diagonal components
    
    if options(1) == false
        % using formula (6) of [1]
        xi_xj = 1.0 ./ xi_xj;
    
        % transform the diagonal elements to 0.0
        xi_xj(linspace(1,(p+1)^2,(p+1))) = 0.0;
    
        % update the diagonal result with the values of formula (6)
        result(linspace(1,(p+1)^2,(p+1))) = sum(xi_xj, 2);
    
    else
        % using formula (9) of [1]
    
        % put all diagonal values equal to 0.0
        result(linspace(1,(p+1)^2,(p+1))) = 0.0;
    
        % compute the diagonal values
        result(linspace(1,(p+1)^2,(p+1))) = -sum(result,2);
    end
    
    if options(2) == true
        % transpose the result
        result = result';
    end
end