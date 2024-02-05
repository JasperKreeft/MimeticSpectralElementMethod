function result = DerivativePoly(x, p, polyType)
%DerivativePoly Returns the derivative of the polynomial interpolant basis
%functions at the points x.
%
%   result = DerivativePoly(x, p, polyType) computes the derivative of the
%   polynomial interpolating basis functions of type polyType, at points x.
%   
%   result(i,j,...,k) :: is the value of dh_{k}(x(j,...,k))
%   
%   The derivative of the Lagrange interpolating basis functions
%   (l_{n}^{p}(x)) are given by:
%   
%   \frac{dl_{n}(x)}{dx} = \sum_{i=1}_{i\neq n}^{p+1}\prod_{j\neq n}_{j\neq i}
%                          \frac{1}{x_{n}-x_{i}}\frac{x-x_{j}}{x_{n}-x_{j}}
%
%   This computation is very slow computed like this. Hence, a clever way
%   is followed.
%   
%       1- compute the derivatives of the basis functions at the nodal
%          points (using DERIVATIVEPOLYNODES)
%       2- compute the basis functions at the x points)
%       3- update the values of the derivatives at the x points
%
%   In this way, instead of having an algorithm of O(p^{2}) one has an
%   algorithm of O(p).
%
%   See also: DERIVATIVEPOLYNODES, LOBATTOPOLY, GAUSSPOLY, EGAUSSPOLY

%   Copyright 2009 Artur Palha
%   $Revision: 1.0 $  $Date: 2009/11/25 13:38:00 $
    
    % check if polyType is a valid one
    if ~TestPolyType(polyType)
        disp(sprintf(':: %s :: is not a valid type of polynomial', polyType));
        return
    end
    
    % shape of x
    sizeOfx = size(x);
    lengthOfx = numel(x);
    
    % transform x in a vector
    x = x(:);
    
    % preallocate memory space for result
    result = zeros(p+1, lengthOfx);
    
    % compute the derivatives at the nodal points
    nodalDerivatives = DerivativePolyNodes(p, polyType, [true false]);
    
    % update the derivatives to the x points
    
    % compute the basis polynomials in the x points
    polyEval = (eval(sprintf('%sPoly(%s, %s)', polyType, 'x', 'p')))';
    
    % compute the derivatives in the x points
    result = (polyEval * nodalDerivatives)';
    
    if min(sizeOfx) ~= 1
        result = reshape(result, [p+1 sizeOfx]);
    end
end