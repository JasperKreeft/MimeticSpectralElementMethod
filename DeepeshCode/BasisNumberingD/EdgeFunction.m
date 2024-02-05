function result = EdgeFunction(x, p, polyType)
%EdgeFunction Returns the edge basis functions at the x points
%associated to polyType.
%
%   result = EdgeFunction(x, p, polyType) computes the edge basis
%   functions [1] at the x points of the associated polyType.
%   
%   result(i,j) :: is the value of edge_{k}(x(j))
%   
%   The edge basis functions are given by, see [1]:
%
%       edge_{i}(x) = -\sum_{j=1}^{i} dh_{j}(x), i=1,...,p
%
%   dh_{j}(x) are computed using DerivativePolyNodes.
%
%   [1] Gerritsma, M.: Edge functions for spectral element methods. 
%       Submitted to the proceedings of ICOSAHOM 2009
%
%   See also: EdgeFunctionNodes, DERIVATIVEPOLYNODES, LOBATTOPOLY,
%             GAUSSPOLY, EGAUSSPOLY

%   Copyright 2009 Artur Palha
%   $Revision: 1.0 $  $Date: 2009/12/04 13:38:00 $
    
    % check if polyType is a valid one
    if ~TestPolyType(polyType)
        disp(sprintf(':: %s :: is not a valid type of polynomial', polyType));
        return
    end
    
    % compute the derivatives at the nodal points
    derivatives = DerivativePoly(x, p, polyType);
    
    % compute the edge basis functions
    result = cumsum(derivatives);
    
    result = -result(1:(end-1),:);
end