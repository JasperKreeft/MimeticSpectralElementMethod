function x = EGaussNodes(p)
%EGaussNodes Return the p+1 Extended Gauss nodes for extended Gauss
%interpolation.
%   x = EGaussNodes(p) gives the p+1 nodes for extended Gauss
%   interpolation. This p+1 nodes consist of the p-1 nodes for Gauss
%   quadrature of order p-2 plus the endpoints -1 and 1.
%   The function GaussQuad(p-2) is used to compute the inner nodes.
%
%   Copyright 2009 Artur Palha
%   $Revision: 1.0 $  $Date: 2009/11/26  $

    x = zeros(p+1,1);
    
    if p > 1,
        x(2:end-1) = GaussQuad(p-2);
    end
    
    x(1) = -1.0;
    x(end) = 1.0;
end