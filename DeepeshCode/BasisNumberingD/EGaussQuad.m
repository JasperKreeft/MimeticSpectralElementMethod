function [x w] = EGaussQuad(p)
%EGaussQuad Return the p+1 Extended Gauss nodes for extended Gauss
%interpolation and the p+1 extend weights.
%   [x w] = EGaussQuad(p) gives the p+1 nodes for extended Gauss
%   interpolation and weights. This p+1 nodes consist of the p-1 nodes for
%   Gauss quadrature of order p-2 plus the endpoints -1 and 1. The weights
%   are the same as the weights for Gauss quadrature of order p-2 but the
%   exterior weights are zero.
%
%   The function GaussQuad(p-2) is used to compute the inner nodes and
%   inner weights.
%
%   Copyright 2009 Artur Palha
%   $Revision: 1.0 $  $Date: 2009/11/26  $

    x = zeros(p+1,1);
    w = zeros(p+1,1);
    
    if p > 1,
        [x(2:end-1) w(2:end-1)] = GaussQuad(p-2);
    end
    
    x(1) = -1.0;
    x(end) = 1.0;
    
    w(1) = 0.0;
    w(end) = 0.0;
end