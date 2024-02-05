function [x weights] = LobattoQuad(p)
%LobattoQuad Return the p+1 Lobatto weights of Gauss-Lobatto quadrature.
%   [x weights] = LobattoWeights(p) gives the p+1 points and weights of the
%   Gauss-Lobatto quadrature of order p. 
%   For the computation of the nodes it uses a Newton method
%   up to machine precision. As initial guess it uses the Chebychev roots.
%   for the roots.
%
%   Copyright 2009 Artur Palha
%   $Revision: 1.0 $  $Date: 2009/11/25 13:38:00 $

    n = p+1;
    roots = cos(pi*(0:n-1)/p)';
    P = zeros(n,n);
    xold = 2.0;
    while max(abs(roots-xold)) > eps,
        xold = roots;
        P(:,1) = 1.0;
        P(:,2) = roots;
        for k=3:n,
            P(:,k) = ((2*(k-1)-1).*roots.*P(:,k-1) - (k-2)*P(:,k-2))/(k-1);
        end
        roots = xold - (roots.*P(:,n) - P(:,n-1))./(n*P(:,n));
    end
    roots = fliplr(roots')';
    weights = 2.0./(p*n*(P(:,n)).^2);
end