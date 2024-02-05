function [x w] = LobattoQuad(p)
%LobattoQuad Return the p+1 Lobatto points and weights of Gauss-Lobatto
%quadrature.
%   [x w] = LobattoQuad(p) gives the p+1 points and weights of the
%   Gauss-Lobatto quadrature of order p. 
%   For the computation of the nodes it uses a Newton method
%   up to machine precision. As initial guess it uses the Chebychev roots.
%   for the roots.
%
%   Copyright 2009 Artur Palha
%   $Revision: 1.0 $  $Date: 2009/11/25 13:38:00 $

    n = p+1;
    x = cos(pi*(0:n-1)/p)';
    P = zeros(n,n);
    xold = 2.0;
    while max(abs(x-xold)) > eps,
        xold = x;
        P(:,1) = 1.0;
        P(:,2) = x;
        for k=3:n,
            P(:,k) = ((2*(k-1)-1).*x.*P(:,k-1) - (k-2)*P(:,k-2))/(k-1);
        end
        x = xold - (x.*P(:,n) - P(:,n-1))./(n*P(:,n));
    end
    x = fliplr(x')';
    w = 2.0./(p*n*(P(:,n)).^2);
end