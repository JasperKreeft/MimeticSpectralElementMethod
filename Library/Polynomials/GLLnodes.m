function [x,w]=GLLnodes(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lglnodes.m
%
% Computes the Gauss-Lobatto-Legendre nodes, weights and the GLL Vandermonde 
% matrix. The GLL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
% integration and spectral methods. 
%
% Reference on LGL nodes and weights: 
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
%
% Changed by Jasper Kreeft - 2009
% Contact: j.j.kreeft@tudelft.nl
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Truncation + 1
N1=N+1;

% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=cos(pi*(0:N)/N)';

% The Legendre Vandermonde Matrix
P=zeros(N1,N1);

% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and 
% update x using the Newton-Raphson method.

xold=2;

while max(abs(x-xold))>eps

    xold=x;

    P(:,1)=1;    P(:,2)=x;
    
    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end

    x=xold-( x.*P(:,N1)-P(:,N) )./( N1*P(:,N1) );

end

w=2./(N*N1*P(:,N1).^2);

x=flipud(x)'; w=flipud(w)';


%% Alternative method
% function [x,w]=GLLnodes(N)
% % zeros and weights for Gauss quadrature
% % 
% % Input N: nr of cells
% 
% [Le,dLe] = LegendrePoly(N);
% 
% innerzeros = roots(dLe(N+1,:))';
% 
% x = sort([-1 innerzeros 1]);
% 
% w = 2./(N*(N+1)*polyval(Le(N+1,:),GLLzeros).^2);
