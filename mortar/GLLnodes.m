function [x w] = GLLnodes(N)
%
%   GLLnodes computes the gauss-lobatto-legendre nodes, weights and the GLL 
%   vandermonde matrix (the GLL nodes are the zeros of (1-x^2)*P'_N(x))
%
%   [x w] = GLLnodes(N)
%
%   input:
%       N       :: is the number of elements 
%
%   output:
%       x       :: GLL nodes
%       w       :: GLL weights
%
%   Copyright 2004 Greg von Winckel (gregvw@chtm.unm.edu)
%   $Revision: 1.0 $  $Date: 04/17/2004 $
%   Revised by: Jasper Kreeft (j.j.kreeft@tudelft.nl)
%   $Revision: 2.0 $  $Date: xx/xx/2009 $
%   Revised by: Peter Kuystermans
%   $Revision: 3.0 $  $Date: 01/10/2011 $  

%-------------------------------------------------------------------------%
% parameters                                                              %
%-------------------------------------------------------------------------%

    N1 = N+1;	% truncation + 1

%-------------------------------------------------------------------------%
% initial guess                                                           %
%-------------------------------------------------------------------------%
    
    x = cos(pi*(0:N)/N)';	% first guess chebyshev-gauss-lobatto nodes

%-------------------------------------------------------------------------%
% storage                                                                 %
%-------------------------------------------------------------------------%
    
    P = zeros(N1,N1);	% legendre vandermonde matrix

%-------------------------------------------------------------------------%
% calculation of weights and nodes                                        %
%-------------------------------------------------------------------------%
    
    %%% compute P_(N) using the recursion relation, compute its first and 
    %%% second derivatives and update x using the Newton-Raphson method
    xold = 2;
    while max(abs(x-xold))>eps

        xold = x;
        P(:,1) = 1;    
        P(:,2) = x;

        for k=2:N
            P(:,k+1) = ((2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1))/k;
        end

        x = xold-(x.*P(:,N1)-P(:,N))./(N1*P(:,N1));

    end

    w = 2./(N*N1*P(:,N1).^2);
    x = flipud(x); 
    w = flipud(w);

end