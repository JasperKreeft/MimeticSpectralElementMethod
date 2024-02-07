function [x,w] = GCnodes(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,w] = GCnodes(N)
%
% Zeros and weights for Chebyshev Gauss quadrature
% 
% Written by Jasper Kreeft - 2010
% Contact: j.j.kreeft@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=N-1;

x = -cos( (2*(0:N)+1)/(2*N+2)*pi );

w = pi/(N+1)*ones(1,N+1);