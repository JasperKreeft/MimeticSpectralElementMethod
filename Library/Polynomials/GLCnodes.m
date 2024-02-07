function [x,w]=GLCnodes(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [x,w]=GLCnodes(N)
% 
% Zeros and weights for Gauss-Lobatto-Chebyshev quadrature
% 
% Written by Jasper Kreeft - 2010
% Contact: j.j.kreeft@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = -cos((0:N)*pi/N);

w = [pi/(2*N) pi/N*ones(1,N-1) pi/(2*N)];