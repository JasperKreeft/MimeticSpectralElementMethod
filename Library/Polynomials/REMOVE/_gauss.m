function [Gzeros,Gweights] = gauss(N)
% zeros and weights for Gauss quadrature
% 
% Input N: nr of cells

[Le,dLe] = LegendrePoly(N);

Gzeros = sort(roots(Le(N+1,:)))';

Gweights = 2./((1-Gzeros.^2).*polyval(dLe(N+1,:),Gzeros).^2);