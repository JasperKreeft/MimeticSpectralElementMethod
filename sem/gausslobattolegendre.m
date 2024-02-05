function [GLLzeros,GLLweights] = gausslobattolegendre(N)
% zeros and weights for Gauss quadrature
% 
% Input N: nr of cells

[Le,dLe] = LegendrePoly(N);

innerzeros = roots(dLe(N+1,:))';

GLLzeros = sort([-1 innerzeros 1]);

GLLweights = 2./(N*(N+1)*polyval(Le(N+1,:),GLLzeros).^2);