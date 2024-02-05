function [A] = lumpedmassmatrix(Ne,Ng,Pe,w,c,h)


A = zeros(Ng);
for l=1:Ne
    p = Pe(l);
    Theta = diag(w(p,1:p+1));
    s = c(l,1):c(l,p+1);
    A(s,s) = A(s,s)+h(l)/2*Theta;
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
% close all
% clc
% 
% % Mass matrix
% 
% p = 4; % order of polynomial expansion
% 
% k = p; % nr of integration points / order of the lobatto quadrature
% 
% 
% % Gauss-Lobatto-Legendre zeros
% P11 = JacobiPoly(p-1,1,1);
% innerzeros = roots(P11(p,end:-1:1));
% xi = sort([-1; innerzeros; 1])';
% 
% % Lobatto integration points
% z = xi;
% 
% % corresponding weights
% Le = JacobiPoly(k,0,0);
% w = zeros(1,k+1);
% for i=1:k+1
%     w(i) = 2/(k*(k+1)*(sum(Le(k+1,:).*z(i).^(0:k)))^2);
% end
% 
% Theta = diag(w);