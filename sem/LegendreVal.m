function [L,dLdx] = LegendreVal(x,N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [L,dLdx] = LegendreVal(x,N)
% This function can be used to calculate the values of the Legendre
% polynomial and its derivative.
% 
% Written by Jasper Kreeft - 2009
% Contact: j.j.kreeft@tudelft.nl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(x,1)>size(x,2)
    x=x';
end

nx = size(x,2);

L      = zeros(N+1,nx);
L(1,:) = ones(1,nx);
if N>0; L(2,:) = x; end

for k = 2:N
    L(k+1,:) = (2*k-1)/k.*x.*L(k,:)-(k-1)/k.*L(k-1,:);
end

if nargout==2

dLdx = zeros(N+1,nx);

dLdx(1,:) = zeros(1,nx);
if N>0; dLdx(2,:) = ones(1,nx); end

for k = 2:N
    dLdx(k+1,:) = dLdx(k-1,:)+(2*k-1)*L(k,:);
end

dLdx = dLdx(N+1,:);
end

L  =  L(N+1,:);