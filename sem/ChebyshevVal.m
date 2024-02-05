function [T,dTdx,d2Tdx2,d3Tdx3,d4Tdx4] = ChebyshevVal(x,N)

% This function can be used to calculate the values of the Chebychev
% polynomial.

if size(x,1)>size(x,2)
    x=x';
end

nx = size(x,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T = zeros(N+1,nx)  ;
T(1,:) = ones(1,nx);
T(2,:) = x         ;

for k=2:N
    T(k+1,:) = 2*x.*T(k,:)-T(k-1,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout>=2

dTdx = zeros(N+1,nx);
dTdx(2,:) = ones(1,nx);

for k=2:N
    dTdx(k+1,:) = 2*T(k,:)+2*x.*dTdx(k,:)-dTdx(k-1,:);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout>=3

d2Tdx2 = zeros(N+1,nx);

for k=2:N
    d2Tdx2(k+1,:) = 4*dTdx(k,:)+2*x.*d2Tdx2(k,:)-d2Tdx2(k-1,:);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout>=4
    
d3Tdx3 = zeros(N+1,nx);

for k=2:N
    d3Tdx3(k+1,:) = 6*d2Tdx2(k,:)+2*x.*d3Tdx3(k,:)-d3Tdx3(k-1,:);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout==5
    
d4Tdx4 = zeros(N+1,nx);
for k=2:N
    d4Tdx4(k+1,:) = 8*d3Tdx3(k,:)+2*x.*d4Tdx4(k,:)-d4Tdx4(k-1,:);
end

end
% d4Tdx4(N+1,:) = 1./((1-x.^2).^3).*(N^4*(1-x.^2).*T(end,:)-N^2*(11*x.^2+4).*T(end,:)+(-6*N^2+9)*x.*(1-x.^2).*dTdx(end,:)+15*x.^3.*dTdx(end,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%