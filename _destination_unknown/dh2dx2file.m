clear all
% close all
clc

n = 1000;
x = linspace(-1,1,n);

N = 6;

xi = -cos((0:N)*pi/N);

[T,dTdx,d2Tdx2,d3Tdx3,d4Tdx4] = ChebyshevVal(x,N);

d2hdx2 = zeros(N+1,n);
for i=1:N+1
    d2hdx2(i,:) = (-1)^(N+i)/N^2/(1+(i==1)+(i==N+1))./(x-xi(i)).^3.* ( (x-xi(i)).^2./(1-x.^2).*(-dTdx(N+1,:)+N^2*(x.*T(N+1,:)-(1-x.^2).*dTdx(N+1,:))) + 2*(N^2*T(N+1,:)+x.*dTdx(N+1,:)).*(x-xi(i))+2*(1-x.^2).*dTdx(N+1,:));
end

figure
for i=1:N+1
    plot(x,d2hdx2(i,:))
    hold on
end
plot(x,d2hdx2(2,:),'r')
grid on


% d2hdx2_xisxi = -1/3./(1-xi.^2).*( 6 + 6*xi./(1-xi.^2) + N^2 - (11*xi.^2+4)./(1-xi.^2) );
d2hdx2_xisxi = -(N^2*(1-xi.^2)+xi.^2+2)./(3*(1-xi.^2).^2);

for i=2:N
    plot(xi(i),d2hdx2_xisxi(i),'xg')
end