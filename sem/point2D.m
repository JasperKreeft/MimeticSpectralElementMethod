clear all
close all
clc

N=4;
kk = 300;
xx=linspace(-1,1,kk);
yy=linspace(-1,1,kk);

[XX,YY] = meshgrid(xx,yy);

xi = gllnodes(N);

ll=MimeticpolyVal(xx,N,1);

ZZ = ll(4,:)'*ll(2,:);

pcolor(XX,YY,ZZ)
shading interp
hold on
plot3(xi(2),xi(4),5,'ok')
axis equal
axis([-1 1 -1 1])
title('asildj')
xlabel('x-axis')
ylabel('y-axis')