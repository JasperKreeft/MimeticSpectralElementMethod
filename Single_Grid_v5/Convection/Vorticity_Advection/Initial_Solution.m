clear all
close all
clc

nn = 1000;
x = linspace(-pi,pi,nn);
X = x'*ones(1,nn);
Y = X';

U = 1;
a = 0.3;


r1 = -0.4;
R = sqrt((X+r1).^2+Y.^2);


W = U/a*(2-R.^2/a^2).*exp(1/2*(1-R.^2/a^2));

r2 = +0.4;

R = sqrt((X+r2).^2+Y.^2);

W = W + U/a*(2-R.^2/a^2).*exp(1/2*(1-R.^2/a^2));

pcolor(X,Y,W)
shading interp
axis equal
axis([-pi pi -pi pi])
axis off