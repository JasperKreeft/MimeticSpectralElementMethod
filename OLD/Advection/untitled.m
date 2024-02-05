close all
clear all
clc
p=1;
e=.1;
x=linspace(-1,1,200);
y=x;
[X,Y]=meshgrid(x,y);
b=2+p*sin(2*pi*(X-Y)/e);
surf(X,Y,1./b)
view([0 0 1])
shading interp
f = -.5*((6*X.^2-1).*(Y.^4-Y.^2)+(6*Y.^2-1).*(X.^4-X.^2));
figure
surf(X,Y,f)
view([0 0 1])
shading interp
