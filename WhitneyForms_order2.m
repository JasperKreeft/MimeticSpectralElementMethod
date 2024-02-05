clear all
close all
clc

WindowSize

xe(1) = 0;
ye(1) = 0;

xe(2) = 1;
ye(2) = 0;

xe(3) = 0.5;
ye(3) = sqrt(3)/2;

xe(4) = 0.5;
ye(4) = 0;

xe(5) = 0.75;
ye(5) = sqrt(3)/4;

xe(6) = 0.25;
ye(6) = sqrt(3)/4;

M = [ ones(1,6) ; xe ; ye ; xe.^2 ; ye.^2 ; xe.*ye ];
% Mi = inv(M);

n = 100;

xi = linspace(0,1,n);
eta = linspace(0,1,n);

[Xi,Eta] = meshgrid(xi,eta);

X = 0.5+(2*Xi-1).*(1-Eta)/2;
Y = sqrt(3)/2*Eta;

for i=1:n
    for j=1:n
        Vec = [ 1 ; X(i,j) ; Y(i,j) ; X(i,j)^2 ; Y(i,j)^2 ; X(i,j)*Y(i,j) ];
        Z = M\Vec;
        Z1(i,j) = Z(1);
        Z2(i,j) = Z(2);
        Z3(i,j) = Z(3);
        Z4(i,j) = Z(4);
        Z5(i,j) = Z(5);
        Z6(i,j) = Z(6);
    end
end

figure(win{2,3}{:})
subplot(2,3,1)
% pcolor(X,Y,Z1)
contour(X,Y,Z1,'ShowText','on')
shading interp
axis equal
colormap('jet')
colorbar
subplot(2,3,2)
pcolor(X,Y,Z2)
shading interp
axis equal
colormap('jet')
colorbar
subplot(2,3,3)
pcolor(X,Y,Z3)
shading interp
axis equal
colormap('jet')
colorbar
subplot(2,3,4)
pcolor(X,Y,Z4)
shading interp
axis equal
colormap('jet')
colorbar
subplot(2,3,5)
pcolor(X,Y,Z5)
shading interp
axis equal
colormap('jet')
colorbar
subplot(2,3,6)
pcolor(X,Y,Z6)
shading interp
axis equal
colormap('jet')
colorbar


figure(win{1,1}{:})
pcolor(X,Y,Z1+Z2+Z3+Z4+Z5+Z6)
shading interp
axis equal
colormap('jet')
colorbar