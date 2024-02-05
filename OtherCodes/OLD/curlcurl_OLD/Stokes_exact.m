% clear all
% close all
% clc


%% Case 1   Proot thesis (pag 72), Gerritsma & Phillips [65]

% x = linspace(0,1,100);
% y = linspace(0,1,100);
% 
% [X,Y] = meshgrid(y,x);
% 
% u_ex  = -sin(2*pi*X).*cos(2*pi*Y);
% v_ex  = cos(2*pi*X).*sin(2*pi*Y);
% p_ex  = sin(pi*X).*sin(pi*Y);
% fx_ex = pi*cos(pi*X).*sin(pi*Y)-8*pi^2*sin(2*pi*X).*cos(2*pi*Y);
% fy_ex = pi*sin(pi*X).*cos(pi*Y)+8*pi^2*cos(2*pi*X).*sin(2*pi*Y);
% F_ex  = sqrt(fx_ex.^2+fy_ex.^2);
% ax = [0 1 0 1];

%% case 2 Adapted Proot

x = linspace(-1,1,100);
y = linspace(-1,1,100);

[X,Y] = meshgrid(y,x);

u_ex  = -sin(pi*X).*cos(pi*Y);
v_ex  = cos(pi*X).*sin(pi*Y);
w_ex = -2*pi*sin(pi*X).*sin(pi*Y);
p_ex  = cos(pi/2*X).*cos(pi/2*Y);
fx_ex = -pi/2*sin(pi/2*X).*cos(pi/2*Y)-2*pi^2*sin(pi*X).*cos(pi*Y);
fy_ex = -pi/2*cos(pi/2*X).*sin(pi/2*Y)+2*pi^2*cos(pi*X).*sin(pi*Y);
F_ex  = sqrt(fx_ex.^2+fy_ex.^2);
ax = [-1 1 -1 1];

%% case 3 Bodard 2008

% x = linspace(-1,1,100);
% y = linspace(-1,1,100);
% 
% [X,Y] = meshgrid(y,x);
% 
% u_ex  =    -cos(pi/2*X).*sin(pi/2*Y);
% v_ex  =     sin(pi/2*X).*cos(pi/2*Y);
% w_ex  =  pi*cos(pi/2*X).*cos(pi/2*Y);
% p_ex  = -pi*sin(pi/2*X).*sin(pi/2*Y);
% fx_ex = -pi^2*cos(pi/2*X).*sin(pi/2*Y);
% fy_ex = 0*X;
% F_ex  = sqrt(fx_ex.^2+fy_ex.^2);
% ax = [-1 1 -1 1];

%% Plotten

figure
subplot(2,3,1)
contourf(X,Y,u_ex,20)
colorbar
axis equal
axis(ax)
title('u')
subplot(2,3,2)
contourf(X,Y,v_ex,20)
colorbar
axis equal
axis(ax)
title('v')
subplot(2,3,3)
contourf(X,Y,w_ex,20)
colorbar
axis equal
axis(ax)
title('\omega')
subplot(2,3,4)
contourf(X,Y,p_ex,20)
colorbar
axis equal
axis(ax)
title('p')
subplot(2,3,5)
contourf(X,Y,fx_ex,20)
colorbar
axis equal
axis(ax)
title('f_x')
subplot(2,3,6)
contourf(X,Y,fy_ex,20)
colorbar
axis equal
axis(ax)
title('f_y')
