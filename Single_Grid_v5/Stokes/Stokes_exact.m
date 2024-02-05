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

% x = linspace(-1,1,100);
% y = linspace(-1,1,100);
% 
% [X,Y] = meshgrid(y,x);
% 
% u_ex  = -sin(pi*X).*cos(pi*Y);
% v_ex  = cos(pi*X).*sin(pi*Y);
% w_ex = -2*pi*sin(pi*X).*sin(pi*Y);
% p_ex  = cos(pi/2*X).*cos(pi/2*Y);
% fx_ex = -pi/2*sin(pi/2*X).*cos(pi/2*Y)-2*pi^2*sin(pi*X).*cos(pi*Y);
% fy_ex = -pi/2*cos(pi/2*X).*sin(pi/2*Y)+2*pi^2*cos(pi*X).*sin(pi*Y);
% F_ex  = sqrt(fx_ex.^2+fy_ex.^2);
% ax = [-1 1 -1 1];

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

%% case AFG 2012

x = linspace(0,1,100);
y = linspace(0,1,100);

[X,Y] = meshgrid(y,x);

u_ex  = -2*X.^2.*(X-1).^2.*Y.*(2*Y-1).*(Y-1);
v_ex  =  2*Y.^2.*(Y-1).^2.*X.*(2*X-1).*(X-1);
w_ex  =  2*Y.^2.*(Y-1).^2.*(2*X-1).*(X-1)+4*Y.^2.*(Y-1).^2.*X.*(X-1)+2*Y.^2.*(Y-1).^2.*X.*(2*X-1)+2*X.^2.*(X-1).^2.*(2*Y-1).*(Y-1)+4*X.^2.*(X-1).^2.*Y.*(Y-1)+2*X.^2.*(X-1).^2.*Y.*(2*Y-1);
p_ex  = (X-1/2).^5+(Y-1/2).^5;
fx_ex = 4*Y.*(Y-1).^2.*(2*X-1).*(X-1)+4*Y.^2.*(Y-1).*(2*X-1).*(X-1)+8*Y.*(Y-1).^2.*X.*(X-1)+...
        8*Y.^2.*(Y-1).*X.*(X-1)+4.*Y.*(Y-1).^2.*X.*(2*X-1)+4*Y.^2.*(Y-1).*X.*(2*X-1)+8*X.^2.*(X-1).^2.*(Y-1)+...
        4*X.^2.*(X-1).^2.*(2*Y-1)+8*X.^2.*(X-1).^2.*Y+5*(X-1/2).^4;
fy_ex = -8*Y.^2.*(Y-1).^2.*(X-1)-4*Y.^2.*(Y-1).^2.*(2*X-1)-8*Y.^2.*(Y-1).^2.*X-4*X.*(X-1).^2.*(2*Y-1).*(Y-1)-...
        4*X.^2.*(X-1).*(2*Y-1).*(Y-1)-8*X.*(X-1).^2.*Y.*(Y-1)-8*X.^2.*(X-1).*Y.*(Y-1)-...
        4*X.*(X-1).^2.*Y.*(2*Y-1)-4*X.^2.*(X-1).*Y.*(2*Y-1)+5*(Y-1/2).^4;
F_ex  = sqrt(fx_ex.^2+fy_ex.^2);
ax = [0 1 0 1];


%% Case 4  Kreeft, Gerritsma SIAM 2012

% x = linspace(0,1,100);
% y = linspace(0,1,100);
% 
% [X,Y] = meshgrid(y,x);
% 
% u_ex  = 2*sin(3/2*pi*X).*cos(3/2*pi*Y);
% v_ex  = cos(3/2*pi*X).*sin(3/2*pi*Y);
% w_ex  = 3/2*pi*sin(3/2*pi*X).*sin(3/2*pi*Y);
% p_ex  = sin(pi*X).*sin(pi*Y);
% fx_ex = pi*cos(pi*X).*sin(pi*Y)+9/4*pi^2*sin(3/2*pi*X).*cos(3/2*pi*Y);
% fy_ex = pi*sin(pi*X).*cos(pi*Y)-9/4*pi^2*cos(3/2*pi*X).*sin(3/2*pi*Y);
% F_ex  = sqrt(fx_ex.^2+fy_ex.^2);
% g_ex  = 9/2*pi*cos(3/2*pi*X).*cos(3/2*pi*Y);
% ax = [0 1 0 1];


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

if exist('g_ex','var')
    figure
    contourf(X,Y,g_ex,20)
    colorbar
    axis equal
    axis(ax)
    title('g')
end