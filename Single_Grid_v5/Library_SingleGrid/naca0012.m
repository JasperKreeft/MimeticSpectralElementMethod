function [y,dydx] = naca0012(x)
% 
% clear all
% clf
% clc
% x=linspace(0,1,10000);

t = .5; %0.12 % thickness
c = 1; % chord-length

x = x+0.5;

xc = x/c;

y = t*c/0.2*( 0.2969*sqrt(xc)-0.1260*xc-0.3516*xc.^2+0.2843*xc.^3-0.1015*xc.^4 ).*(xc~=1);

% plot([x x(end:-1:1)]-0.5,[y -y(end:-1:1)])
% axis equal
% axis([-0.5 0.5 -0.1 0.1])

dydx = t/0.2*( -0.2969./(2*sqrt(xc+eps))-0.1260-2*0.3516*xc+3*0.2843*xc.^2-4*0.1015*xc.^3 );
% dydx(1) = max(-10,2*dydx(2)-dydx(3));
dydx(xc==0) = -10;


% hold on
% plot([x x(end:-1:1)]-0.5,[dydx -dydx(end:-1:1)],'r')
% axis equal
% axis([-0.5 0.5 -0.1 0.1])