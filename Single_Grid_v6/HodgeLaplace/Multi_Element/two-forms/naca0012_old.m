function y = naca0012(x)

% clear all
% clf
% clc

t = 0.12; % thickness
c = 1; % chord-length

% x=linspace(0,1,1000);
x = x+0.5;

xc = x/c;

y = t*c/0.2*( 0.2969*sqrt(xc)-0.1260*xc-0.3516*xc.^2+0.2843*xc.^3-0.1015*xc.^4 );

plot([x x(end:-1:1)]-0.5,[y -y(end:-1:1)])
axis equal
axis([-0.5 0.5 -0.1 0.1])