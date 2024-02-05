function [F] = force_zeroform(x,y)

global m

F = -2*m^2*pi^2.*sin(m*pi*x).*sin(m*pi*y);