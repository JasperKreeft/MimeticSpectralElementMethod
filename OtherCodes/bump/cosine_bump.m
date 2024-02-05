clear all
close all
clc

a = 0.2;

nn = 30;

xi  = linspace(-1,1,nn)'*ones(1,nn);
eta = xi';

x = xi;
y = ((1-a)-a*cos(pi*xi)).*(1-eta)/2;

surf(x,y,zeros(nn))
view([0 0 1])