function [F]=force_v2(m,x,y)

Nx = length(x)-1;
Ny = length(y)-1;

F = zeros(Nx,Ny);
for j=1:Ny
    for i=1:Nx
        F(i,j) = -2*(cos(m*pi*x(i+1))-cos(m*pi*x(i)))*(cos(m*pi*y(j+1))-cos(m*pi*y(j)));
%         F(i,j) = 2/3*(x(i+1)^3-x(i)^3)*(y(j+1)-y(j))+2/3*(x(i+1)-x(i))*(y(j+1)^3-y(j)^3)-4*(x(i+1)-x(i))*(y(j+1)-y(j));
    end
end