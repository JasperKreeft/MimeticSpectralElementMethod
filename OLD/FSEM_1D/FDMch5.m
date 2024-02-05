clear
clc

N = 5;
h = 1/N;
x = 0+(0:N)/N;
dx = 1/N;

bc_R = 0;
bc_L = 0;

A = zeros(N+1);
for i=2:N
    A(i,i-1) = -1/dx^2;
    A(i,i) = 2/dx^2+1;
    A(i,i+1) = -1/dx^2;
end
A(1,1) = 2/dx^2+1;
A(1,2) = -1/dx^2;
A(end,end-1) = -1/dx^2;
A(end,end) = 2/dx^2+1;

f = x;
f(2) = f(2)-bc_L*A(2,1);
f(end-1) = f(end-1)-bc_R*A(end-1,end);

u = inv(A(2:end-1,2:end-1))*f(2:end-1)';

u = [bc_L; u; bc_R];

plot(x,u,'-ro','linewidth',2)
grid on