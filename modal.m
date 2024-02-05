%%
clear all
close all
clc

N = 22;

n=1:2:N;

sgn = zeros(size(n)); sgn(1) = 1;
for i=2:length(n)
    sgn(i) = -sgn(i-1);
end

f = zeros(N,1);
f(n) = sgn.*(pi).^n./factorial(n);

% plot(x,f)

xx = linspace(0,1,1000);

FF = f(1)*xx;
for i=2:N
    FF = FF + f(i)*xx.^i;
end

plot(xx,FF)
% break



A = zeros(N);
B = zeros(N);
for i=1:N
    for j=1:N
        A(i,j) = i*j/(i+j-1);
        B(i,j) = 1/(i+j+1);
    end
end

C = B\A;

% c = cond(C)

coeff = (A\B)*f;
% coeff = C\f;

xx = linspace(0,1,1000);

FFh = coeff(1)*xx;
for i=2:N
    FFh = FFh + coeff(i)*xx.^i;
end
% break

% FF = sin(2*pi*XX);

hold on
plot(xx,FFh,'r')