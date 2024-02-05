clear all
close all
clc



N = 14;

xi = GLLnodes(N);
eta = xi;

Xi  = xi'*ones(1,N+1);
Eta = ones(N+1,1)*eta;

subplot(1,3,1)
plot3(Xi,Eta,zeros(N+1),'bx')
axis equal
view([0 0 1])

a = 0.5;
b = 0.5;
R = b*(Xi+1)/2+a;
k = 1;
T = (Eta+1)/2*pi/2+k*pi/2;

subplot(1,3,2)
plot3(R,T,zeros(N+1),'bx')
axis equal
view([0 0 1])

X = R.*cos(T);
Y = R.*sin(T);

subplot(1,3,3)
plot3(X,Y,zeros(N+1),'bx')
axis equal
view([0 0 1])

%%

figure
for k = 0:3;
T = (Eta+1)/2*pi/2+k*pi/2;

X = R.*cos(T);
Y = R.*sin(T);

plot3(X,Y,zeros(N+1),'bx')
hold on
axis equal
view([0 0 1])
end