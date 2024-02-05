clear all
clf%ose all
clc

N = 8;
x = linspace(0,1,N+1);
X = x'*ones(1,N+1);
Y = X';

F0 = sin(pi*X).*sin(pi*Y);
F1x = -1/pi*cos(pi*X).*sin(pi*Y);
F1y = -1/pi*sin(pi*X).*cos(pi*Y);
F2 = 1/pi^2*cos(pi*X).*cos(pi*Y);


XX = []; YY = []; FF = [];
for i=1:N+1
    XX = [ XX ; X(1:N,i) X(2:N+1,i) ];
    YY = [ YY ; Y(1:N,i) Y(2:N+1,i) ];
    FF = [ FF ; F1x(2:N+1,i)-F1x(1:N,i) F1x(2:N+1,i)-F1x(1:N,i) ];
end

for i=1:N+1
    XX = [ XX ; X(i,1:N)' X(i,2:N+1)' ];
    YY = [ YY ; Y(i,1:N)' Y(i,2:N+1)' ];
    FF = [ FF ; (F1y(i,2:N+1)-F1y(i,1:N))' (F1y(i,2:N+1)-F1y(i,1:N))' ];
end

% figure
subplot(2,2,2)
dotcolorplot(F0,X,Y,0,1)
view([0 0 1])
grid off
axis equal
axis tight

subplot(2,2,3)
linecolorplot(FF,XX,YY)
view([0 0 1])
grid off
axis equal
axis tight

subplot(2,2,4)
fmin = 0; fmax = 0;
for i=1:N
    for j=1:N
        XXX = [ X(i,j) X(i+1,j) ; X(i,j+1) X(i+1,j+1) ];
        YYY = [ Y(i,j) Y(i+1,j) ; Y(i,j+1) Y(i+1,j+1) ];
        FFF = (F2(i+1,j+1)-F2(i,j+1)-F2(i+1,j)+F2(i,j))*ones(2);
        surf(XXX,YYY,FFF)
        hold on
        fmin = min(fmin,FFF(1,1));
        fmax = max(fmax,FFF(1,1));
    end
end
colorbar('yticklabel',roundn(linspace(fmin,fmax,6),-3))
view([0 0 1])
grid off
axis equal
axis tight


xx = linspace(0,1,1000)'*ones(1,1000);
yy = xx';
ff = sin(pi*xx).*sin(pi*yy);
subplot(2,2,1)
surf(xx,yy,ff)
shading interp
colorbar('yticklabel',linspace(0,1,6))
view([0 0 1])
grid off
axis equal
axis tight