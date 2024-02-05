clear all
close all
clc


P = 2;
Nx = 5;
% zeta = linspace(0,1,Nx);
zeta = [0 .25 .5 .5 .75 1]
x = linspace(zeta(1),zeta(end),1000);
xk = GrevilleAbscissa(zeta,P);
X = xk;
size(zeta)
size(xk)
% break

B = Bspline(X,zeta,P);

iB = inv(B);

[b,e] = Bspline(x,zeta,P);

SI0 = iB*b;

% F = cos(pi*X);
a = 0-eps;
F = (X>a);

f = F*SI0;

K = F*iB;

figure%(1)
plot(x,x>a)
hold on
plot(x,f,'r')
plot(X,F,'or','markerface','r')
plot(xk,K,'-sk','linewidth',1,'markerface','k')
% break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% X = linspace(zeta(1),zeta(end),length(zeta)+P-2);
% 
% [~,E] = Bspline(X,zeta,P);
% 
% iE = inv(E);
% 
% SI1 = iE*e;

figure%(2)
% subplot(2,2,1)
plot(x,b')
hold on
plot(xk,zeros(size(xk)),'sk')
plot(zeta,zeros(size(zeta)),'ok','markerface','k')

% subplot(2,2,2)
figure
plot(x,SI0')
hold on
plot(xk,zeros(size(xk)),'sk')
% plot(zeta,zeros(size(zeta)),'ok','markerface','k')
% subplot(2,2,3)
% plot(x,e')
% subplot(2,2,4)
% plot(x,SI1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xgl,wgl] = GLLnodes(20);

B_L2 = Bspline(xgl,zeta,P);

M0l = zeros(P+Nx-1);
for i=1:P+Nx-1
    for j=1:P+Nx-1
        M0l(i,j) = sum(wgl.*B_L2(i,:).*B_L2(j,:));
    end
end


M0r = zeros(P+Nx-1);
for j=1:P+Nx-1
    M0r(j,j) = sum(wgl.*B_L2(j,:));
end


C_L2 = inv(M0l)*M0r*F'


cc = C_L2'*b;