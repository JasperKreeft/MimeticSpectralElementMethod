clear all
close all
clc

P = 2;
Nx = 5;
xi = linspace(0,1,Nx)
x = linspace(xi(1),xi(end),1000);
xk = GrevilleAbscissa(xi,P);
X = xk;


[xw,ww] = GLLnodes(10);
D = zeros(Nx);
for i=1:length(xk)-1
    
    xww = (xk(i)+xk(i+1))/2+(xk(i+1)-xk(i))/2*xw;  
    
[~,dw] = Bspline(xww,xi,P);

for j=1:size(dw,1)
    D(i,j) = sum(ww.*dw(j,:));
end

end
% keyboard


iD = inv(D);

[~,d] = Bspline(x,xi,P);

SI1 = iD*d;

% % F = cos(pi*X);
% a = 0-eps;
% F = (X>a);
% 
% f = F*SI1;
% 
% K = F*iD;
% 
% figure%(1)
% plot(x,x>a)
% hold on
% plot(x,f,'r')
% plot(X,F,'or','markerface','r')
% plot(xk,K,'-sk','linewidth',1,'markerface','k')
% % break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% X = linspace(xi(1),xi(end),length(xi)+P-2);
% 
% [~,E] = Bspline(X,xi,P);
% 
% iE = inv(E);
% 
% SI1 = iE*e;

figure%(2)
% subplot(2,2,1)
plot(x,d')
hold on
plot(xk,zeros(size(xk)),'sk')
plot(xi,zeros(size(xi)),'ok','markerface','k')

% subplot(2,2,2)
figure
plot(x,SI1')
hold on
plot(xk,zeros(size(xk)),'sk')
plot(xi,zeros(size(xi)),'ok','markerface','k')
% subplot(2,2,3)
% plot(x,e')
% subplot(2,2,4)
% plot(x,SI1')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [xgl,wgl] = GLLnodes(20);
% 
% B_L2 = Bspline(xgl,xi,P);
% 
% M0l = zeros(P+Nx-1);
% for i=1:P+Nx-1
%     for j=1:P+Nx-1
%         M0l(i,j) = sum(wgl.*B_L2(i,:).*B_L2(j,:));
%     end
% end
% 
% 
% M0r = zeros(P+Nx-1);
% for j=1:P+Nx-1
%     M0r(j,j) = sum(wgl.*B_L2(j,:));
% end
% 
% 
% C_L2 = inv(M0l)*M0r*F'
% 
% 
% cc = C_L2'*b;