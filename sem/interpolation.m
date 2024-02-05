clear all
% clf%ose all
clc

x = linspace(-1,1,200);
X = x'*ones(1,200);
Y = X';

N = 8;

xG = Gnodes(N);
xEG = [-1 xG 1];
XEG = xEG'*ones(1,N+2);
YEG = XEG';

hEG = LagrangeVal(x,N,3);
hG = LagrangeVal(x,N,2);

% phi = ones(N+2);

phi = sin(pi*3/2*XEG).*sin(pi*3/2*YEG);

phi_ex = sin(pi*3/2*X).*sin(pi*3/2*Y);



% % case 1
% Z = hEG'*phi*hEG;

% % case 2
% Z = hG'*phi(2:N+1,2:N+1)*hG;

% % case 3
% Z = hEG(2:N+1,:)'*phi(2:N+1,2:N+1)*hEG(2:N+1,:);
% Z = Z+hG'*phi(2:N+1,1)*hEG(1,:);
% Z = Z+hG'*phi(2:N+1,N+2)*hEG(N+2,:);
% Z = Z+hEG(1,:)'*phi(1,2:N+1)*hG;
% Z = Z+hEG(N+2,:)'*phi(N+2,2:N+1)*hG;

% case 4
Z = hEG(2:N+1,:)'*phi(2:N+1,2:N+1)*hEG(2:N+1,:);
Z = Z+(hG+hEG(2:N+1,:))'/2*phi(2:N+1,1)*hEG(1,:);
Z = Z+(hG+hEG(2:N+1,:))'/2*phi(2:N+1,N+2)*hEG(N+2,:);
Z = Z+hEG(1,:)'*phi(1,2:N+1)*(hG+hEG(2:N+1,:))/2;
Z = Z+hEG(N+2,:)'*phi(N+2,2:N+1)*(hG+hEG(2:N+1,:))/2;


error = (phi_ex-Z);

figure(1)
surf(X,Y,error)%phi_ex
% zlim([0.8 1.2])
% set(gca,'clim',zlim)
shading interp
colorbar

figure(2)
plot(X(:,1),error(:,1))