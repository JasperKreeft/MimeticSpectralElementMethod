% JACOBIAN'S AND OTHER TRANSFORMATION MATRICES MUST BE INCLUDED

clear all
close all
clc

if ispc
addpath('\\AMSDC1-NA-V511\Jasper.Kreeft$\Cached\My Documents\MATLAB\MSEM\OneD\Convection_1D\Library\')
elseif unix
path(path,'/media/My Passport/MSEM/MSEM_codes/OneD/Convection_1D/Library')
end


if ispc; figure('windowstyle','docked'); else figure; end

N = 1;
H = 11;
v = 1;
T = 10;
L = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

dxdxi = L/(2*H);

x = [];
for h=1:H
    xe = (2*h-1)*dxdxi+xi*dxdxi;
    x  = [ x(1:end-1) xe ]; 
end
plot(x,zeros(1,N*H+1),'sk')
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = spdiags([-ones(N*H+1,1) ones(N*H+1,1)],[0 1],N*H,N*H+1);

E1 = zeros(N*H,N*H+1);
for h=1:H
    E1((h-1)*N+(1:N),(h-1)*N+(1:N+1)) = e;
end
% E1(:,1) = E1(:,1)+E1(:,N*H+1);
% E1(:,N*H+1) = [];

% Central
for h=1:H
    E1(:,(h-1)*N+1) = E1(:,(h-1)*N+1)/2;
end
% Upwind
% for h=1:H-1
%     E1(h*N+(1:N),h*N+1) = zeros(N,1);
% end
% % Downwind
% for h=1:H-1
%     E1((h-1)*N+(1:N),h*N+1) = zeros(N,1);
% end
% % 1/2 Upwind, 1/2 central
% a = 2;      % a=2 is central 
% E1(1:N,1)     = E1(1:N,1)*(a-1)/a;
% E1(N+(1:N),1) = E1(N+(1:N),1)/a;
% for h=1:H-1
%     E1((h-1)*N+(1:N),h*N+1) = E1((h-1)*N+(1:N),h*N+1)*(a-1)/a;
%     E1(h*N+(1:N),h*N+1) = E1(h*N+(1:N),h*N+1)/a;
% end

We = zeros(N*H+1,1);
for h=1:H
    ind2 = (h-1)*N+(1:N+1);
    We(ind2,1) = We(ind2,1)+w';
end
% We(1) = We(1)+We(N*H+1);
% We(N*H+1) = [];
W = spdiags(We,0,N*H+1,N*H+1);

A = W*E1'*D;


M1 = (v/dxdxi)*(A-A')/2;
% M1 = (v/dxdxi)*A;


f = cos(pi*x)*pi;
% plot(x,f)

M1(:,1) = []; M1(1,:) = []; W(:,1) = []; W(1,:) = []; f(1) = [];

y = inv(M1)*W*f';

y = [ 0 ; y];

plot(x,y)

% nn = 100;
% xixi = linspace(-1,1,nn+1);
% hh = MimeticpolyVal(xixi,N,1);
% 
% xx = linspace(x(1),x(end),H*nn+1);
% 
% 
% %     phi_post = phi';
%     phi_post = [phi' phi(1)];
%     for h=1:H
%     ind = (h-1)*N+(1:N+1);
%     pphi = phi_post(ind)*hh;
% 
%     plot(x(ind),phi_post(ind),'o')
%     ind = (h-1)*nn+(1:nn+1);
%     plot(xx(ind),pphi)
%     end
% %     ylim([-.2 1.2])
% %     ylim([.98 1.02])
% %     ylim([-.01 .02])
%     pause(0.1)
%     hold off

if ispc
rmpath('\\AMSDC1-NA-V511\Jasper.Kreeft$\Cached\My Documents\MATLAB\MSEM\OneD\Convection_1D\Library\')
elseif unix
rmpath('/media/My Passport/MSEM/MSEM_codes/OneD/Convection_1D/Library')
end