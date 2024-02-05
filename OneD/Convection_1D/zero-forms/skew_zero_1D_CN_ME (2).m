clear all
close all
clc

if ispc
path(path,'C:\Users\Jasper\Documents\MATLAB/MSEM_codes/OneD/Convection_1D/Library')
elseif isunix
addpath('/home/jjkreeft/Convection/OneD/Convection_1D/Library')
end

if ispc; figure('windowstyle','docked'); else figure; end

N = 8;
H = 6;
v = 1;
T = 40;
L = 2;
dx = L/H;
dt = 0.05; CFL = v*dt/dx
% CFL = 1.2; dt = CFL*dx/v


InitInfo.shape = 'cosine1'; %'step'; %'Gosse1'; %
InitInfo.X0    = 1;%[0.25 0.75];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

nn = 100;
[xixi,ww] = GLLnodes(nn-1);
hh = MimeticpolyVal(xixi,N,1);

% % X = GLLnodes(H)+1;
% X = linspace(0,2,H+1);
% % X = sign((X-1)).*((X-1).^2)+1;
% x = [];
% dxdxi = zeros(1,H);
% for h=1:H
%     dxdxi(h) = (X(h+1)-X(h))/2;
%     xe = (X(h)+X(h+1))/2+dxdxi(h)*xi;
%     x = [ x(1:end-1) xe ];
% end
% plot(x,zeros(1,N*H+1),'sk')
% 
% 
% nn = 100;
% [xixi,ww] = GLLnodes(nn-1);
% hh = MimeticpolyVal(xixi,N,1);
% 
% xx = zeros(1,H*(nn-1)+1);
% for h=1:H
%     ind = (h-1)*(nn-1)+(1:nn);
%     xx(ind) = (X(h)+X(h+1))/2 + dxdxi(h)*xixi;
% end

% dxdxi = 


Xi = []; XiXi = [];
for h=1:H
    Xi = [ Xi(1:end-1) (2*h-1)*L/(2*H)+L/(2*H)*xi ];
    XiXi = [ XiXi(1:end-1) (2*h-1)*L/(2*H)+L/(2*H)*xixi ];
end
% X = sign((Xi-1)).*((Xi-1).^2)+1;
% dXdXi = 2*sign((Xi-1)).*(Xi-1);

% x = Xi;
% dxdxi = ones(size(Xi))/H;
% xx = XiXi;
% dxxdxi = ones(size(XiXi))/H;

x = ((Xi-1).^5+Xi+1)/2;
dxdxi = (5*(Xi-1).^4+1)/2/H;
xx = ((XiXi-1).^5+XiXi+1)/2;
dxxdxi = (5*(XiXi-1).^4+1)/2/H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = initialfunction(x(1:N*H),InitInfo);
pphi_ex = initialfunction(xx,InitInfo)';
phi_post = [phi' phi(1)];

plot(x,phi_post,'-o')
grid on

L2error = zeros(1,ceil(T/dt)+1);
L2norm = 0;
Energy = zeros(1,ceil(T/dt)+1);
for h=1:H
    ind = (h-1)*N+(1:N+1);
    ind2 = (h-1)*(nn-1)+(1:nn);
    pphi = phi_post(ind)*hh;
    L2error(1) = L2error(1) + sqrt(sum((ww.*dxxdxi(ind2)).*(pphi-pphi_ex(ind2)).^2));
    L2norm = L2norm + sqrt(sum((ww.*dxxdxi(ind2)).*pphi_ex(ind2).^2));
    Energy(1) = Energy(1) + 1/2*sum(pphi.*(ww.*dxxdxi(ind2)).*pphi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = spdiags([-ones(N*H,1) ones(N*H,1)],[0 1],N*H,N*H); D(N*H,1) = 1;

E1 = zeros(N*H,N*H+1);
for h=1:H
    ind = (h-1)*N+(1:N+1);
    E1((h-1)*N+(1:N),(h-1)*N+(1:N+1)) = e*spdiags(1./dxdxi(ind)',0,N+1,N+1);
end
E1(:,1) = E1(:,1)+E1(:,N*H+1);
E1(:,N*H+1) = [];

% Central
for h=1:H
    E1(:,(h-1)*N+1) = E1(:,(h-1)*N+1)/2;
end
% % Upwind
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
%     E1((h-1)*N+(1:N),h*N+1) = E1((h-1)*N+(1:N),h*N+1)*(a-1)/a * dxdxi(h+1)/dxdxi(h); % Upw
%     E1(h*N+(1:N),h*N+1) = E1(h*N+(1:N),h*N+1)/a * dxdxi(h)/dxdxi(h+1); % Dwnw
% end


We = zeros(N*H+1,1);
for h=1:H
    ind = (h-1)*N+(1:N+1);
    We(ind,1) = We(ind,1)+(w./dxdxi(ind))';
end
We(1) = We(1)+We(N*H+1);
We(N*H+1) = [];
W = spdiags(We,0,N*H,N*H);

A = W*E1'*D;


% M1 = v*(A-A')/2;
M1 = v*A;
M2 = W;



t=dt; j=1;
while t<T

    t = t+dt;
    j = j+1;

    %%%%%%%%%%%%%%%%%%%%
    % nodally exact, interpolated solution
    s = x-v*(j-1)*dt;
    z1 = (s<0); z2 = (s>L);
    s = s - L*fix(s/L-z1+z2).*(z1+z2);
    u_interp = initialfunction(s,InitInfo)';
    
    ss = xx-v*(j-1)*dt;
    z1 = (ss<0); z2 = (ss>L);
    ss = ss - L*fix(ss/L-z1+z2).*(z1+z2);
    pphi_ex = initialfunction(ss,InitInfo)';
    
    

    %%%%%%%%%%%%%%%%%%%%%
%     if j==2
%         phinew = timemarching(M2,M1,dt,phi,'BE');
%     else
%         phinew = timemarching(M2,M1,dt,[phi phiold],'BDF2');
%     end
%     phiold = phi;

phinew = timemarching(M2,M1,dt,phi,'gauss4');%
% phinew = mimetictimemarching(M2,M1,dt,phi,3);

    phi    = phinew;

%     phi_post = phi';
    phi_post = [phi' phi(1)];
    for h=1:H
        
        ind = (h-1)*N+(1:N+1);
        ind2 = (h-1)*(nn-1)+(1:nn);
        
        uu_interp = u_interp(ind)*hh;
        pphi = phi_post(ind)*hh;

        Energy(j) = Energy(j) + 1/2*sum(pphi.*(ww.*dxxdxi(ind2)).*pphi);

        L2error(j) = L2error(j) + sqrt(sum((ww.*dxxdxi(ind2)).*(pphi-pphi_ex(ind2)).^2));

    %%%%%%%%%%%%%%%%%%%
    
        plot(x(ind),u_interp(ind),'xr')
        hold on
        plot(xx(ind2),uu_interp,'r')
        plot(x(ind),phi_post(ind),'o')
        plot(xx(ind2),pphi)
        ylim([-.2 1.2])

    end
    pause(0.05)
    hold off

end

taxis = linspace(0,t,j);
figure
subplot(1,2,1)
plot(taxis,Energy(1:j),'-o')
subplot(1,2,2)
plot(taxis,Energy(1:j)-Energy(1),'-o')

figure
plot(taxis,L2error(1:j)/L2norm,'-o')


if ispc
rmpath('C:\Users\Jasper\Documents\MATLAB/MSEM_codes/OneD/Convection_1D/Library')
elseif isunix
rmpath('/home/jjkreeft/Convection/OneD/Convection_1D/Library')
end