clear all
close all
clc

% if ispc
% path(path,'C:\Users\Jasper\Documents\MATLAB/MSEM_codes/OneD/Convection_1D/Library')
% else isunix
% path(path,'/media/My Passport/MSEM/MSEM_codes/OneD/Convection_1D/Library')
% end

if ispc; figure('windowstyle','docked'); else figure; end

dx = 2;
v = 1;
dt = 0.1; CFL = v*dt/dx
% CFL = 0.1; dt = CFL*dx/v

N = 20;
T = 10;  % Veelfout van 2

InitInfo.shape = 'cosine1'; %'step'; %
InitInfo.X0    = 0; %0; %[-0.5 0.5]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);
[h,e] = MimeticpolyVal(xi,N,1);

nn = 1000;
[xx,ww] = GLLnodes(nn-1);
hh = MimeticpolyVal(xx,N,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = initialfunction(xi(1:N),InitInfo);
pphi_ex = initialfunction(xx,InitInfo)';

pphi = [phi' phi(1)]*hh;

L2error(1) = sqrt(sum(ww.*(pphi-pphi_ex).^2));
L2norm = sqrt(sum(ww.*pphi_ex.^2));
Energy(1) = 1/2*sum(pphi.*ww.*pphi);
Mass(1) = sum(ww.*pphi);

plot(xi,[phi ; phi(1)],'-o')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N); D(N,1) = 1;

E1 = e;
E1(:,1) = E1(:,1)+E1(:,N+1);
E1(:,N+1) = [];

% Central
E1(:,1) = E1(:,1)/2;

We = zeros(N+1,1);
ind2 = (1:N+1);
We = w';
We(1) = We(1)+We(N+1);
We(N+1) = [];
W = spdiags(We,0,N,N);

A = W*E1'*D;

%%%


M = v*(A-A')/2;
% M = v*A;

t=0; j=1;
while t<T

    t = t+dt
    j = j+1;
    
    %%%%%%%%%%%%%%%%%%%%
    % nodally exact, interpolated solution
    s = xi-v*(j-1)*dt;
    z1 = (s<-1); z2 = (s>1);
    s = s - 2*fix((s-z1+z2)/2).*(z1+z2);
    u_interp = initialfunction(s,InitInfo);
    uu_interp = u_interp'*hh;
    
    ss = xx-v*(j-1)*dt;
    z1 = (ss<-1); z2 = (ss>1);
    ss = ss - 2*fix((ss-z1+z2)/2).*(z1+z2);
    pphi_ex = initialfunction(ss,InitInfo)';
    
    %%%%%%%%%%%%%%%%%%%%%

    phinew = timemarching(W,M,dt,phi,'gauss4');%ESDIRK4BE
%         phinew = mimetictimemarching(W,M,dt,phi,5);

%     phiold = phi;
    phi    = phinew;
    
    pphi = [phi' phi(1)]*hh;

    Energy(j) = 1/2*sum(pphi.*ww.*pphi);
    Mass(j)    = sum(ww.*pphi);
    L2error(j) = sqrt(sum(ww.*(pphi-pphi_ex).^2));
    
    %%%%%%%%%%%%%%%%%%%%%

    plot(xi,u_interp,'xr')
    hold on
    plot(xx,uu_interp,'r')
    plot(xx,pphi_ex,'g')
    plot(xi,[phi ; phi(1)],'o')
    plot(xx,pphi)
    ylim([-.2 1.2])
    pause(0.05)
    hold off

end

taxis = linspace(0,t,j);
figure
subplot(2,2,1)
plot(taxis,Energy,'-o')
subplot(2,2,2)
plot(taxis,Energy-Energy(1),'-o')
subplot(2,2,3)
plot(taxis,Mass,'-o')
subplot(2,2,4)
plot(taxis,Mass-Mass(1),'-o')

figure
% plot(taxis,L2error/L2norm,'-o')
loglog(taxis,L2error/L2norm,'-o')

% rmpath('C:\Users\Jasper\Documents\MATLAB/MSEM_codes/OneD/Convection_1D/Library')