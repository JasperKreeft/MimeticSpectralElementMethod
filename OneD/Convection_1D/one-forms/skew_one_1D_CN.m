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

N = 20;
T = 1;  % Veelfout van 2

InitInfo.shape = 'cosine1'; %'step'; %
InitInfo.X0    = 0; %0; %[-0.5 0.5]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);
[h,e] = MimeticpolyVal(xi,N,1);

nn = 1000;
[xx,ww] = GLLnodes(nn-1);
[hh,ee] = MimeticpolyVal(xx,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = zeros(N,1);
for i=1:N
    if xi(i+1)<=0
        phi(i,1) = 1/2*(xi(i+1)-xi(i))-1/(4*pi)*(sin(2*pi*xi(i+1))-sin(2*pi*xi(i)));
    elseif xi(i+1)>0 && xi(i)<0
        phi(i,1) = 1/2*(0-xi(i))-1/(4*pi)*(0-sin(2*pi*xi(i)));
    else
        phi(i,1) = 0;
    end
end

pphi = phi'*ee;

phidx = phi./diff(xi)';
plot([xi(2:N+1) ; xi(1:N)],[phidx phidx]','g')
hold on
grid on

pphi_ex = initialfunction(xx,InitInfo)';

L2error(1) = sqrt(sum(ww.*(pphi-pphi_ex).^2));
L2norm = sqrt(sum(ww.*pphi_ex.^2));
Energy(1) = 1/2*sum(pphi.*ww.*pphi);
Mass(1) = sum(ww.*pphi);

plot(xx,pphi)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N); D(N,1) = 1;
% E = e(:,1:N);
% E(:,1) = E(:,1)+e(:,N+1);
% wp = w(1:N)';
% W = spdiags(wp,0,N,N);
% A = E*W*E'*D*E';
% K2 = E*W*E';
% K1 = M2*D*E';


D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);
E = e;
W = spdiags(w',0,N+1,N+1);
M1 = E*W*E';
K1 = M1*D*E';
K2 = M1;


t=0; j=1;
while t<T

    t = t+dt;
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
    
    plot(xi,u_interp,'xr')
    hold on
    plot(xx,uu_interp,'r')
    ylim([-1 1])
    %%%%%%%%%%%%%%%%%%%%%

    phinew = timemarching(K2,K1,dt,phi,'gauss4');
    
    pphinew = phinew'*ee;
    
    Energy(j)  = 1/2*sum(pphinew.*ww.*pphinew);
    Mass(j)    = sum(ww.*pphinew);
    L2error(j) = sqrt(sum(ww.*(pphinew-pphi_ex).^2));

    phidx = phinew./diff(xi)';
    plot([xi(2:N+1) ; xi(1:N)],[phidx phidx]','g')
    plot(xx,pphinew)
%     ylim([-.2 1])
    pause(0.05)
    hold off
    phi = phinew;

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
plot(taxis,L2error/L2norm,'-o')

% rmpath('C:\Users\Jasper\Documents\MATLAB/MSEM_codes/OneD/Convection_1D/Library')