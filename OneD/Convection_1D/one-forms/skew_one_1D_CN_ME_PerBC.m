% JACOBIAN'S AND OTHER TRANSFORMATION MATRICES MUST BE INCLUDED

clear all
close all
clc

% if ispc
% path(path,'C:\Users\Jasper\Documents\MATLAB/MSEM_codes/OneD/Convection_1D/Library')
% else isunix
% path(path,'/media/My Passport/MSEM/MSEM_codes/OneD/Convection_1D/Library')
% end

if ispc
    fig = figure('windowstyle','docked');
else
    fig = figure('Units','normalized','Position',[0, 0, 0.5, 1]);
end

N = 4;
H = 12;
v = 1;
T = 10;  T = T-1e-10;
L = 1;
dx = L/H;
dt = 0.05; CFL = v*dt/dx
figplot = 1;

InitInfo.shape = 'cosine3'; %'step'; %
InitInfo.X0    = [0.25 0.75]; %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

nn = 1000;
[xixi,ww] = GLLnodes(nn-1);
[hh,ee] = MimeticpolyVal(xixi,N,1);

Xi = []; XiXi = [];
for h=1:H
    Xi = [ Xi(1:end-1) (2*h-1)/H+1/H*xi-1 ];
    XiXi = [ XiXi(1:end-1) (2*h-1)/H+1/H*xixi-1 ];
end

x = (Xi+1)/2;
dxdxi = 1/2/H+0*Xi;
xx = (XiXi+1)/2;
dxxdxi = 1/2/H+0*XiXi;

% x = (Xi.^3+Xi+2)/4;
% dxdxi = (3*Xi.^2+1)/4/H;
% xx = (XiXi.^3+XiXi+2)/4;
% dxxdxi = (3*XiXi.^2+1)/4/H;

% x = (Xi.^5+Xi+2)/4;
% dxdxi = (5*Xi.^4+1)/4/H;
% xx = (XiXi.^5+XiXi+2)/4;
% dxxdxi = (5*XiXi.^4+1)/4/H;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = initialfunction(x,InitInfo);
pphi_ex =1/2+1/2*cos(4*pi*xx.*((xx>=0.25).*(xx<=0.75)))-((xx<0.25)+(xx>0.75));

if figplot
phidx = phi./diff(x)';
plot([x(2:N*H+1) ; x(1:N*H)],[phidx phidx]','g')
hold on
grid on
end

L2error = zeros(1,ceil(T/dt)+1);
L2norm = 0;
Energy = zeros(1,ceil(T/dt)+1);
Mass = zeros(1,ceil(T/dt)+1);
for h=1:H
    ind1 = (h-1)*N+(1:N);
    ind2 = (h-1)*(nn-1)+(1:nn);
    pphi(ind2) = (phi(ind1)'*ee)./dxxdxi(ind2);
    L2error(1) = L2error(1) + sqrt(sum((ww.*dxxdxi(ind2)).*(pphi(ind2)-pphi_ex(ind2)).^2));
    L2norm     = L2norm + sqrt(sum((ww.*dxxdxi(ind2)).*pphi_ex(ind2).^2));
    Energy(1)  = Energy(1) + 1/2*sum(pphi(ind2).*(ww.*dxxdxi(ind2)).*pphi(ind2));
    Mass(1)    = Mass(1) + sum(ww.*dxxdxi(ind2).*pphi(ind2));
    if figplot
        plot(xx(ind2),pphi(ind2),'r')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = spdiags([-ones(N*H,1) ones(N*H,1)],[0 1],N*H,N*H); D(N*H,1) = 1;

E1 = zeros(N*H,N*H+1);
for h=1:H
    ind = (h-1)*N+(1:N+1);
    E1((h-1)*N+(1:N),(h-1)*N+(1:N+1)) = e/spdiags(dxdxi(ind)',0,N+1,N+1);
end
% Periodic BC
E1(:,1) = E1(:,1)+E1(:,N*H+1);
E1(:,N*H+1) = [];

% Central
for h=1:H
    E1(:,(h-1)*N+1) = E1(:,(h-1)*N+1)/2;
end

I1 = [];
for h=1:H
    I1_ = zeros(N);
    for i=1:N
        for j=1:N
%             ind = (h-1)*N+(1:N+1);
%             I1_(i,j) = sum(w.*e(i,:).*e(j,:)./dxdxi(ind));
            ind = (h-1)*(nn-1)+(1:nn);
            I1_(i,j) = sum(ww.*ee(i,:).*ee(j,:)./dxxdxi(ind));
        end
    end
    I1 = blkdiag(I1,I1_);
end

A = I1*D*E1';

M1 = v*0.5*(A-A');
% M1 = v*A;
M2 = I1;


t=0; j=1;
ref.shape = 'cosine2'; ref.X0 = InitInfo.X0;
while t<T

    t = t+dt
    j = j+1;
    
    %%%%%%%%%%%%%%%%%%%%
    % nodally exact, interpolated solution

    ss = xx-v*(j-1)*dt;
    z1 = (ss<0); z2 = (ss>L);
    ss = ss - L*fix(ss/L-z1+z2).*(z1+z2);
    pphi_ex = initialfunction(ss,ref)';
    
    %%%%%%%%%%%%%%%%%%%%%

    phinew = timemarching(M2,M1,dt,phi,'ESDIRK4'); %gauss4
    
%     if ~exist('phiold','var')
%         phinew = timemarching(M2,M1,dt,phi,'CN');
%     else
%         phinew = timemarching(M2,M1,dt,[phi phiold],'BDF2');
%     end

    for h=1:H
        ind1 = (h-1)*N+(1:N);
        ind2 = (h-1)*(nn-1)+(1:nn);
        ind3 = (h-1)*N+(1:N+1);
        pphinew(ind2) = (phinew(ind1)'*ee)./dxxdxi(ind2);

        if figplot
            phidx = phinew(ind1)./diff(x(ind3))';
            plot(xx(ind2),pphi_ex(ind2),'r','linewidth',4)
            hold on
            plot([x(ind1+1) ; x(ind1)],[phidx phidx]','g','linewidth',4)
            plot(xx(ind2),pphinew(ind2))
            ylim([-.2 1.2])
        end
        
        Energy(j) = Energy(j) + 1/2*sum(pphinew(ind2).*(ww.*dxxdxi(ind2)).*pphinew(ind2));
        Mass(j)    = Mass(j) + sum(ww.*dxxdxi(ind2).*pphinew(ind2));
        L2error(j) = L2error(j) + sqrt(sum((ww.*dxxdxi(ind2)).*(pphinew(ind2)-pphi_ex(ind2)).^2));
        
    end
    if figplot
        pause(0.15)
        hold off
    end
    
    phiold = phi;
    phi    = phinew;

end

    for h=1:H
        ind1 = (h-1)*N+(1:N);
        ind2 = (h-1)*(nn-1)+(1:nn);
        ind3 = (h-1)*N+(1:N+1);
        pphinew(ind2) = (phinew(ind1)'*ee)./dxxdxi(ind2);
        phidx = phinew(ind1)./diff(x(ind3))';
        plot(xx(ind2),pphi_ex(ind2),'r','linewidth',4)
        hold on
        plot([x(ind1+1) ; x(ind1)],[phidx phidx]','g','linewidth',4)
        plot(xx(ind2),pphinew(ind2))
        ylim([-.2 1.2])
    end

taxis = linspace(0,t,j);
figure
subplot(2,2,1)
plot(taxis,Energy(1:j),'-o')
subplot(2,2,2)
plot(taxis,Energy(1:j)-Energy(1),'-o')
subplot(2,2,3)
plot(taxis,Mass(1:j),'-o')
subplot(2,2,4)
plot(taxis,Mass(1:j)-Mass(1),'-o')

figure
plot(taxis,L2error(1:j)/L2norm,'-o')

L2error(end)