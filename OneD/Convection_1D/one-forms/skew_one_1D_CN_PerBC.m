% JACOBIAN'S AND OTHER TRANSFORMATION MATRICES MUST BE INCLUDED

clear all
close all
clc

if ispc
path(path,'C:\Users\Jasper\Documents\MATLAB/MSEM_codes/OneD/Convection_1D/Library')
else isunix
path(path,'/media/My Passport/MSEM/MSEM_codes/OneD/Convection_1D/Library')
end

if ispc
    fig = figure('windowstyle','docked');
else
    fig = figure('Units','normalized','Position',[0, 0, 0.5, 1]);
end

dx = 2;
v = 1;
dt = 0.01; CFL = v*dt/dx

N = 30;
T = 100;  % Veelfout van 2

% InitInfo.shape = 'cosine1'; %'step'; %
% InitInfo.X0    = 0; %0; %[-0.5 0.5]

InitInfo.shape = 'cosine3'; %'step'; %
InitInfo.X0    = [0.25 0.75]; %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);
[h,e] = MimeticpolyVal(xi,N,1);

nn = 1000;
[xx,ww] = GLLnodes(nn-1);
[hh,ee] = MimeticpolyVal(xx,N,1);

x = xi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = initialfunction(x,InitInfo);

phidx = phi./diff(x)';
pphi = phi'*ee;

% pphi_ex = initialfunction(xx,InitInfo)';
pphi_ex =1/2+1/2*cos(4*pi*xx.*((xx>=0.25).*(xx<=0.75)))-((xx<0.25)+(xx>0.75));

L2error(1) = sqrt(sum(ww.*(pphi-pphi_ex).^2));
L2norm = sqrt(sum(ww.*pphi_ex.^2));
Energy(1) = 1/2*sum(pphi.*ww.*pphi);
Mass(1) = sum(ww.*pphi);

plot([x(2:N+1) ; x(1:N)],[phidx phidx]','g')
hold on
grid on
plot(xx,pphi,'r')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N); D(N,1) = 1;

E1 = e;
% Periodic BC
E1(:,1) = E1(:,1)+E1(:,N+1);
E1(:,N+1) = [];

% Central
E1(:,1) = E1(:,1)/2;

I1 = zeros(N);
for i=1:N
    for j=1:N
        I1(i,j) = sum(w.*e(i,:).*e(j,:));
    end
end

A = I1*D*E1';

M1 = v*0.5*(A-A');
% M1 = v*A;
M2 = I1;



t=0; j=1;
ref.shape = 'cosine2'; ref.X0 = InitInfo.X0;
while t<T

    t = t+dt;
    j = j+1;
    
    %%%%%%%%%%%%%%%%%%%%
    % nodally exact, interpolated solution
    s = xi-v*(j-1)*dt;
    z1 = (s<-1); z2 = (s>1);
    s = s - 2*fix((s-z1+z2)/2).*(z1+z2);
    u_interp = initialfunction(s,InitInfo);
    uu_interp = u_interp'*ee;
    
    ss = xx-v*(j-1)*dt;
    z1 = (ss<-1); z2 = (ss>1);
    ss = ss - 2*fix((ss-z1+z2)/2).*(z1+z2);
%     pphi_ex = initialfunction(ss,InitInfo)';
    pphi_ex =1/2+1/2*cos(4*pi*xx.*((ss>=0.25).*(ss<=0.75)))-((ss<0.25)+(ss>0.75));
    
    
   
    %%%%%%%%%%%%%%%%%%%%%

    phinew = timemarching(M2,M1,dt,phi,'ERK4b');
%     phinew = mimetictimemarching(M2,M1,dt,phi,4);
    
%     if ~exist('phiold','var')
%         phinew = timemarching(M2,M1,dt,phi,'BE');
%     else
%         phinew = timemarching(M2,M1,dt,[phi phiold],'BDF2');
%     end


    pphinew = phinew'*ee;
    

    
    Energy(j)  = 1/2*sum(pphinew.*ww.*pphinew);
    Mass(j)    = sum(ww.*pphinew);
    L2error(j) = sqrt(sum(ww.*(pphinew-pphi_ex).^2));

    phidx = phinew./diff(x)';
    
    
    plot(xx,uu_interp,'r','linewidth',4)
    hold on
    plot([x(2:N+1) ; x(1:N)],[phidx phidx]','g','linewidth',4)
    plot(xx,pphinew)
%         ylim([-.2 1.2])
    hold off
    pause(0.05)
    
%     avimovie('test',fig,j,t>=T);
    
    
    phiold = phi;
    phi    = phinew;
    

end

taxis = linspace(0,t,j);
figure
subplot(2,2,1)
plot(taxis,Energy,'-')
subplot(2,2,2)
plot(taxis,Energy-Energy(1),'-')
subplot(2,2,3)
plot(taxis,Mass,'-')
subplot(2,2,4)
plot(taxis,Mass-Mass(1),'-')

figure
% plot(taxis,L2error/L2norm,'-o')
plot(taxis,L2error,'-o')

rmpath('C:\Users\Jasper\Documents\MATLAB/MSEM_codes/OneD/Convection_1D/Library')


