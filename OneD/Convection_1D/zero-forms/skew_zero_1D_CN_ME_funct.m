function antwoord=skew_zero_1D_CN_ME_funct(N,H,T,dt)

% addpath('C:\Users\Jasper.Kreeft\Data\MATLAB\Codes\MimeticSpectralElementMethod\OneD\Convection_1D\Library')

figure
% if ispc; figure('windowstyle','docked'); else figure; end

% N = 4;
% H = 12;
v = 1;
% T = 10;
L = 1;
dx = L/H;
% dt = 0.07; CFL = v*dt/dx
% CFL = 1.0; dt = CFL*dx/v
figplot = 1;

InitInfo.shape = 'cosine2'; %'step'; %'Gosse1'; %
InitInfo.X0    = [0.25 0.75];

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

% x = (Xi+1)/2*L;
% dxdxi = 1/2/H*L+0*Xi;
% xx = (XiXi+1)/2*L;
% dxxdxi = 1/2/H*L+0*XiXi;

P = 1;
% if even(P); disp('Warning: P must be odd!! Calculation stopped!!'); break; end
x = (Xi.^P+Xi+2)/4*L;
dxdxi = (P*Xi.^(P-1)+1)/4/H*L;
xx = (XiXi.^P+XiXi+2)/4*L;
dxxdxi = (P*XiXi.^(P-1)+1)/4/H*L;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = initialfunction(x(1:N*H),InitInfo);
pphi_ex = initialfunction(xx,InitInfo)';
phi_post = [phi' phi(1)];

plot(x,phi_post,'-o')
grid on
% break
L2error = zeros(1,ceil(T/dt)+1);
L2norm = 0;
Energy = zeros(1,ceil(T/dt)+1);
Mass = zeros(1,ceil(T/dt)+1);
for h=1:H
    ind = (h-1)*N+(1:N+1);
    ind2 = (h-1)*(nn-1)+(1:nn);
    pphi = phi_post(ind)*hh;
    L2error(1) = L2error(1) + sqrt(sum((ww.*dxxdxi(ind2)).*(pphi-pphi_ex(ind2)).^2));
    L2norm = L2norm + sqrt(sum((ww.*dxxdxi(ind2)).*pphi_ex(ind2).^2));
    Energy(1) = Energy(1) + 1/2*sum(pphi.*(ww.*dxxdxi(ind2)).*pphi);
    Mass(1)    = Mass(1) + sum(ww.*dxxdxi(ind2).*pphi);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = spdiags([-ones(N*H,1) ones(N*H,1)],[0 1],N*H,N*H); D(N*H,1) = 1;

E1 = zeros(N*H,N*H+1);
for h=1:H
    ind = (h-1)*N+(1:N+1);
    E1((h-1)*N+(1:N),(h-1)*N+(1:N+1)) = e/spdiags(dxdxi(ind)',0,N+1,N+1);
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


M1 = v*(A-A')/2;
% M1 = v*A;
M2 = W;



t=dt; j=1;
while t<T

    t = t+dt
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
        Mass(j)    = Mass(j) + sum(ww.*dxxdxi(ind2).*pphi);
        L2error(j) = L2error(j) + sqrt(sum((ww.*dxxdxi(ind2)).*(pphi-pphi_ex(ind2)).^2));

    %%%%%%%%%%%%%%%%%%%
        if figplot
        plot(x(ind),u_interp(ind),'xr')
        hold on
        plot(xx(ind2),uu_interp,'r')
        plot(x(ind),phi_post(ind),'o')
        plot(xx(ind2),pphi)
        ylim([-.2 1.2])
        end
    end
    if figplot
    pause(0.05)
    hold off
    end

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


% rmpath('C:\Users\Jasper.Kreeft\Data\MATLAB\Codes\MimeticSpectralElementMethod\OneD\Convection_1D\Library')

antwoord = L2error(j);

keyboard