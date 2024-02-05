% JACOBIAN'S AND OTHER TRANSFORMATION MATRICES MUST BE INCLUDED

clear all
close all
clc

if ispc; figure('windowstyle','docked'); else figure; end

N = 7;
H = 20;
v = 1;
T = 5;
L = 4;
dt = .05;

InitInfo.shape = 'step'; %'cosine0'; %
InitInfo.X0    = [0.25 0.75]; %

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
% Initial condition

phi = initialfunction(x(1:N*H),InitInfo);

plot(x,[phi ; phi(1)],'o')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = spdiags([-ones(N*H,1) ones(N*H,1)],[0 1],N*H,N*H); D(N*H,1) = 1;

E1 = zeros(N*H,N*H+1);
for h=1:H
    E1((h-1)*N+(1:N),(h-1)*N+(1:N+1)) = e;
end
E1(:,1) = E1(:,1)+E1(:,N*H+1);
E1(:,N*H+1) = [];

% % Central
% for h=1:H
%     E1(:,(h-1)*N+1) = E1(:,(h-1)*N+1)/2;
% end
% % Upwind
% for h=1:H-1
%     E1(h*N+(1:N),h*N+1) = zeros(N,1);
% end
% % Downwind
% for h=1:H-1
%     E1((h-1)*N+(1:N),h*N+1) = zeros(N,1);
% end
% 1/2 Upwind, 1/2 central
a = 2;      % a=2 is central 
E1(1:N,1)     = E1(1:N,1)*(a-1)/a;
E1(N+(1:N),1) = E1(N+(1:N),1)/a;
for h=1:H-1
    E1((h-1)*N+(1:N),h*N+1) = E1((h-1)*N+(1:N),h*N+1)*(a-1)/a;
    E1(h*N+(1:N),h*N+1) = E1(h*N+(1:N),h*N+1)/a;
end

We = zeros(N*H+1,1);
for h=1:H
    ind2 = (h-1)*N+(1:N+1);
    We(ind2,1) = We(ind2,1)+w';
end
We(1) = We(1)+We(N*H+1);
We(N*H+1) = [];
W = spdiags(We,0,N*H,N*H);

A = W*E1'*D/dxdxi;


% M1 = v*0.5*(A-A');
M1 = v*A;
M2 = W;



nn = 100;
xixi = linspace(-1,1,nn+1);
hh = MimeticpolyVal(xixi,N,1);

xx = linspace(x(1),x(end),H*nn+1);

t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;

    %%%%%%%%%%%%%%%%%%%%
    % nodally exact, interpolated solution

    u_interp = zeros(1,N*H);
    for i=1:N*H
        s = x(i)-v*j*dt;
        while s<0 || s>L
            if s>1
                s = s-L;
            elseif s<0
                s = s+L;
            end
        end

        u_interp(i) = initialfunction(s,InitInfo);
    end

    u_interp = [u_interp u_interp(1)];
    for h=1:H
        ind = (h-1)*N+(1:N+1);
        uu_interp = u_interp(ind)*hh;
    
        plot(x(ind),u_interp(ind),'xr')
        hold on
        ind = (h-1)*nn+(1:nn+1);
        plot(xx(ind),uu_interp,'r')
    end
    %%%%%%%%%%%%%%%%%%%%%

    phinew = timemarching(M2,M1,dt,phi,'ESDIRK5');

%     phiold = phi;
    phi    = phinew;

%     phi_post = phi';
    phi_post = [phi' phi(1)];
    for h=1:H
    ind = (h-1)*N+(1:N+1);
    pphi = phi_post(ind)*hh;

    plot(x(ind),phi_post(ind),'o')
    ind = (h-1)*nn+(1:nn+1);
    plot(xx(ind),pphi)
    end
    ylim([-.2 1.2])
    pause(0.01)
    hold off

end