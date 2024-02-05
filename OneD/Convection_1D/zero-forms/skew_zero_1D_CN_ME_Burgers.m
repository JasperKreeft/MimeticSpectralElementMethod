% JACOBIAN'S AND OTHER TRANSFORMATION MATRICES MUST BE INCLUDED

clear all
close all
clc

if ispc; figure('windowstyle','docked'); else figure; end

N = 1;
H = 52;
T = 0.9;
dt = .01;

InitInfo.shape = 'SineBurgers';
InitInfo.X0    = [0.5 2.5]; %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

dxdxi = 2/(2*H);

x = [];
for h=1:H
    xe = (2*h-1)*dxdxi+xi*dxdxi-1;
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

% Central
for h=1:H
    E1(:,(h-1)*N+1) = E1(:,(h-1)*N+1)/2;
end

We = zeros(N*H+1,1);
for h=1:H
    ind2 = (h-1)*N+(1:N+1);
    We(ind2,1) = We(ind2,1)+w';
end
We(1) = We(1)+We(N*H+1);
We(N*H+1) = [];
W = spdiags(We,0,N*H,N*H);


M2 = W;



nn = 1000;
xixi = linspace(-1,1,nn+1);
hh = MimeticpolyVal(xixi,N,1);

xx = linspace(x(1),x(end),H*nn+1);

t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;


    phiP = phi;
    phiPold = zeros(N*H,1);
    k=0;
    while max(abs(phiP-phiPold))>1e-4
        k=k+1;

        phiPold = phiP;

        PHI = spdiags(phiP,0,N*H,N*H);

        A = W*PHI*E1'*D/dxdxi;

        M1 = 0.5*(A-A');

        phiP = timemarching(M2,M1,dt,phi,'ESDIRK5');

    end
    phi = phiP;
    
    phi_post = [phi' phi(1)];
    for h=1:H
    ind = (h-1)*N+(1:N+1);
    pphi = phi_post(ind)*hh;

    plot(x(ind),phi_post(ind),'o')
    hold on
    ind = (h-1)*nn+(1:nn+1);
    plot(xx(ind),pphi)
    end
    ylim([-0.5 0.5])
    pause(0.05)
    hold off

end