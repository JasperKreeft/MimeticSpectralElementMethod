clear all
close all
clc

if ispc; figure('windowstyle','docked'); else figure; end

v = 2;
dt = .02;

N = 30;
T = 30;

InitInfo.shape = 'step'; %'cosine'; %
InitInfo.X0    = [-0.5 0.5]; %0; %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = initialfunction(xi(1:N),InitInfo);

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


M = v*0.5*(A-A');
% M = v*A;

nn = 1000;
xx = linspace(-1,1,nn);
hh = MimeticpolyVal(xx,N,1);

t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;
    
    %%%%%%%%%%%%%%%%%%%%
    % nodally exact, interpolated solution

    u_interp = zeros(1,N);
    for i=1:N
        s = xi(i)-v*j*dt;
        while s<-1 || s>1
            if s>1
                s = s-2;
            elseif s<0
                s = s+2;
            end
        end

        u_interp(i) = initialfunction(s,InitInfo);
    end
    
    uu_interp = [u_interp u_interp(1)]*hh;
    
    plot(xi,[u_interp u_interp(1)],'xr')
    hold on
    plot(xx,uu_interp,'r')
    %%%%%%%%%%%%%%%%%%%%%


    phinew = timemarching(W,M,dt,phi,'ESDIRK3');

%     phiold = phi;
    phi    = phinew;
    
    pphi = [phi' phi(1)]*hh;

    plot(xi,[phi ; phi(1)],'o')
    plot(xx,pphi)
    ylim([-.2 1.2])
    pause(0.05)
    hold off

end