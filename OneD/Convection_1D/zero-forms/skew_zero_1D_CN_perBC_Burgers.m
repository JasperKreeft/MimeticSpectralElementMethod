clear all
close all
clc

% path(path,'C:\Users\Jasper\Documents\MATLAB/MSEM_codes/OneD/Convection_1D/Library')

if ispc; figure('windowstyle','docked'); else figure; end

dt = .02;

N = 40;
T = 0.8;

InitInfo.shape = 'SineBurgers';
InitInfo.X0    = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = initialfunction(xi(1:N),InitInfo);

plot(xi,[phi ; phi(1)],'-o')
grid on
ylim([-0.5 0.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N); D(N,1) = 1;

E1 = e;
E1(:,1) = E1(:,1)+E1(:,N+1);
E1(:,N+1) = [];

% Central
E1(:,1) = E1(:,1)/2;

We = w';
We(1) = We(1)+We(N+1);
We(N+1) = [];
W = spdiags(We,0,N,N);



%%%


nn = 1000;
xx = linspace(-1,1,nn);
hh = MimeticpolyVal(xx,N,1);

t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;

    phiP = phi;
    phiPold = zeros(N,1);
    k=0;
    while max(abs(phiP-phiPold))>1e-4
        k=k+1;

        phiPold = phiP;

        PHI = spdiags(phiP,0,N,N);

        A = W*PHI*E1'*D;

        M = 0.5*(A-A');
%         M = A;

        phiP = timemarching(W,M,dt,phi,'ESDIRK5');

    end
    phi = phiP;

    pphi = [phi' phi(1)]*hh;

    plot(xi,[phi ; phi(1)],'o')
    hold on
    plot(xx,pphi)
    grid on
    ylim([-0.5 0.5])
    pause(0.05)
    hold off

end

% rmpath('C:\Users\Jasper\Documents\MATLAB/MSEM_codes/OneD/Convection_1D/Library')