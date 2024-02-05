clear all
close all
clc

if ispc; figure('windowstyle','docked'); else figure; end

dt = .02;

N = 5;
T = 0.8;

InitInfo.shape = 'SineBurgers';
InitInfo.X0    = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = initialfunction(xi,InitInfo);

plot(xi,phi,'-o')
grid on
ylim([-0.5 0.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);

W = spdiags(w',0,N+1,N+1);

P = ceil(3/2*N);
[xw,ww] = GLLnodes(P);
[Hw,Ew]  = MimeticpolyVal(xw,N,1);

LieMatrix = zeros(N+1,N);
for k=1:N+1
    for i=1:N
        for p=1:P+1
            LieMatrix(k,i) = LieMatrix(k,i) + ww(p)*Ew(i,p)*Hw(k,p)*sum(phi.*Ew(:,p))
            
            


%%%

I = ones(N+1);

nn = 1000;
xx = linspace(-1,1,nn);
hh = MimeticpolyVal(xx,N,1);

t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;

    phiP = phi;
    phiPold = zeros(N+1,1);
    k=0;
    while max(abs(phiP-phiPold))>1e-4
        k=k+1;

        phiPold = phiP;

        PHI = spdiags(phiP,0,N+1,N+1);

        A = M0*((Hk'*PHI*I).*(Ek'*D));

        M = 0.5*(A-A');
%         M = A;

%         phiP = timemarching(W,M,dt,phi,'ESDIRK5');

    end
    phi = phiP;

    pphi = phi'*hh;

    plot(xi,phi,'o')
    hold on
    plot(xx,pphi)
    grid on
    ylim([-0.5 0.5])
    pause(0.05)
    hold off

end