% With intermediate higher-order interpolation

clear all
close all
clc

if ispc; figure('windowstyle','docked'); else figure; end

dt = .02;

N = 20;
T = 0.7;

InitInfo.shape = 'SineBurgers';
InitInfo.X0    = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

nn = 1000;
xx = linspace(-1,1,nn);
[hh,ee] = MimeticpolyVal(xx,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial condition

phi = zeros(N,1);
for i=1:N
    phi(i,1) = 1/(4*pi)*(cos(pi*(xi(i+1)))-cos(pi*xi(i)));
end

phidx = phi./diff(xi)';
plot([xi(2:N+1) ; xi(1:N)],[phidx phidx]','g')
hold on
grid on

pphi = phi'*ee;
plot(xx,pphi,'r')
ylim([-0.5 0.5])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nd = 2*N-2;
D = spdiags([-ones(Nd,1) ones(Nd,1)],[0 1],Nd,Nd+1);

xd = GLLnodes(Nd);
[dummy,Ed] = MimeticpolyVal(xd,N,1);

P = ceil(3/2*(N-1));
[xw,ww] = GLLnodes(P);
[dummy,Ehw] = MimeticpolyVal(xw,Nd,1);
[dummy,Ew]  = MimeticpolyVal(xw,N ,1);

M1 = zeros(N,Nd);
for i=1:N
    for j=1:Nd
        M1(i,j) = sum(ww.*Ew(i,:).*Ehw(j,:));
    end
end

K2 = zeros(N);
for i=1:N
    for j=1:N
        K2(i,j) = sum(w.*e(i,:).*e(j,:));
    end
end

I = ones(N);


t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;

    phiP = phi;
    phiPold = zeros(N,1);
    kk=0;
    % Picard Iterations
    while max(abs(phiP-phiPold))>1e-6
        kk=kk+1;

        phiPold = phiP;

        PHI = spdiags(phiP,0,N,N); % NOTE: PHI*I = repmat(phiP,1,N);

        A = M1*D*((Ed'*PHI*I).*Ed');

%         K1 = 0.5*(A-A');
        K1 = A;

        phiP = timemarching(K2,K1,dt,phi,'CN');


    end
    phi = phiP;

    pphi = phi'*ee;

    phidx = phi./diff(xi)';
    plot([xi(2:N+1) ; xi(1:N)],[phidx phidx]','g')
    hold on
    plot(xx,pphi)
    ylim([-0.5 0.5])
    pause(0.2)
    hold off

end