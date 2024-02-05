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

D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);

M0 = spdiags(w',0,N+1,N+1);

M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(w.*e(i,:).*e(j,:));
    end
end

P = ceil(3/2*N);
[xw,ww] = GLLnodes(P);
[Hw,Ew]  = MimeticpolyVal(xw,N,1);

t=0; j=0;
while t<T

    t = t+dt;
    j = j+1;

    phiP = phi;
    phiPold = zeros(N,1);
    kk=0; error=1;
    % Picard Iterations
    while error>1e-6
        kk=kk+1;

        phiPold = phiP;

        LieMatrix = zeros(N+1,N);
        for k=1:N+1
            for i=1:N
                for p=1:P+1
                    LieMatrix(k,i) = LieMatrix(k,i) + ww(p)*Ew(i,p)*Hw(k,p)*sum(phi.*Ew(:,p));
                end
            end
        end

        Matrix = [ M1/dt M1*D ; LieMatrix -M0 ];
        RHS = [ M1/dt*phi ; zeros(N+1,1) ];

        Phi_Q = Matrix\RHS;

        phiP = Phi_Q(1:N);

        error = max(abs(phiP-phiPold));


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