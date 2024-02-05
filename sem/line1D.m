clear all
close all
clc

N=4;
kk = 801;
xx=linspace(-1,1,kk);


% Gauss-Lobatto-Legendre zeros
xi = gausslobattolegendre(N);
% xi = linspace(-1,1,N+1);

dLa = zeros(N+1,kk);

for i=1:N+1
    xj     = xi;
    xj(i)  = [];
    ci(i)  = prod(xi(i)-xj);
end

La = zeros(N+1,kk);
for r=1:N+1
    xj    = xi;
    xj(r) = [];

    R=1:N+1;
    R(r) = [];
    for s=1:N
        S=R;
        S(s) = [];
        for k=1:kk
            La(r,k) = 1/ci(r)*prod(xx(k)-xj);
            dLa(r,k) = dLa(r,k)+1/ci(r)*prod(xx(k)-xi(S));
        end
    end
end


psi = zeros(N,kk);
for k=1:N
    for j=1:kk
        psi(k,j) = -sum(dLa(1:k,j));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(xx,La)
grid
legend(num2str((0:N)'))

figure
plot(xx,dLa)
grid
legend(num2str((0:N)'))

figure
plot(xx,psi)
legend(num2str((1:N)'))
hold on
for i=1:N+1
    plot([xi(i) xi(i)],[floor(min(min(psi))) ceil(max(max(psi)))],'--k')
end
for i=1:N
    plot([xi(i) xi(i+1)],[1/(xi(i+1)-xi(i)) 1/(xi(i+1)-xi(i))],'k')
end
plot([-1 1],[0 0],'k')
plot(xi,zeros(1,N+1),'sk')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
