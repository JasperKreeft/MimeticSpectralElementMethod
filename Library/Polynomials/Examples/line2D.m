clear all
% close all
clc

N=4;
kk = 300;
xx=linspace(-1,1,kk);
yy=linspace(-1,1,kk);

[XX,YY] = meshgrid(xx,yy);

% Gauss-Lobatto-Legendre zeros
% P11 = JacobiPoly(N,1,1);
% innerzeros = roots(P11(N,N:-1:1))';
% xi = sort([-1 innerzeros 1]);
xi = GLLnodes(N);
% xi=linspace(-1,1,N+1);
yi = xi;


for i=1:N+1
    xj     = xi;
    xj(i)  = [];
    ci(i)  = prod(xi(i)-xj);
end

La = zeros(N+1,kk);
for i=1:N+1
    xj    = xi;
    xj(i) = [];
    for k=1:kk
        La(i,k) = 1/ci(i)*prod(xx(k)-xj);
        if N==1
            dLa(i,k) = 1/ci(i);
        elseif N==2
            dLa(i,k) = 1/ci(i)*sum([2 -sum(xj)].*xx(k).^(N-1:-1:0));
        elseif N==3
            dLa(i,k) = 1/ci(i)*sum([3 -2*sum(xj) (xj(1)*xj(2)+xj(1)*xj(3)+xj(2)*xj(3))].*xx(k).^(N-1:-1:0));
        elseif N==4
            dLa(i,k) = 1/ci(i)*sum([4 -3*sum(xj) 2*(xj(1)*xj(2)+xj(1)*xj(3)+xj(1)*xj(4)+xj(2)*xj(3)+xj(2)*xj(4)+xj(3)*xj(4)) -(xj(1)*xj(2)*xj(3)+xj(1)*xj(2)*xj(4)+xj(1)*xj(3)*xj(4)+xj(2)*xj(3)*xj(4))].*xx(k).^(N-1:-1:0));
        end
    end
end

for k=1:N
    for j=1:kk
%         dpsi(k,j) = 2/N/(N+1)*(-(N+1-k)*sum(dLa(1:k,j))+k*sum(dLa(k+1:N+1,j)));
        dpsi(k,j) = -2/N*sum(dLa(1:k,j));
    end
end

figure
for i=1:N
    for j=1:N+1
        for k1=1:kk
            for k2=1:kk
                Z11(k1,k2) = La(j,k1)*dpsi(i,k2);
            end
        end
        surf(XX,YY,Z11')
%         pcolor(XX,YY,Z11)
        shading interp
        hold on
        plot3([xi(j) xi(j)],[yi(i) yi(i+1)],[2/N/(xi(i+1)-xi(i)) 2/N/(xi(i+1)-xi(i))]+1,'k','linewidth',4)
%         plot([xi(i) xi(i+1)],[yi(j) yi(j)],'k','linewidth',4)
        axis([-1 1 -1 1 -1 2.5])
        title(['Gradient in x-direction:   i=',num2str(i),' j=',num2str(j)])
        xlabel('x-axis')
        ylabel('y-axis')
        view([0 0 1])
        axis equal
        keyboard
        hold off
    end
end