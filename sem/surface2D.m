clear all
close all
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
% xi=linspace(-1,1,N+1);
xi = GLLnodes(N);
yi = xi;

for i=1:N+1
    xj    = xi;
    xj(i) = [];
    ci(i)  = prod(xi(i)-xj);
    for k=1:kk
%       N=3;
%         dLa(i,k) = 1/ci(i)*sum([3 -2*sum(xj) (xj(1)*xj(2)+xj(1)*xj(3)+xj(2)*xj(3))].*xx(k).^(N-1:-1:0));
%       N=4
        dLa(i,k) = 1/ci(i)*sum([4 -3*sum(xj) 2*(xj(1)*xj(2)+xj(1)*xj(3)+xj(1)*xj(4)+xj(2)*xj(3)+xj(2)*xj(4)+xj(3)*xj(4)) -(xj(1)*xj(2)*xj(3)+xj(1)*xj(2)*xj(4)+xj(1)*xj(3)*xj(4)+xj(2)*xj(3)*xj(4))].*xx(k).^(N-1:-1:0));
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
    for j=1:N
        for k1=1:kk
            for k2=1:kk
                Z11(k1,k2) = dpsi(j,k1)*dpsi(i,k2);
            end
        end
        surf(XX,YY,Z11)
%         pcolor(XX,YY,Z11)
        shading interp
        hold on
        plot3([xi(i) xi(i+1) xi(i+1) xi(i) xi(i)],[yi(j) yi(j) yi(j+1) yi(j+1) yi(j)],[2/N/(xi(i+1)-xi(i))/N/(yi(j+1)-yi(j)) 2/N/(xi(i+1)-xi(i))/N/(yi(j+1)-yi(j)) 2/N/(xi(i+1)-xi(i))/N/(yi(j+1)-yi(j)) 2/N/(xi(i+1)-xi(i))/N/(yi(j+1)-yi(j)) 2/N/(xi(i+1)-xi(i))/N/(yi(j+1)-yi(j))]+5,'k','linewidth',4)
%         plot([xi(i) xi(i+1)],[yi(j) yi(j)],'k','linewidth',4)
%         axis([-1 1 -1 1 -1 2.5])
        title(['Surface in x-direction:   i=',num2str(i),' j=',num2str(j)])
        xlabel('x-axis')
        ylabel('y-axis')
        N^2*sum(sum(Z11))/kk^2
        axis equal
        view([0 0 1])
        keyboard
        hold off

    end
end