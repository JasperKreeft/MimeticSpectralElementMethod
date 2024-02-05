function [EGzeros,EGweights] = extendedgauss(N)
% zeros and weights for extended Gauss quadrature
% 
% Input N: nr of cells

[Gzeros,Gweights] = gauss(N);

EGzeros = sort([-1 Gzeros 1]);
EGweights = [0 Gweights 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% La = LagrangePoly(EGzeros,2);
% kk = 1000;
% x = linspace(min(EGzeros),max(EGzeros),kk);
% f = zeros(N,kk);
% for i=1:N+2
%     xj    = EGzeros;
%     xj(i) = [];
%     for k=1:kk
%         f(i,k) = La(i,1)*prod(x(k)-xj);
%     end
% end
% 
% figure
% plot(x,f)
% grid
% xlim([-1 1])
% legend(num2str((1:N+2)'))
% hold on
% plot(EGzeros,zeros(N+2),'o','MarkerFaceColor','k','MarkerEdgeColor','k')