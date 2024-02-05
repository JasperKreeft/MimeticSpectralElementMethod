%%% MESH

clear
close all
clc

N = 2;

% Gauss-Lobatto-Legendre zeros
P11 = JacobiPoly(N,1,1);
innerzeros = roots(P11(N,N:-1:1))';
gausslobattozeros = sort([-1 innerzeros 1]);


% Gauss zeros
P00 = JacobiPoly(N+1,0,0);
gausszeros = sort(roots(P00(N+1,N+1:-1:1)))';

% Bouwman zeros
bouwmanzeros = sort([-1 gausszeros 1]);

figure
hold on
% dual grid
for i=2:N+1
plot([bouwmanzeros(i) bouwmanzeros(i)],[-1 1],'--g')
end
for j=2:N+1
plot([-1 1],[bouwmanzeros(j) bouwmanzeros(j)],'--g')
end

% primal grid
for i=1:N+1
plot([gausslobattozeros(i) gausslobattozeros(i)],[-1 1])
end
for j=1:N+1
plot([-1 1],[gausslobattozeros(j) gausslobattozeros(j)])
end
axis equal
axis([-1.1 1.1 -1.1 1.1])

% for i=1:N+1
% for j=1:N+1
% plot(gausslobattozeros(i),gausslobattozeros(j),'o','LineWidth',2,'MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',6)
% end
% end

for i=1:N+2
for j=1:N+2
plot(bouwmanzeros(i),bouwmanzeros(j),'o','MarkerFaceColor','r','MarkerEdgeColor','r','MarkerSize',6)
end
end