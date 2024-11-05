function plot_eigenvalues(E)

global N

figure
Idx = min(length(E),10);
ymax = ceil(E(Idx)/10)*10;
Marker(plot(E(1:Idx),'o'));
grid on
set(gca,'xtick',1:10,'ytick',0:ymax)
axis equal
axis([0 10 0 ymax])
xlabel('number of eigenvalue')
ylabel('value of eigenvalue')
title("Result for N="+N)



% figure(1)
% plot(E(1:10),'o','markerface','b')
% grid on
% set(gca,'xtick',1:10,'ytick',0:2:20)
% % axis equal
% axis([0 10 0 20])
% xlabel('number of eigenvalue')
% ylabel('value of eigenvalue')
% title(['Result for N=' num2str(NrCellRange(end))])