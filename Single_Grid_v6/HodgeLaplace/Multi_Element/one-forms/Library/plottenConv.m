figure
plot(E(1:10),'o','markerface','b')
grid on
set(gca,'xtick',1:10,'ytick',0:10)
axis equal
axis([0 10 0 10])
xlabel('number of eigenvalue')
ylabel('value of eigenvalue')
title(['Result for N=' num2str(NrCellRange(end)) ', c=' num2str(cc)])

if length(NrCellRange)>=9
figure
handle(9) = semilogy(NrCellRange,error(9,:)',':dk','markerface','k');
hold on
handle(8) = semilogy(NrCellRange,error(8,:)','--sk','markerface','k');
handle(7) = semilogy(NrCellRange,error(7,:)','-ok','markerface','k');
handle(6) = semilogy(NrCellRange,error(6,:)','-oy','markerface','y');
handle(5) = semilogy(NrCellRange,error(5,:)','-oc','markerface','c');
handle(4) = semilogy(NrCellRange,error(4,:)','-om','markerface','m');
handle(3) = semilogy(NrCellRange,error(3,:)','-or','markerface','r');
handle(2) = semilogy(NrCellRange,error(2,:)','-og','markerface','g');
handle(1) = semilogy(NrCellRange,error(1,:)','-ob','markerface','b');
grid on
legend(handle,'1','2','3','4','5','6','7','8','9',4,'orientation','horizontal')
axis([0 N 1e-10 1e1])
xlabel('N')
ylabel('error eigenvalues')
title(['Convergence of first nine non-zero eigenvalues for c=' num2str(cc)])
end

nrH = length(HconvRange);
if nrH > 1
figure
if cc==0
handle(10) =loglog(pi./HconvRange,error(10,1:nrH),'-ok','markerface','k');
hold on
handle(9) = loglog(pi./HconvRange,error(9,1:nrH),'-ok','markerface','k');
handle(8) = loglog(pi./HconvRange,error(8,1:nrH),'-oc','markerface','c');
handle(7) = loglog(pi./HconvRange,error(7,1:nrH),'-om','markerface','m');
handle(6) = loglog(pi./HconvRange,error(6,1:nrH),'-om','markerface','m');
handle(5) = loglog(pi./HconvRange,error(5,1:nrH),'-og','markerface','g');
handle(4) = loglog(pi./HconvRange,error(4,1:nrH),'-og','markerface','g');
handle(3) = loglog(pi./HconvRange,error(3,1:nrH),'-or','markerface','r');
handle(2) = loglog(pi./HconvRange,error(2,1:nrH),'-ob','markerface','b');
handle(1) = loglog(pi./HconvRange,error(1,1:nrH),'-ob','markerface','b');
% grid on
legend(handle,'1','2','3','4','5','6','7','8','9','10',4)%,'orientation','horizontal')
ylim([1e-4 10])
xlabel('h')
ylabel('error eigenvalues')
title(['h-Convergence of first ten non-zero eigenvalues for N=' num2str(N) ' and c=' num2str(cc)])
set(gca,'ytick',10.^(-4:1))
else
handle(9) = loglog(pi./HconvRange,error(9,1:nrH),':dk','markerface','k');
hold on
handle(8) = loglog(pi./HconvRange,error(8,1:nrH),'--sk','markerface','k');
handle(7) = loglog(pi./HconvRange,error(7,1:nrH),'-ok','markerface','k');
handle(6) = loglog(pi./HconvRange,error(6,1:nrH),'-oy','markerface','y');
handle(5) = loglog(pi./HconvRange,error(5,1:nrH),'-oc','markerface','c');
handle(4) = loglog(pi./HconvRange,error(4,1:nrH),'-om','markerface','m');
handle(3) = loglog(pi./HconvRange,error(3,1:nrH),'-or','markerface','r');
handle(2) = loglog(pi./HconvRange,error(2,1:nrH),'-og','markerface','g');
handle(1) = loglog(pi./HconvRange,error(1,1:nrH),'-ob','markerface','b');
grid on
legend(handle,'1','2','3','4','5','6','7','8','9',4,'orientation','horizontal')
% axis([0 N 1e-10 1e2])
xlabel('h')
ylabel('error eigenvalues')
title(['h-Convergence of first nine non-zero eigenvalues for N=' num2str(N) ' and c=' num2str(cc)])
end
end

figure
plot(EE,'.')
grid on