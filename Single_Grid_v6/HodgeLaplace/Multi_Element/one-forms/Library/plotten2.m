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
axis([0 N 1e-10 1e2])
xlabel('N')
ylabel('error eigenvalues')
title(['Convergence of first nine non-zero eigenvalues for c=' num2str(cc)])