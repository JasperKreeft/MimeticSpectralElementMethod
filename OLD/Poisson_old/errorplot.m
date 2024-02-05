% figure
% semilogy(NrCellRange,errorL1,color(m))
% hold on
% xlim([NrCellRange(1) z])
% xlabel('N')
% ylabel('L^1-error')

figure
semilogy(NrCellRange,errorL2(NrCellRange))%,color(m))
hold on
semilogy(NrCellRange,errorL2_interp(NrCellRange),'--r')%['--' color(m)])
xlim([NrCellRange(1) N])
xlabel('N')
ylabel('L^2-error')
title('potential \phi')

legend('numerical','interpolation')

figure
semilogy(NrCellRange,errorL2_interp_u(NrCellRange),'--r')
hold on
semilogy(NrCellRange,errorL2_uv(NrCellRange),'b')
xlim([NrCellRange(1) N])
xlabel('N')
ylabel('L^2-error')
title('velocity u')