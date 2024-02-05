figure
semilogy(NrCellRange,L2_u(NrCellRange),'-o')
hold on
semilogy(NrCellRange,L2_v(NrCellRange),'-^g')
semilogy(NrCellRange,L2_w(NrCellRange),'-xr')
semilogy(NrCellRange,L2_p(NrCellRange),'-sc')
legend('u','v','\omega','p')