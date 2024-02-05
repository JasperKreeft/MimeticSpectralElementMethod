figure(10)
subplot(2,2,1)
semilogy(NrCellRange,errorL2(NrCellRange),'-^b')%,color(m))
hold on
semilogy(NrCellRange,errorL2_interp(NrCellRange),'--or')%['--' color(m)])
xlim([NrCellRange(1) N])
xlabel('N')
ylabel('L^2-error')
title('potential \phi')

legend('numerical','interpolation')

figure(10)
subplot(2,2,2)
if c==0.0
    semilogy(NrCellRange,errorL2_interp_u(NrCellRange),'--m')
    hold on
end
semilogy(NrCellRange,errorL2_uv(NrCellRange),'-sg')
hold on
xlim([NrCellRange(1) N])
xlabel('N')
ylabel('L^2-error')
title('velocity u')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)
subplot(2,2,3)
semilogy(NrCellRange,ConditionNumber(NrCellRange),'-')
xlabel('N')
ylabel('Condition Number')
grid
xlim([0 2*ceil((N+1)/2)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(10)
subplot(2,2,4)
semilogy(NrCellRange,Linv_errorDiv(NrCellRange),'o')
hold on
semilogy(NrCellRange,L1_errorDiv(NrCellRange),'^')
H=legend('$||Dq-f||_{L^\infty(\Omega)}$','$||Dq-f||_{L^1(\Omega)}$',2);
semilogy([0 2*ceil(NrCellRange(end)/2)],[eps eps],'-g')
xlabel('N')
title(['c = ',num2str(c)])
ylabel('error')


% set(H,'Interpreter','Latex','fontSize',12)