figure(10)
subplot(2,2,1)
semilogy(NrCellRange,errorL2(NrCellRange),'-^b','markerface','b')
hold on
% semilogy(NrCellRange,errorL2_interp(NrCellRange),'--or','markerface','r')
xlim([NrCellRange(1) N])
xlabel('N')
ylabel('L^2-error')
title('potential \phi')

legend('numerical','interpolation')

subplot(2,2,2)
% if c==0.0
%     semilogy(NrCellRange,errorL2_interp_q(NrCellRange),'--m')
%     hold on
% end
% semilogy(NrCellRange,errorL2_q(NrCellRange),'-sg')
% hold on
% xlim([NrCellRange(1) N])
% xlabel('N')
% ylabel('L^2-error')
% title('velocity u')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3)
loglog(NrCellRange,ConditionNumber(NrCellRange),'-')
hold on
loglog(NrCellRange,CondRef(NrCellRange),'-r')
xlabel('N')
ylabel('Condition Number')
grid on
legend('numerical',['ref: c=' c_str],4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,4)
if exist('D','var')
semilogy(NrCellRange,Linv_errorDiv(NrCellRange),'o')
hold on
semilogy(NrCellRange,L1_errorDiv(NrCellRange),'^')
H=legend('$||Dq-f||_{L^\infty(\Omega)}$','$||Dq-f||_{L^1(\Omega)}$',2);
semilogy([0 2*ceil(NrCellRange(end)/2)],[eps eps],'-g')
xlabel('N')
if strcmp(DomInfo,'SinDeformGrid'); title(['c = ',num2str(DomInfo)]); end
ylabel('error')
end

% set(H,'Interpreter','Latex','fontSize',12)