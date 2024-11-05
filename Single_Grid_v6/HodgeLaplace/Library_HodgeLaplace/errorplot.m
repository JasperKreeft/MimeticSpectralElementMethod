
figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,1)
if length(NrCellRange)>1
    if exist('errorL2_interp','var')
        Marker(semilogy(NrCellRange,errorL2_interp,'-o'));
        hold on
    end
    Marker(semilogy(NrCellRange,errorL2,'-o'));
    hold on
    Marker(semilogy(NrCellRange,abs(errorL2-errorL2_interp),'-o'));
    xlim([NrCellRange(1) N])
    xlabel('N')
elseif length(HconvRange)>1
    if exist('errorL2_interp','var')
        Marker(loglog(2./(HconvRange),errorL2_interp,'-o'));
        hold on
    end
    Marker(loglog(2./(HconvRange),errorL2,'-o'));
    hold on
%     xlim([])
    xlabel('h')
end
grid on
ylabel('L^2-error')
title('potential \phi')
legend('interpolation','numerical','|num-int|','location','northeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,2)
if length(NrCellRange)>1
    if exist('errorL2_q_interp','var')
        Marker(semilogy(NrCellRange,errorL2_q_interp,'-o'));
        hold on
    end
    Marker(semilogy(NrCellRange,errorL2_q,'-o'));
    hold on
    xlim([NrCellRange(1) N])
    xlabel('N')
elseif length(HconvRange)>1
    if exist('errorL2_q_interp','var')
        Marker(loglog(2./(HconvRange),errorL2_q_interp,'-o'));
        hold on
    end
    Marker(loglog(2./(HconvRange),errorL2_q,'-o'));
    hold on
%     xlim([])
    xlabel('h')
end
grid on
ylabel('L^2-error')
title('flux q')
legend('interpolation','numerical','location','northeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,3)
if length(NrCellRange)>1
    loglog(NrCellRange,ConditionNumber,'-')
    hold on
    loglog(NrCellRange,CondRef,'-r')
    xlabel('N')
%     legend('numerical',['ref: c=' c_str],'location','northeast')
elseif length(HconvRange)>1
    loglog(2./HconvRange,ConditionNumber,'-')
    hold on
    loglog(2./HconvRange,CondRef,'-r')
    xlabel('h')
%     legend('numerical',['ref: c=' c_str],'location','northeast')
end
ylabel('Condition Number')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,2,4)
if exist('D','var')
    if length(NrCellRange)>1
        semilogy(NrCellRange,Linv_errorDiv,'o')
        hold on
        semilogy(NrCellRange,L1_errorDiv,'^')
        H=legend('$||Dq-f||_{L^\infty(\Omega)}$','$||Dq-f||_{L^1(\Omega)}$','location','northeast');
        semilogy([0 2*ceil(NrCellRange(end)/2)],[eps eps],'-g')
        xlabel('N')
    elseif length(HconvRange)>1
        loglog(2./HconvRange,Linv_errorDiv,'o')
        hold on
        loglog(2./HconvRange,L1_errorDiv,'^')
        H=legend('$||Dq-f||_{L^\infty(\Omega)}$','$||Dq-f||_{L^1(\Omega)}$','location','northeast');
        loglog(xlim,[eps eps],'-g')
        xlabel('h')
    end
    ylabel('error')
    if strcmp(DomInfo,'SinDeformGrid')
        title(['c = ',num2str(DomInfo)]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['errorL2          = ' num2str(errorL2,4)])

if exist('phi_interp','var')
disp(['errorL2_interp   = ' num2str(errorL2_interp,4)])
end

disp(['errorL2_q        = ' num2str(errorL2_q,4)])

if exist('qx_interp','var')
disp(['errorL2_q_interp = ' num2str(errorL2_q_interp,4)])
end

if exist('ConditionNumber','var')
disp(['ConditionNumber  = ' num2str(ConditionNumber,4)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(H,'Interpreter','Latex','fontSize',12)