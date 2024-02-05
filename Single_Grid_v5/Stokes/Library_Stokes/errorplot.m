
% figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(NrCellRange)>1
%     if exist('errorL2_interp','var')
%     semilogy(NrCellRange,errorL2_interp,'--or','markerface','r')
%     hold on
%     end
    semilogy(NrCellRange,L2_u,'-vc','markerface','c')%'-ob','markerface','b')
    hold on
%     semilogy(NrCellRange,L2_v,'-^g','markerface','g')
%     semilogy(NrCellRange,L2_w,'-xr','markerface','r')
%     semilogy(NrCellRange,L2_p,'-sc','markerface','c')
    xlim([NrCellRange(1) N])
    xlabel('N')
elseif length(HconvRange)>1
%     if exist('errorL2_interp','var')
%     loglog(2./(HconvRange),errorL2_interp,'--or','markerface','r')
%     hold on
%     end
    ax = 2./(HconvRange);
    loglog(ax,L2_u,'-ob','markerface','b')
    hold on
    loglog(ax,L2_v,'-^g','markerface','g')
    loglog(ax,L2_w,'-xr','markerface','r')
    loglog(ax,L2_p,'-sc','markerface','c')

a = (log(L2_u(2))-log(L2_u(4))) / (log(ax(2))-log(ax(4)));
loglog(ax(2:4),L2_u(2:4)/(1+N/10),'k')
text(ax(3),L2_u(3)/(1+(N+1)/10),['u=' num2str(a,2)],'color','k','fontsize',18)

a = (log(L2_w(2))-log(L2_w(4))) / (log(ax(2))-log(ax(4)));
loglog(ax(2:4),L2_w(2:4)*(1+N/10),'k')
text(ax(3),L2_w(3)*(1+(N+1)/10),['w=' num2str(a,2)],'color','k','fontsize',18)

a = (log(L2_p(2))-log(L2_p(4))) / (log(ax(2))-log(ax(4)));
loglog(ax(2:4),L2_p(2:4)/(1+N/10),'k')
text(ax(3),L2_p(3)/(1+(N+1)/10),['p=' num2str(a,2)],'color','k','fontsize',18)


    grid on
%     xlim([])
    xlabel('h')
%     set(gca,'xtick',roundn(get(gca,'xtick'),-2))
end
ylabel('L^2-error')
title('error Stokes')
legend('u','v','\omega','p',0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subplot(2,2,2)
% if length(NrCellRange)>1
%     if exist('errorL2_q_interp','var')
%     semilogy(NrCellRange,errorL2_q_interp,'--dm','markerface','m')
%     hold on
%     end
%     semilogy(NrCellRange,errorL2_q,'-sg','markerface','g')
%     hold on
%     xlim([NrCellRange(1) N])
%     xlabel('N')
% elseif length(HconvRange)>1
%     if exist('errorL2_q_interp','var')
%     loglog(2./(HconvRange),errorL2_interp,'--dm','markerface','m')
%     hold on
%     end
%     loglog(2./(HconvRange),errorL2_q,'-sg','markerface','g')
%     hold on
% %     xlim([])
%     xlabel('h')
% end
% ylabel('L^2-error')
% title('flux q')
% legend('interpolation','numerical',0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subplot(2,2,3)
% if length(NrCellRange)>1
%     loglog(NrCellRange,ConditionNumber,'-')
%     hold on
%     loglog(NrCellRange,CondRef,'-r')
%     xlabel('N')
%     legend('numerical',['ref: c=' c_str],0)
% elseif length(HconvRange)>1
%     loglog(2./HconvRange,ConditionNumber,'-')
%     hold on
%     loglog(2./HconvRange,CondRef,'-r')
%     xlabel('h')
%     legend('numerical',['ref: c=' c_str],0)
% end
% ylabel('Condition Number')
% grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subplot(2,2,4)
% if exist('D','var')
% if length(NrCellRange)>1
%     semilogy(NrCellRange,Linv_errorDiv,'o')
%     hold on
%     semilogy(NrCellRange,L1_errorDiv,'^')
%     H=legend('$||Dq-f||_{L^\infty(\Omega)}$','$||Dq-f||_{L^1(\Omega)}$',0);
%     semilogy([0 2*ceil(NrCellRange(end)/2)],[eps eps],'-g')
%     xlabel('N')
% elseif length(HconvRange)>1
%     loglog(2./HconvRange,Linv_errorDiv,'o')
%     hold on
%     loglog(2./HconvRange,L1_errorDiv,'^')
%     H=legend('$||Dq-f||_{L^\infty(\Omega)}$','$||Dq-f||_{L^1(\Omega)}$',0);
%     loglog(xlim,[eps eps],'-g')
%     xlabel('h')
% end
% ylabel('error')
% if strcmp(DomInfo,'SinDeformGrid'); title(['c = ',num2str(DomInfo)]); end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% disp(['errorL2          = ' num2str(errorL2,4)])
% 
% if exist('phi_interp','var')
% disp(['errorL2_interp   = ' num2str(errorL2_interp,4)])
% end
% 
% disp(['errorL2_q        = ' num2str(errorL2_q,4)])
% 
% if exist('qx_interp','var')
% disp(['errorL2_q_interp = ' num2str(errorL2_interp,4)])
% end
% 
% if exist('ConditionNumber','var')
% disp(['ConditionNumber  = ' num2str(ConditionNumber,4)])
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(H,'Interpreter','Latex','fontSize',12)