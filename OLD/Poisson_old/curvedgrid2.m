% clear all; close all; clc; N=6; c=0.2;

[xiGLL,wGLL] = GLLnodes(N);    etaGLL = xiGLL;
[xiG,wG]   = Gnodes(N);        etaG   = xiG;
xiEG     = [-1 xiG 1];         etaEG  = xiEG;

XiGLLGLL = xiGLL'*ones(1,N+1); EtaGLLGLL = XiGLLGLL';
XiGLLG  = xiGLL'*ones(1,N);    EtaGLLG   = ones(N+1,1)*etaG;
XiGGLL  = xiG'*ones(1,N+1);    EtaGGLL   = ones(N,1)*etaGLL;
XiEGEG  = xiEG'*ones(1,N+2);   EtaEGEG   = XiEGEG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XGLLGLL = XiGLLGLL + c*sin(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
YGLLGLL = EtaGLLGLL+ c*sin(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);

XGLLG = XiGLLG + c*sin(pi*XiGLLG).*sin(pi*EtaGLLG);
YGLLG = EtaGLLG+ c*sin(pi*XiGLLG).*sin(pi*EtaGLLG);

XGGLL = XiGGLL + c*sin(pi*XiGGLL).*sin(pi*EtaGGLL);
YGGLL = EtaGGLL+ c*sin(pi*XiGGLL).*sin(pi*EtaGGLL);

XEGEG = XiEGEG + c*sin(pi*XiEGEG).*sin(pi*EtaEGEG);
YEGEG = EtaEGEG+ c*sin(pi*XiEGEG).*sin(pi*EtaEGEG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dXdXiGLLGLL  = 1+pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
dXdEtaGLLGLL = pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);
dYdXiGLLGLL  = pi*c*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
dYdEtaGLLGLL = 1+pi*c*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);

dXdXiGLLG  = 1+pi*c*cos(pi*XiGLLG).*sin(pi*EtaGLLG);
dXdEtaGLLG = pi*c*sin(pi*XiGLLG).*cos(pi*EtaGLLG);
dYdXiGLLG  = pi*c*cos(pi*XiGLLG).*sin(pi*EtaGLLG);
dYdEtaGLLG = 1+pi*c*sin(pi*XiGLLG).*cos(pi*EtaGLLG);

dXdXiGGLL  = 1+pi*c*cos(pi*XiGGLL).*sin(pi*EtaGGLL);
dXdEtaGGLL = pi*c*sin(pi*XiGGLL).*cos(pi*EtaGGLL);
dYdXiGGLL  = pi*c*cos(pi*XiGGLL).*sin(pi*EtaGGLL);
dYdEtaGGLL = 1+pi*c*sin(pi*XiGGLL).*cos(pi*EtaGGLL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JGLLGLL = dXdXiGLLGLL.*dYdEtaGLLGLL-dXdEtaGLLGLL.*dYdXiGLLGLL;
JGLLG   = dXdXiGLLG.*dYdEtaGLLG-dXdEtaGLLG.*dYdXiGLLG;
JGGLL   = dXdXiGGLL.*dYdEtaGGLL-dXdEtaGGLL.*dYdXiGGLL;

S11GLLG = dXdXiGLLG.^2+dYdXiGLLG.^2;
S12GLLGLL = dXdXiGLLGLL.*dXdEtaGLLGLL+dYdXiGLLGLL.*dYdEtaGLLGLL;
S21GLLGLL = S12GLLGLL;
S22GGLL = dXdEtaGGLL.^2+dYdEtaGGLL.^2;

S11J = S11GLLG./JGLLG;
S12J = S12GLLGLL./JGLLGLL;
S21J = S12J;
S22J = S11J';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cxi  = ?????
% Ceta = ?????

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% l = linspace(-1,1,200);
% x_xi = zeros(N+1,200); y_xi = x_xi; x_eta = x_xi; y_eta = x_xi;
% for i=1:N+1
%     % xi const
%     x_xi(i,:)  = xiGLL(i)  + c*sin(pi*xiGLL(i))*sin(pi*l);
%     y_xi(i,:)  = l + c*sin(pi*xiGLL(i))*sin(pi*l);
%     % eta const
%     x_eta(i,:) = l + c*sin(pi*l)*sin(pi*etaGLL(i));
%     y_eta(i,:) = etaGLL(i) + c*sin(pi*l)*sin(pi*etaGLL(i));
% end
% 
% figure
% subplot(1,2,1)
% hold on
% % eta const
% for i=1:N
%     for j=1:N+1
%         plot([XGLLGLL(i,j) XGLLGLL(i+1,j)],[YGLLGLL(i,j) YGLLGLL(i+1,j)],'.-')
%     end
% end
% % xi const
% for i=1:N+1
%     for j=1:N
%         plot([XGLLGLL(i,j) XGLLGLL(i,j+1)],[YGLLGLL(i,j) YGLLGLL(i,j+1)],'.-')
%     end
% end
% axis('square')
% title(['N = ' num2str(N)])
% 
% subplot(1,2,2)
% hold on
% % eta const
% for i=1:N+1
%     plot(x_xi(i,:),y_xi(i,:))
%     plot(x_eta(i,:),y_eta(i,:))
% end
% axis('square')
% title(['N = ' num2str(N)])
