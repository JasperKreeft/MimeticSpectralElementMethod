% clear all; close all; clc; N=6; c=0.2; gridtype = 'sinecurve';

[xiGLL,wGLL] = GLLnodes(N);    etaGLL = xiGLL;  % Gauss-Lobotto-Legendre
[xiG,wG]   = Gnodes(N);        etaG   = xiG;    % Gauss
xiEG     = [-1 xiG 1];         etaEG  = xiEG;   % Extended Gauss

XiGLLGLL = xiGLL'*ones(1,N+1); EtaGLLGLL = XiGLLGLL';
XiGLLG  = xiGLL'*ones(1,N);    EtaGLLG   = ones(N+1,1)*etaG;
XiGGLL  = xiG'*ones(1,N+1);    EtaGGLL   = ones(N,1)*etaGLL;
XiEGEG  = xiEG'*ones(1,N+2);   EtaEGEG   = XiEGEG';
XiGG    = xiG'*ones(1,N);      EtaGG     = XiGG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

XGLLGLL = XiGLLGLL + cc*sin(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
YGLLGLL = EtaGLLGLL+ cc*sin(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);

XGLLG = XiGLLG + cc*sin(pi*XiGLLG).*sin(pi*EtaGLLG);
YGLLG = EtaGLLG+ cc*sin(pi*XiGLLG).*sin(pi*EtaGLLG);

XGGLL = XiGGLL + cc*sin(pi*XiGGLL).*sin(pi*EtaGGLL);
YGGLL = EtaGGLL+ cc*sin(pi*XiGGLL).*sin(pi*EtaGGLL);

XEGEG = XiEGEG + cc*sin(pi*XiEGEG).*sin(pi*EtaEGEG);
YEGEG = EtaEGEG+ cc*sin(pi*XiEGEG).*sin(pi*EtaEGEG);

XGG = XiGG + cc*sin(pi*XiGG).*sin(pi*EtaGG);
YGG = EtaGG+ cc*sin(pi*XiGG).*sin(pi*EtaGG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dXdXiGLLGLL  = 1+pi*cc*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
dXdEtaGLLGLL = pi*cc*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);
dYdXiGLLGLL  = pi*cc*cos(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
dYdEtaGLLGLL = 1+pi*cc*sin(pi*XiGLLGLL).*cos(pi*EtaGLLGLL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JGLLGLL = dXdXiGLLGLL.*dYdEtaGLLGLL-dXdEtaGLLGLL.*dYdXiGLLGLL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = reshape(JGLLGLL,1,(N+1)^2);

Jacobian = spdiags(kron(j,[1 1])',0,2*(N+1)^2,2*(N+1)^2);

qinv11 = kron(reshape(( dXdXiGLLGLL./JGLLGLL),1,(N+1)^2),[1 0])';
qinv22 = kron(reshape((dYdEtaGLLGLL./JGLLGLL),1,(N+1)^2),[0 1])';
qinv12 = kron(reshape((dXdEtaGLLGLL./JGLLGLL),1,(N+1)^2),[0 1])';
qinv21 = kron(reshape(( dYdXiGLLGLL./JGLLGLL),1,(N+1)^2),[1 0])';

Qinv = spdiags([qinv21 qinv11+qinv22 qinv12],-1:1,2*(N+1)^2,2*(N+1)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
