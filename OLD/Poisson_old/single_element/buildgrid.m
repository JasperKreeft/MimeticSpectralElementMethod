% clear all; close all; clc; N=6; c=0.2;

[xiGLL,wGLL] = GLLnodes(N);    etaGLL = xiGLL;  % Gauss-Lobotto-Legendre
[xiG,wG]   = Gnodes(N);        etaG   = xiG;    % Gauss
xiEG     = [-1 xiG 1];         etaEG  = xiEG;   % Extended Gauss

XiGLLGLL = xiGLL'*ones(1,N+1); EtaGLLGLL = XiGLLGLL';
XiGLLG  = xiGLL'*ones(1,N);    EtaGLLG   = ones(N+1,1)*etaG;
XiGGLL  = xiG'*ones(1,N+1);    EtaGGLL   = ones(N,1)*etaGLL;
XiEGEG  = xiEG'*ones(1,N+2);   EtaEGEG   = XiEGEG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch gridtype
    case 'standard'
        XGLLGLL = XiGLLGLL;
        YGLLGLL = EtaGLLGLL;

        XGLLG = XiGLLG;
        YGLLG = EtaGLLG;

        XGGLL = XiGGLL;
        YGGLL = EtaGGLL;

        XEGEG = XiEGEG;
        YEGEG = EtaEGEG;

    case 'sinecurve'
        XGLLGLL = XiGLLGLL + c*sin(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);
        YGLLGLL = EtaGLLGLL+ c*sin(pi*XiGLLGLL).*sin(pi*EtaGLLGLL);

        XGLLG = XiGLLG + c*sin(pi*XiGLLG).*sin(pi*EtaGLLG);
        YGLLG = EtaGLLG+ c*sin(pi*XiGLLG).*sin(pi*EtaGLLG);

        XGGLL = XiGGLL + c*sin(pi*XiGGLL).*sin(pi*EtaGGLL);
        YGGLL = EtaGGLL+ c*sin(pi*XiGGLL).*sin(pi*EtaGGLL);

        XEGEG = XiEGEG + c*sin(pi*XiEGEG).*sin(pi*EtaEGEG);
        YEGEG = EtaEGEG+ c*sin(pi*XiEGEG).*sin(pi*EtaEGEG);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
