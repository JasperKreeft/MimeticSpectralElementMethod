% clear all; close all; clc; N=6; cc=0.2; numColumns = 4; numRows = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% standard element mesh and weights                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xiGLL,wGLL] = GLLnodes(N);    etaGLL = xiGLL;  % Gauss-Lobotto-Legendre
[xiG,wG]   = Gnodes(N);        etaG   = xiG;    % Gauss
xiEG     = [-1 xiG 1];         etaEG  = xiEG;   % Extended Gauss

XiGLLGLL = xiGLL'*ones(1,N+1); EtaGLLGLL = XiGLLGLL';
XiGLLG  = xiGLL'*ones(1,N);    EtaGLLG   = ones(N+1,1)*etaG;
XiGGLL  = xiG'*ones(1,N+1);    EtaGGLL   = ones(N,1)*etaGLL;
XiEGEG  = xiEG'*ones(1,N+2);   EtaEGEG   = XiEGEG';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]
% x = linspace(-1,1,numColumns+1);
% y = linspace(-1,1,numRows+1);

xibLR  = linspace(-1,1,numColumns+1);
etabAB = linspace(-1,1,numRows+1)   ;

% Global mesh
XGLLGLL = zeros(N*numColumns+1,N*numRows+1);
YGLLGLL = zeros(N*numColumns+1,N*numRows+1);
XGG   = zeros(N*numColumns  ,N*numRows  );
YGG   = zeros(N*numColumns  ,N*numRows  );
for i=1:numRows
    for j=1:numColumns

        XibGLLGLL  = ( (xibLR(i)+xibLR(i+1))/2+(xibLR(i+1)-xibLR(i))/2*xiGLL )'*ones(1,N+1);
        EtabGLLGLL = ones(N+1,1)*( (etabAB(j)+etabAB(j+1))/2+(etabAB(j+1)-etabAB(j))/2*etaGLL );

        XibGG  = ( (xibLR(i)+xibLR(i+1))/2+(xibLR(i+1)-xibLR(i))/2*xiG )'*ones(1,N);
        EtabGG = ones(N,1)*( (etabAB(j)+etabAB(j+1))/2+(etabAB(j+1)-etabAB(j))/2*etaG );

        % gridtype = 'sinecurve'
        XeGLLGLL = XibGLLGLL + cc*sin(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);
        YeGLLGLL = EtabGLLGLL+ cc*sin(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);
        
        XeGG = XibGG + cc*sin(pi*XibGG).*sin(pi*EtabGG);
        YeGG = EtabGG+ cc*sin(pi*XibGG).*sin(pi*EtabGG);

        if i<numColumns && j<numRows
            XGLLGLL((i-1)*N+(1:N),(j-1)*N+(1:N)) = XeGLLGLL(1:N,1:N);
            YGLLGLL((i-1)*N+(1:N),(j-1)*N+(1:N)) = YeGLLGLL(1:N,1:N);
        elseif i<numColumns && j==numRows
            XGLLGLL((i-1)*N+(1:N),(j-1)*N+(1:N+1)) = XeGLLGLL(1:N,1:N+1);
            YGLLGLL((i-1)*N+(1:N),(j-1)*N+(1:N+1)) = YeGLLGLL(1:N,1:N+1);
        elseif i==numColumns && j<numRows
            XGLLGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N)) = XeGLLGLL(1:N+1,1:N);
            YGLLGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N)) = YeGLLGLL(1:N+1,1:N);
        elseif i==numColumns && j==numRows
            XGLLGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N+1)) = XeGLLGLL;
            YGLLGLL((i-1)*N+(1:N+1),(j-1)*N+(1:N+1)) = YeGLLGLL;
        end
        XGG((i-1)*N+(1:N),(j-1)*N+(1:N)) = XeGG;
        YGG((i-1)*N+(1:N),(j-1)*N+(1:N)) = YeGG;


    end
end

% surf(XGLLGLL,YGLLGLL,zeros(size(XGLLGLL))); view([0 0 1])
