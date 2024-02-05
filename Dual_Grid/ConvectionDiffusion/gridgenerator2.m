function [XGLLGLL,YGLLGLL,XGG,YGG,QGLLGLL,JGLLGLL,JGG] = gridgenerator2()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global cc

% grid option
cc = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N N2 numColumns numRows
global xiGLL xiG xiEG
global wGLL wG

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

xibLR  = linspace(-1,1,numColumns+1);
etabAB = linspace(-1,1,numRows+1)   ;

% XI  = xibLR'*ones(1,numRows+1);
% ETA = ones(numColumns+1,1)*etabAB;
% 
% Xmain = XI + cc*sin(pi*XI).*sin(pi*ETA);
% Ymain = ETA+ cc*sin(pi*XI).*sin(pi*ETA);

% Global mesh
XibGLLGLL = zeros(N*numColumns+1,N*numRows+1);
EtabGLLGLL = zeros(N*numColumns+1,N*numRows+1);
XGLLGLL = zeros(N*numColumns+1,N*numRows+1);
YGLLGLL = zeros(N*numColumns+1,N*numRows+1);
XGG   = zeros(N*numColumns  ,N*numRows  );
YGG   = zeros(N*numColumns  ,N*numRows  );
JGLLGLL = zeros(2*(N+1)^2,numRows*numColumns);
JGG     = zeros(2*N2,numRows*numColumns);
QGLLGLL = zeros(2*(N+1)^2,3*numRows*numColumns);

for r=1:numRows
    for c=1:numColumns
        
        rc = c+(r-1)*numColumns;

        XibGLLGLL  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xiGLL )'*ones(1,N+1);
        EtabGLLGLL = ones(N+1,1)*( (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*etaGLL );

        XibGG  = ( (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*xiG )'*ones(1,N);
        EtabGG = ones(N,1)*( (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*etaG );
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gridtype = 'sinecurve' % This might be more advanced !!!!!
        XeGLLGLL = XibGLLGLL + cc*sin(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);
        YeGLLGLL = EtabGLLGLL+ cc*sin(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);

        XeGG = XibGG + cc*sin(pi*XibGG).*sin(pi*EtabGG);
        YeGG = EtabGG+ cc*sin(pi*XibGG).*sin(pi*EtabGG);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if c<numColumns && r<numRows
            XGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N)) = XeGLLGLL(1:N,1:N);
            YGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N)) = YeGLLGLL(1:N,1:N);
        elseif c<numColumns && r==numRows
            XGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = XeGLLGLL(1:N,1:N+1);
            YGLLGLL((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = YeGLLGLL(1:N,1:N+1);
        elseif c==numColumns && r<numRows
            XGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = XeGLLGLL(1:N+1,1:N);
            YGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = YeGLLGLL(1:N+1,1:N);
        elseif c==numColumns && r==numRows
            XGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = XeGLLGLL;
            YGLLGLL((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = YeGLLGLL;
        end
        XGG((c-1)*N+(1:N),(r-1)*N+(1:N)) = XeGG;
        YGG((c-1)*N+(1:N),(r-1)*N+(1:N)) = YeGG;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dXibdXi   = 1/numColumns;
        dXibdEta  = 0;
        dEtabdXi  = 0;
        dEtabdEta = 1/numRows;
        
        % gridtype = 'sinecurve'
        dXdXibGLLGLL  = 1+pi*cc*cos(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);
        dXdEtabGLLGLL = pi*cc*sin(pi*XibGLLGLL).*cos(pi*EtabGLLGLL);
        dYdXibGLLGLL  = pi*cc*cos(pi*XibGLLGLL).*sin(pi*EtabGLLGLL);
        dYdEtabGLLGLL = 1+pi*cc*sin(pi*XibGLLGLL).*cos(pi*EtabGLLGLL);

        dXdXiGLLGLL  = dXdXibGLLGLL*dXibdXi +dXdEtabGLLGLL*dEtabdXi ;
        dXdEtaGLLGLL = dXdXibGLLGLL*dXibdEta+dXdEtabGLLGLL*dEtabdEta;
        dYdXiGLLGLL  = dYdXibGLLGLL*dXibdXi +dYdEtabGLLGLL*dEtabdXi ;
        dYdEtaGLLGLL = dYdXibGLLGLL*dXibdEta+dYdEtabGLLGLL*dEtabdEta;
        
        dXdXibGG  = 1+pi*cc*cos(pi*XibGG).*sin(pi*EtabGG);
        dXdEtabGG = pi*cc*sin(pi*XibGG).*cos(pi*EtabGG);
        dYdXibGG  = pi*cc*cos(pi*XibGG).*sin(pi*EtabGG);
        dYdEtabGG = 1+pi*cc*sin(pi*XibGG).*cos(pi*EtabGG);

        dXdXiGG  = dXdXibGG*dXibdXi +dXdEtabGG*dEtabdXi ;
        dXdEtaGG = dXdXibGG*dXibdEta+dXdEtabGG*dEtabdEta;
        dYdXiGG  = dYdXibGG*dXibdXi +dYdEtabGG*dEtabdXi ;
        dYdEtaGG = dYdXibGG*dXibdEta+dYdEtabGG*dEtabdEta;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        JacobianGLLGLL =  dXdXiGLLGLL.*dYdEtaGLLGLL-dXdEtaGLLGLL.*dYdXiGLLGLL;

        JacobianGG =  dXdXiGG.*dYdEtaGG-dXdEtaGG.*dYdXiGG;

        JGLLGLL(:,rc) = kron(reshape(JacobianGLLGLL,1,(N+1)^2),ones(1,2))';

        JGG(:,rc) = kron(reshape(JacobianGG,1,N2),ones(1,2))';

        % JGLLGLL = spdiags(kron(reshape(JacobianGLLGLL,1,(N+1)^2),ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

        % JGG = spdiags(kron(reshape(JacobianGG,1,N2),ones(1,2))',0,2*N2,2*N2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        q11 = kron(reshape(( dXdXiGLLGLL./JacobianGLLGLL),1,(N+1)^2),[1 0])';
        q22 = kron(reshape((dYdEtaGLLGLL./JacobianGLLGLL),1,(N+1)^2),[0 1])';
        q12 = kron(reshape((dXdEtaGLLGLL./JacobianGLLGLL),1,(N+1)^2),[0 1])';
        q21 = kron(reshape(( dYdXiGLLGLL./JacobianGLLGLL),1,(N+1)^2),[1 0])';

        QGLLGLL(:,3*(rc-1)+(1:3)) = [q21 q11+q22 q12];

        % QGLLGLL = spdiags([q21 q11+q22 q12],-1:1,2*(N+1)^2,2*(N+1)^2);

    end
end