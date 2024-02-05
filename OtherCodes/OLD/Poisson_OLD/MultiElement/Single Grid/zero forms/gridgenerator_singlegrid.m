function [X,Y,Qinv,J] = gridgenerator_singlegrid()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N numColumns numRows
global xi w
global cc


[xi,w] = GLLnodes(N);  % Gauss-Lobotto-Legendre

Xi = xi'*ones(1,N+1); Eta = Xi';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = linspace(-1,1,numColumns+1);
etabAB = linspace(-1,1,numRows+1)   ;


% Global mesh
X    = zeros(N*numColumns+1,N*numRows+1);
Y    = zeros(N*numColumns+1,N*numRows+1);
% J    = zeros(2*(N+1)^2,numRows*numColumns);
J    = zeros((N+1)^2,numRows*numColumns);
Qinv = zeros(2*(N+1)^2,3*numRows*numColumns);

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        
        % Everything is in GLL-GLL nodes
        Xib  = (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*Xi;
        Etab = (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*Eta;

        % gridtype = 'sinecurve' % This might be more advanced !!!!!
        Xe = Xib + cc*sin(pi*Xib).*sin(pi*Etab);
        Ye = Etab+ cc*sin(pi*Xib).*sin(pi*Etab);

        if c<numColumns && r<numRows
            X((c-1)*N+(1:N),(r-1)*N+(1:N)) = Xe(1:N,1:N);
            Y((c-1)*N+(1:N),(r-1)*N+(1:N)) = Ye(1:N,1:N);
        elseif c<numColumns && r==numRows
            X((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = Xe(1:N,1:N+1);
            Y((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = Ye(1:N,1:N+1);
        elseif c==numColumns && r<numRows
            X((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = Xe(1:N+1,1:N);
            Y((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = Ye(1:N+1,1:N);
        elseif c==numColumns && r==numRows
            X((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = Xe;
            Y((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = Ye;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         dXibdXi   = 1/numColumns;
%         dXibdEta  = 0;
%         dEtabdXi  = 0;
%         dEtabdEta = 1/numRows;
        
            dXibdXi   = (Xib(N+1,1)-Xib(1,1))/2*ones(N+1);
            dXibdEta  = zeros(N+1);
            dEtabdXi  = zeros(N+1);
            dEtabdEta = (Etab(1,N+1)-Etab(1,1))/2*ones(N+1);

%             dXibdXi   = 2/(Xib(N+1,1)-Xib(1,1))*ones(N+1);
%             dXibdEta  = zeros(N+1);
%             dEtabdXi  = zeros(N+1);
%             dEtabdEta = 2/(Etab(1,N+1)-Etab(1,1))*ones(N+1);
% keyboard
        % gridtype = 'sinecurve' [-1,1] -> [0,pi]
        dXdXib  = 1+pi*cc*cos(pi*Xib).*sin(pi*Etab);
        dXdEtab = pi*cc*sin(pi*Xib).*cos(pi*Etab);
        dYdXib  = pi*cc*cos(pi*Xib).*sin(pi*Etab);
        dYdEtab = 1+pi*cc*sin(pi*Xib).*cos(pi*Etab);

        dXdXi  = dXdXib.*dXibdXi +dXdEtab.*dEtabdXi ;
        dXdEta = dXdXib.*dXibdEta+dXdEtab.*dEtabdEta;
        dYdXi  = dYdXib.*dXibdXi +dYdEtab.*dEtabdXi ;
        dYdEta = dYdXib.*dXibdEta+dYdEtab.*dEtabdEta;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Jacobian =  dXdXi.*dYdEta-dXdEta.*dYdXi;

        J(:,rc) = reshape(Jacobian,1,(N+1)^2)';
%         J(:,rc) = kron(reshape(Jacobian,1,(N+1)^2),ones(1,2))';
%         J(:,:,rc) = Jacobian;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        qinv11 = kron(reshape(( dXdXi./Jacobian),1,(N+1)^2),[1 0])';
        qinv22 = kron(reshape((dYdEta./Jacobian),1,(N+1)^2),[0 1])';
        qinv12 = kron(reshape((dXdEta./Jacobian),1,(N+1)^2),[0 1])';
        qinv21 = kron(reshape(( dYdXi./Jacobian),1,(N+1)^2),[1 0])';

        Qinv(:,3*(rc-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];

%         Qinv = spdiags([qinv21 qinv11+qinv22 qinv12],-1:1,2*(N+1)^2,2*(N+1)^2);
% keyboard
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%