function Mesh = meshgenerator(Domain,DomInfo)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N numColumns numRows
global xi w

if isempty(numColumns) || isempty(numRows)
    numColumns = 1;
    numRows    = 1;
end
    

[xi,w] = GLLnodes(N);  % Gauss-Lobotto-Legendre

Xi = xi'*ones(1,N+1); Eta = Xi';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = linspace(-1,1,numColumns+1);
etabAB = linspace(-1,1,numRows+1)   ;


% Global Mesh
Mesh.X    = zeros(N*numColumns+1,N*numRows+1);
Mesh.Y    = zeros(N*numColumns+1,N*numRows+1);
Mesh.J    = zeros((N+1)^2,numRows*numColumns);
Mesh.Qinv = zeros(2*(N+1)^2,3*numRows*numColumns);
Mesh.dXdXi  = zeros((N+1)^2,numRows*numColumns);
Mesh.dXdEta = zeros((N+1)^2,numRows*numColumns);
Mesh.dYdXi  = zeros((N+1)^2,numRows*numColumns);
Mesh.dYdEta = zeros((N+1)^2,numRows*numColumns);

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        
        % Everything is in GLL-GLL nodes
        Xib  = (xibLR(c)+xibLR(c+1))/2+(xibLR(c+1)-xibLR(c))/2*Xi;
        Etab = (etabAB(r)+etabAB(r+1))/2+(etabAB(r+1)-etabAB(r))/2*Eta;

        [Xe,Ye,dXdXib,dXdEtab,dYdXib,dYdEtab] = coordinatemap(Xib,Etab,Domain,DomInfo);

        if c<numColumns && r<numRows
            Mesh.X((c-1)*N+(1:N),(r-1)*N+(1:N)) = Xe(1:N,1:N);
            Mesh.Y((c-1)*N+(1:N),(r-1)*N+(1:N)) = Ye(1:N,1:N);
        elseif c<numColumns && r==numRows
            Mesh.X((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = Xe(1:N,1:N+1);
            Mesh.Y((c-1)*N+(1:N),(r-1)*N+(1:N+1)) = Ye(1:N,1:N+1);
        elseif c==numColumns && r<numRows
            Mesh.X((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = Xe(1:N+1,1:N);
            Mesh.Y((c-1)*N+(1:N+1),(r-1)*N+(1:N)) = Ye(1:N+1,1:N);
        elseif c==numColumns && r==numRows
            Mesh.X((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = Xe;
            Mesh.Y((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)) = Ye;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dXibdXi   = (Xib(N+1,1)-Xib(1,1))/2*ones(N+1);
        dXibdEta  = zeros(N+1);
        dEtabdXi  = zeros(N+1);
        dEtabdEta = (Etab(1,N+1)-Etab(1,1))/2*ones(N+1);

        dXdXi  = dXdXib.*dXibdXi +dXdEtab.*dEtabdXi ;
        dXdEta = dXdXib.*dXibdEta+dXdEtab.*dEtabdEta;
        dYdXi  = dYdXib.*dXibdXi +dYdEtab.*dEtabdXi ;
        dYdEta = dYdXib.*dXibdEta+dYdEtab.*dEtabdEta;
        
        Mesh.dXdXi(:,rc)  = reshape(dXdXi ,(N+1)^2,1);
        Mesh.dXdEta(:,rc) = reshape(dXdEta,(N+1)^2,1);
        Mesh.dYdXi(:,rc)  = reshape(dYdXi ,(N+1)^2,1);
        Mesh.dYdEta(:,rc) = reshape(dYdEta,(N+1)^2,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Mesh.J(:,rc) = Mesh.dXdXi(:,rc).*Mesh.dYdEta(:,rc)-Mesh.dXdEta(:,rc).*Mesh.dYdXi(:,rc);

        qinv11 = kron(Mesh.dXdXi(:,rc)./Mesh.J(:,rc) ,[1 ; 0]);
        qinv22 = kron(Mesh.dYdEta(:,rc)./Mesh.J(:,rc),[0 ; 1]);
        qinv12 = kron(Mesh.dXdEta(:,rc)./Mesh.J(:,rc),[0 ; 1]);
        qinv21 = kron(Mesh.dYdXi(:,rc)./Mesh.J(:,rc) ,[1 ; 0]);

        Mesh.Qinv(:,3*(rc-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration variables for force and boundary conditions

global n xibar Gw

n = 8;

[xibar,Gw] = Gnodes(n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%