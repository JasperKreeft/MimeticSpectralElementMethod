function [Meshx,Meshy] = meshgenerator_square_covector(Domain,DomInfo,corners)
% Gecontroleerd en klopt !!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N numElements numRows numColumns
global xi w

if isempty(numElements)
    numElements = 1;
end
if isempty(numRows)
    numRows = sqrt(numElements);
    numColumns = numRows;
end


[xi,w] = GLLnodes(N);  % Gauss-Lobotto-Legendre
eta = xi;
[xieg,weg] = EGnodes(N);
etaeg = xieg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = kron(ones(1,numRows+1),linspace(-1,1,numColumns+1));
etabAB = kron(linspace(-1,1,numRows+1),ones(1,numColumns+1));

% Global Mesh
Mesh.X      = zeros((N+1)*(N+2),numElements);
Mesh.Y      = zeros((N+1)*(N+2),numElements);
Mesh.dXdXi  = zeros((N+1)*(N+2),numElements);
Mesh.dXdEta = zeros((N+1)*(N+2),numElements);
Mesh.dYdXi  = zeros((N+1)*(N+2),numElements);
Mesh.dYdEta = zeros((N+1)*(N+2),numElements);
Mesh.J      = zeros((N+1)*(N+2),numElements);
Mesh.Qinv   = zeros(2*(N+1)*(N+2),3*numElements);

Xi  = kron(ones(1,N+1),xieg);
Eta = kron(eta,ones(1,N+2));

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

        % Everything is in GLL-GLL nodes
        ind = c+(r-1)*(numColumns+1);
        Xib  = (xibLR(ind)+xibLR(ind+1))/2+(xibLR(ind+1)-xibLR(ind))/2*Xi;
        Etab = (etabAB(ind)+etabAB(ind+numColumns+1))/2+(etabAB(ind+numColumns+1)-etabAB(ind))/2*Eta;

        [Xe,Ye,dXdXib,dXdEtab,dYdXib,dYdEtab] = coordinatemap(Xib,Etab,Domain,DomInfo);

        Mesh.X(:,rc) = Xe;
        Mesh.Y(:,rc) = Ye;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dXibdXi   = (Xib(1,(N+1)*(N+2))-Xib(1,1))/2*ones(1,(N+1)*(N+2));
        dXibdEta  = zeros(1,(N+1)*(N+2));
        dEtabdXi  = zeros(1,(N+1)*(N+2));
        dEtabdEta = (Etab(1,(N+1)*(N+2))-Etab(1,1))/2*ones(1,(N+1)*(N+2));

        Mesh.dXdXi(:,rc)  = dXdXib.*dXibdXi +dXdEtab.*dEtabdXi ;
        Mesh.dXdEta(:,rc) = dXdXib.*dXibdEta+dXdEtab.*dEtabdEta;
        Mesh.dYdXi(:,rc)  = dYdXib.*dXibdXi +dYdEtab.*dEtabdXi ;
        Mesh.dYdEta(:,rc) = dYdXib.*dXibdEta+dYdEtab.*dEtabdEta;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Mesh.J(:,rc) = Mesh.dXdXi(:,rc).*Mesh.dYdEta(:,rc)-Mesh.dXdEta(:,rc).*Mesh.dYdXi(:,rc);

        qinv11 = kron(Mesh.dXdXi(:,rc)./Mesh.J(:,rc) ,[1 ; 0]);
        qinv22 = kron(Mesh.dYdEta(:,rc)./Mesh.J(:,rc),[0 ; 1]);
        qinv12 = kron(Mesh.dXdEta(:,rc)./Mesh.J(:,rc),[0 ; 1]);
        qinv21 = kron(Mesh.dYdXi(:,rc)./Mesh.J(:,rc) ,[1 ; 0]);

        Mesh.Qinv(:,3*(rc-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];

    end
end

Meshx = Mesh;



% Global Mesh
Mesh.X      = zeros((N+1)*(N+2),numElements);
Mesh.Y      = zeros((N+1)*(N+2),numElements);
Mesh.dXdXi  = zeros((N+1)*(N+2),numElements);
Mesh.dXdEta = zeros((N+1)*(N+2),numElements);
Mesh.dYdXi  = zeros((N+1)*(N+2),numElements);
Mesh.dYdEta = zeros((N+1)*(N+2),numElements);
Mesh.J      = zeros((N+1)*(N+2),numElements);
Mesh.Qinv   = zeros(2*(N+1)*(N+2),3*numElements);

Xi  = kron(ones(1,N+2),xi);
Eta = kron(etaeg,ones(1,N+1));

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

        % Everything is in GLL-GLL nodes
        ind = c+(r-1)*(numColumns+1);
        Xib  = (xibLR(ind)+xibLR(ind+1))/2+(xibLR(ind+1)-xibLR(ind))/2*Xi;
        Etab = (etabAB(ind)+etabAB(ind+numColumns+1))/2+(etabAB(ind+numColumns+1)-etabAB(ind))/2*Eta;

        [Xe,Ye,dXdXib,dXdEtab,dYdXib,dYdEtab] = coordinatemap(Xib,Etab,Domain,DomInfo);

        Mesh.X(:,rc) = Xe;
        Mesh.Y(:,rc) = Ye;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dXibdXi   = (Xib(1,(N+1)*(N+2))-Xib(1,1))/2*ones(1,(N+1)*(N+2));
        dXibdEta  = zeros(1,(N+1)*(N+2));
        dEtabdXi  = zeros(1,(N+1)*(N+2));
        dEtabdEta = (Etab(1,(N+1)*(N+2))-Etab(1,1))/2*ones(1,(N+1)*(N+2));

        Mesh.dXdXi(:,rc)  = dXdXib.*dXibdXi +dXdEtab.*dEtabdXi ;
        Mesh.dXdEta(:,rc) = dXdXib.*dXibdEta+dXdEtab.*dEtabdEta;
        Mesh.dYdXi(:,rc)  = dYdXib.*dXibdXi +dYdEtab.*dEtabdXi ;
        Mesh.dYdEta(:,rc) = dYdXib.*dXibdEta+dYdEtab.*dEtabdEta;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Mesh.J(:,rc) = Mesh.dXdXi(:,rc).*Mesh.dYdEta(:,rc)-Mesh.dXdEta(:,rc).*Mesh.dYdXi(:,rc);

        qinv11 = kron(Mesh.dXdXi(:,rc)./Mesh.J(:,rc) ,[1 ; 0]);
        qinv22 = kron(Mesh.dYdEta(:,rc)./Mesh.J(:,rc),[0 ; 1]);
        qinv12 = kron(Mesh.dXdEta(:,rc)./Mesh.J(:,rc),[0 ; 1]);
        qinv21 = kron(Mesh.dYdXi(:,rc)./Mesh.J(:,rc) ,[1 ; 0]);

        Mesh.Qinv(:,3*(rc-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];

    end
end

Meshy = Mesh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration variables for force and boundary conditions

global n xibar Gw

n = 8;

[xibar,Gw] = Gnodes(n);
xibar = xibar'; Gw = Gw';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%