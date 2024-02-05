function Mesh = meshgenerator_adaptive(Domain,DomInfo,corners)
% Gecontroleerd en klopt !!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Nadaptive numElements numRows numColumns
global xi w

if isempty(numElements)
    numElements = 1;
end
if isempty(numRows)
    numRows = sqrt(numElements);
    numColumns = numRows;
end

Nmax = max(Nadaptive);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = kron(ones(1,numRows+1),linspace(-1,1,numColumns+1));
etabAB = kron(linspace(-1,1,numRows+1),ones(1,numColumns+1));

% xibLR  = kron(ones(1,numRows+1),GLLnodes(numColumns));
% etabAB = kron(GLLnodes(numRows),ones(1,numColumns+1));


% Global Mesh
Mesh.X      = zeros((Nmax+1)^2,numElements);
Mesh.Y      = zeros((Nmax+1)^2,numElements);
Mesh.dXdXi  = zeros((Nmax+1)^2,numElements);
Mesh.dXdEta = zeros((Nmax+1)^2,numElements);
Mesh.dYdXi  = zeros((Nmax+1)^2,numElements);
Mesh.dYdEta = zeros((Nmax+1)^2,numElements);
Mesh.J      = zeros((Nmax+1)^2,numElements);
Mesh.Qinv   = zeros(2*(Nmax+1)^2,3*numElements);



for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        Nel = Nadaptive(rc);

        [xi,w] = GLLnodes(Nel);  % Gauss-Lobotto-Legendre
        eta = xi;
        Xi  = kron(ones(1,Nel+1),xi);
        Eta = kron(eta,ones(1,Nel+1));
        
        % Everything is in GLL-GLL nodes
        ind = c+(r-1)*(numColumns+1);
        Xib  = (xibLR(ind)+xibLR(ind+1))/2+(xibLR(ind+1)-xibLR(ind))/2*Xi;
        Etab = (etabAB(ind)+etabAB(ind+numColumns+1))/2+(etabAB(ind+numColumns+1)-etabAB(ind))/2*Eta;

        [Xe,Ye,dXdXib,dXdEtab,dYdXib,dYdEtab] = coordinatemap(Xib,Etab,Domain,DomInfo);

        ind = 1:(Nel+1)^2;
        Mesh.X(ind,rc) = Xe;
        Mesh.Y(ind,rc) = Ye;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        dXibdXi   = (Xib(1,(Nel+1)^2)-Xib(1,1))/2*ones(1,(Nel+1)^2);
        dXibdEta  = zeros(1,(Nel+1)^2);
        dEtabdXi  = zeros(1,(Nel+1)^2);
        dEtabdEta = (Etab(1,(Nel+1)^2)-Etab(1,1))/2*ones(1,(Nel+1)^2);

        Mesh.dXdXi(ind,rc)  = dXdXib.*dXibdXi +dXdEtab.*dEtabdXi ;
        Mesh.dXdEta(ind,rc) = dXdXib.*dXibdEta+dXdEtab.*dEtabdEta;
        Mesh.dYdXi(ind,rc)  = dYdXib.*dXibdXi +dYdEtab.*dEtabdXi ;
        Mesh.dYdEta(ind,rc) = dYdXib.*dXibdEta+dYdEtab.*dEtabdEta;

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
xibar = xibar'; Gw = Gw';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%