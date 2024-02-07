function Mesh = meshgenerator_square(Domain,DomInfo,corners)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Main Mesh creation. Here for multi-element on unit square [0,1]x[0,1] or
% standard square [-1,1]x[-1,1]

xibLR  = kron(ones(1,numRows+1),linspace(-1,1,numColumns+1));
etabAB = kron(linspace(-1,1,numRows+1),ones(1,numColumns+1));

% xibLR  = kron(ones(1,numRows+1),GLLnodes(numColumns));
% etabAB = kron(GLLnodes(numRows),ones(1,numColumns+1));

% Global Mesh
Mesh.X      = zeros((N+1)^2,numElements);
Mesh.Y      = zeros((N+1)^2,numElements);
Mesh.dXdXi  = zeros((N+1)^2,numElements);
Mesh.dXdEta = zeros((N+1)^2,numElements);
Mesh.dYdXi  = zeros((N+1)^2,numElements);
Mesh.dYdEta = zeros((N+1)^2,numElements);
Mesh.J      = zeros((N+1)^2,numElements);
Mesh.Qinv   = zeros(2*(N+1)^2,3*numElements);

% for i=1:numElements
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Transfinite elements only !!!
% 
% [Mesh.X(:,i),Mesh.dXdXi(:,i),Mesh.dXdEta(:,i)] = ...
%                               transfinitemapping_v2(xibLR(corners(i,:)));
% [Mesh.Y(:,i),Mesh.dYdXi(:,i),Mesh.dYdEta(:,i)] = ...
%                               transfinitemapping_v2(etabAB(corners(i,:)));
% 
% pcolor(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),(i-1)/2*ones(N+1))
% hold on
% axis equal
% axis([-1.2 1.2 -1.2 1.2])
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     Mesh.J(:,i) = Mesh.dXdXi(:,i).*Mesh.dYdEta(:,i)-Mesh.dXdEta(:,i).*Mesh.dYdXi(:,i);
% 
%     qinv11 = kron(Mesh.dXdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);
%     qinv22 = kron(Mesh.dYdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
%     qinv12 = kron(Mesh.dXdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
%     qinv21 = kron(Mesh.dYdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);
% 
%     Mesh.Qinv(:,3*(i-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];
% 
% end


Xi  = kron(ones(1,N+1),xi);
Eta = kron(eta,ones(1,N+1));

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

        dXibdXi   = (Xib(1,(N+1)^2)-Xib(1,1))/2*ones(1,(N+1)^2);
        dXibdEta  = zeros(1,(N+1)^2);
        dEtabdXi  = zeros(1,(N+1)^2);
        dEtabdEta = (Etab(1,(N+1)^2)-Etab(1,1))/2*ones(1,(N+1)^2);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration variables for force and boundary conditions

global n xibar Gw

n = 8;

[xibar,Gw] = GLnodes(n);
xibar = xibar'; Gw = Gw';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%