function Mesh = meshgenerator_square_smooth(Domain,DomInfo)

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
if strcmp(Domain,'SinDeformGrid'); L=2; end
if strcmp(Domain,'SinDeformGrid_01'); L=1; end
if strcmp(Domain,'Vorticity'); L=DomInfo; end

z = [];
for h=1:numColumns
    z = [ z  (2*h-1)*L/(2*numColumns)-L/2+L/(2*numColumns)*xi' ];
end
Xi = repmat(z,N+1,numRows);
z = [];
for h=1:numRows
    z = [ z  (2*h-1)*L/(2*numRows)-L/2+L/(2*numRows)*eta' ];
end
Eta = kron(z,ones(N+1,numColumns));

[Mesh.X,Mesh.dXdXi]  = map(Xi,L,numColumns);
[Mesh.Y,Mesh.dYdEta] = map(Eta,L,numRows);

Mesh.dXdEta = zeros((N+1)^2,numElements);
Mesh.dYdXi  = zeros((N+1)^2,numElements);


for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

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