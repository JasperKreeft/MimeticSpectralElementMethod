function Mesh = meshgenerator_triangle
% Gecontroleerd en klopt !!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N
global xi w


[xi,w] = GLLnodes(N);  % Gauss-Lobotto-Legendre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = sqrt(3)/3;

xglobal = [ -1  0  1 1/2 0  -1/2 0 ];
yglobal = [ -a -a -a a/2 2*a a/2 0 ];


corners = [  1 2 7 6
             2 3 4 7
             7 4 5 6 ];


% Global Mesh
Mesh.X      = zeros((N+1)^2,3);
Mesh.Y      = zeros((N+1)^2,3);
Mesh.J      = zeros((N+1)^2,3);
Mesh.Qinv   = zeros(2*(N+1)^2,3*3);
Mesh.dXdXi  = zeros((N+1)^2,3);
Mesh.dXdEta = zeros((N+1)^2,3);
Mesh.dYdXi  = zeros((N+1)^2,3);
Mesh.dYdEta = zeros((N+1)^2,3);

for i=1:3

[Mesh.X(:,i),Mesh.dXdXi(:,i),Mesh.dXdEta(:,i)] = ...
                              transfinitemapping_v2(xglobal(corners(i,:)));
[Mesh.Y(:,i),Mesh.dYdXi(:,i),Mesh.dYdEta(:,i)] = ...
                              transfinitemapping_v2(yglobal(corners(i,:)));

% pcolor(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),(i-1)/2*ones(N+1))
% hold on
% axis equal
% axis([-1.2 1.2 -.6 1.2])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Mesh.J(:,i) = Mesh.dXdXi(:,i).*Mesh.dYdEta(:,i)-Mesh.dXdEta(:,i).*Mesh.dYdXi(:,i);

    qinv11 = kron(Mesh.dXdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);
    qinv22 = kron(Mesh.dYdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv12 = kron(Mesh.dXdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv21 = kron(Mesh.dYdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);

    Mesh.Qinv(:,3*(i-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration variables for force and boundary conditions

global n xibar Gw

n = 8;

[xibar,Gw] = Gnodes(n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%