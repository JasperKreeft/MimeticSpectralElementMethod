function Mesh = meshgenerator_triangle_E9
% Gecontroleerd en klopt !!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N
global xi w


[xi,w] = GLLnodes(N);  % Gauss-Lobotto-Legendre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = sqrt(3)/3;
b = 1/10;

xglobal = [ -1 -1+b  -1+b   -1+b/2       0 0 1-b  1-b     1 1-b/2       -1/2 1/2   -b/2          0       b/2         0 ];
yglobal = [ -a  -a  a*(b-1) a*(3/2*b-1) -a 0 -a  a*(b-1) -a a*(3/2*b-1)  a/2 a/2 a*(2-3/2*b) 2*a*(1-b) a*(2-3/2*b) 2*a ];

corners = [  1  2  3  4
             2  5  6  3
             5  7  8  6
             7  9 10  8
             3  6 11  4
             8 10 12  6
            11  6 14 13
             6 12 15 14
            14 15 16 13 ];


% Global Mesh
Mesh.X      = zeros((N+1)^2,9);
Mesh.Y      = zeros((N+1)^2,9);
Mesh.J      = zeros((N+1)^2,9);
Mesh.Qinv   = zeros(2*(N+1)^2,3*9);
Mesh.dXdXi  = zeros((N+1)^2,9);
Mesh.dXdEta = zeros((N+1)^2,9);
Mesh.dYdXi  = zeros((N+1)^2,9);
Mesh.dYdEta = zeros((N+1)^2,9);

for i=1:9

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