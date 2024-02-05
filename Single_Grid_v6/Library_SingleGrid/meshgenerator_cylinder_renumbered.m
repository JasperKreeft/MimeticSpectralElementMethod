function Mesh = meshgenerator_cylinder_renumbered(R)

global N numElements
global xi w

[xi,w] = GLLnodes(N);  % Gauss-Lobotto-Legendre


xglobal = [ -1.5 -0.75 0 0.75 3 -1.5 -0.75 0.75 3 -1.5 -0.75 0 0.75 3 ...
                  -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 0 -sqrt(2)/4 ];
yglobal = [ -0.75 -0.75 -0.75 -0.75 -0.75 0 0 0 0 0.75 0.75 0.75 0.75 0.75 ...
                  0 -sqrt(2)/4 -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 ];


% plot(xglobal,yglobal,'o','markerface','b')
% hold on
% axis equal
% axis([-1.5 3 -0.75 0.75])


corners = [  6  1  2  7
            10  6  7 11
             7  2 16 15
            11  7 15 22
             2  3 17 16
             3  4 18 17
            12 11 22 21
            13 12 21 20
             4  8 19 18
             8 13 20 19
             5  9  8  4
             9 14 13  8 ];


% Global Mesh
Mesh.X      = zeros((N+1)^2,numElements);
Mesh.Y      = zeros((N+1)^2,numElements);
Mesh.dXdXi  = zeros((N+1)^2,numElements);
Mesh.dXdEta = zeros((N+1)^2,numElements);
Mesh.dYdXi  = zeros((N+1)^2,numElements);
Mesh.dYdEta = zeros((N+1)^2,numElements);
Mesh.J      = zeros((N+1)^2,numElements);
Mesh.Qinv   = zeros(2*(N+1)^2,3*numElements);

for i=[1 2 11 12]

[ Mesh.X(:,i),Mesh.dXdXi(:,i),Mesh.dXdEta(:,i) ] = transfinitemapping(xglobal(corners(i,:)),N+1,xi);
[ Mesh.Y(:,i),Mesh.dYdXi(:,i),Mesh.dYdEta(:,i) ] = transfinitemapping(yglobal(corners(i,:)),N+1,xi);
                          
%     pcolor(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),(i-1)/11*ones(N+1))

    Mesh.J(:,i) = Mesh.dXdXi(:,i).*Mesh.dYdEta(:,i)-Mesh.dXdEta(:,i).*Mesh.dYdXi(:,i);

    qinv11 = kron(Mesh.dXdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);
    qinv22 = kron(Mesh.dYdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv12 = kron(Mesh.dXdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv21 = kron(Mesh.dYdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);

    Mesh.Qinv(:,3*(i-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];

end



tblr = zeros(8,4);
% tblr(1,4) = 1; tblr(2,1) = 1; tblr(3,1) = 1; tblr(4,3) = 1;
% tblr(5,4) = 1; tblr(6,2) = 1; tblr(7,2) = 1; tblr(8,3) = 1;
tblr(:,1) = 1;

theta = [     -pi -3/4*pi
          -3/4*pi   -pi/2 
            -pi/2   -pi/4
            -pi/4       0
                0    pi/4
             pi/4    pi/2
             pi/2  3/4*pi
           3/4*pi      pi ];

k=0;
for i=[ 3 5 6 9 10 8 7 4 ]

    k=k+1;
    cor = [ xglobal(corners(i,:))
            yglobal(corners(i,:)) ];
        
[Mesh.X(:,i),Mesh.Y(:,i),Mesh.dXdXi(:,i),Mesh.dXdEta(:,i),Mesh.dYdXi(:,i),Mesh.dYdEta(:,i)] = ...
                transfinitemapping_cylinder(cor,R,theta(k,:),tblr(k,:),N+1,xi);

%     pcolor(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),(i-1)/11*ones(N+1))

    Mesh.J(:,i) = Mesh.dXdXi(:,i).*Mesh.dYdEta(:,i)-Mesh.dXdEta(:,i).*Mesh.dYdXi(:,i);

    qinv11 = kron(Mesh.dXdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);
    qinv22 = kron(Mesh.dYdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv12 = kron(Mesh.dXdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv21 = kron(Mesh.dYdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);

    Mesh.Qinv(:,3*(i-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cilinder

% Rc = R;
% xc = [linspace(-R,R,100) linspace(R,-R,100)]; % coordinates of cilinder
% yc = sqrt(Rc^2-xc.^2); yc(101:200) = -yc(101:200);
% 
% plot(xc,yc,'k','linewidth',4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration variables for force and boundary conditions

global n xibar Gw

n = 8;

[xibar,Gw] = GLnodes(n);
xibar = xibar'; Gw = Gw';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% meshplot