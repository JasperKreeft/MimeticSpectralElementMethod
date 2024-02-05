function Mesh = meshgenerator_naca0012_E16(R)

global N
global xi w
% close all
% R = 0.5; N=10;

[xi,w] = GLLnodes(N);  % Gauss-Lobotto-Legendre


xglobal = [ -2 0 2 2 2 0 -2 -2 ...
            -1 0 1 1 1 0 -1 -1 ...
            -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4 0 -sqrt(2)/4 -0.5 ];
yglobal = [ -2 -2 -2 0 2 2 2 0 ...
            -1 -1 -1 0 1 1 1 0 ...
            -sqrt(2)/4 -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 sqrt(2)/4  0 ];
         

% plot(xglobal,yglobal,'o','markerface','b')
hold on
axis equal
axis([-2 2 -2 2])


corners = [  1  2 10  9 % outer-ring
             2  3 11 10
             3  4 12 11
             4  5 13 12
             5  6 14 13
             6  7 15 14
             7  8 16 15
             8  1  9 16
             9 10 18 17 % inner-ring
            10 11 19 18
            11 12 20 19
            12 13 21 20
            13 14 22 21
            14 15 23 22
            15 16 24 23
            16  9 17 24 ];


% Global Mesh
Mesh.X      = zeros((N+1)^2,16);
Mesh.Y      = zeros((N+1)^2,16);
Mesh.dXdXi  = zeros((N+1)^2,16);
Mesh.dXdEta = zeros((N+1)^2,16);
Mesh.dYdXi  = zeros((N+1)^2,16);
Mesh.dYdEta = zeros((N+1)^2,16);
Mesh.J      = zeros((N+1)^2,16);
Mesh.Qinv   = zeros(2*(N+1)^2,3*16);

for i=1:8

[ Mesh.X(:,i),Mesh.dXdXi(:,i),Mesh.dXdEta(:,i) ] = transfinitemapping(xglobal(corners(i,:)),N+1,xi);
[ Mesh.Y(:,i),Mesh.dYdXi(:,i),Mesh.dYdEta(:,i) ] = transfinitemapping(yglobal(corners(i,:)),N+1,xi);
                          
    pcolor(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),(i-1)/15*ones(N+1))

    Mesh.J(:,i) = Mesh.dXdXi(:,i).*Mesh.dYdEta(:,i)-Mesh.dXdEta(:,i).*Mesh.dYdXi(:,i);

    qinv11 = kron(Mesh.dXdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);
    qinv22 = kron(Mesh.dYdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv12 = kron(Mesh.dXdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv21 = kron(Mesh.dYdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);

    Mesh.Qinv(:,3*(i-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];

end



tblr = zeros(8,4);
tblr(1,1) = 1; tblr(2,1) = 1; tblr(3,1) = 1; tblr(4,1) = 1;
tblr(5,1) = 1; tblr(6,1) = 1; tblr(7,1) = 1; tblr(8,1) = 1;

theta = [ -3/4*pi   -pi/2 
            -pi/2   -pi/4
            -pi/4       0
                0  1/4*pi
           1/4*pi    pi/2
             pi/2    3*pi/4
           3*pi/4    pi
              -pi   -3/4*pi ];

k=0;
for i=9:16
i
    k=k+1;
    cor = [ xglobal(corners(i,:))
            yglobal(corners(i,:)) ];

[Mesh.X(:,i),Mesh.Y(:,i),Mesh.dXdXi(:,i),Mesh.dXdEta(:,i),Mesh.dYdXi(:,i),Mesh.dYdEta(:,i)] = ...
                transfinitemapping_naca0012(cor,R,theta(k,:),tblr(k,:),N+1,xi);

    pcolor(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),(i-1)/15*ones(N+1))

    Mesh.J(:,i) = Mesh.dXdXi(:,i).*Mesh.dYdEta(:,i)-Mesh.dXdEta(:,i).*Mesh.dYdXi(:,i);

    qinv11 = kron(Mesh.dXdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);
    qinv22 = kron(Mesh.dYdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv12 = kron(Mesh.dXdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv21 = kron(Mesh.dYdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);

    Mesh.Qinv(:,3*(i-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];
% if i==11; keyboard; end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cilinder

% Rc = R;
% xc = [linspace(-R,R,100) linspace(R,-R,100)]; % coordinates of cilinder
% yc = sqrt(Rc^2-xc.^2); yc(101:200) = -yc(101:200);
% 
% plot(xc,yc,'k','linewidth',4)
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration variables for force and boundary conditions

global n xibar Gw

n = 8;

[xibar,Gw] = Gnodes(n);
xibar = xibar'; Gw = Gw';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%