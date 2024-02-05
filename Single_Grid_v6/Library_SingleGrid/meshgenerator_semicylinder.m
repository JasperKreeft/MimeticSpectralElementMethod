function Mesh = meshgenerator_semicylinder(R)

global N numElements
global xi w

[xi,w] = GLLnodes(N);  % Gauss-Lobotto-Legendre


xglobal = [ -1.5 -0.75 0 0.75 3 -1.5 -0.75 0.75 3 -0.5 -sqrt(2)/4 0 sqrt(2)/4 0.5 ];
yglobal = [ -0.75 -0.75 -0.75 -0.75 -0.75 0 0 0 0 0 -sqrt(2)/4 -0.5 -sqrt(2)/4 0 ];


% plot(xglobal,yglobal,'o','markerface','b')
% hold on
% axis equal
% axis([-1.5 3 -0.75 0])


corners = [  1  2  7  6
             2 11 10  7
             2  3 12 11
             3  4 13 12
            13  4  8 14
             4  5  9  8 ];


% Global Mesh
Mesh.X      = zeros((N+1)^2,numElements);
Mesh.Y      = zeros((N+1)^2,numElements);
Mesh.dXdXi  = zeros((N+1)^2,numElements);
Mesh.dXdEta = zeros((N+1)^2,numElements);
Mesh.dYdXi  = zeros((N+1)^2,numElements);
Mesh.dYdEta = zeros((N+1)^2,numElements);
Mesh.J      = zeros((N+1)^2,numElements);
Mesh.Qinv   = zeros(2*(N+1)^2,3*numElements);

for i=[1 6]

[ Mesh.X(:,i),Mesh.dXdXi(:,i),Mesh.dXdEta(:,i) ] = transfinitemapping(xglobal(corners(i,:)),N+1,xi);
[ Mesh.Y(:,i),Mesh.dYdXi(:,i),Mesh.dYdEta(:,i) ] = transfinitemapping(yglobal(corners(i,:)),N+1,xi);
                          
%     pcolor(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),(i-1)/5*ones(N+1))

    Mesh.J(:,i) = Mesh.dXdXi(:,i).*Mesh.dYdEta(:,i)-Mesh.dXdEta(:,i).*Mesh.dYdXi(:,i);

    qinv11 = kron(Mesh.dXdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);
    qinv22 = kron(Mesh.dYdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv12 = kron(Mesh.dXdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv21 = kron(Mesh.dYdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);

    Mesh.Qinv(:,3*(i-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];

%     surf(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),reshape(Mesh.dYdXi(:,i),N+1,N+1))
%     hold on
end


tblr = zeros(4,4);
tblr(1,4) = 1; tblr(2,1) = 1; tblr(3,1) = 1; tblr(4,3) = 1;

theta = [ -3/4*pi     -pi
          -3/4*pi   -pi/2 
            -pi/2   -pi/4
            -pi/4       0 ];

k=0;
for i=2:5

    k=k+1;
    cor = [ xglobal(corners(i,:))
            yglobal(corners(i,:)) ];

[Mesh.X(:,i),Mesh.Y(:,i),Mesh.dXdXi(:,i),Mesh.dXdEta(:,i),Mesh.dYdXi(:,i),Mesh.dYdEta(:,i)] = ...
                transfinitemapping_cylinder(cor,R,theta(k,:),tblr(k,:),N+1,xi);
            
%     pcolor(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),(i-1)/5*ones(N+1))


    Mesh.J(:,i) = Mesh.dXdXi(:,i).*Mesh.dYdEta(:,i)-Mesh.dXdEta(:,i).*Mesh.dYdXi(:,i);

    qinv11 = kron(Mesh.dXdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);
    qinv22 = kron(Mesh.dYdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv12 = kron(Mesh.dXdEta(:,i)./Mesh.J(:,i),[0 ; 1]);
    qinv21 = kron(Mesh.dYdXi(:,i)./Mesh.J(:,i) ,[1 ; 0]);

    Mesh.Qinv(:,3*(i-1)+(1:3)) = [qinv21 qinv11+qinv22 qinv12];
    
%     surf(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),reshape(Mesh.dYdXi(:,i),N+1,N+1))
%     hold on

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cylinder

% Rc = R;
% xc = linspace(-R,R,100); % coordinates of cilinder
% yc = -sqrt(Rc^2-xc.^2);
% 
% plot(xc,yc,'k','linewidth',4)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration variables for force and boundary conditions

global n xibar Gw

n = 8;

[xibar,Gw] = GLnodes(n);
xibar = xibar'; Gw = Gw';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%