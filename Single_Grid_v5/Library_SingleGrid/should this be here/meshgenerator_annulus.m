function Mesh = meshgenerator_annulus(R)

global N numElements
global xi w

[xi,w] = GLLnodes(N);  % Gauss-Lobotto-Legendre


% xglobal = [ 0 1 0 -1 0 R 0 -R ];
% yglobal = [ -1 0 1 0 -R 0 R 0 ];

a = 1/2*sqrt(2);
xglobal = [ a a -a -a R*a R*a -R*a -R*a ];
yglobal = [ -a a a -a -R*a R*a R*a -R*a ];


% plot(xglobal,yglobal,'o','markerface','b')
% hold on
% axis equal
% axis([-1.2 1.2 -1.2 1.2])


corners = [ 1 2 6 5
            2 3 7 6
            3 4 8 7
            4 1 5 8 ];


% Global Mesh
Mesh.X      = zeros((N+1)^2,numElements);
Mesh.Y      = zeros((N+1)^2,numElements);
Mesh.dXdXi  = zeros((N+1)^2,numElements);
Mesh.dXdEta = zeros((N+1)^2,numElements);
Mesh.dYdXi  = zeros((N+1)^2,numElements);
Mesh.dYdEta = zeros((N+1)^2,numElements);
Mesh.J      = zeros((N+1)^2,numElements);
Mesh.Qinv   = zeros(2*(N+1)^2,3*numElements);

% theta = [ -pi/2      0   
%               0   pi/2
%            pi/2     pi
%             -pi -1/2*pi ];
        
theta = [ -pi/4    pi/4   
           pi/4    3*pi/4
         3*pi/4    5*pi/4
         -3*pi/4   -pi/4 ];

for i=1:4
    cor = [ xglobal(corners(i,:))
            yglobal(corners(i,:)) ];

[Mesh.X(:,i),Mesh.Y(:,i),Mesh.dXdXi(:,i),Mesh.dXdEta(:,i),Mesh.dYdXi(:,i),Mesh.dYdEta(:,i)] = ...
                transfinitemapping_annulus(cor,[1 R],theta(i,:));

%     pcolor(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),(i-1)/3*ones(N+1))

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

[xibar,Gw] = Gnodes(n);
xibar = xibar'; Gw = Gw';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%