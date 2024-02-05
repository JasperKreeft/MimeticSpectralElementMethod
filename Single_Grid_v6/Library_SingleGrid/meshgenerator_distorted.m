function Mesh = meshgenerator_distorted(R)
% R = distortion number
global N numElements
global xi w

[xi,w] = GLLnodes(N);  % Gauss-Lobotto-Legendre

H = sqrt(numElements);

x = linspace(0,pi,H+1);

X = x'*ones(1,H+1);

Y = X';

Y(1:2:H+1,2:2:H+1) = Y(1:2:H+1,2:2:H+1)-R*pi/H;
Y(2:2:H+1,2:2:H+1) = Y(2:2:H+1,2:2:H+1)+R*pi/H;

xglobal = reshape(X,1,[]);
yglobal = reshape(Y,1,[]);


% plot(xglobal,yglobal,'o','markerface','b')
% hold on
% axis equal
% axis([ -0.5 3.5 -0.5 3.5])

number_corners = reshape(1:(H+1)^2,H+1,H+1);
corners = zeros(numElements,4);
for i=1:H
    for j=1:H
        ij = i+(j-1)*H;
corners(ij,:) = [ number_corners([i i+1],j) ; number_corners([i+1 i],j+1) ];
    end
end


% Global Mesh
Mesh.X      = zeros((N+1)^2,numElements);
Mesh.Y      = zeros((N+1)^2,numElements);
Mesh.dXdXi  = zeros((N+1)^2,numElements);
Mesh.dXdEta = zeros((N+1)^2,numElements);
Mesh.dYdXi  = zeros((N+1)^2,numElements);
Mesh.dYdEta = zeros((N+1)^2,numElements);
Mesh.J      = zeros((N+1)^2,numElements);
Mesh.Qinv   = zeros(2*(N+1)^2,3*numElements);

for i=1:numElements;

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Integration variables for force and boundary conditions

global n xibar Gw

n = 8;

[xibar,Gw] = GLnodes(n);
xibar = xibar'; Gw = Gw';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%