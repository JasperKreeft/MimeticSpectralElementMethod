function [U,V] = VelocityField(velofield,Mesh)

global N numElements
global globalnr_1v globalnr_1h

U = zeros(size(globalnr_1v));
V = zeros(size(globalnr_1h));

switch velofield

%% uniform velocity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'uniform'

for i=1:numElements
    U(:,i) = kron(diff(Mesh.Y(1:N+1:(N+1)^2,i)),ones(N+1,1));
    V(:,i) = kron(diff(Mesh.X(1:N+1,i)),ones(N+1,1));
end

%% circular velocity field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'circular'

for i=1:numElements
    U(:,i) = kron(+pi*(Mesh.Y(N+2:N+1:(N+1)^2,i)'.^2-Mesh.Y(1:N+1:N*(N+1),i)'.^2),ones(1,N+1))';
    V(:,i) = kron(-pi*(Mesh.X(2:N+1,i)'.^2-Mesh.X(1:N,i)'.^2),ones(1,N+1))';
end

%% Rudmann vortex %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

case 'rudmann'

% U = -cos(pi*X).*sin(pi*Y);
% V =  sin(pi*X).*cos(pi*Y);

for i=1:numElements
    U(:,i) = kron( 2*(cos(pi/2*Mesh.Y(N+2:N+1:(N+1)^2,i))-cos(pi/2*Mesh.Y(1:N+1:N*(N+1),i)))/pi , cos(pi/2*Mesh.X(1:N+1,i)) );
    V(:,i) = kron( -2*( cos(pi/2*Mesh.X(2:N+1,i)) - cos(pi/2*Mesh.X(1:N,i)) )/pi , cos(pi/2*Mesh.Y(1:N+1:(N+1)^2,i))  );
end




end