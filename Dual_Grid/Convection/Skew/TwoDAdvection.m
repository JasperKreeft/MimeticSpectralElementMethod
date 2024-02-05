clear all
close all
clc

if ispc; figure('windowstyle','docked'); else figure; end

global N w

N = 20;
T = 10;
dt = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

x = xi; y = x;
X = x'*ones(1,N+1);
Y = ones(N+1,1)*y;

subplot(1,2,1)
for i=1:N+1
    plot([ X(i,1) X(i,N+1) ],[Y(i,1) Y(i,N+1)],'linewidth',1)
    hold on
    plot([ X(1,i) X(N+1,i) ],[Y(1,i) Y(N+1,i)],'linewidth',1)
end
axis equal
axis([min(min(X)) max(max(X)) min(min(Y)) max(max(Y))])

nn = 100;
xx = linspace(-1,1,nn); yy = xx;
[hh,ee] = MimeticpolyVal(xx,N,1);
XX = xx'*ones(1,nn);
YY = ones(nn,1)*yy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

lambda = 1/8;
x0 = -1/2;
y0 = 1/2;

% x0 = 0;
% y0 = -1/2;


% Gaussian function
% u0 = exp(-((X-x0).^2+(Y-y0).^2)/(2*lambda^2));

U0 = zeros(N);
for i=1:N
    for j=1:N
        U0(i,j) = pi/2*lambda^2*( ...
( erf( (x(i+1)-x0)/(sqrt(2)*lambda) ) - erf( (x(i)-x0)/(sqrt(2)*lambda) ) ) * ...
( erf( (y(j+1)-y0)/(sqrt(2)*lambda) ) - erf( (y(j)-y0)/(sqrt(2)*lambda) ) ) );
    end
end

UU0 = ee'*U0*ee;

subplot(1,2,2)
% surf(XX,YY,UU0)
pcolor(XX,YY,UU0)
shading interp
% colorbar
axis equal
% axis([-1 1 -1 1 0 1])
axis([-1 1 -1 1])
% set(gca,'Ztick',[0 1])
set(gca,'clim',[0 1])
pause(0.05)

% surf(X,Y,U0)
% shading interp
% colorbar
% axis equal
% axis([-1 1 -1 1 0 1])
% set(gca,'Ztick',[0 1])
% set(gca,'clim',[0 1])
% pause(0.05)


% Ax =  pi*(Y(:,2:N+1).^2-Y(:,1:N).^2);
% Ay = -pi*(X(2:N+1,:).^2-X(1:N,:).^2);

% Ax = pi*(Y(:,2:N+1)-Y(:,1:N));
% Ay = zeros(N,N+1);

% Ax =  zeros(N+1,N);
% Ay = pi*(X(2:N+1,:)-X(1:N,:));

Ax = ones(N+1,N);
Ay = ones(N,N+1);

% Ex operator
Axe = zeros(N*(N+1),N*N);
for n=1:N
    AX((1:N+1)+(n-1)*(N+1),(1:N)+(n-1)*N) = repmat(Ax((1:N+1),n),1,N);
end

% Ey operator
AYe = zeros(N*(N+1),N);
for n=1:N
    AYe((n-1)*(N+1)+(1:N+1),:) = [ zeros(N+1,n-1) Ay(n,1:N+1)' zeros(N+1,N-n) ];
end
AY = repmat(AYe,1,N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization


Ex = kron(eye(N),e');

Ey = zeros(N*(N+1),N*N);
for i=1:N
    Ey(1:N*(N+1),(i-1)*N+(1:N)) = kron(eye(N),e(i,:)');
end

vE = [ AX.*Ex; AY.*Ey];

D = div(N);

% Periodic BC
for i=2*N:-1:1
    ind1 = (i-1)*(N+1)+1;
    ind2 = i*(N+1);
    
    vE(ind1,:) = (vE(ind1,:)+vE(ind2,:))/2;  % central discretization
    vE(ind2,:) = [];

    D(:,ind1) = D(:,ind1)+D(:,ind2);
    D(:,ind2) = [];
    
end


% Inner product two-forms
J = ones((N+1)^2,1);
W = kron(w,w)';
JW = spdiags(J.*W,0,(N+1)^2,(N+1)^2);
I = kron(e,e)';
M2 = I'*JW*I;



A = M2*D*vE;

K1 = (A-A')/2;  % Skew-symmetric
% K1 = A;         % Standard
K2 = M2;



U = reshape(U0,N*N,1);

t = 0;
while t<T
    
    t = t+dt;
    
    U = timemarching(K2,K1,dt,U,'ESDIRK5');
    
    UU = ee'*reshape(U,N,N)*ee;
    
%     surf(XX,YY,UU)
    pcolor(XX,YY,UU)
    shading interp
    colorbar
    axis equal
%     axis([-1 1 -1 1 -0.2 1.2])
    axis([-1 1 -1 1])
%     set(gca,'Ztick',[0 1])
    set(gca,'clim',[0 1])
    pause(0.05)


end
