clear all
close all
clc

if ispc; figure('windowstyle','docked'); else figure; end

global N numRows numColumns
global xi
global h e
global nr_1 nr_2
global globalnr_1v globalnr_1h globalnr_2

N = 5;
H = 6;
T = 100;
dt = 0.5;

numRows    = H;
numColumns = H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

numbering('multi');

Mesh = meshgenerator('SinDeformGrid',0);
X = Mesh.X;
Y = Mesh.Y;

[Meshp,hp,ep] = postproces_grid('SinDeformGrid',0);
Xp = Meshp.X;
Yp = Meshp.Y;

subplot(1,2,1)
for i=1:N*H+1
    plot([ X(i,1) X(i,N*H+1) ],[Y(i,1) Y(i,N*H+1)],'linewidth',1)
    hold on
    plot([ X(1,i) X(N*H+1,i) ],[Y(1,i) Y(N*H+1,i)],'linewidth',1)
end
axis equal
axis([min(min(X)) max(max(X)) min(min(Y)) max(max(Y))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

lambda = 1/8;
x0 = -1/2;
y0 = 1/2;

% x0 = 0;
% y0 = +1/2;

% x0 = 0;
% y0 = 0;



% Gaussian function
% u0 = exp(-((X-x0).^2+(Y-y0).^2)/(2*lambda^2));

U0 = zeros(N*numColumns,N*numRows);
for i=1:N*numColumns
    for j=1:N*numRows
        U0(i,j) = pi/2*lambda^2*( ...
( erf( (X(i+1,j)-x0)/(sqrt(2)*lambda) ) - erf( (X(i,j)-x0)/(sqrt(2)*lambda) ) ) * ...
( erf( (Y(i,j+1)-y0)/(sqrt(2)*lambda) ) - erf( (Y(i,j)-y0)/(sqrt(2)*lambda) ) ) );
    end
end

UU0 = reconstruct_twoforms(U0,ep,1./Meshp.J);

subplot(1,2,2)
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
%         surf(Xp(:,:,rc),Yp(:,:,rc),UU0(:,:,rc))
        pcolor(Xp(:,:,rc),Yp(:,:,rc),UU0(:,:,rc))
        hold on
    end
end
shading interp
% colorbar('EastOutside')
axis equal
% axis([-1 1 -1 1 0 1])
axis([-1 1 -1 1])
% set(gca,'Ztick',[0 1])
set(gca,'clim',[0 1])
pause(0.05)



% Ax =  pi*(Y(:,2:N*H+1).^2-Y(:,1:N*H).^2);
% Ay = -pi*(X(2:N*H+1,:).^2-X(1:N*H,:).^2);

Ax = ones(N*H+1,N*H);
Ay = ones(N*H,N*H+1);

% Ax = pi*(Y(:,2:N*H+1)-Y(:,1:N*H));
% Ay = ones(N*H,N*H+1);

% Ax = zeros(N*H+1,N*H);
% Ay = pi*(X(2:N*H+1,:)-X(1:N*H,:));

% Ax = ones(N*H+1,N*H);
% Ay = ones(N*H,N*H+1);

% Ax = pi*(Y(:,2:N*H+1)-Y(:,1:N*H));
% Ay = pi*(X(2:N*H+1,:)-X(1:N*H,:));

% AX = reshape(Ax ,N*H*(N*H+1),1)*ones(1,(N*H)^2);
% AY = reshape(Ay',N*H*(N*H+1),1)*ones(1,(N*H)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization



De  = div(N);
Exe = kron(eye(N),e');
Eye = zeros(N*(N+1),N*N);
for i=1:N
    Eye(1:N*(N+1),(i-1)*N+(1:N)) = kron(eye(N),e(i,:)');
end

D   = zeros(nr_2,nr_1);
M2  = zeros(nr_2);
vEx  = zeros(nr_1,nr_2);
vEy  = zeros(nr_1,nr_2);

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

ind1 = [ reshape(globalnr_1v((c-1)*N+(1:N+1),(r-1)*N+(1:N)),N*(N+1),1)
         reshape(globalnr_1h((c-1)*N+(1:N),(r-1)*N+(1:N+1))',N*(N+1),1) ];
ind2 = globalnr_2((c-1)*N+(1:N),(r-1)*N+(1:N));

% Divergence operator
D(ind2,ind1) = De;

% two-forms
M2e = innerproduct(2,Mesh.J(:,rc));
M2(ind2,ind2) = M2e;

% Ex operator
Axe = zeros(N*(N+1),N*N);
for n=1:N
    Axe((1:N+1)+(n-1)*(N+1),(1:N)+(n-1)*N) = repmat(Ax((1:N+1)+(c-1)*N,n+(r-1)*N),1,N);
end

ind = reshape(globalnr_1v((c-1)*N+(1:N+1),(r-1)*N+(1:N)),N*(N+1),1);
vEx(ind,ind2) = vEx(ind,ind2) + Axe.*Exe;


% Ey operator
Ayee = zeros(N*(N+1),N);
for n=1:N
    Ayee((n-1)*(N+1)+(1:N+1),:) = [ zeros(N+1,n-1) Ay(n+(r-1)*N,(1:N+1)+(c-1)*N)' zeros(N+1,N-n) ];
end
Aye = repmat(Ayee,1,N);

ind = reshape(globalnr_1h((c-1)*N+(1:N),(r-1)*N+(1:N+1))',N*(N+1),1);
vEy(ind,ind2) = vEy(ind,ind2) + Aye.*Eye;

    end
end

vE = vEx+vEy;


% Periodic BC
ind1 = [ globalnr_1v(    1,:) globalnr_1h(:,1)'     ];
ind2 = [ globalnr_1v(N*H+1,:) globalnr_1h(:,N*H+1)' ];

D(:,ind1) = D(:,ind1)+D(:,ind2);
D(:,ind2) = [];

vE(ind1,:) = ( vE(ind1,:)+vE(ind2,:) )/2; % central discretization
vE(ind2,:) = [];

% Vertical cell-interfaces
ind = 0;
for r=1:numRows
    for c=1:numColumns-1
        rc = c+(r-1)*numColumns;

        ind = ind + (c==1)*( (r>1)*(N*(N-1)+N*N) + (r==2)*N );
        
        ind1 = ind + ((N+(c==1)):(N+(c==1)):N*(N+(c==1)));

        vE(ind1,:) = vE(ind1,:)/2; % central discretization

        ind = ind+(N*N+N*(N-1+(r<numRows)))+(r==1)*N+(c==1)*N;

    end
end


% Horizontal cell-interfaces
ind = 0;
for r=1:numRows-1
    for c=1:numColumns

        ind = ind + N*(N-1+(c<numColumns))+(c==1)*N;
        
        ind1 = ind + ((N+(r==1)):(N+(r==1)):N*(N+(r==1)));

        vE(ind1,:) = vE(ind1,:)/2; % central discretization

        ind = ind+(N*(N-1+(r<numRows)))+(r==1)*N;

    end
end


D  = sparse(D);
M2 = sparse(M2);
vE = sparse(vE);


A = M2*D*vE;

K1 = (A-A')/2;  % Skew-symmetric
% K1 = A;         % Standard
K2 = M2;



U = zeros(N*numColumns*N*numRows,1);
U(globalnr_2) = U0;

t = 0;
while t<T
    
    t = t+dt

    U = timemarching(K2,K1,dt,U,'ESDIRK5');

UU = reconstruct_twoforms(U(globalnr_2),ep,1./Meshp.J);

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
%         surf(Xp(:,:,rc),Yp(:,:,rc),UU(:,:,rc))
        pcolor(Xp(:,:,rc),Yp(:,:,rc),UU(:,:,rc))
        hold on
    end
end
shading interp
colorbar('EastOutside')
axis equal
% axis([-1 1 -1 1 0 1])
axis([-1 1 -1 1])
% set(gca,'Ztick',[0 1])
set(gca,'clim',[0 1])
pause(0.1)
hold off




end
