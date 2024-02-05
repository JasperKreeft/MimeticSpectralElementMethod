clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N nn numElements numRows numColumns
global xi
global h e
global nr_1 nr_2
global globalnr_1v globalnr_1h globalnr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'SinDeformGrid';
DomInfo = 0.0;

N = 12;
H = 1;
T = 100;
dt = 0.5;

numRows    = H;
numColumns = H;
numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N); eta = xi;

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

numbering('square')

Mesh = meshgenerator_square(Domain,DomInfo);

[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

subplot(1,2,1)
meshplot
axis([min(min(Mesh.X)) max(max(Mesh.X)) min(min(Mesh.Y)) max(max(Mesh.Y))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

U0 = zeros(N*N,numElements);
for c=1:numColumns
    for r=1:numRows
        rc = c+(r-1)*numColumns;
        U0(:,rc) = initialforms(2,r,c,'gaussian',Domain,DomInfo);
    end
end

Up = reconstruct(2,U0,ep,Meshp);

subplot(1,2,2)
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
        
%         surf(reshape(Meshp.X(:,rc),nn,nn),reshape(Meshp.Y(:,rc),nn,nn),reshape(Up(:,rc),nn,nn))
        pcolor(reshape(Meshp.X(:,rc),nn,nn),reshape(Meshp.Y(:,rc),nn,nn),reshape(Up(:,rc),nn,nn))
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

x0 = -1/2; y0 = 0; lambda = 1/8; % MOET ZELFDE ZIJN ALS IN initialforms
Up_ex = exp(-((Meshp.X-x0).^2+(Meshp.Y-y0).^2)/(2*lambda^2));

L2error = zeros(1,ceil(T/dt)+1);
L2norm = 0;
Energy = zeros(1,ceil(T/dt)+1);
for i=1:numElements
    L2error(1) = L2error(1) + sqrt(sum(Meshp.W.*Meshp.J(:,i).*(Up(:,i)-Up_ex(:,i)).^2));
    L2norm = L2norm + sqrt(sum(Meshp.W.*Meshp.J(:,i).*Up_ex(:,i).^2));
    Energy(1) = Energy(1) + 1/2*sum(Up(:,i).*Meshp.W.*Meshp.J(:,i).*Up(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization

De  = div(N);
Exe = kron(ones(N),e');
Eye = zeros(N*(N+1),N*N);
for i=1:N
    Eye(1:N*(N+1),(i-1)*N+(1:N)) = kron(ones(N),e(i,:)');
end

D   = zeros(nr_2,nr_1);
M2  = zeros(nr_2);
vEx  = zeros(nr_1,nr_2);
vEy  = zeros(nr_1,nr_2);

for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

ind1 = [ globalnr_1v(:,rc) ; globalnr_1h(:,rc) ];
ind2 = globalnr_2(:,rc);

% Divergence operator
D(ind2,ind1) = De;

% two-forms
M2e = innerproduct(2,Mesh.J(:,rc));
M2(ind2,ind2) = M2e;

% Ex operator
% Axe = zeros(N*(N+1),N*N);
% for n=1:N
%     Axe((1:N+1)+(n-1)*(N+1),(1:N)+(n-1)*N) = repmat(Ax((1:N+1)+(c-1)*N,n+(r-1)*N),1,N);
% end
Axe = 1;

ind = globalnr_1v(:,rc);
vEx(ind,ind2) = vEx(ind,ind2) + Axe.*Exe;


% Ey operator
% Ayee = zeros(N*(N+1),N);
% for n=1:N
%     Ayee((n-1)*(N+1)+(1:N+1),:) = [ zeros(N+1,n-1) Ay(n+(r-1)*N,(1:N+1)+(c-1)*N)' zeros(N+1,N-n) ];
% end
% Aye = repmat(Ayee,1,N);
Aye = 1;

ind = globalnr_1h(:,rc);
vEy(ind,ind2) = vEy(ind,ind2) + Aye.*Eye;

    end
end

vE = vEx+vEy;






% NIET AF !!!!!!!!!










WPe = zeros(N);
for i=1:N
    xp = (xi(i+1)+xi(i))/2+(xi(i+1)-xi(i))/2*xi;
    [~,ep] = MimeticpolyVal(xp,N,1);
    u = -2*pi*xp;
    for k=1:N
        WPe(i,k) = sum(w.*u.*ep(k,:))*(xi(i+1)-xi(i))/2;
    end
end
WPee = kron(WPe,ones(N+1,1));
WP = repmat(WPee,1,N);

WQe = zeros(N);
for j=1:N
    yq = (eta(j+1)+eta(j))/2+(eta(j+1)-eta(j))/2*eta;
    [~,eq] = MimeticpolyVal(yq,N,1);
    v = 2*pi*yq;
    for l=1:N
        WQe(j,l) = sum(w.*v.*eq(l,:))*(eta(j+1)-eta(j))/2;
    end
end
WQ = kron(WQe,ones(N+1,N));

V = [ Ex.*WQ ; Ey.*WP ];


D = div(N);

% Periodic BC
for i=2*N:-1:1
    ind1 = (i-1)*(N+1)+1;
    ind2 = i*(N+1);

    V(ind1,:) = (V(ind1,:)+V(ind2,:))/2;  % central discretization
    V(ind2,:) = [];

    D(:,ind1) = D(:,ind1)+D(:,ind2);
    D(:,ind2) = [];

end


% Inner product two-forms
J = ones((N+1)^2,1);
W = kron(w,w)';
JW = spdiags(J.*W,0,(N+1)^2,(N+1)^2);
I = kron(e,e)';
M2 = I'*JW*I;



A = M2*D*V;

% K1 = (A-A')/2;  % Skew-symmetric
K1 = A;         % Standard
K2 = M2;



U = reshape(U0,N*N,1);

t = 0;
while t<T
    
    t = t+dt;
    
    U = timemarching(K2,K1,dt,U,'ESDIRK5');
    
    Up = ee'*reshape(U,N,N)*ee;
    
%     surf(XX,YY,Up)
    pcolor(XX,YY,Up)
    shading interp
    colorbar
    axis equal
%     axis([-1 1 -1 1 -0.2 1.2])
    axis([-1 1 -1 1])
%     set(gca,'Ztick',[0 1])
    set(gca,'clim',[0 1])
    pause(0.05)

    Xh = Meshp.X-x0*cos(2*pi*t)-y0*sin(2*pi*t);
    Yh = Meshp.Y+x0*sin(2*pi*t)-y0*cos(2*pi*t);
    Up_ex = exp(-(Xh.^2+Yh.^2)/(2*lambda^2));
    for i=1:numElements
        L2error(j) = L2error(j) + sqrt(sum(Meshp.W.*Meshp.J(:,i).*(Up(:,i)-Up_ex(:,i)).^2));
        Energy(j) = Energy(j) + 1/2*sum(Up(:,i).*Meshp.W.*Meshp.J(:,i).*Up(:,i));
    end
    
end

taxis = linspace(0,t,j);
figure
subplot(1,2,1)
plot(taxis,Energy(1:j),'-o')
subplot(1,2,2)
plot(taxis,Energy(1:j)-Energy(1),'-o')

figure
plot(taxis,L2error(1:j)/L2norm,'-o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%