% Zero-forms, Outer-oriented, constant vector-field

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_Convection/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N nn numElements numRows numColumns
global globalnr_0 globalnr_1h globalnr_1v
global nr_0 nr_1
global w h e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'Vorticity';%'CurlCurl';%'SinDeformGrid';%
DomInfo = 3;%1.2;

InitVortex = 'taylor';%'gaussian';%'constant';%

velofield = 'circular';%'rudmann';%'uniform';%

plot_fig  = 0;
avimatlab = 0;
Tecplot   = 0;

N = 7;
H = 9;
T = 0.2;
dt = .05;

filename = ['Zero_Vort_ME_N' num2str(N) '_dt' num2str(dt)];

numRows    = H;
numColumns = H;
numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

numbering('square');

Mesh = meshgenerator_square(Domain,DomInfo);

[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize figure and plot Mesh

if ispc && (plot_fig || avimatlab)
    fig = figure('windowstyle','docked');
elseif  plot_fig || avimatlab
    fig = figure('Units','normalized','Position',[0, 0, 0.5, 1]);
end
if plot_fig
    set(gcf,'visible','on')
else
    set(gcf,'visible','off')
end

if plot_fig || avimatlab
    subplot(1,2,1)
    meshplot
    axis([min(min(Mesh.X)) max(max(Mesh.X)) min(min(Mesh.Y)) max(max(Mesh.Y))])
    subplot(1,2,2)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

t = 0; j=1;

disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])


W0 = InitialVortex(InitVortex,Mesh);

WW0(globalnr_0) = W0; WW0 = WW0';

Wp = reconstruct(0,W0,hp);

Clim = [min(min(Wp)) max(max(Wp))];
posten_vorticity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

Exe = zeros((N+1)^2,N*(N+1));
for i=1:N
    Exe(:,(i-1)*(N+1)+(1:N+1)) = kron(eye(N+1),e(i,:)');
end
Exe = sparse(Exe);
Eye = kron(e',speye(N+1));


% Matrix assembly
NG = normalgradient_assembly();
M0 = innerproduct_assembly(0,Mesh);
M1 = innerproduct_assembly(1,Mesh);

PoissonMatrix = -NG'*M1*NG;

[PM_L,PM_U] = lu(PoissonMatrix);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary stuff

ind0B = unique(reshape(globalnr_0(1:N+1,1:numColumns),[],1));
ind0T = unique(reshape(globalnr_0(N*(N+1)+(1:N+1),(numRows-1)*numColumns+1:numElements),[],1));
ind0L = unique(reshape(globalnr_0(1:N+1:(N+1)^2,1:numColumns:numElements),[],1));
ind0R = unique(reshape(globalnr_0(N+1:N+1:(N+1)^2,numColumns:numColumns:numElements),[],1));

boundary_points = unique([ ind0B ind0T ind0L ind0R ]);
interior_points = 1:nr_0; interior_points(boundary_points) = [];

ind1L = unique(reshape(globalnr_1v(1:(N+1):N*(N+1),1:numColumns:numElements),[],1));
ind1R = unique(reshape(globalnr_1v(N+1:(N+1):N*(N+1),numColumns:numColumns:numElements),[],1));
ind1B = unique(reshape(globalnr_1h(1:N+1:N*(N+1),1:numColumns),[],1));
ind1T = unique(reshape(globalnr_1h(N+1:N+1:N*(N+1),(numRows-1)*numColumns+1:numElements),[],1));

% Boundary condition for stream function
PSIbc = zeros(length(boundary_points),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Velocity field

% [U,V] = VelocityField(velofield,Mesh);

% [U,V,PSI] = VorticityInducedVelocityField(WW0,PSIbc,PoissonMatrix,NG,M0,boundary_points,interior_points);

disp('LU decomposition Poisson matrix')
[PM_L,PM_U] = lu(PoissonMatrix);

clear PoissonMatrix

[U,V,PSI] = VorticityInducedVelocityField(WW0,M0,PM_L,PM_U,NG);%,M0,boundary_points,interior_points);

[~,~,velo] = reconstruct(1, U,V,hp,ep,Meshp);

Vorticity = zeros(ceil(T/dt)+1,1);
Energy    = zeros(ceil(T/dt)+1,1);
Enstrophy = zeros(ceil(T/dt)+1,1);
for i=1:numElements
Vorticity(1) = Vorticity(1) + sum(Wp(:,i).*Meshp.J(:,i).*Meshp.W);
Energy(1)    = Energy(1) + 1/2*sum(velo(:,i).*Meshp.J(:,i).*Meshp.W.*velo(:,i));
Enstrophy(1) = Enstrophy(1) + 1/2*sum(Wp(:,i).*Meshp.J(:,i).*Meshp.W.*Wp(:,i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interior Product matrix

he = kron(e,h)';
eh = zeros((N+1)^2,N*(N+1));
for i=1:N
    for j=1:N+1
        ij = j+(i-1)*(N+1);
        eh(:,ij) = kron(h(j,:),e(i,:));
    end
end

Exe = zeros((N+1)^2,N*(N+1));
for i=1:N
    Exe(:,(i-1)*(N+1)+(1:N+1)) = kron(eye(N+1),e(i,:)');
end
Exe = sparse(Exe);
Eye = kron(e',speye(N+1));

vE01 = interiorProduct(U,V,Exe,Eye,he,eh,Mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Periodic Boundaries

M0d = diag(M0);
NGp = NG;

% Periodic Boundary conditions
NGp([ind1L;ind1B],:) = NGp([ind1L;ind1B],:)+NGp([ind1R;ind1T],:);
NGp([ind1R;ind1T],:) = [];

NGp(:,ind0B) = NGp(:,ind0B)+NGp(:,ind0T);
NGp(:,ind0L) = NGp(:,ind0L)+NGp(:,ind0R);
NGp(:,[ind0R;ind0T]) = [];

vE01(:,[ind1L;ind1B]) = (vE01(:,[ind1L;ind1B])+vE01(:,[ind1R;ind1T]))/2; % central discretization
vE01(:,[ind1R;ind1T]) = [];

vE01(ind0B,:) = (vE01(ind0B,:)+vE01(ind0T,:));
vE01(ind0L,:) = (vE01(ind0L,:)+vE01(ind0R,:));
vE01([ind0R;ind0T],:) = [];

% Inner product zero-forms
M0d(ind0B) = M0d(ind0B) + M0d(ind0T);
M0d(ind0L) = M0d(ind0L) + M0d(ind0R);
M0d([ind0R;ind0T]) = [];

ind = length(M0d);
M0p = spdiags(M0d,0,ind,ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convection matrices

A = M0p*vE01*NGp;

K1 = (A-A')/2;  % Skew-symmetric
K2 = M0p;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indw3 = 1:nr_0; indw3([ind0R;ind0T]) = [];

W = zeros(nr_0,1);
W(globalnr_0) = W0;
W([ind0R;ind0T]) = [];

t=0; j=1;
WW = zeros(nr_0,1); Wold = W;
while t<T
    j = j+1;
    t = t+dt;

    disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

    %%%
error = 1; er = 0;
while error>1e-5 && er<30 %0.9
 er = er+1;
 if er==20; disp('er=20 !!!'); end
 
    W_ = timemarching(K2,K1,dt,W,'gauss4');

    WW(indw3) = W_;
    WW(ind0R) = WW(ind0L);
    WW(ind0T) = WW(ind0B);

    %%%

% [U,V,PSI] = VorticityInducedVelocityField(WW,PSIbc,PoissonMatrix,NG,M0,boundary_points,interior_points);
[U,V,PSI] = VorticityInducedVelocityField(WW,M0,PM_L,PM_U,NG);
    %%%

    vE01 = interiorProduct(U,V,Exe,Eye,he,eh,Mesh);

    % Periodic Boundary conditions
    vE01(:,[ind1L;ind1B]) = (vE01(:,[ind1L;ind1B])+vE01(:,[ind1R;ind1T]))/2; % central discretization
    vE01(:,[ind1R;ind1T]) = [];

    vE01(ind0B,:) = (vE01(ind0B,:)+vE01(ind0T,:));
    vE01(ind0L,:) = (vE01(ind0L,:)+vE01(ind0R,:));
    vE01([ind0R;ind0T],:) = [];

    %%%

    A = M0p*vE01*NGp;

    K1 = (A-A')/2;  % Skew-symmetric

    %%%

    error = max(abs(W_-Wold));
    Wold = W_;

    disp(['er = ' num2str(er) ', error = ' num2str(error)])
end

W = W_;

    %%%
    Wp = reconstruct(0,WW(globalnr_0),hp);
    [~,~,velo] = reconstruct(1,U,V,hp,ep,Meshp);
    
    for i=1:numElements
    Vorticity(j) = Vorticity(j) + sum(Wp(:,i).*Meshp.J(:,i).*Meshp.W);
    Energy(j)    = Energy(j) + 1/2*sum(velo(:,i).*Meshp.J(:,i).*Meshp.W.*velo(:,i));
    Enstrophy(j) = Enstrophy(j) + 1/2*sum(Wp(:,i).*Meshp.J(:,i).*Meshp.W.*Wp(:,i));
    end

    posten_vorticity

end

plot(Energy,'-o')
title('Energy')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
run Library_Convection/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%