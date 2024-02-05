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
global Mesh

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'Vorticity';%'SinDeformGrid';%'CurlCurl';%
DomInfo = 2%2*pi;%0%1.2;

InitVortex = 'taylor';%'gaussian';%'constant';%

velofield = 'circular';%'rudmann';%'uniform';%

plot_fig  = 0;
avimatlab = 0;
Tecplot   = 0;

N = 10; %9
H = 3;  %13
T = 0.1;
dt = 0.01
DiffCoeff = dt^2;

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

Mesh = meshgenerator_square_smooth(Domain,DomInfo);
[Meshp,hp,ep] = postproces_grid_square_smooth(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize figure and plot Mesh

if plot_fig || avimatlab
    meshplot
    % axis([-1 1 -1 1])
    pause(0.2)
end

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

t = 0; j=1;

disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

W0 = InitialVortex(InitVortex,Mesh);

WW0(globalnr_0) = W0; WW0 = WW0';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% Matrix assembly
NG = normalgradient_assembly();
M0 = innerproduct_assembly(0,Mesh);
M1 = innerproduct_assembly(1,Mesh);

disp('Assemble Poisson matrix')
PoissonMatrix = -NG'*M1*NG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Velocity field

disp('LU decomposition Poisson matrix')
[PM_L,PM_U] = lu(PoissonMatrix+0*speye(size(PoissonMatrix)));

% W0 = W0+1e-10*rand(size(W0));
% WW0(globalnr_0,1) = forcefunction(0,Mesh.X,Mesh.Y,'sine');
[U,V] = VorticityInducedVelocityField(W0,M0,PM_L,PM_U,NG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wp = reconstruct(0,W0,hp);
[Up,Vp,Velo] = reconstruct(1,U,V,hp,ep,Meshp);

Clim = [min(min(Wp)) max(max(Wp))];
posten_vorticity

Enstrophy = zeros(ceil(T/dt)+1,1);
for i=1:numElements
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
%% Convection matrices

disp('Assemble convection matrix')
A = M0*vE01*NG;

K1 = (A-A')/2;  % Skew-symmetric
K2 = M0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0; j=1;

W = zeros(nr_0,1);
W(globalnr_0) = W0;

WW = zeros(nr_0,1); Wold = W;
while t<T
    j = j+1;
    t = t+dt;

    disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

    %%%
    error = 1; er = 0;
while error>1e-8 && er<10
    er = er+1;
 
    if er==20; disp('er=10 !!!'); end
 
    WW = timemarching(K2,K1,dt,W,'gauss4');

    [U,V] = VorticityInducedVelocityField(WW,M0,PM_L,PM_U,NG);

    vE01 = interiorProduct(U,V,Exe,Eye,he,eh,Mesh);
    
    disp('Assemble convection matrix')
    A = M0*vE01*NG;

    K1 = (A-A')/2;  % Skew-symmetric
% K1 = (A-A')/2 + -DiffCoeff*spdiags(((WW>-0.1)&(WW<0.1)),0,length(WW),length(WW))*MatrixFull;
% K1 = (A-A')/2 + -DiffCoeff*spdiags((WW<0.1),0,length(WW),length(WW))*PoissonMatrix;

    %%%

    error = max(abs(WW-Wold));
    
    disp(['er = ' num2str(er) ', error = ' num2str(error)])
    
    Wold = WW;

end

    W = WW;

    disp('Postprocessing')
    Wp = reconstruct(0,WW(globalnr_0),hp);
    [Up,Vp,Velo] = reconstruct(1,U,V,hp,ep,Meshp);
    
    posten_vorticity
    
    for i=1:numElements
    Enstrophy(j) = Enstrophy(j) + 1/2*sum(Wp(:,i).*Meshp.J(:,i).*Meshp.W.*Wp(:,i));
    end

end

plot_fig=1;
posten_vorticity

figure
plot(Enstrophy,'-o')
title('Enstrophy')

figure
plot(Energy-Energy(1),'-or')
title('Energy')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
run Library_Convection/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%