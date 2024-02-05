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

global N
global w e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'Vorticity';%'SinDeformGrid';%'CurlCurl';%
DomInfo = 2;%0%1.2;

InitVortex = 'taylor';%'gaussian';%'constant';%

velofield = 'circular';%'rudmann';%'uniform';%

plot_fig  = 1;
avimatlab = 0;
Tecplot   = 0;

N = 30;
T = 2;
dt = 0.05;

filename = ['Zero_Vort_SE_N' num2str(N) '_dt' num2str(dt)];

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

Wp = reconstruct(0,W0,hp);

Enstrophy = zeros(ceil(T/dt)+1,1);
Enstrophy(1) = 1/2*sum(Wp.*Meshp.J.*Meshp.W.*Wp);

Clim = [min(min(W0)) max(max(W0))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Velocity field

M0 = innerproduct(0,Mesh.J);
M1 = innerproduct(1,Mesh.J,Mesh.Qinv);
NG = normalgrad(N);

disp('Assemble Poisson matrix')
PoissonMatrix = -NG'*M1*NG;

disp('LU decomposition Poisson matrix')
[PM_L,PM_U] = lu(PoissonMatrix+1e-10*speye(size(PoissonMatrix)));

[U,V] = VorticityInducedVelocityField(W0,M0,PM_L,PM_U,NG);

[uu,vv,Velo] = reconstruct(1,U,V,hp,ep,Meshp);

posten_vorticity

if plot_fig
subplot(1,2,1)
hold on
quiver(Meshp.X,Meshp.Y,uu,vv,'w')
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

W = W0;

Wold = W;
while t<T
    j = j+1;
    t = t+dt;

    disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

    %%%
    error = 1; er = 0;
    while error>1e-3 && er<10
        er = er+1
        
        if er==10; disp('er=10 !!!'); end
        
        WW = timemarching(K2,K1,dt,W,'gauss4');

        %%%

        [U,V] = VorticityInducedVelocityField(WW,M0,PM_L,PM_U,NG);

        vE01 = interiorProduct(U,V,Exe,Eye,he,eh,Mesh);
        %%%

        disp('Assemble convection matrix')
        A = M0*vE01*NG;

        K1 = (A-A')/2;  % Skew-symmetric
        % K1 = A;         % Standard
        
        error = max(abs(WW-Wold))
        Wold = WW;
        
    end

    %%%

    W = WW;
    
    disp('Postprocessing')
    Wp = reconstruct(0,W,hp);
    [~,~,Velo] = reconstruct(1,U,V,hp,ep,Meshp);
    
    posten_vorticity
    
    Enstrophy(j) = 1/2*sum(Wp.*Meshp.J.*Meshp.W.*Wp);

end

figure
plot(Enstrophy,'-o')
title('Enstrophy')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
run Library_Convection/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%