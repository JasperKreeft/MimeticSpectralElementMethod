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

global N nr_1 nr_2
global w e
global Re

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'SinDeformGrid_01';
DomInfo = 0.0;

plot_fig  = 0;
avimatlab = 0;
Tecplot   = 0;

Re = 1;
N = 30;
T = 2000;
dt = 1;

numColumns = 1;
numRows = 1;
numElements = 1;

filename = ['Zero_Vort_SE_N' num2str(N) '_dt' num2str(dt)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

numbering('square');

Mesh = meshgenerator_square(Domain,DomInfo);

[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

[psi_ex,u_ex,v_ex,w_ex,p_ex,Fy_ex] = ReguralizedLDC(Re,Mesh.X,Mesh.Y);
Energy_exact = 3.676492819329155e-02;
Enstrophy_exact = 2.01433106575958;

Clim = [min(min(w_ex)) max(max(w_ex))];

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

[PSIbc,Wbc,boundary_points,interior_points] = boundaryconditions_0_square_NS(Re,Mesh);

t = 0; j=1;

disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

W0 = Wbc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Velocity field

M0 = innerproduct(0,Mesh.J);
M1 = innerproduct(1,Mesh.J,Mesh.Qinv);
M2 = innerproduct(2,Mesh.J);
NG = normalgrad(N);
D = div(N);

disp('Assemble Poisson matrix')
MatrixFull = -NG'*M1*NG;

Matrixrhs = MatrixFull(:,boundary_points);

PoissonMatrix = MatrixFull;
PoissonMatrix(:,boundary_points) = [];
PoissonMatrix(boundary_points,:) = [];

disp('LU decomposition Poisson matrix')
[PM_L,PM_U] = lu(PoissonMatrix);

[U,V,PSI] = VorticityInducedVelocityField_old(W0,PSIbc,M0,Matrixrhs,PM_L,PM_U,NG,boundary_points,interior_points);

Wp = reconstruct(0,W0,hp);
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
%% Forcing term

[fx fy] = forcefunction(1,1,1,'reguralizedLDC',Domain,DomInfo);

RHS = -NG'*M1*[ fx ; fy ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convection matrices

disp('Assemble convection matrix')
A = M0*vE01*NG;

K1 = M0;
K2 = (A-A')/2 - 1/Re*MatrixFull;

K1_ = K1;
K1_(boundary_points,:) = [];
K1_(:,boundary_points) = [];
K2_ = K2;
K2_(boundary_points,:) = [];
Mbc = K2_;
K2_(:,boundary_points) = [];

RHS(boundary_points) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 0; j=1;

W_ = W0;
W_(boundary_points) = [];

Wold = W0;
error = 1;
while t<T && error>1e-12
    j = j+1;
    t = t+dt;

    disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

    %%%
    error = 1; er = 0;
%     while error>1e-10 && er<100
        er = er+1;
        
%         if er==100; disp('er=100 !!!'); end
        
        W_ = timemarching_bc(K1_,K2_,dt,W_,Mbc,Wbc,RHS,'BE');

        WW = Wbc;
        WW(interior_points) = W_;

        %%%

[U,V] = VorticityInducedVelocityField_old(WW,PSIbc,M0,Matrixrhs,PM_L,PM_U,NG,boundary_points,interior_points);

        vE01 = interiorProduct(U,V,Exe,Eye,he,eh,Mesh);
        %%%

        disp('Assemble convection matrix')
        A = M0*vE01*NG;

        K2 = (A-A')/2 - 1/Re*MatrixFull;  % Skew-symmetric
        K2_ = K2;
        K2_(boundary_points,:) = [];
        Mbc = K2_;
        K2_(:,boundary_points) = [];
        
        error = max(abs(WW-Wold));
        
        disp(['er = ' num2str(er) ', error = ' num2str(error)])
        
        Wold = WW;
        
%     end

    %%%

    W = WW;
    
    disp('Postprocessing')
    Wp = reconstruct(0,W,hp);
    [Up,Vp,Velo] = reconstruct(1,U,V,hp,ep,Meshp);
    posten_vorticity
    
    Energy(j) = 1/2*( sum(Up.*Meshp.J.*Meshp.W.*Up) ...
                    + sum(Vp.*Meshp.J.*Meshp.W.*Vp) );
    Enstrophy(j) = 1/2*sum(Wp.*Meshp.J.*Meshp.W.*Wp);

end

figure
% plot(Enstrophy,'-o')
% hold on
% plot(Energy,'-r^')
semilogy(abs(Enstrophy-Enstrophy_exact),'-o')
hold on
semilogy(abs(Energy-Energy_exact),'-^r')
title('Kinetic Energy & Enstrophy')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pressure

[Fx Fy] = forcefunction(1,1,1,'reguralizedLDC_p',Domain,DomInfo);
[~,~,boundary_flux,interior_flux] = boundaryconditions_1_square('reguralizedLDC',Domain,DomInfo,[1 0 0 0]);

RHS = [zeros(nr_1,1) ; M2*D*[ Fx ; Fy ]];

Matrix = [ -M1 D'*M2
           M2*D spalloc(nr_2,nr_2,0) ];

% Matrix = [ -M1 D'*M2
%            M2*D speye(nr_2) ];

Matrix(:,boundary_flux) = [];
Matrix(boundary_flux,:) = [];

RHS(boundary_flux) = [];

PW = Matrix\RHS;

ind = nr_1-length(boundary_flux);
P = PW(ind+1:end);


pp = reconstruct(2,P,ep,Meshp) + Energy(end);

figure
Xp = reshape(Meshp.X(:,1),nn,nn);
Yp = reshape(Meshp.Y(:,1),nn,nn);
pcolor(Xp,Yp,reshape(pp,nn,nn))
shading interp
axis equal
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

% in = 'finish';                                                  %#ok<NASGU>
% run Library_Convection/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%