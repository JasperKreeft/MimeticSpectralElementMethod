% Zero-forms, Outer-oriented, constant vector-field

clear all
% close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_Convection/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N nn numElements numRows numColumns
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2
global w h e
global Mesh
global Re

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'SinDeformGrid_01';
DomInfo = 0.0;

plot_fig  = 0;
avimatlab = 0;
Tecplot   = 0;

NrCellRange = 6;
HconvRange = 4%[ 4 8 16];
T = 2000;
dt = 1
Re = 100;

filename = ['Zero_Vort_ME_N' num2str(N) '_dt' num2str(dt)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start H-convergence loop

conv_i = 0;
L2_u = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_v = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_w = zeros(1,max(length(NrCellRange),length(HconvRange)));
% L2_p = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_du = zeros(1,max(length(NrCellRange),length(HconvRange)));
% L2_cw1 = zeros(1,max(length(NrCellRange),length(HconvRange)));
% L2_cw2 = zeros(1,max(length(NrCellRange),length(HconvRange)));
% L1_du = zeros(1,max(length(NrCellRange),length(HconvRange)));
% Linf_du = zeros(1,max(length(NrCellRange),length(HconvRange)));

for Hconv = HconvRange

numRows    = Hconv;
numColumns = Hconv;
numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

for N=NrCellRange

disp('  ')
disp(['N = ' num2str(N) ', H = ' num2str(Hconv)])
disp('  ')

conv_i = conv_i+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

numbering('square');

Mesh = meshgenerator_square(Domain,DomInfo);

[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

[psi_ex,u_ex,v_ex,w_ex,p_ex,Fy_ex] = ReguralizedLDC(Re,Meshp.X,Meshp.Y);
Energy_exact    = 3.676492819329155e-02;
Enstrophy_exact = 2.01433106575958;
divu_ex = 0*u_ex;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize figure and plot Mesh

if plot_fig || avimatlab
    meshplot
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

% NIEUW !!!!!
[PSIbc,Wbc,boundary_points,interior_points] = boundaryconditions_0_square_NS(Re,Mesh);

t = 0; j=1;

disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

W0 = Wbc;

WW0 = W0(globalnr_0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

% Matrix assembly
NG = normalgradient_assembly();
D  = divergence_assembly();
M0 = innerproduct_assembly(0,Mesh);
M1 = innerproduct_assembly(1,Mesh);
M2 = innerproduct_assembly(2,Mesh);

disp('Assemble Poisson matrix')
MatrixFull = -NG'*M1*NG;

Matrixrhs = MatrixFull(:,boundary_points);

PoissonMatrix = MatrixFull;
PoissonMatrix(:,boundary_points) = [];
PoissonMatrix(boundary_points,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Velocity field

disp('LU decomposition Poisson matrix')
[PM_L,PM_U] = lu(PoissonMatrix);

clear PoissonMatrix

[U,V,PSI] = VorticityInducedVelocityField_old(W0,PSIbc,M0,Matrixrhs,PM_L,PM_U,NG,boundary_points,interior_points);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Wp = reconstruct(0,WW0,hp);
[Up,Vp,Velo] = reconstruct(1,U,V,hp,ep,Meshp);

Clim = [min(min(Wp)) max(max(Wp))];
posten_vorticity

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

disp('creation force function')
f  = zeros(nr_1,1);
Fx = 0*globalnr_1v;
Fy = 0*globalnr_1h;

for i=1:numElements

r = ceil(i/numColumns);
c = i-(r-1)*numColumns;

[fx,fy] = forcefunction(1,r,c,'reguralizedLDC',Domain,DomInfo);

Fx(:,i) = fx;
Fy(:,i) = fy;

end

f(globalnr_1v) = Fx;
f(globalnr_1h) = Fy;

RHS = -NG'*M1*f;

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

WW = zeros(nr_0,1); Wold = W0;
error = 1; er = 0;
while t<T && error>1e-13
    j = j+1;
    t = t+dt;
    er = er+1;

    disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

    if er==100; disp('er=100 !!!'); end

    W_ = timemarching_bc(K1_,K2_,dt,W_,Mbc,Wbc,RHS,'BE');

    WW = Wbc;
    WW(interior_points) = W_;

[U,V] = VorticityInducedVelocityField_old(WW,PSIbc,M0,Matrixrhs,PM_L,PM_U,NG,boundary_points,interior_points);

    vE01 = interiorProduct(U,V,Exe,Eye,he,eh,Mesh);

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

end

W = WW;

disp('Postprocessing')
Wp = reconstruct(0,WW(globalnr_0),hp);
[Up,Vp,Velo] = reconstruct(1,U,V,hp,ep,Meshp);
posten_vorticity

UV(globalnr_1v) = U;
UV(globalnr_1h) = V;
    
%     curlW = NG*W;
%     curlW1 = curlW(globalnr_1v);
%     curlW2 = curlW(globalnr_1h);

D  = divergence_assembly();
divU = D*UV';
divU = divU(globalnr_2);
divu = reconstruct(2,divU,ep,Meshp);

Energy(conv_i) = 0;
for i=1:numElements
Energy(conv_i) = Energy(conv_i) + 1/2*( sum(Up(:,i).*Meshp.J(:,i).*Meshp.W.*Up(:,i)) ...
                            + sum(Vp(:,i).*Meshp.J(:,i).*Meshp.W.*Vp(:,i)) );
end

Enstrophy(conv_i,1) = 0;
for i=1:numElements
Enstrophy(conv_i,1) = Enstrophy(conv_i,1) + 1/2*sum(Wp(:,i).*Meshp.J(:,i).*Meshp.W.*Wp(:,i));
end

% figure
% plot(Enstrophy,'-o')
% title('Enstrophy & Energy')
% hold on
% plot(Energy,'-^r')

% if error_figures
    disp('calculation of errors')
    fouten;
% end

end % for N
end % for H

% for i=1:length(HconvRange)-1
%    ax = 1./(HconvRange);
%    rate_w(i) = (log(L2_w(i+1))-log(L2_w(i))) / (log(ax(i+1))-log(ax(i)));
%    rate_U(i) = (log(L2_U(i+1))-log(L2_U(i))) / (log(ax(i+1))-log(ax(i)));
% %    rate_p(i) = (log(L2_p(i+1))-log(L2_p(i))) / (log(ax(i+1))-log(ax(i)));
%    rate_du(i) = (log(L2_du(i+1))-log(L2_du(i))) / (log(ax(i+1))-log(ax(i)));
% %    rate_cw(i) = (log(L2_cw(i+1))-log(L2_cw(i))) / (log(ax(i+1))-log(ax(i)));
% 
% end
% 
% 
% ax = 1./(HconvRange);
% axisXY = [ 0.05 3 1e-10 1];
% str2 = 'c';
% str3 = ['-o' str2];
% figure(11)
% loglog(ax,L2_U,str3,'markerface',str2)
% % axis(axisXY)
% title('u error Stokes')
% figure(12)
% loglog(ax,L2_w,str3,'markerface',str2)
% % axis(axisXY)
% title('w error Stokes')
% 
% 
% filename = [ 'LDC_Hconv_N' num2str(N) '_c' num2str(10*DomInfo) ];
% save(filename,'NrCellRange','HconvRange','L2_w','L2_u','L2_v','L2_U','L2_du','Energy','Energy_exact','Enstrophy','Enstrophy_exact')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pressure

% [Fx Fy] = forcefunction(1,1,1,'reguralizedLDC_p',Domain,DomInfo);

disp('creation force function')
F  = zeros(nr_1,1);
Fx = 0*globalnr_1v;
Fy = 0*globalnr_1h;
for i=1:numElements

r = ceil(i/numColumns);
c = i-(r-1)*numColumns;

[fx fy] = forcefunction(1,r,c,'reguralizedLDC_p',Domain,DomInfo);

Fx(:,i) = fx;
Fy(:,i) = fy;

end
F(globalnr_1v) = Fx;
F(globalnr_1h) = Fy;

[~,~,boundary_uv,interior_uv] = boundaryconditions_1_square('reguralizedLDC',Domain,DomInfo,[1 0 1 0]);

RHS = [zeros(nr_1,1) ; M2*D*F];

Matrix = [ -M1 D'*M2
           M2*D spalloc(nr_2,nr_2,0) ];

Matrix(:,boundary_uv) = [];
Matrix(boundary_uv,:) = [];

RHS(boundary_uv) = [];

PW = Matrix\RHS;

ind = nr_1-length(boundary_uv);
P = PW(ind+1:end);
P = P(globalnr_2);


pp = reconstruct(2,P,ep,Meshp) + Energy(end);

figure
for i=1:numElements
Xp = reshape(Meshp.X(:,i),nn,nn);
Yp = reshape(Meshp.Y(:,i),nn,nn);
pcolor(Xp,Yp,reshape(pp(:,i),nn,nn))
hold on
end
shading interp
axis equal
axis tight


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
run Library_Convection/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%