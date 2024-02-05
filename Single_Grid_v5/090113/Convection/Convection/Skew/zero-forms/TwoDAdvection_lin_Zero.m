% Zero-forms, Outer-oriented, constant vector-field

clear all
close all
clc
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N nn globalnr_0
global w e

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'SinDeformGrid';
DomInfo = 0.0;

plot_fig  = 1;
avimatlab = 0;
Tecplot   = 0;

N = 20;
T = 5;
dt = 0.05;

filename = ['Zero_lin_SE_N' num2str(N) '_dt' num2str(dt)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

numbering('square')

Mesh = meshgenerator_square(Domain,DomInfo);

[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

Xp = reshape(Meshp.X,nn,nn);
Yp = reshape(Meshp.Y,nn,nn);

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
    axis([min(Mesh.X) max(Mesh.X) min(Mesh.Y) max(Mesh.Y)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

t = 0; j=1;

disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

lambda = 1/8;
x0 = -1/2;
y0 =  1/2;

% Gaussian function
U0 = exp(-((Mesh.X-x0).^2+(Mesh.Y-y0).^2)/(2*lambda^2));

Up = reconstruct(0,U0,hp);

if plot_fig || avimatlab
    subplot(1,2,2)
end
if Tecplot
    dataXY = [ Mesh.X Mesh.Y ];
end

posten_v2

Up_ex = exp(-((Meshp.X-x0).^2+(Meshp.Y-y0).^2)/(2*lambda^2));

L2error(1) = sqrt(sum(Meshp.W.*(Up-Up_ex).^2));
L2norm = sqrt(sum(Meshp.W.*Up_ex.^2));
Energy(1) = 1/2*sum(Up.*Meshp.W.*Up);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization

Ax = 1;
Ay = 1;

NG = normalgrad(N);

Ex = zeros((N+1)^2,N*(N+1));
for i=1:N
    Ex(1:(N+1)^2,(i-1)*(N+1)+(1:N+1)) = kron(eye(N+1),e(i,:)');
end
Ey = kron(e',eye(N+1));

vE01 = [ -Ay*Ey Ax*Ex ];

JW = Mesh.J.*kron(w,w)';

% Periodic BC
for i=2*N:-1:1
    ind1 = (i-1)*(N+1)+1;
    ind2 = i*(N+1);
    
    NG(ind1,:) = NG(ind1,:)+NG(ind2,:);
    NG(ind2,:) = [];
    
    vE01(:,ind1) = (vE01(:,ind1)+vE01(:,ind2))/2;  % central discretization
    vE01(:,ind2) = [];

end

ind1 = 1:N+1;
ind2 = N*(N+1)+(1:N+1);
NG(:,ind1) = NG(:,ind1)+NG(:,ind2);
NG(:,ind2) = [];
vE01(ind1,:) = (vE01(ind1,:)+vE01(ind2,:))/2;  % central discretization
vE01(ind2,:) = [];
JW(ind1) = JW(ind1) + JW(ind2);
JW(ind2) = [];

ind1 = 1:(N+1):N*(N+1);
ind2 = (N+1):(N+1):N*(N+1);
NG(:,ind1) = NG(:,ind1)+NG(:,ind2);
NG(:,ind2) = [];
vE01(ind1,:) = (vE01(ind1,:)+vE01(ind2,:))/2;  % central discretization
vE01(ind2,:) = [];
JW(ind1) = JW(ind1) + JW(ind2);
JW(ind2) = [];

% Inner product zero-forms
M0 = spdiags(JW,0,N^2,N^2);

A = M0*vE01*NG;

K1 = sparse( (A-A')/2 );  % Skew-symmetric
% K1 = sparse(A);           % Standard
K2 = M0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind1 = N*(N+1)+(1:N+1);
ind2 = N+1:N+1:N*(N+1);
ind3 = globalnr_0; ind3(ind1) = []; ind3(ind2) = [];

U = U0; U(ind1) = []; U(ind2) = [];

Xh = Meshp.X-x0; Yh = Meshp.Y-y0;
UU = zeros((N+1)^2,1);
while t<T
    j = j+1;
    t = t+dt;

    disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

%     U = timemarching(K2,K1,dt,U,'gauss4');
    U = mimetictimemarching(K2,K1,dt,U,5);

    UU(ind3) = U;
    UU(ind1) = U(1:N+1);
    UU(ind2) = U(1:N:N^2);

    Up = reconstruct(0,UU,hp);

    posten_v2
    
    Xh = Xh - Ax*dt + 2*((Xh<-1)-(Xh>1));
    Yh = Yh - Ay*dt + 2*((Yh<-1)-(Yh>1));
    Up_ex = exp(-(Xh.^2+Yh.^2)/(2*lambda^2));
    
    L2error(j) = sqrt(sum(Meshp.W.*(Up-Up_ex).^2));

    Energy(j) = 1/2*sum(Up.*Meshp.W.*Up);

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
toc