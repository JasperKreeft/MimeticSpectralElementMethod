% Zero-forms, Outer-oriented, constant vector-field

clear all
close all
% clc
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N nn globalnr_0 globalnr_1v globalnr_1h
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
U0 = exp(-((Mesh.X-x0).^2+(Mesh.Y-y0).^2)/(2*lambda^2)); % U0 = ones(size(Mesh.X));

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
% TOT HIER HETZELFDE ALS V1

% Inner product zero-forms
M0 = spdiags(Mesh.J.*kron(w,w)',0,(N+1)^2,(N+1)^2);

% Periodic BC
indL = unique(reshape(globalnr_1v(1:(N+1):N*(N+1),1),[],1));
indR = unique(reshape(globalnr_1v(N+1:(N+1):N*(N+1),1),[],1));
indB = unique(reshape(globalnr_1h(1:N+1:N*(N+1),1),[],1));
indT = unique(reshape(globalnr_1h(N+1:N+1:N*(N+1),1),[],1));

NG([indL;indB],:) = NG([indL;indB],:)+NG([indR;indT],:);
NG([indR;indT],:) = [];

vE01(:,[indL;indB]) = (vE01(:,[indL;indB])+vE01(:,[indR;indT]))/2; % central discretization
vE01(:,[indR;indT]) = [];

indB = unique(reshape(globalnr_0(1:N+1,1),[],1));
indT = unique(reshape(globalnr_0(N*(N+1)+(1:N+1),1),[],1));
indL = unique(reshape(globalnr_0(1:N+1:(N+1)^2,1),[],1)); indL(end) = [];
indR = unique(reshape(globalnr_0(N+1:N+1:(N+1)^2,1),[],1)); indR(end) = [];

NG(:,indB) = NG(:,indB)+NG(:,indT);
NG(:,indL) = NG(:,indL)+NG(:,indR);
NG(:,[indR;indT]) = [];

vE01(indB,:) = (vE01(indB,:)+vE01(indT,:))/2;
vE01(indL,:) = (vE01(indL,:)+vE01(indR,:))/2;
vE01([indR;indT],:) = [];

M0(:,indB) = M0(:,indB)+M0(:,indT);
M0(indB,:) = M0(indB,:)+M0(indT,:);
M0(:,indL) = M0(:,indL)+M0(:,indR);
M0(indL,:) = M0(indL,:)+M0(indR,:);
M0(:,[indR;indT]) = [];
M0([indR;indT],:) = [];

A = M0*vE01*NG;

K1 = (A-A')/2;  % Skew-symmetric
% K1 = A;         % Standard
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

    U = timemarching(K2,K1,dt,U,'gauss4');
%     U = mimetictimemarching(K2,K1,dt,U,5);

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