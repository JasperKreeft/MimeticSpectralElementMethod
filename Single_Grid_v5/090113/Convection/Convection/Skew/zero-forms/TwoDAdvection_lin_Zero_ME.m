clear all
close all
clc

system('rm -f *.dat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
if ispc
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'
elseif isunix
run '/home/jjkreeft/Convection/Convection/Skew/Library_Convection/GetLibrary.m'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns
global xi
global h e
global nr_0 nr_1
global globalnr_0 globalnr_1v globalnr_1h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'SinDeformGrid';
DomInfo = 0.0;

plot_fig  = 1;
avimatlab = 0;
Tecplot   = 0;

filename = 'tecplottest';

N = 5;
H = 8;
T = 15;
dt = 0.02;

numRows    = H;
numColumns = H;
numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

numbering('square');

% Mesh = meshgenerator_square(Domain,DomInfo);
% [Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);
Mesh = meshgenerator_square_smooth(Domain);
[Meshp,hp,ep] = postproces_grid_square_smooth(Domain);

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

t = 0; j=1;

disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

lambda = 1/8;
x0 = -1/2;
y0 = 0;


% Gaussian function
U0 = exp(-((Mesh.X-x0).^2+(Mesh.Y-y0).^2)/(2*lambda^2));
% U0 = ones(size(Mesh.X));

Up = reconstruct(0,U0,hp);

UU=U0;
if plot_fig || avimatlab; subplot(1,2,2); end

posten_v3

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

NGe = normalgrad(N);

Exe = zeros((N+1)^2,N*(N+1));
for i=1:N
    Exe(:,(i-1)*(N+1)+(1:N+1)) = kron(eye(N+1),e(i,:)');
end
Eye = kron(e',eye(N+1));

NG  = zeros(nr_1,nr_0);
M0  = zeros(nr_0);
vEx  = zeros(nr_0,nr_1);
vEy  = zeros(nr_0,nr_1);

for i=1:numElements

ind0 = globalnr_0(:,i);
ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];

% normal gradient operator
NG(ind1,ind0) = NGe;

% zero-forms
M0e = innerproduct(0,Mesh.J(:,i));
M0(ind0,ind0) = M0(ind0,ind0) + M0e;


% Axe = 1;
% Aye = 1;
Axe = (Mesh.dYdEta(:,i) - Mesh.dXdEta(:,i))./Mesh.J(:,i);
Aye = (Mesh.dXdXi(:,i)  - Mesh.dYdXi(:,i) )./Mesh.J(:,i);

% keyboard
ind1_h = globalnr_1h(:,i);
vEx(ind0,ind1_h) = vEx(ind0,ind1_h) + diag(Axe)*Exe;

ind1_v = globalnr_1v(:,i);
vEy(ind0,ind1_v) = vEy(ind0,ind1_v) + diag(Aye)*Eye;

end

vE = -vEy+vEx;

% Periodic BC
indL = unique(reshape(globalnr_1v(1:(N+1):N*(N+1),1:numColumns:numElements),[],1));
indR = unique(reshape(globalnr_1v(N+1:(N+1):N*(N+1),numColumns:numColumns:numElements),[],1));
indB = unique(reshape(globalnr_1h(1:N+1:N*(N+1),1:numColumns),[],1));
indT = unique(reshape(globalnr_1h(N+1:N+1:N*(N+1),(numRows-1)*numColumns+1:numElements),[],1));

NG([indL;indB],:) = NG([indL;indB],:)+NG([indR;indT],:);
NG([indR;indT],:) = [];

vE(:,[indL;indB]) = (vE(:,[indL;indB])+vE(:,[indR;indT]))/2; % central discretization
vE(:,[indR;indT]) = [];

indB = unique(reshape(globalnr_0(1:N+1,1:numColumns),[],1));
indT = unique(reshape(globalnr_0(N*(N+1)+(1:N+1),(numRows-1)*numColumns+1:numElements),[],1));
indL = unique(reshape(globalnr_0(1:N+1:(N+1)^2,1:numColumns:numElements),[],1)); indL(end) = [];
indR = unique(reshape(globalnr_0(N+1:N+1:(N+1)^2,numColumns:numColumns:numElements),[],1)); indR(end) = [];

NG(:,indB) = NG(:,indB)+NG(:,indT);
NG(:,indL) = NG(:,indL)+NG(:,indR);
NG(:,[indR;indT]) = [];

ind_element_boundary = [1:N+1 N+2:N+1:N*(N+1) 2*(N+1):N+1:N*(N+1) N*(N+1)+1:(N+1)^2];
ind_all_el_bcs = unique(reshape(globalnr_0(ind_element_boundary,:),[],1));
vE(ind_all_el_bcs,:) = vE(ind_all_el_bcs,:)/2; % central discretization

ind = unique(globalnr_0([1 N+1 N*(N+1)+1 (N+1)^2],:));
vE(ind,:) = vE(ind,:)/2; % central discretization

vE(indB,:) = (vE(indB,:)+vE(indT,:));
vE(indL,:) = (vE(indL,:)+vE(indR,:));
vE([indR;indT],:) = [];

M0(:,indB) = M0(:,indB)+M0(:,indT);
M0(indB,:) = M0(indB,:)+M0(indT,:);
M0(:,indL) = M0(:,indL)+M0(:,indR);
M0(indL,:) = M0(indL,:)+M0(indR,:);
M0(:,[indR;indT]) = [];
M0([indR;indT],:) = [];


NG = sparse(NG);
M0 = sparse(M0);
vE = sparse(vE);

A = M0*vE*NG;

K1 = (A-A')/2;  % Skew-symmetric
% K1 = A;         % Standard
K2 = M0;

ind3 = 1:nr_0; ind3([indR;indT]) = [];

U = zeros(nr_0,1);
U(globalnr_0) = U0;
U([indR;indT]) = [];

Xh = Meshp.X-x0; Yh = Meshp.Y-y0;
UUU = zeros(nr_0,1);
while t<T
    j = j+1;
    t = t+dt;

    disp(['nr = ' num2str(j) ', t = ' num2str(t,4)])

    U = timemarching(K2,K1,dt,U,'gauss4');%
%     U = mimetictimemarching(K2,K1,dt,U,5);
    
    UUU(ind3) = U;
    UUU(indR) = UUU(indL);
    UUU(indT) = UUU(indB);
    
    UU = UUU(globalnr_0);

    Up = reconstruct(0,UU,hp);

    posten_v3

    Ax = 1; Ay = 1;
    Xh = Xh - Ax/H*dt + 2*((Xh<-1)-(Xh>1));
    Yh = Yh - Ay/H*dt + 2*((Yh<-1)-(Yh>1));
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

in = 'finish';
if ispc
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'
elseif isunix
run '/home/jjkreeft/Convection/Convection/Skew/Library_Convection/GetLibrary.m'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#ok<*UNRCH>