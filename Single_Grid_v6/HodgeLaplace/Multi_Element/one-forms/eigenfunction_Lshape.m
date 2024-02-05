%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% L-shape domain curl-curl problem                                        %
%                                                                         %
% written by Jasper Kreeft (2011)                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library/GetLibrary.m

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns
global xi w
global h e
global nr_0 nr_1
global globalnr_0 globalnr_1v globalnr_1h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'CurlCurl';
Domain       = 'Lshape';
DomInfo      = 0.;

N = 7;
numElements = 3;
numRows = 2;
numColumns = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square')

nr_0 = 3*(N+1)^2-2*(N+1);
nr_1 = 3*2*N*(N+1)-2*N;

[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NGe = normalgrad(N);

NG = zeros(nr_1,nr_0);
M0 = zeros(nr_0);
M1 = zeros(nr_1);

R = [1 1 2];
C = [1 2 1];
for i = 1:3;
    r = R(i);
    c = C(i);

ind0 = globalnr_0(:,i);
ind1 = [ globalnr_1v(:,i) globalnr_1h(:,i) ];

% Gradient operator
NG(ind1,ind0) = NGe;
     
% zero-forms
M0e = innerproduct(0,1);
M0(ind0,ind0) = M0(ind0,ind0) + M0e;

% one-forms
M1e = innerproduct(1,1,1);
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

end

M0 = sparse(M0);
M1 = sparse(M1);
NG = sparse(NG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

load NonZeroEigenvalues_Lshape.mat

% Matrix = M1*NG/M0*NG'*M1 - E(2)*M1;
% RHS = zeros(nr_1,1);

% U = Matrix\RHS;

Matrix = M1*NG/M0*NG'*M1;
RHS = M1;

[Veig,Deig] = eig(full(Matrix),full(RHS));

EE = 4*sort(diag(Deig));

E = EE(abs(EE)>.2);

E(1:min(length(E),20))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

global nn

for i1=1:9

ind = (4*diag(Deig)==E(i1));

U = Veig(:,ind);

U = [U ; zeros(2*N*N,1) ];
UU = U(globalnr_1v);
VV = U(globalnr_1h);

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square('SinDeformGrid',0);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hp,ep,Meshp);

Xp = Meshp.X;
Yp = Meshp.Y;

XYlim = [min(min(min(Xp))) max(max(max(Xp))) min(min(min(Yp))) max(max(max(Yp)))];


for i=1:3

Xp = reshape(Meshp.X(:,i),nn,nn);
Yp = reshape(Meshp.Y(:,i),nn,nn);
Vp = reshape(velo(:,i),nn,nn);
figure(1)
subplot(3,3,i1)
colormap jet%bone
pcolor(Xp,-Yp,Vp)
hold on
% quiver(Xp(:,:,rc),Yp(:,:,rc),uu(:,:,rc),vv(:,:,rc),'w')
shading interp
% colorbar
% set(gca,'clim',[0 14.5])
text(0.2,-0.6,['eigenvalue ' num2str(i1)])
axis equal
axis(XYlim)
axis off
box off

% figure(2)
% subplot(3,3,i1)
% quiver(Xp(:,:,rc),Yp(:,:,rc),uu(:,:,rc),vv(:,:,rc))
% hold on
% text(0.2,+0.6,['eigenvalue ' num2str(i1)])
% axis equal
% axis(XYlim)
% axis off
% box off

% figure(3)
% surf(Xp(:,:,rc),Yp(:,:,rc),velo(:,:,rc))
% hold on
% shading interp

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%