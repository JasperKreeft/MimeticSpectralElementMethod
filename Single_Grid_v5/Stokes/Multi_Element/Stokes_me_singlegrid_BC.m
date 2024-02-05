%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Multi-element Stokes problem                                            %
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns
global xi
global h e
global globalnr_0 globalnr_1v globalnr_1h globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'sine2';
Domain       = 'SinCosDeformGrid';%'MuST';%'SinDeformGrid_01';
DomInfo      = 0.1;

% bc = [ 1 0 0 0 ]; % 1 = Dirichlet, 0 = Neumann

NrCellRange = 6;
HconvRange  = 6;%[ 1 2 4 8 16 32 ];


plot_figures  = 1;
error_figures = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start H-convergence loop

er = 0;
L2_u = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_v = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_w = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_p = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_du = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_cw1 = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_cw2 = zeros(1,max(length(NrCellRange),length(HconvRange)));
L1_du = zeros(1,max(length(NrCellRange),length(HconvRange)));
Linf_du = zeros(1,max(length(NrCellRange),length(HconvRange)));

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square');

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NGe = normalgrad(N);

De = div(N);

NG = zeros(nr_1,nr_0);
D  = zeros(nr_2,nr_1);
M0 = spalloc(nr_0,nr_0,nr_0);
M1 = zeros(nr_1);
M2 = zeros(nr_2);
f  = zeros(nr_1,1);
Fx = 0*globalnr_1v;
Fy = 0*globalnr_1h;
Wbc0 = zeros(nr_0);
Wbc1 = zeros(nr_1);

for i=1:numElements

ind0 = globalnr_0(:,i);
ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
ind2 = globalnr_2(:,i);

% Normal Gradient operator
NG(ind1,ind0) = NGe;

% Divergence operator
D(ind2,ind1) = De;

% zero-forms
M0e = innerproduct(0,Mesh.J(:,i));
M0(ind0,ind0) = M0(ind0,ind0) + M0e;

% one-forms
M1e = innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)));
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

% two-forms
M2e = innerproduct(2,Mesh.J(:,i));
M2(ind2,ind2) = M2e;

% boundary matrix zero form
Wbc0(ind0,ind0) = Wbc0(ind0,ind0) + boundaryIntegral([0 1]);

% boundary matrix one form
Wbc1(ind1,ind1) = Wbc1(ind1,ind1) + boundaryIntegral([1 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

[Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_square(Mesh,FunctionType,[0 0 0 1]);

[Pbc,UVbc,boundary_uv,interior_uv] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,[0 0 0 1]);

% [fx fy] = force_oneform_square(i,FunctionType,Domain,DomInfo);
r = ceil(i/numColumns);
c = i-(r-1)*numColumns;
[fx fy] = forcefunction(1,r,c,FunctionType,Domain,DomInfo);

Fx(:,i) = fx;
Fy(:,i) = fy;

end

D = sparse(D);
NG = sparse(NG);
M0 = sparse(M0);
M1 = sparse(M1);
M2 = sparse(M2);
Wbc0 = sparse(Wbc0);
Wbc1 = sparse(Wbc1);

f(globalnr_1v) = Fx;
f(globalnr_1h) = Fy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

% Matrix = [       M0           NG'*M1        spalloc(nr_0,nr_2,0)
%                M1*NG   spalloc(nr_1,nr_1,0)      D'*M2
%          spalloc(nr_2,nr_0,0)  M2*D         spalloc(nr_2,nr_2,0) ];

Matrix = [       M0           NG'*M1   spalloc(nr_0,nr_2,0)
               M1*NG   1e-10*speye(nr_1)      D'*M2
         spalloc(nr_2,nr_0,0)  M2*D      1e-10*speye(nr_2) ];

RHS = [ zeros(nr_0,1)
        -M1*f
        zeros(nr_2,1) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

% 0
RHS(1:nr_0) = RHS(1:nr_0) - Wbc0*TangentialVelocity_bc;
if ~isempty(boundary_w)
RHS = RHS - Matrix(:,boundary_w)*Vorticity_bc;
end

% 1
RHS(nr_0+1:nr_0+nr_1) = RHS(nr_0+1:nr_0+nr_1) + Wbc1*Pbc;
if ~isempty(boundary_uv)
RHS = RHS - Matrix(:,boundary_uv+nr_0)*UVbc;
end

ind = sort([boundary_w ; boundary_uv+nr_0]);

Matrix(:,ind) = [];
Matrix(ind,:) = [];

RHS(ind,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

% WUVP = minres(Matrix,RHS,1e-10,2000);
WUVP = Matrix\RHS;

ind1 = nr_0-length(boundary_w);
W_in = WUVP(1:ind1);
ind2 = ind1+nr_1-length(boundary_uv);
UV_in = WUVP(ind1+1:ind2);
ind3 = ind2+nr_2;
P_in = WUVP(ind2+1:ind3);

UV = zeros(nr_1,1);
W  = zeros(nr_0,1);
P  = zeros(nr_2,1);

W(interior_w)   = W_in;
W(boundary_w)   = Vorticity_bc;
UV(boundary_uv) = UVbc;
UV(interior_uv) = UV_in;
P = P_in;

UU = UV(globalnr_1v);
VV = UV(globalnr_1h);
WW = W(globalnr_0);
PP = P(globalnr_2);
if N==1
    PP = PP';
end

divU = D*UV;
divU = divU(globalnr_2);

curlW = NG*W;
curlW1 = curlW(globalnr_1v);
curlW2 = curlW(globalnr_1h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hp,ep,Meshp);
ww           = reconstruct(0,WW,hp);
pp           = reconstruct(2,PP,ep,Meshp);
pp = pp-pp(1,1)+exact_solution(Meshp.X(1,1),Meshp.Y(1,1),FunctionType,'two');

[ffx,ffy,ff] = reconstruct(1,Fx,Fy,hp,ep,Meshp);

[curlw1,curlw2] = reconstruct(1,curlW1,curlW2,hp,ep,Meshp);
divu = reconstruct(2,divU,ep,Meshp);

% Exact
w_ex = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[u_ex,v_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');
p_ex = exact_solution(Meshp.X,Meshp.Y,FunctionType,'two');

divu_ex = zeros(size(Meshp.X));

curlw1_ex =  8*pi^2*sin(2*pi*Meshp.X).*cos(2*pi*Meshp.Y);
curlw2_ex = -8*pi^2*cos(2*pi*Meshp.X).*sin(2*pi*Meshp.Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotten

if plot_figures
    meshplot
    plotten_me
end

%% Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if error_figures
    fout_Stokes;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end % for N
end % for H

if error_figures
%     errorplot

% ax = NrCellRange;
% axisXY = [ 1 17 1e-12 1];
% str1 = '-vc';
% str2 = 'c';

% figure(1)
% semilogy(ax,Hd_u,str1,'markerface',str2)
% hold on
% axis(axisXY)
% title('u error Stokes')
% figure(2)
% semilogy(ax,Hd_w,str1,'markerface',str2)
% hold on
% title('w error Stokes')
% axis(axisXY)
% figure(3)
% semilogy(ax,L2_p,str1,'markerface',str2)
% hold on
% title('p error Stokes')
% axis(axisXY)

ax = 1./(HconvRange);
axisXY = [ 0.02 2 1e-10 1];
str1 = '-ob';
str2 = 'b';
figure(1)
loglog(ax,Hd_u,str1,'markerface',str2)
hold on
axis(axisXY)
title('u error Stokes')
figure(2)
loglog(ax,Hd_w,str1,'markerface',str2)
hold on
axis(axisXY)
title('w error Stokes')
figure(3)
loglog(ax,L2_p,str1,'markerface',str2)
hold on
axis(axisXY)
title('p error Stokes')

% figure(1)
%     loglog(ax,L2_u,str1,'markerface',str2)
%     hold on
% figure(2)
%     loglog(ax,L2_v,str1,'markerface',str2)
%     hold on
% figure(3)
%     loglog(ax,L2_w,str1,'markerface',str2)
%     hold on
% figure(4)
%     loglog(ax,L2_p,str1,'markerface',str2)
%     hold on
% figure(5)
%     loglog(ax,L2_du,str1,'markerface',str2)
%     hold on
% figure(6)
%     loglog(ax,Hd_u,str1,'markerface',str2)
%     hold on
% %     loglog(ax,L2_U,'--sr','markerface','r')
% figure(7)
%     loglog(ax,L2_w,str1,'markerface',str2)
%     hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%