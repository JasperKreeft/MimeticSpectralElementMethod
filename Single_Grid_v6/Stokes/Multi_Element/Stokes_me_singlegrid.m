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

FunctionType = 'sine1';
Domain       = 'SinDeformGrid_01';%'MuST';%
DomInfo      = 0.0;

NrCellRange = 2%:2:16;
HconvRange  = [ 2 4 8 16 ];%32


plot_figures  = 0;
error_figures = 1;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Stokes Problem ')
disp(['Testfunction is ' FunctionType])
disp(['Computational Domain is ' Domain ' with settings ' num2str(DomInfo)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


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

NG = normalgradient_assembly();

D = divergence_assembly();

M0 = innerproduct_assembly(0,Mesh);

M1 = innerproduct_assembly(1,Mesh);

M2 = innerproduct_assembly(2,Mesh);

Wbc0 = boundaryIntegral_assembly();

Wbc1 = boundaryIntegral_assembly_1();

disp('creation force function')
f  = zeros(nr_1,1);
Fx = 0*globalnr_1v;
Fy = 0*globalnr_1h;
G  = 0*globalnr_2;

for i=1:numElements

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

% [fx fy] = force_oneform_square(i,FunctionType,Domain,DomInfo);
r = ceil(i/numColumns);
c = i-(r-1)*numColumns;

[fx fy] = forcefunction(1,r,c,FunctionType,Domain,DomInfo);
g = zeros(N^2,1); % g       = massfunction(r,c,FunctionType,Domain,DomInfo);

Fx(:,i) = fx;
Fy(:,i) = fy;

G(:,i) = g;

end

f(globalnr_1v) = Fx;
f(globalnr_1h) = Fy;

g(globalnr_2) = G;


disp('creation of boundary conditions')

[Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_square(Mesh,FunctionType,[0 0 0 0]);

[Pbc,UVbc,boundary_uv,interior_uv] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,[1 1 1 1]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

disp('assembly system matrix and righthandside')

% if Hconv>16
if N>8
Matrix = [ spalloc(nr_1,nr_1,0)      M1*NG               D'*M2
                NG'*M1                 M0          spalloc(nr_0,nr_2,0)
                 M2*D         spalloc(nr_2,nr_0,0) spalloc(nr_2,nr_2,0) ];

else
Matrix = [ 1e-12*speye(nr_1)        M1*NG               D'*M2
                NG'*M1              M0          spalloc(nr_0,nr_2,0)
                 M2*D      spalloc(nr_2,nr_0,0)   1e-12*speye(nr_2)   ];
end
% Matrix = [ spalloc(nr_1,nr_1,0)      M1*NG               D'*M2
%                 NG'*M1                 M0          spalloc(nr_0,nr_2,0)
%                  D         spalloc(nr_2,nr_0,0) spalloc(nr_2,nr_2,0) ];

RHS = [    -M1*f
        zeros(nr_0,1)
           M2*g ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

% 0
RHS(nr_1+1:nr_1+nr_0) = RHS(nr_1+1:nr_1+nr_0) - Wbc0*TangentialVelocity_bc;
if ~isempty(boundary_w)
RHS = RHS - Matrix(:,boundary_w+nr_1)*Vorticity_bc;
end

% 1
RHS(1:nr_1) = RHS(1:nr_1) + Wbc1*Pbc;
if ~isempty(boundary_uv)
RHS = RHS - Matrix(:,boundary_uv)*UVbc;
end

ind = sort([boundary_uv ; boundary_w+nr_1]);

Matrix(:,ind) = [];
Matrix(ind,:) = [];

RHS(ind,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

disp('Solve system')

UVWP = Matrix\RHS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

disp('Postprocessen')

ind1 = nr_1-length(boundary_uv);
UV_in = UVWP(1:ind1);
ind2 = ind1+nr_0-length(boundary_w);
W_in = UVWP(ind1+1:ind2);
ind3 = ind2+nr_2;
P_in = UVWP(ind2+1:ind3);

W  = zeros(nr_0,1);
UV = zeros(nr_1,1);
P  = zeros(nr_2,1);

W(interior_w)   = W_in;
W(boundary_w)   = Vorticity_bc;
UV(interior_uv) = UV_in;
UV(boundary_uv) = UVbc;
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



% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hp,ep,Meshp);
ww           = reconstruct(0,WW,hp);
pp           = reconstruct(2,PP,ep,Meshp);
% pp = pp-pp(1,1)+exact_solution(Meshp.X(1,1),Meshp.Y(1,1),FunctionType,'two');
pp = pp-pp(end,end)+exact_solution(Meshp.X(end,end),Meshp.Y(end,end),FunctionType,'two');

[ffx,ffy,ff] = reconstruct(1,Fx,Fy,hp,ep,Meshp);
gg           = reconstruct(2,G,ep,Meshp);

[curlw1,curlw2] = reconstruct(1,curlW1,curlW2,hp,ep,Meshp);
divu = reconstruct(2,divU,ep,Meshp);

% Exact
w_ex = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[u_ex,v_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');
p_ex = exact_solution(Meshp.X,Meshp.Y,FunctionType,'two');
[fx_ex fy_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'force');
divu_ex = exact_solution(Meshp.X,Meshp.Y,FunctionType,'mass');

% curlw1_ex = -9/4*pi^2*sin(3/2*pi*Meshp.X).*cos(3/2*pi*Meshp.Y);
% curlw2_ex = 9/4*pi^2*cos(3/2*pi*Meshp.X).*sin(3/2*pi*Meshp.Y);
curlw1_ex = 8*pi^2*sin(2*pi*Meshp.X).*cos(2*pi*Meshp.Y);
curlw2_ex = -8*pi^2*cos(2*pi*Meshp.X).*sin(2*pi*Meshp.Y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotten

if plot_figures
    disp('creation colorplots')
    meshplot
    plotten_me
end

%% Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if error_figures
    disp('calculation of errors')
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
axisXY = [ 0.05 3 1e-10 1];
str1 = '-^g';
str2 = 'g';
str3 = '--g';
figure(1)
loglog(ax,Hd_u,str1,'markerface',str2)
hold on
loglog(ax,L2_U,str3,'markerface',str2)
% axis(axisXY)
title('u error Stokes')
figure(2)
loglog(ax,Hd_w,str1,'markerface',str2)
hold on
loglog(ax,L2_w,str3,'markerface',str2)
% axis(axisXY)
title('w error Stokes')
figure(3)
loglog(ax,L2_p,str1,'markerface',str2)
hold on
% axis(axisXY)
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


for i=1:length(HconvRange)-1
    
   rate_w(i) = (log(L2_w(i+1))-log(L2_w(i))) / (log(ax(i+1))-log(ax(i)));
   rate_U(i) = (log(L2_U(i+1))-log(L2_U(i))) / (log(ax(i+1))-log(ax(i)));
   rate_p(i) = (log(L2_p(i+1))-log(L2_p(i))) / (log(ax(i+1))-log(ax(i)));
   rate_du(i) = (log(L2_du(i+1))-log(L2_du(i))) / (log(ax(i+1))-log(ax(i)));
   rate_cw(i) = (log(L2_cw(i+1))-log(L2_cw(i))) / (log(ax(i+1))-log(ax(i)));
   
end


end

% filename = [ 'LDC_Hconv_N' num2str(Hconv) '_c' num2str(10*DomInfo) ];
% save(filename,'N','HconvRange','L2_w','L2_cw1','L2_cw2','L2_cw','Hd_w','L2_u','L2_v','L2_U','Linf_du','L1_du','L2_du','Hd_u','L2_p','rate_w','rate_cw','rate_U','rate_du','rate_p')
filename = [ 'LDC_Pconv_H' num2str(Hconv) '_c' num2str(10*DomInfo) ];
save(filename,'NrCellRange','Hconv','L2_w','L2_cw1','L2_cw2','L2_cw','Hd_w','L2_u','L2_v','L2_U','Linf_du','L1_du','L2_du','Hd_u','L2_p')

%%%%%%%%%%%%%%%%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%