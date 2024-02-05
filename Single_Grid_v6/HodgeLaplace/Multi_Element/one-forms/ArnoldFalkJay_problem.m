%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Multi-element Arnold-Falk-Jay problem                                   %
%                                                                         %
% written by Jasper Kreeft (2012)                                         %
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
global nr_0 nr_1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

BC = 2;

Domain = 'SinDeformGrid_01';
if BC==1
    FunctionType = 'arnoldfalkjay1';
elseif BC==2
    FunctionType = 'arnoldfalkjay2';
end
DomInfo = 0.0;

NrCellRange = 3;
HconvRange  = 4;%2.^(5:7);

plot_figures  = 1;
error_figures = 0;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Hodge-Laplace for one-forms ')
disp('The Arnold-Falk-Gopalakrishnan (2012) testcases')
disp(['Testfunction is ' FunctionType])
disp(['Computational Domain is ' Domain ' with settings ' num2str(DomInfo)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start H-convergence loop

er = 0;
L2_s = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_u = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_v = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_du = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_ngx = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_ngy = zeros(1,max(length(NrCellRange),length(HconvRange)));

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

f  = zeros(nr_1,1);
Fx = 0*globalnr_1v;
Fy = 0*globalnr_1h;

for i=1:numElements

% [fx fy] = force_oneform_square(i,FunctionType,Domain,DomInfo);
r = ceil(i/numColumns);
c = i-(r-1)*numColumns;
[fx fy] = forcefunction(1,r,c,FunctionType,Domain,DomInfo);

Fx(:,i) = fx;
Fy(:,i) = fy;

end

f(globalnr_1v) = Fx;
f(globalnr_1h) = Fy;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

disp('assembly system matrix and righthandside')

Matrix = -[   M0    NG'*M1
            M1*NG -D'*M2*D ];

RHS = [ zeros(nr_0,1)
        M1*f ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

if BC==1
    boundary = [];
    interior = 1:nr_1;
elseif BC==2
[~,~,boundary,interior] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,[ 0 0 0 0 ]);
boundary = sort(boundary);
end

Matrix(:,nr_0+boundary) = [];
Matrix(nr_0+boundary,:) = [];

RHS(nr_0+boundary,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

disp('Solve system')

SigmaUV = Matrix\RHS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

disp('Postprocessen')

Sigma = SigmaUV(1:nr_0);
UV_in = SigmaUV(nr_0+1:end);

UV = zeros(nr_1,1);
UV(interior) = UV_in;
UV(boundary) = zeros(length(boundary),1);

SSigma = Sigma(globalnr_0);
UU = UV(globalnr_1v);
VV = UV(globalnr_1h);

DivUV = D*UV;
DivUV = DivUV(globalnr_2);

NG_sigma = NG*Sigma;
NG_sigma_x = NG_sigma(globalnr_1v);
NG_sigma_y = NG_sigma(globalnr_1h);



% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
ssigma        = reconstruct(0,SSigma,hp);
[uu,vv,velo]  = reconstruct(1,UU,VV,hp,ep,Meshp);
[ffx,ffy]     = reconstruct(1,Fx,Fy,hp,ep,Meshp);
[ng_sx,ng_sy] = reconstruct(1,NG_sigma_x,NG_sigma_y,hp,ep,Meshp);
div_uv        = reconstruct(2,DivUV,ep,Meshp);

% Exact
sigma_ex          = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[ u_ex v_ex ]     = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');
[fx_ex fy_ex]     = exact_solution(Meshp.X,Meshp.Y,FunctionType,'force');
[ng_x_ex,ng_y_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'d_zero');
div_ex            = exact_solution(Meshp.X,Meshp.Y,FunctionType,'d_one');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotten

if plot_figures
    disp('creation colorplots')
%     meshplot
    plotten_AFG
end

%% Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if error_figures
    disp('calculation of errors')
    fout_AFG;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keyboard

end % for N
end % for H

if error_figures
%     errorplot

% ax = NrCellRange;
% axisXY = [ 1 17 1e-12 1];
% str1 = '-vc';
% str2 = 'c';
% 
% figure(1)
% % semilogy(ax,Hd_u,str1,'markerface',str2)
% semilogy(ax,L2_U,str1,'markerface',str2)
% hold on
% % axis(axisXY)
% title('u error')
% figure(2)
% % semilogy(ax,Hd_s,str1,'markerface',str2)
% semilogy(ax,L2_s,str1,'markerface',str2)
% hold on
% title('\sigma error')
% % axis(axisXY)
% figure(3)
% semilogy(ax,L2_du,str1,'markerface',str2)
% hold on
% title('div u error')
% % axis(axisXY)

ax = 1./(HconvRange);
axisXY = [ 0.05 3 1e-10 1];
str1 = '-^r';
str2 = 'r';
figure(1)
% loglog(ax,Hd_u,str1,'markerface',str2)
loglog(ax,L2_U,str1,'markerface',str2)
hold on
% axis(axisXY)
title('u error')
figure(2)
% loglog(ax,Hd_s,str1,'markerface',str2)
loglog(ax,L2_s,str1,'markerface',str2)
hold on
% axis(axisXY)
title('\sigma error')
figure(3)
loglog(ax,L2_du,str1,'markerface',str2)
hold on
% axis(axisXY)
title('du error')
figure(4)
loglog(ax,L2_ngs,str1,'markerface',str2)
hold on
% axis(axisXY)
title('ngs error')

for i=1:length(HconvRange)-1

   rate_s(i) = (log(L2_s(i+1))-log(L2_s(i))) / (log(ax(i+1))-log(ax(i)));
   rate_U(i) = (log(L2_U(i+1))-log(L2_U(i))) / (log(ax(i+1))-log(ax(i)));
   rate_du(i) = (log(L2_du(i+1))-log(L2_du(i))) / (log(ax(i+1))-log(ax(i)));
   rate_ngs(i) = (log(L2_ngs(i+1))-log(L2_ngs(i))) / (log(ax(i+1))-log(ax(i)));

end

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

disp('Calculations finished')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%