% JJK: 15-3-2011, check dirichlet boundary conditions, i.e. bc = 1 !!!
% JJK: 30-6-2011, maximaal een boundary bc = 1 !!

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N
global xi w
global h e
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'AFG';
Domain       = 'SinDeformGrid_01';
DomInfo      = 0.0;

NrCellRange = 8%1:8;

plot_figures  = 1;
error_figures = 0;
casenr        = FunctionType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

for N=NrCellRange

disp(['N = ' num2str(N)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square')

[xi,w] = GLLnodes(N);     eta = xi;  % Gauss-Lobotto-Legendre

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

NG = normalgrad(N);

D  = div(N);

[h,e] = MimeticpolyVal(xi,N,1);

M0 = innerproduct(0,Mesh.J);

M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

M2 = innerproduct(2,Mesh.J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

[~,~,boundary,interior] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,[ 0 0 0 0 ]);
boundary = sort(boundary);

[fx fy] = forcefunction(1,1,1,FunctionType,Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

% Matrix = [ spalloc(nr_1,nr_1,0)    M1*NG                 D'*M2
%                NG'*M1               M0           spalloc(nr_0,nr_2,0)
%                 M2*D        spalloc(nr_2,nr_0,0) spalloc(nr_2,nr_2,0) ];

Matrix = [ eps*speye(nr_1)    M1*NG                 D'*M2
               NG'*M1               M0           spalloc(nr_0,nr_2,0)
                M2*D        spalloc(nr_2,nr_0,0) eps*speye(nr_2) ];


RHS = [ -M1*[ fx ; fy ]
         zeros(nr_0,1)
         zeros(nr_2,1) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

Matrix(:,boundary) = [];
Matrix(boundary,:) = [];

RHS(boundary,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

UVWP = Matrix\RHS;

ind1 = nr_1-length(boundary);
UV_in = UVWP(1:ind1);
ind2 = ind1+nr_0;
W_in = UVWP(ind1+1:ind2);
ind3 = ind2+nr_2;
P_in = UVWP(ind2+1:ind3);

UV = zeros(nr_1,1);
UV(interior) = UV_in;
% UV(boundary) = UVbc;

W  = W_in;
P = P_in;

UU = UV(globalnr_1v);
VV = UV(globalnr_1h);
WW = W(globalnr_0);
PP = P(globalnr_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hGLp,eGLp] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hGLp,eGLp,Meshp);
ww           = reconstruct(0,WW,hGLp);
pp           = reconstruct(2,PP,eGLp,Meshp); pp = pp-mean(mean(pp));
[ffx,ffy]    = reconstruct(1,fx,fy,hGLp,eGLp,Meshp);

% Interpolated solution

% Exact Solution
w_ex          = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[ u_ex v_ex ] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');
p_ex          = exact_solution(Meshp.X,Meshp.Y,FunctionType,'two'); p_ex = p_ex-mean(mean(p_ex));
[fx_ex fy_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'force');

if plot_figures
    plotten
end

if error_figures
    fout
end


end % for N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convergence plots

if error_figures
    foutplot
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%