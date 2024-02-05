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

FunctionType = 'sine4';
Domain       = 'SinCosDeformGrid';%'SinDeformGrid_01';%'MuST';%'parallellogram';%'MuST2';%
DomInfo      = 0.0;

% bc = [ 0 0 0 0 ]; % 1 = Dirichlet, 0 = Neumann

NrCellRange = 16;%8:2:12;

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

Wbc0 = boundaryIntegral([0 1]);

Wbc1 = boundaryIntegral([1 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

% [Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_square(Mesh,FunctionType,[ 0 0 0 0 ]);
% 
% [Pbc,UVbc,boundary_uv,interior_uv] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,[ 1 1 1 1 ]);
% 
% [fx fy] = forcefunction(1,1,1,FunctionType,Domain,DomInfo);
% 
% g = massfunction(1,1,FunctionType,Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Sys = [


% Matrix = [ spalloc(nr_1,nr_1,0)    M1*NG                 D'*M2
%                NG'*M1               M0           spalloc(nr_0,nr_2,0)
%                 M2*D        spalloc(nr_2,nr_0,0) spalloc(nr_2,nr_2,0) ];
% 
% 
% RHS = [ -M1*[ fx ; fy ]
%          zeros(nr_0,1)
%          M2*g ];

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

RHS(ind) = [];

% % 0
% RHS(1:nr_0) = RHS(1:nr_0) - Wbc0*TangentialVelocity_bc;
% if ~isempty(boundary_w)
% RHS = RHS - Matrix(:,boundary_w)*Vorticity_bc;
% end
% 
% % 1
% RHS(nr_0+1:nr_0+nr_1) = RHS(nr_0+1:nr_0+nr_1) + Wbc1*Pbc;
% if ~isempty(boundary_uv)
% RHS = RHS - Matrix(:,boundary_uv+nr_0)*UVbc;
% end
% 
% ind = sort([boundary_w ; boundary_uv+nr_0]);
% 
% Matrix(:,ind) = [];
% Matrix(ind,:) = [];
% 
% RHS(ind,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

UVWP = Matrix\RHS;

ind1 = nr_1-length(boundary_uv);
UV_in = UVWP(1:ind1);
ind2 = ind1+nr_0-length(boundary_w);
W_in = UVWP(ind1+1:ind2);
ind3 = ind2+nr_2;
P_in = UVWP(ind2+1:ind3);

% WUVP = Matrix\RHS;
% 
% ind1 = nr_0-length(boundary_w);
% W_in = WUVP(1:ind1);
% ind2 = ind1+nr_1-length(boundary_uv);
% UV_in = WUVP(ind1+1:ind2);
% ind3 = ind2+nr_2;
% P_in = WUVP(ind2+1:ind3);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hGLp,eGLp] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hGLp,eGLp,Meshp);
ww           = reconstruct(0,WW,hGLp);
pp           = reconstruct(2,PP,eGLp,Meshp); % pp = pp-mean(mean(pp));
[ffx,ffy]    = reconstruct(1,fx,fy,hGLp,eGLp,Meshp);
gg           = reconstruct(2,g,eGLp,Meshp);

% Interpolated solution

% Exact Solution
w_ex          = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[ u_ex v_ex ] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');
p_ex          = exact_solution(Meshp.X,Meshp.Y,FunctionType,'two'); p_ex = p_ex-mean(mean(p_ex));
[fx_ex fy_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'force');
g_ex = exact_solution(Meshp.X,Meshp.Y,FunctionType,'mass');

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