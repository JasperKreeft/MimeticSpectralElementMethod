clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_ZeroForms/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements
global xi w
global h e
global globalnr_0 globalnr_1h globalnr_1v
global nr_0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'sine';
Domain       = 'SinDeformGrid';%'SinCosDeformGrid';%'MuST';%'LavalNozzle';%'parallellogram';%
DomInfo      = 0.2;

bc = [ 1 1 1 1 ]; % 1 = Dirichlet, 0 = Neumann

NrCellRange = 18;%3:2:25;
numElements = 1;

plot_figures  = 1;
error_figures = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

er = 0;
errorL2          = zeros(size(NrCellRange));
errorL2_interp   = zeros(size(NrCellRange));
errorL2_q        = zeros(size(NrCellRange));
errorL2_q_interp = zeros(size(NrCellRange));
ConditionNumber  = zeros(size(NrCellRange));

for N=NrCellRange

disp(['N = ' num2str(N)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square');

[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

M0 = innerproduct(0,Mesh.J);

M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

NG = normalgrad(N);

Wbc = boundaryIntegral([0 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

F(globalnr_0,1) = forcefunction(0,Mesh.X,Mesh.Y,FunctionType);

boundary_points = unique([N+1:N+1:(N+1)^2 N*(N+1)+(1:N+1)]);
periodic_points = unique([1:N+1 1:N+1:(N+1)^2]);
interior_points = (1:nr_0)';
interior_points(boundary_points) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = -NG'*M1*NG;

RHS = M0*F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

Matrix(:,periodic_points) = Matrix(:,periodic_points) + Matrix(:,boundary_points);
Matrix(periodic_points,:) = Matrix(periodic_points,:) + Matrix(boundary_points,:);

Matrix(:,boundary_points) = [];
Matrix(boundary_points,:) = [];

RHS(periodic_points) = RHS(periodic_points) + RHS(boundary_points);
RHS(boundary_points) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

PHI_in = Matrix\RHS;

PHI = zeros(nr_0,1);
PHI(interior_points) = PHI_in;
PHI(boundary_points) = PHI(periodic_points);

% Gradient values
Q = NG*PHI;

PHI = PHI(globalnr_0);
Qxi  = Q(globalnr_1v);
Qeta = Q(globalnr_1h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hGLp,eGLp] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
phi          = reconstruct(0,PHI,hGLp);
[qy,qx,qMag] = reconstruct(1,Qxi,Qeta,hGLp,eGLp,Meshp);
qy = -qy; % WHY ???

% Exact Solution
phi_ex          = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[ qx_ex qy_ex ] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');

% phi_exact
PHI_exact    = exact_solution(Mesh.X,Mesh.Y,FunctionType,'zero');  % Nodal exact
phi_interp   = reconstruct(0,PHI_exact,hGLp);
[Qxi_interp,Qeta_interp] = fluxValue(FunctionType,Domain,DomInfo);
[qx_interp,qy_interp,q_interp] = reconstruct(1,Qxi_interp,Qeta_interp,hGLp,eGLp,Meshp);

% Plotten
if plot_figures
    meshplot
    plotten
end

% error
if error_figures
    fout_Poisson
end


end % for N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for convergence plots

if length(NrCellRange)>1 && error_figures
    CondRef = 0.05*NrCellRange.^3;  %2.95
    c_str = '0.05N^{3}';
    errorplot
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%