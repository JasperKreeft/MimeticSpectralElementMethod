% EXACT INTERPOLATED FLUX SOLUTION DOES NOT EXIST !!!

clear all
close all
clc
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_ZeroForms/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements
global xi
global h e
global nr_0
global globalnr_0 globalnr_1v globalnr_1h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'sine';
Domain       = 'SinDeformGrid';%'LavalNozzle';%'parallellogram';%
DomInfo      = 0.0;

bc = [ 1 1 1 1 ]; % 1 = Dirichlet, 0 = Neumann

NrCellRange = 3;
HconvRange  = 5;%[3 5 7 9 11 13 15]; %[ 2 4 8 16 ];
numElements = HconvRange^2;

plot_figures  = 1;
error_figures = 0;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Hodge-Laplace for zero-forms ')
disp(['Testfunction is ' FunctionType])
disp(['Computational Domain is ' Domain ' with settings ' num2str(DomInfo)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start H-convergence loop

er = 0;
errorL2        = zeros(1,max(length(NrCellRange),length(HconvRange)));
errorL2_interp = zeros(1,max(length(NrCellRange),length(HconvRange)));
errorL2_q      = zeros(1,max(length(NrCellRange),length(HconvRange)));
errorL2_q_interp = zeros(1,max(length(NrCellRange),length(HconvRange)));

for Hconv = HconvRange

numRows     = Hconv;
numColumns  = Hconv;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

for N = NrCellRange

disp(['Hconv = ' num2str(Hconv) ', N = ' num2str(N)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square');

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NG = normalgradient_assembly();

M0 = innerproduct_assembly(0,Mesh);

M1 = innerproduct_assembly(1,Mesh);

Wbc = boundaryIntegral_assembly();

f = forcefunction_assembly(0,Mesh,FunctionType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary conditions

[PHIbc,Ubc,boundary_points,interior_points] = boundaryconditions_0_square(Mesh,FunctionType,bc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

disp('assembly system matrix and righthandside')

Matrix = -NG'*M1*NG;
RHS = M0*f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove boundary Conditions

RHS = RHS - Wbc*Ubc;

if ~isempty(boundary_points)
RHS = RHS - Matrix(:,boundary_points)*PHIbc;
end

Matrix(:,boundary_points) = [];
Matrix(boundary_points,:) = [];

RHS(boundary_points) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

disp('Solve system')

PHI_in = Matrix\RHS;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

disp('Postprocessen')

PPHI = zeros(nr_0,1);
PPHI(interior_points) = PHI_in;
PPHI(boundary_points) = PHIbc;

% Gradient values
QQ = NG*PPHI;

PHI  = PPHI(globalnr_0);
Qxi  = QQ(globalnr_1v);
Qeta = QQ(globalnr_1h);

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
phi          = reconstruct(0,PHI,hp);
[qy,qx,qMag] = reconstruct(1,Qxi,Qeta,hp,ep,Meshp);
qy = -qy;

% Interpolated solution
PHI_exact = exact_solution(Mesh.X,Mesh.Y,FunctionType,'zero');   % Nodal exact
phi_interp = reconstruct(0,PHI_exact,hp);
[Qxi_interp,Qeta_interp] = fluxValue(FunctionType,Domain,DomInfo);

% Exact Solution
phi_ex        = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[qx_ex,qy_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');

% Plotten
if plot_figures
    disp('creation colorplots')
    meshplot
    plotten
end

% error
if error_figures
    disp('calculation of errors')
    fout_Poisson
end


end % for N
end % for H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for convergence plots

if error_figures
    if length(NrCellRange)>1
        CondRef = NrCellRange.^3;  %2.95
        c_str = 'N^{3}';
    elseif length(HconvRange)>1
        CondRef = 0*HconvRange;
        c_str = '???';
    end
   errorplot
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

disp('Calculations finished')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc