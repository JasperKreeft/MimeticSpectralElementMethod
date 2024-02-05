clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_TwoForms/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N
global xi w
global h e
global globalnr_1h globalnr_1v globalnr_2
global nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'sine';
Domain       = 'SinDeformGrid';%'LavalNozzle';%'parallellogram';
DomInfo      = 0.2;

bc = [ 1 0 0 0 ]; % 1 = Dirichlet, 0 = Neumann

NrCellRange = 24%2:2:20;%4:4:40;

plot_figures  = 1;
error_figures = 1;

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

M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

M2 = innerproduct(2,Mesh.J);

D = div(N);

Wbc = boundaryIntegral([1 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

[PHIbc,Qbc,boundary_flux,interior_flux] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,bc);

F = forcefunction(2,1,1,FunctionType,Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = [   M1          D'*M2
            M2*D  spalloc(nr_2,nr_2,0) ];

RHS = [zeros(nr_1,1) ; M2*F];

A = M2*D*inv(M1)*D'*M2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

RHS(1:nr_1) = RHS(1:nr_1) + Wbc*PHIbc;
if ~isempty(boundary_flux)
RHS = RHS - Matrix(:,boundary_flux)*Qbc;

Matrix(:,boundary_flux) = [];
Matrix(boundary_flux,:) = [];

RHS(boundary_flux) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

QPHI = Matrix\RHS;

ind = nr_1-length(boundary_flux);
Q_in = QPHI(1:ind);
PHI  = QPHI(ind+1:end);

Q = zeros(nr_1,1);
Q(interior_flux) = Q_in;
Q(boundary_flux) = Qbc;

PHI = PHI(globalnr_2);
Qxi  = Q(globalnr_1v);
Qeta = Q(globalnr_1h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hGLp,eGLp] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
phi          = reconstruct(2,PHI,eGLp,Meshp);
[qx,qy,qMag] = reconstruct(1,Qxi,Qeta,hGLp,eGLp,Meshp);

% Interpolated solution
PHI_exact    = potentialValue(xi,eta,FunctionType,Domain,DomInfo); % Nodal exact
phi_interp   = reconstruct(2,PHI_exact,eGLp,Meshp);
[Qxi_interp Qeta_interp] = fluxValue(FunctionType,Domain,DomInfo);
[qx_interp,qy_interp,q_interp] = reconstruct(1,Qxi_interp,Qeta_interp,hGLp,eGLp,Meshp);

% Exact Solution
phi_ex        = exact_solution(Meshp.X,Meshp.Y,FunctionType,'two');
[qx_ex,qy_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');

% Plotten
if plot_figures
    % Plot grid points
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
    CondRef = NrCellRange.^4;
    c_str = 'N^{4}';
    errorplot
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%