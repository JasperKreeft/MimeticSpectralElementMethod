% EXACT INTERPOLATED FLUX SOLUTION DOES NOT EXIST !!!
% NEUMANN BOUNDARY CONDITIONS NOT WORKING
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
global xi w
global h e
global nr_0 nr_1
global globalnr_0 globalnr_1v globalnr_1h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'sine';
Domain       = 'SinDeformGrid';%'LavalNozzle';%'parallellogram';%
DomInfo      = 0.0;

bc = [ 1 1 1 0 ]; % 1 = Dirichlet, 0 = Neumann

NrCellRange = 6;
HconvRange  = 1;%[3 5 7 9 11 13 15]; %[ 2 4 8 16 ];
numElements = HconvRange^2;

plot_figures  = 1;
error_figures = 0;

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

[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NGe = normalgrad(N);

Wbce = boundaryIntegral([0 1]);

NG = zeros(nr_1,nr_0);
M0 = zeros(nr_0);
M1 = zeros(nr_1);
Wbc = zeros(nr_0);
f  = zeros(nr_0,1);

for i=1:numElements

ind0 = globalnr_0(:,i);
ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];

% Gradient operator
NG(ind1,ind0) = NGe;
     
% zero-forms
M0e = innerproduct(0,Mesh.J(:,i));
M0(ind0,ind0) = M0(ind0,ind0) + M0e;

% one-forms
M1e = innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)));
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

% boundary matrix
Wbc(ind0,ind0) = Wbc(ind0,ind0) + Wbce;

M0 = sparse(M0);
M1 = sparse(M1);
NG = sparse(NG);
Wbc = sparse(Wbc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

F = forcefunction(0,Mesh.X(:,i),Mesh.Y(:,i),FunctionType);
f(ind0) = F;

end

[PHIbc,Ubc,boundary_points,interior_points] = boundaryconditions_0_square(Mesh,FunctionType,bc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

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

PHI_in = Matrix\RHS;

PPHI = zeros(nr_0,1);
PPHI(interior_points) = PHI_in;
PPHI(boundary_points) = PHIbc;

% Gradient values
QQ = NG*PPHI;

PHI  = PPHI(globalnr_0);
Qxi  = QQ(globalnr_1v);
Qeta = QQ(globalnr_1h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
phi          = reconstruct(0,PHI,hp);
[qy,qx,qMag] = reconstruct(1,Qxi,Qeta,hp,ep,Meshp);
qy = -qy;

% Interpolated solution
PHI_exact = exact_solution(Mesh.X,Mesh.Y,FunctionType,'zero');   % Nodal exact
phi_interp = reconstruct(0,PHI_exact,hp);

% Exact Solution
phi_ex        = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[qx_ex,qy_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc