clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_ZeroForms/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns
global xi w
global h e
global nr_0 nr_1
global globalnr_0 globalnr_1v globalnr_1h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'nozzle';
Domain       = 'LavalNozzle';
DomInfo      = 0.25;

bc = [ 1 0 0 0 ]; % 1 = Dirichlet, 0 = Neumann

N = 2;

numRows     = 2;
numColumns  = 2;

numElements = numRows*numColumns;

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

ind1 = globalnr_0(:,i);
ind2 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];

% Gradient operator
NG(ind2,ind1) = NGe;
     
% zero-forms
M0e = innerproduct(0,Mesh.J(:,i));
M0(ind1,ind1) = M0(ind1,ind1) + M0e;

% one-forms
M1e = innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)));
M1(ind2,ind2) = M1(ind2,ind2) + M1e;

% boundary matrix
Wbc(ind1,ind1) = Wbc(ind1,ind1) + Wbce;

M0 = sparse(M0);
M1 = sparse(M1);
NG = sparse(NG);
Wbc = sparse(Wbc);

end

[PHIbc,Ubc,boundary_points,interior_points] = boundaryconditions_0_square(Mesh,FunctionType,bc);

Ubc = zeros(nr_0,1);
ind1 = N+1:N+1:(N+1)^2;
ind2 = numColumns:numColumns:numElements;
Ubc(globalnr_0(ind1,ind2)) = 1/numRows;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = -NG'*M1*NG;
RHS = zeros(nr_0,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove boundary Conditions

RHS = RHS - Wbc*Ubc;

RHS = RHS - Matrix(:,boundary_points)*PHIbc;

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
phi_ex          = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[ qx_ex qy_ex ] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');


meshplot
plotten


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%