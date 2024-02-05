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
global xi
global h e
global nr_0 nr_1
global globalnr_0 globalnr_1v globalnr_1h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'cylinder';
Domain       = 'cylinder';

bc = [ 1 0 0 0 ]; % 1 = Dirichlet, 0 = Neumann

N = 8;

numElements = 12;

plot_figures  = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('cylinder')

Mesh = meshgenerator_cylinder(0.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NGe = normalgrad(N);

Wbce = boundaryIntegral([0 1]);

NG = zeros(nr_1,nr_0);
M0 = zeros(nr_0);
M1 = zeros(nr_1);
Wbc = zeros(nr_0);

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

end

[PHIbc,Ubc,boundary_points,interior_points] = boundaryconditions_0_cylinder();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = -NG'*M1*NG;
RHS = zeros(nr_0,1);

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
[Meshp,hp,ep] = postproces_grid_cylinder();

% Reconstruction
phi          = reconstruct(0,PHI,hp);
[qy,qx,qMag] = reconstruct(1,Qxi,Qeta,hp,ep,Meshp);
qy = -qy;

% Plotten
if plot_figures
    meshplot
    plotten
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%