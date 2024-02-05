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

FunctionType = 'semicylinder';
Domain       = 'semicylinder';

bc = [ 1 0 0 0 ]; % 1 = Dirichlet, 0 = Neumann

N = 12;

numElements = 6;

plot_figures  = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('semicylinder')

R = 0.5;
Mesh = meshgenerator_semicylinder(R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NG = normalgradient_assembly();

M0 = innerproduct_assembly(0,Mesh);

M1 = innerproduct_assembly(1,Mesh);

Wbc = boundaryIntegral_assembly();

[PHIbc,Ubc,boundary_points,interior_points] = boundaryconditions_0_semicylinder();

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
[Meshp,hp,ep] = postproces_grid_semicylinder(R);

% Reconstruction
phi          = reconstruct(0,PHI,hp);
[qy,qx,qMag] = reconstruct(1,Qxi,Qeta,hp,ep,Meshp);
qy = -qy;

% Plotten
if plot_figures
    meshplot
    plotten
end

figure
for i=1:numElements
    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
    Qp = reshape(qMag(:,i),nn,nn);
    pcolor(Xp,Yp,Qp)
    axis equal
    axis(XYaxis)
    hold on
end
shading interp
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

% in = 'finish';
% GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%