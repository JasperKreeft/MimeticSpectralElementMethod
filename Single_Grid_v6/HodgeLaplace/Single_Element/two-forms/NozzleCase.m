% WERKT NIET, OA RANDVOORWAARDEN KLOPPEN NIET !!!
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

N = 12;

FunctionType = 'nozzle';
Domain       = 'LavalNozzle';
DomInfo      = 0.35;

bc = [ 1 0 0 0 ]; % 1 = Dirichlet, 0 = Neumann

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
%% Construct system matrix and righthandside

Matrix = [   M1          D'*M2
            M2*D  spalloc(nr_2,nr_2,0) ];

RHS = zeros(nr_1+nr_2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Nozzle
% Randvoorwaarden kloppen niet !!!
[PHIbc,Qbc,boundary_flux,interior_flux] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,bc);

UbcL = boundary_oneforms(1,1,FunctionType,Domain,DomInfo,'left','zero');

boundary_points = globalnr_1v(1,:);
keyboard
Ubc = zeros(nr_1,1);
Ubc(boundary_points) = UbcL;

QbcR = UbcL;
QbcB = zeros(N,1);
QbcA = zeros(N,1);

Qbc = [ QbcR' ; QbcB ; QbcA ];

interior_points = sort([ reshape(globalnr_1v(1:end-1,:),1,[]) ...
                         reshape(globalnr_1h(:,2:end-1),1,[]) ]);
boundary_points = [ globalnr_1v(end,:) globalnr_1h(:,1)' globalnr_1h(:,end)'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

RHS(1:nr_1) = RHS(1:nr_1) + Wbc*Ubc;
if ~isempty(boundary_points)
RHS = RHS - Matrix(:,boundary_points)*Qbc;

Matrix(:,boundary_points) = [];
Matrix(boundary_points,:) = [];

RHS(boundary_points) = [];
end

qphi = Matrix\RHS;

ind = nr_1-length(boundary_points);
Q_in = qphi(1:ind);
PHI  = qphi(ind+1:end);

Q = zeros(nr_1,1);
Q(interior_points) = Q_in;
Q(boundary_points) = Qbc;

PHI = PHI(globalnr_2);
Qxi  = Q(globalnr_1v);
Qeta = Q(globalnr_1h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hGLp,eGLp] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
phi          = reconstruct(2,PHI,eGLp,Meshp.J);
[qx,qy,qMag] = reconstruct(1,Qxi,Qeta,hGLp,eGLp,Meshp);

% Exact Solution
% exact_solution
phi_ex        = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zerotwo');
[qx_ex,qy_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'flux');

% phi_exact = ???
% phi_interp = reconstruct_twoforms(phi_exact,eGL,Jp);
[Qxi_interp Qeta_interp] = fluxValue(FunctionType,Domain,DomInfo);
Q_interp = zeros(nr_2,1);
[qx_int,qy_int,q_interp] = reconstruct(1,Qxi_interp,Qeta_interp,hGLp,eGLp,Meshp);

% Plotten
% Plot grid points
meshplot
plotten

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%