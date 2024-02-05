
clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_TwoForms/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns
global xi w
global h e
global nr_1 nr_2
global globalnr_1v globalnr_1h globalnr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'nozzle';
Domain       = 'LavalNozzle';
DomInfo      = 0.35;

N = 8;

numRows     = 2;
numColumns  = 2;
numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square')

[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

De = div(N);

Wbce = boundaryIntegral([1 2]);

D   = zeros(nr_2,nr_1);
M1  = zeros(nr_1);
M2  = zeros(nr_2);
Wbc = zeros(nr_1);
f   = zeros(nr_2,1);

for i=1:numElements

ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
ind2 = globalnr_2(:,i);

% Divergence operator
D(ind2,ind1) = De;
     
% one-forms
M1e = innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)));
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

% two-forms
M2e = innerproduct(2,Mesh.J(:,i));
M2(ind2,ind2) = M2e;

% boundary matrix
Wbc(ind1,ind1) = Wbc(ind1,ind1) + Wbce;

end

M1  = sparse(M1);
M2  = sparse(M2);
D   = sparse(D);
Wbc = sparse(Wbc);



% Nozzle boundary conditions
ind1 = 1:N+1:N*(N+1);
ind2 = 1:numColumns:numElements;
boundary_points = globalnr_1v(ind1,ind2);
ind1 = 1:N+1:N*(N+1);
ind2 = N+1:N+1:N*(N+1);
ind3 = numColumns:numColumns:numElements;
ind4 = 1:numColumns;
ind5 = (numRows-1)*numColumns+(1:numColumns);
boundary_flux = [ reshape(globalnr_1v(ind2,ind3),1,[]) ...
                  reshape(globalnr_1h(ind1,ind4),1,[]) ...
                  reshape(globalnr_1h(ind2,ind5),1,[]) ];

interior_flux = 1:nr_1; interior_flux(boundary_flux) = [];

PHIbcL = [];
for r=1:numRows
PHIbcL_ = (1/numRows)*boundary_oneforms(r,1,FunctionType,Domain,DomInfo,'left','two');
PHIbcL = [ PHIbcL ; PHIbcL_];
end
PHIbc = zeros(nr_1,1);
PHIbc(boundary_points) = PHIbcL;

QbcR = PHIbcL;
QbcB = zeros(N*numColumns,1);
QbcA = zeros(N*numColumns,1);

Qbc = [ QbcR ; QbcB ; QbcA ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = [   M1        D'*M2
           M2*D spalloc(nr_2,nr_2,0) ];

RHS = zeros(nr_1+nr_2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Element- and Domain boundary Conditions

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
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
phi          = reconstruct(2,PHI,ep,Meshp);
[qx,qy,qMag] = reconstruct(1,Qxi,Qeta,hp,ep,Meshp);

% Interpolated solution
% PHI_exact    = potentialValue(xi,eta,FunctionType,Domain,DomInfo); % Does not exist
% phi_interp   = reconstruct(2,PHI_exact,eGLp,Jp);
global nn; 
phi_interp = zeros(nn^2,numElements);

% Exact Solution
phi_ex          = exact_solution(Meshp.X,Meshp.Y,FunctionType,'two');
[ qx_ex qy_ex ] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');

% Plotten
meshplot
plotten

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%