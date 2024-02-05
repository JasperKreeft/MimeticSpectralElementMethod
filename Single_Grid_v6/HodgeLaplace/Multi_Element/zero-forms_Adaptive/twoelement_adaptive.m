clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

% in = 'start';                                                   %#ok<NASGU>
% run Library_ZeroForms/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N Nmax Nadaptive numElements numRows numColumns
global xi w
global h e
global nr_0 nr_1
global globalnr_0 globalnr_1v globalnr_1h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'adaptive1';%'cosine';%'sine';%
Domain       = 'SinDeformGrid_x';%'CosDeformGrid';%'SinDeformGrid';
DomInfo      = 0.2;

bc = [ 1 1 1 1 ]; % 1 = Dirichlet, 0 = Neumann

Nadaptive = [4 8];

Nmax = max(Nadaptive);
Nmin = min(Nadaptive);

plot_figures  = 1;
error_figures = 0;

numRows     = 1;
numColumns  = 2;

numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering_adaptive;

Mesh = meshgenerator_adaptive(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

NG = zeros(nr_1,nr_0);
M0 = zeros(nr_0);
M1 = zeros(nr_1);
f  = zeros(nr_0,1);

for i=1:numElements

N = Nadaptive(i);
    
[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre
[h,e] = MimeticpolyVal(xi,N,1);

indN = 1:(N+1)^2;
ind0 = globalnr_0(indN,i);
ind1 = [ globalnr_1v(1:N*(N+1),i) ; globalnr_1h(1:N*(N+1),i) ];

% Gradient operator
NG(ind1,ind0) = normalgrad(N);
     
% zero-forms
M0e = innerproduct(0,Mesh.J(indN,i));
M0(ind0,ind0) = M0(ind0,ind0) + M0e;

% one-forms
M1e = innerproduct(1,Mesh.J(indN,i),Mesh.Qinv(1:2*(N+1)^2,3*(i-1)+(1:3)));
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

F = force_zeroform(Mesh.X(indN,i),Mesh.Y(indN,i),FunctionType);
f(ind0) = F;

end

M0 = sparse(M0);
M1 = sparse(M1);
NG = sparse(NG);

ME = MortarElement();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = ME'*(-NG'*M1*NG)*ME;
RHS = ME'*(M0*f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove boundary Conditions

Na1 = Nadaptive(1);
Na2 = Nadaptive(2);

boundary_points = [ 1:Na1 ...
                    Na1+1:Na1:Na1^2 ...
                    Na1*Na1+1:Na1*(Na1+1) ...
                    Na1*(Na1+1)+(1:Na2) ...
                    Na1*(Na1+1)+(2*Na2:Na2:Na2^2) ...
                    Na1*(Na1+1)+(Na2^2+1:Na2*(Na2+1)) ...
                    Na1*(Na1+1)+Na2*(Na2+1)+[1 Nmin+1] ];

interior_points = 1:(Na1*(Na1+1)+Na2*(Na2+1)+Nmin+1);
interior_points(boundary_points) = [];

Matrix(:,boundary_points) = [];
Matrix(boundary_points,:) = [];

RHS(boundary_points) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

PHI_in = Matrix\RHS;

PPHI = zeros(Na1*(Na1+1)+Na2*(Na2+1)+Nmin+1,1);
PPHI(interior_points) = PHI_in;
PPHI(boundary_points) = zeros(size(boundary_points));

PPHI = ME*PPHI;

PHI = zeros((Nmax+1)^2,2);
for i=1:2
    ind = 1:(Nadaptive(i)+1)^2;
    PHI(ind,i)  = PPHI(globalnr_0(ind,i));
end
figure
dotcolorplot(PHI,Mesh.X,Mesh.Y)

% Gradient values
QQ = NG*PPHI;

Qxi  = zeros(Nmax*(Nmax+1),2);
Qeta = zeros(Nmax*(Nmax+1),2);
for i=1:2
    N = Nadaptive(i);
    Qxi(1:N*(N+1),i)  = QQ(globalnr_1v(1:N*(N+1),i));
    Qeta(1:N*(N+1),i) = QQ(globalnr_1h(1:N*(N+1),i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square_adaptive(Domain,DomInfo);

% Reconstruction
phi          = reconstruct_adaptive(0,PHI,hp);

[qy,qx,qMag] = reconstruct_adaptive(1,Qxi,Qeta,hp,ep,Meshp);
qy = -qy;

% Interpolated solution
PHI_exact = exact_solution(Mesh.X,Mesh.Y,FunctionType,'zero');   % Nodal exact
phi_interp = reconstruct_adaptive(0,PHI_exact,hp);
[Qxi_interp,Qeta_interp] = fluxValue(FunctionType,Domain,DomInfo);
[qy_interp,qx_interp,qMag_interp] = reconstruct_adaptive(1,Qxi_interp,Qeta_interp,hp,ep,Meshp);
qy_interp = -qy_interp;

% Exact Solution
phi_ex        = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[qx_ex,qy_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');

% Plotten
if plot_figures
    meshplot_adaptive
    plotten
end

% error
if error_figures
    fout_Poisson
end

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

% in = 'finish';
% GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
