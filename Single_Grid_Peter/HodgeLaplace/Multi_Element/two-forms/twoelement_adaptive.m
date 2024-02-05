clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_TwoForms/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N Nmax Nadaptive numElements numRows numColumns
global xi w
global h e
global nr_1 nr_2
global globalnr_1v globalnr_1h globalnr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'adaptive1';%'cosine';%'sine';%
Domain       = 'SinDeformGrid';
DomInfo      = 0.0;

bc = [ 1 1 1 1 ]; % 1 = Dirichlet, 0 = Neumann

Nadaptive = [6 12];

if Nadaptive(1)>Nadaptive(2)
    warning('Na1 > Na2')
    break
end

Nmax = max(Nadaptive);

plot_figures  = 1;
error_figures = 0;

numRows     = 1;
numColumns  = 2;

numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering_adaptive;

Mesh = meshgenerator_adaptive(Domain,DomInfo);

% meshplot_adaptive

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

D   = zeros(nr_2,nr_1);
M1  = zeros(nr_1);
M2  = zeros(nr_2);

for i=1:numElements

N = Nadaptive(i);
    
[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre
[h,e] = MimeticpolyVal(xi,N,1);

indN = 1:(N+1)^2;
ind1 = [ globalnr_1v(1:N*(N+1),i) ; globalnr_1h(1:N*(N+1),i) ];
ind2 = globalnr_2(1:N^2,i);

% Divergence operator
D(ind2,ind1) = div(N);
     
% one-forms
M1e = innerproduct(1,Mesh.J(indN,i),Mesh.Qinv(1:2*(N+1)^2,3*(i-1)+(1:3)));
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

% two-forms
M2e = innerproduct(2,Mesh.J(indN,i));
M2(ind2,ind2) = M2e;

end

M1  = sparse(M1);
M2  = sparse(M2);
D   = sparse(D);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

F   = zeros(nr_2,1);
for r=1:numRows
    for c=1:numColumns
        i = c+(r-1)*numColumns;
        N = Nadaptive(i);
        [xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre
%         [h,e] = MimeticpolyVal(xi,N,1);
        f = forcefunction(2,r,c,FunctionType,Domain,DomInfo);
        F(globalnr_2(1:N^2,i)) = f;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mortar

ME = MortarElement();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = [   M1        D'*M2
           M2*D spalloc(nr_2,nr_2,0) ];

RHS = [  zeros(nr_1,1)
            M2*F      ];
        
Matrix = ME'*Matrix*ME;
RHS    = ME'*RHS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove boundary Conditions

% RHS(1:nr_1) = RHS(1:nr_1) + Wbc*PHIbc;
% if ~isempty(boundary_flux)
% RHS = RHS - Matrix(:,boundary_flux)*Qbc;
% 
% Matrix(:,boundary_flux) = [];
% Matrix(boundary_flux,:) = [];
% 
% RHS(boundary_flux) = [];
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

QPHI = Matrix\RHS;

QPHI = ME*QPHI;

q = QPHI(1:nr_1);
phi = QPHI(nr_1+1:end);

PHI  = zeros(Nmax^2,numElements);
Qxi  = zeros(Nmax*(Nmax+1),numElements);
Qeta = zeros(Nmax*(Nmax+1),numElements);
for i=1:numElements
    N = Nadaptive(i);
    PHI(1:N^2,i) = phi(globalnr_2(1:N^2,i));
    Qxi(1:N*(N+1),i)  = q(globalnr_1v(1:N*(N+1),i));
    Qeta(1:N*(N+1),i) = q(globalnr_1h(1:N*(N+1),i));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

figure
hold on
for nel=1:numElements
    N = Nadaptive(nel);
for j=1:N
    for i=1:N
        k = i+(j-1)*(N+1);
        l = i+(j-1)*N;
        Phi = PHI(l,nel)/((Mesh.X(k+1,nel)-Mesh.X(k,nel))*(Mesh.Y(k+N+1,nel)-Mesh.Y(k,nel)));
        surf([Mesh.X(k,nel) Mesh.X(k+1,nel)],[Mesh.Y(k,nel) Mesh.Y(k+N+1,nel)],[Phi Phi ; Phi Phi])
    end
end
end


z = 0;
for nel=1:numElements
    N = Nadaptive(nel);
for j=1:N
    for i=1:N+1
        k = i+(j-1)*(N+1);
        z = z+1;
        qxi(z,1:2) = Qxi(k,nel)/(Mesh.Y(k+N+1,nel)-Mesh.Y(k,nel))*[1 1];
        x(z,1:2) = [Mesh.X(k,nel) Mesh.X(k,nel)];
        y(z,1:2) = [Mesh.Y(k,nel) Mesh.Y(k+N+1,nel)];
    end
end
end
figure
linecolorplot(qxi,x,y,min(qxi(:,1)),max(qxi(:,1)))

break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
phi          = reconstruct(2,PHI,ep,Meshp.J);
[qx,qy,qMag] = reconstruct(1,Qxi,Qeta,hp,ep,Meshp);

% Interpolated solution
% PHI_exact    = potentialValue(xi,eta,FunctionType,Domain,DomInfo); % Does not exist
% phi_interp   = reconstruct(2,PHI_exact,eGLp,Jp);
[Qxi_interp,Qeta_interp] = fluxValue(FunctionType,Domain,DomInfo);
[qx_interp,qy_interp,q_interp] = reconstruct(1,Qxi_interp,Qeta_interp,hp,ep,Meshp);

% Exact Solution
phi_ex        = exact_solution(Meshp.X,Meshp.Y,FunctionType,'two');
[qx_ex,qy_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');

% Plotten
if plot_figures
    meshplot
    plotten
end

% error
if error_figures
%     fout_singlegrid
    fout_Poisson
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for convergence plots

if error_figures
    if length(NrCellRange)>1
        CondRef = NrCellRange.^4;
        c_str = 'N^{4}';
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