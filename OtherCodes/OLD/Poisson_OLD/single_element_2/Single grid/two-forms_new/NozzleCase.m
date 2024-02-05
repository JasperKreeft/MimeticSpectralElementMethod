clear all
close all
clc

if ispc
path(path,'O:\msem\MSEM_codes\Poisson\single_element_2\library_PoissonSingleElement')
path(path,'O:\msem\MSEM_codes\Poisson\single_element_2\Single grid\Library_SingleGrid')
elseif isunix
path(path,'/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/library_PoissonSingleElement')
path(path,'/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/Single grid/Library_SingleGrid')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N
global xi w
global globalnr_1h globalnr_1v
global nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instellingen

N = 24;

FunctionType = 'nozzle';
Domain       = 'LavalNozzle';   %'SinDeformGrid';
DomInfo      = 0.35;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numbering

[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre

[Xi,Eta,X,Y,J,Qinv,dXdXi,dXdEta,dYdXi,dYdEta] = buildgrid(xi,N,Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xi,N,1);

M1 = innerproduct_oneforms(e,J,Qinv);

M2 = innerproduct_twoforms(e,J);

D = topology(N); D = D(1:nr_2,:);

F = reshape(force_twoform(xi,eta,FunctionType,Domain,DomInfo),nr_2,1);

Wbc = boundaryIntegral(e,[1 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Element- and Domain boundary Conditions                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Matrix = [   M1          D'*M2
            M2*D  spalloc(nr_2,nr_2,0) ];

RHS = [zeros(nr_1,1) ; M2*F];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nozzle
UbcL = boundary_oneforms(FunctionType,Domain,DomInfo,'left','potential');

boundary_points = globalnr_1v(1,:);

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSEN                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid, basis-functions and weights for post-processen
postproces_grid

% Reconstruction
phi          = reconstruct_twoforms(PHI,eGLp,Jp);
[qx,qy,velo] = reconstruct_oneforms(Q,hGLp,eGLp,Jp,dXdXip,dXdEtap,dYdXip,dYdEtap);

% Exact Solution
% exact_solution
phi_ex        = exact_solution(Xp,Yp,FunctionType,'potential');
[qx_ex,qy_ex] = exact_solution(Xp,Yp,FunctionType,'flux');

% phi_exact
% phi_interp = reconstruct_twoforms(phi_exact,eGL,Jp);
[Qxi Qeta] = fluxValue(FunctionType,Domain,DomInfo);
Q_interp = zeros(nr_2,1);
Q_interp(globalnr_1v) = Qxi;
Q_interp(globalnr_1h) = Qeta;
[qx_int,qy_int,q_interp] = reconstruct_oneforms(Q_interp,hGLp,eGLp,Jp,dXdXip,dXdEtap,dYdXip,dYdEtap);

% Plotten
% Plot grid points
figure; mesh(X,Y,zeros(size(X)),'EdgeColor','black'); view([0 0 1]);

plotten

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
rmpath('O:\msem\MSEM_codes\Poisson\single_element_2\library_PoissonSingleElement')
rmpath('O:\msem\MSEM_codes\Poisson\single_element_2\Single grid\Library_SingleGrid')
elseif isunix
rmpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/library_PoissonSingleElement')
rmpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/Single grid/Library_SingleGrid')
end
