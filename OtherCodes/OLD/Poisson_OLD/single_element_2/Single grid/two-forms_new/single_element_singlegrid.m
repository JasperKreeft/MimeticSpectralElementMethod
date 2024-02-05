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
global h e
global globalnr_1h globalnr_1v
global nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instellingen

plot_figures  = 1;
error_figures = 0;

NrCellRange = 12;

FunctionType = 'sine';
Domain       = 'SinDeformGrid';
DomInfo      = 0.0;

bc = [ 1 1 1 1 ]; % 1 = Dirichlet, 0 = Neumann

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

errorL2 = zeros(size(NrCellRange));
ConditionNumber = zeros(size(NrCellRange));

for N=NrCellRange

disp(['N = ' num2str(N)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid generation                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numbering

[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre

[Xi,Eta,X,Y,J,Qinv,dXdXi,dXdEta,dYdXi,dYdEta] = buildgrid(xi,N,Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Topology and Metric matrices                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xi,N,1);

M1 = innerproduct(1,J,Qinv);

M2 = innerproduct(2,J);

D = topology(N); D = D(1:nr_2,:);

Wbc = boundaryIntegral(e,[1 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boundary conditions and Force vector                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[PHIbc,Qbc,boundary_flux,interior_flux] = boundaryconditions(FunctionType,Domain,DomInfo,bc);

F = reshape(force_twoform(xi,eta,FunctionType,Domain,DomInfo),nr_2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrix and righthand-side construction                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Matrix = [   M1          D'*M2
            M2*D  spalloc(nr_2,nr_2,0) ];

RHS = [zeros(nr_1,1) ; M2*F];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Element- and Domain boundary Conditions                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RHS(1:nr_1) = RHS(1:nr_1) + Wbc*PHIbc;
if ~isempty(boundary_flux)
RHS = RHS - Matrix(:,boundary_flux)*Qbc;

Matrix(:,boundary_flux) = [];
Matrix(boundary_flux,:) = [];

RHS(boundary_flux) = [];
end

qphi = Matrix\RHS;

ind = nr_1-length(boundary_flux);
Q_in = qphi(1:ind);
PHI  = qphi(ind+1:end);

Q = zeros(nr_1,1);
Q(interior_flux) = Q_in;
Q(boundary_flux) = Qbc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSEN                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Grid, basis-functions and weights for post-processen
postproces_grid

% Reconstruction
phi          = reconstruct_twoforms(PHI,eGLp,Jp);
[qx,qy,velo] = reconstruct_oneforms(Q,hGLp,eGLp,Jp,dXdXip,dXdEtap,dYdXip,dYdEtap);

% Exact Solution
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
if plot_figures
    % Plot grid points
    figure; mesh(X,Y,zeros(size(X)),'EdgeColor','black'); view([0 0 1]);
    
    plotten
end

% error
if error_figures
    fout_singlegrid
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end % for N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(NrCellRange)>1 && error_figures
    CondRef(NrCellRange) = NrCellRange.^4;
    c_str = 'N^{4}';
    errorplot
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ispc
rmpath('O:\msem\MSEM_codes\Poisson\single_element_2\library_PoissonSingleElement')
rmpath('O:\msem\MSEM_codes\Poisson\single_element_2\Single grid\Library_SingleGrid')
elseif isunix
rmpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/library_PoissonSingleElement')
rmpath('/media/FREECOM HDD/msem/MSEM_codes/Poisson/single_element_2/Single grid/Library_SingleGrid')
end
