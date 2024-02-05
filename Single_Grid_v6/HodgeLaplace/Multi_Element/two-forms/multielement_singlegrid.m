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
global xi
global h e
global nr_1 nr_2
global globalnr_1v globalnr_1h globalnr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'sine';
Domain       = 'SinDeformGrid';
DomInfo      = 0;

bc = [ 0 0 0 1 ]; % 1 = Dirichlet, 0 = Neumann

NrCellRange = 3;
HconvRange  = [ 2 4 8 16 32 ];

plot_figures  = 0;
error_figures = 1;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Hodge-Laplace for two-forms ')
disp(['Testfunction is ' FunctionType])
disp(['Computational Domain is ' Domain ' with settings ' num2str(DomInfo)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start H-convergence loop

er = 0;
errorL2          = zeros(1,max(length(NrCellRange),length(HconvRange)));
errorL2_interp   = zeros(1,max(length(NrCellRange),length(HconvRange)));
errorL2_q        = zeros(1,max(length(NrCellRange),length(HconvRange)));
errorL2_q_interp = zeros(1,max(length(NrCellRange),length(HconvRange)));
ConditionNumber  = zeros(1,max(length(NrCellRange),length(HconvRange)));

for Hconv = HconvRange

numRows     = Hconv;
numColumns  = Hconv;
numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

for N = NrCellRange

disp(['Hconv = ' num2str(Hconv) ', N = ' num2str(N)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square')

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

M1 = innerproduct_assembly(1,Mesh);

M2 = innerproduct_assembly(2,Mesh);

D = divergence_assembly();

Wbce = boundaryIntegral([1 2]);

Wbc = zeros(nr_1);

for i=1:numElements

ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];

% boundary matrix
Wbc(ind1,ind1) = Wbc(ind1,ind1) + Wbce;

end

Wbc = sparse(Wbc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

F   = zeros(nr_2,1);
for r=1:numRows
    for c=1:numColumns
        f = forcefunction(2,r,c,FunctionType,Domain,DomInfo);
        i = c+(r-1)*numColumns;
        F(globalnr_2(:,i)) = f;
    end
end

[PHIbc,Qbc,boundary_flux,interior_flux] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,bc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

disp('assembly system matrix and righthandside')

Matrix = [   M1        D'*M2
           M2*D spalloc(nr_2,nr_2,0) ];

RHS = [  zeros(nr_1,1)
            M2*F      ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove boundary Conditions

RHS(1:nr_1) = RHS(1:nr_1) + Wbc*PHIbc;
if ~isempty(boundary_flux)
RHS = RHS - Matrix(:,boundary_flux)*Qbc;

Matrix(:,boundary_flux) = [];
Matrix(boundary_flux,:) = [];

RHS(boundary_flux) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

disp('Solve system')

QPHI = Matrix\RHS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

disp('Postprocessen')

ind = nr_1-length(boundary_flux);
Q_in = QPHI(1:ind);
PHI  = QPHI(ind+1:end);

Q = zeros(nr_1,1);
Q(interior_flux) = Q_in;
Q(boundary_flux) = Qbc;

PHI = PHI(globalnr_2);
Qxi  = Q(globalnr_1v);
Qeta = Q(globalnr_1h);

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
phi          = reconstruct(2,PHI,ep,Meshp);
[qx,qy,qMag] = reconstruct(1,Qxi,Qeta,hp,ep,Meshp);

% Interpolated solution
PHI_exact    = potentialValue(FunctionType,Domain,DomInfo);
phi_interp   = reconstruct(2,PHI_exact,ep,Meshp);
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


end % for N
end % for H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for convergence plots

if error_figures
    if length(NrCellRange)>1
        CondRef = NrCellRange.^4;
        c_str = 'N^{4}';
        filename = ['Pconv_H' num2str(Hconv) '_c' num2str(DomInfo) '.mat'];
    elseif length(HconvRange)>1
        CondRef = 0*HconvRange;
        c_str = '???';
        filename = ['Hconv_N' num2str(N) '_c' num2str(DomInfo) '.mat'];
    end
    errorplot
end

% D = [];
% save(filename,'DomInfo','c_str','N','NrCellRange','errorL2','errorL2_interp','errorL2_q','errorL2_q_interp','ConditionNumber','Linv_errorDiv','L1_errorDiv','CondRef','D','HconvRange');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%