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
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Instellingen

plot_figures  = 0;
error_figures = 1;

NrCellRange = 2:2:20;

FunctionType = 'exp';
Domain       = 'SinDeformGrid';
DomInfo      = 0;

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

figure; mesh(X,Y,zeros(size(X)),'EdgeColor','black'); view([0 0 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xi,N,1);

M0 = innerproduct_zeroforms(J);

M1 = innerproduct_oneforms(e,J,Qinv);

NG = normalgrad(N);

F(globalnr_0) = force_zeroform(X,Y,FunctionType);               %#ok<AGROW>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove Element- and Domain boundary Conditions                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bc = 'phi';

if strcmp(bc,'phi')

    Matrix = -NG'*M1*NG;

    RHS = M0*F';

    interior_points = globalnr_0(2:N,2:N);
    boundary_points = sort([ globalnr_0(:,1)' globalnr_0(:,end)' globalnr_0(1,2:end-1) globalnr_0(end,2:end-1) ]);

    PHIbc = exact_solution(X(boundary_points),Y(boundary_points),FunctionType,'potential')';

    RHS = RHS - Matrix(:,boundary_points)*PHIbc;

    Matrix(:,boundary_points) = [];
    Matrix(boundary_points,:) = [];

    RHS(boundary_points) = [];

elseif strcmp(bc,'curl')

    % Only 'curl', i.e. tangential velocity boundary conditions possible
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

PHI_in = Matrix\RHS;

toc
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POST-PROCESSEN                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(bc,'phi')
    PHI = zeros(nr_0,1);
    PHI(interior_points) = PHI_in;
    PHI(boundary_points) = PHIbc;
elseif strcmp(bc,'curl')
    PHI = PHI_in;
end

% Gradient values
Q = NG*PHI;

% Grid, basis-functions and weights for post-processen
postproces_grid

% Reconstruction
phi          = reconstruct_zeroforms(PHI,hGLp);
[qx,qy,velo] = reconstruct_oneforms(Q,hGLp,eGLp,Jp,dXdXip,dXdEtap,dYdXip,dYdEtap);

% Exact Solution
phi_ex          = exact_solution(Xp,Yp,FunctionType,'potential');
[ qx_ex qy_ex ] = exact_solution(Xp,Yp,FunctionType,'flux');

% phi_exact
% phi_interp = reconstruct_zeroforms(phi_exact,hGLp);

% Plotten
if plot_figures
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
    CondRef(NrCellRange) = 0.05*NrCellRange.^3;  %2.95
    c_str = '0.05N^{3}';
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