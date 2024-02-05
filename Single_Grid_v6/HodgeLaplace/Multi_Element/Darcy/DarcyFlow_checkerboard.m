clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_Darcy/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns nn
global xi
global h e
global nr_1 nr_2
global globalnr_1v globalnr_1h globalnr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'checkerboard';%'wheelerxueyotov2012b';%'darcy_1';
Domain       = 'SinDeformGrid_01';%'wheelerxueyotoy2012b';%'SinCosDeformGrid';%
DomInfo      = 0.0;

bc = [ 0 1 0 0 ];

NrCellRange = 8;
HconvRange = 15;%[ 4 8 16 32 64];

plot_figures  = 0;
error_figures = 0;
Tecplot = 1;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Darcy for two-forms ')
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

M1 = innerproduct_Darcy_assembly(1,Mesh);

M2 = innerproduct_assembly(2,Mesh);

D = divergence_assembly();

Wbc = boundaryIntegral_assembly_1();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

F   = zeros(nr_2,1);
for r=1:numRows
    for c=1:numColumns
        f = forcefunction(2,r,c,FunctionType,Domain,DomInfo,Mesh);
        i = c+(r-1)*numColumns;
        F(globalnr_2(:,i)) = f;
    end
end

[PHIbc,Qbc,boundary_flux,interior_flux] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,bc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

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
PHI_exact    = potentialValue(FunctionType,Domain,DomInfo);
% PHI_exact    = potentialValue(xi,eta,FunctionType,Domain,DomInfo); % Does not exist
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

% % error
% if error_figures
%     fout_Poisson
% end


end % for N
end % for H

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Code for convergence plots
% 
% if error_figures
%     if length(NrCellRange)>1
%         CondRef = NrCellRange.^4;
%         c_str = 'N^{4}';
%         filename = ['Pconv_H' num2str(Hconv) '.mat'];
%     elseif length(HconvRange)>1
%         CondRef = 0*HconvRange;
%         c_str = '???';
%         filename = ['Hconv_N' num2str(N) '.mat'];
%     end
%     errorplot
% end
% 
% D = [];
% save(filename,'DomInfo','c_str','N','NrCellRange','errorL2','errorL2_interp','errorL2_q','errorL2_q_interp','ConditionNumber','Linv_errorDiv','L1_errorDiv','CondRef','D','HconvRange');


% plotten











if Tecplot

    data = zeros(nn^2*numElements,5);
    corners = zeros((nn-1)^2*numElements,4);
    for nel=1:numElements

        for i=1:nn-1
            for k=1:nn-1

                iknel = i+(k-1)*(nn-1)+(nel-1)*(nn-1)^2;

                ll = i+(k-1)*nn; % lowerleft corner

                corners(iknel,:) = [ ll ll+1 ll+nn+1 ll+nn ] + (nel-1)*nn^2;

            end
        end


        ind = (1:nn^2)+(nel-1)*nn^2;
        data(ind,:) = [ Meshp.X(:,nel) Meshp.Y(:,nel) qx(:,nel) qy(:,nel) qMag(:,nel) ];

        name = ['checkerboard_N' num2str(N)];

    end

    MatlabToTecplot('FE',name,name,'"X" "Y" "U" "V" "Velo"',[numElements*nn^2 numElements*(nn-1)^2],data,2,corners);

end

% keyboard
% data = [ reshape(Meshp.X,[],1) reshape(Meshp.Y,[],1) reshape(qx,[],1) reshape(qy,[],1) reshape(qMag,[],1) ];
% MatlabToTecplot('checkerboard','checkerboard','"X" "Y" "U" "V" "Velo"',[nn nn],data,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Close libraries
% 
% in = 'finish';
% GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%