%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Fichera Corner curl-curl problem                                        %
%                                                                         %
% written by Jasper Kreeft (2011)                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

if ispc
    path(path,'O:\MSEM\MSEM_codes\ThreeDim\Library');
else
    path(path,'/media/My Passport/MSEM/MSEM_codes/ThreeDim/Library');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements
global xi w
global h e
global globalnr_1x globalnr_1y globalnr_1z
global globalnr_2x globalnr_2y globalnr_2z
global nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

NrCellRange = 14;
% N = 2;

numElements = 7;
% numRows = 2;
% numColumns = 2;

% error_figures = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

% er = 0;
% error = zeros(10);

for N = NrCellRange

disp(['N = ' num2str(N)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering_FicheraCorner;

[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

C = curl_out_3D_assembly();

M1 = innerproduct_3D_assembly(1);
M2 = innerproduct_3D_assembly(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

disp('Matrix assembly')

Matrix = C'*M2*C;
RHS = M1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

disp('Solve matrix system')

% [Veig,Deig] = eig(full(Matrix),full(RHS));
[Veig,Deig] = eigs(Matrix,RHS,10,2);%20,3);

Deig = real(diag(Deig));

EE = 4*sort(Deig);

E = EE(abs(EE)>.001);

E(1:min(length(E),20))


save(['FicheraCorner_eigVal_eigVect_N' num2str(N) '.mat'],'N','Veig','Deig','E','globalnr_1x','globalnr_1y','globalnr_1z','nr_1');


% er = er+1;
% error(er,1:min(length(E),10)) = E(1:min(length(E),10));

end % for N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

if ispc
    rmpath('O:\MSEM\MSEM_codes\ThreeDim\Library');
else
    rmpath('/media/My Passport/MSEM/MSEM_codes/ThreeDim/Library');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%