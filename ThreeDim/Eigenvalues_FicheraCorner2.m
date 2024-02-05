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

NrCellRange = 9
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

Ce = curl_out_3D(N);

C = zeros(nr_2,nr_1);
M1 = zeros(nr_1);
M2 = zeros(nr_2);

for i = 1:numElements

ind1 = [ globalnr_1z(:,i) ; globalnr_1y(:,i) ; globalnr_1x(:,i) ];
ind2 = [ globalnr_2x(:,i) ; globalnr_2y(:,i) ; globalnr_2z(:,i) ];

% Gradient operator
C(ind2,ind1) = Ce;

% one-forms
M1e = innerproduct_3D(1);
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

% two-forms
M2e = innerproduct_3D(2);
M2(ind2,ind2) = M2(ind2,ind2) + M2e;

end

M1 = sparse(M1);
M2 = sparse(M2);
C  = sparse(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = M2*C/M1*C'*M2;
RHS = M2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

% [Veig,Deig] = eig(full(Matrix),full(RHS));
[Veig,Deig] = eigs(full(Matrix),full(RHS),10,2);

Deig = real(diag(Deig));

EE = 4*sort(Deig);

E = EE(abs(EE)>.001);

E(1:min(length(E),20))


save(['FicheraCorner2_eigVal_eigVect_N' num2str(N) '.mat'],'N','Veig','Deig','E','globalnr_2x','globalnr_2y','globalnr_2z','nr_2');


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