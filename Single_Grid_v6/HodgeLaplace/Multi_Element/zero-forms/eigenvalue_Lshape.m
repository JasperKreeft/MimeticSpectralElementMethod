clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_ZeroForms/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns
global xi w
global h e
global nr_0 nr_1
global globalnr_0 globalnr_1v globalnr_1h

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'CurlCurl';
Domain       = 'Lshape';
DomInfo      = 0.;

NrCellRange = 5; %2:2:12;
numElements = 3;
numRows = 2;
numColumns = 2;

error_figures = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

er = 0;
error = zeros(10);

for N = NrCellRange

disp(['N = ' num2str(N)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square');

nr_0 = 3*(N+1)^2-2*(N+1);
nr_1 = 3*2*N*(N+1)-2*N;

[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NGe = normalgrad(N);

NG = zeros(nr_1,nr_0);
M0 = zeros(nr_0);
M1 = zeros(nr_1);

R = [1 1 2];
C = [1 2 1];
for i = 1:3;
    r = R(i);
    c = C(i);

ind0 = globalnr_0(:,i);
ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];

% Gradient operator
NG(ind1,ind0) = NGe;
     
% zero-forms
M0e = innerproduct(0,1);
M0(ind0,ind0) = M0(ind0,ind0) + M0e;

% one-forms
M1e = innerproduct(1,1,1);
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

M0 = sparse(M0);
M1 = sparse(M1);
NG = sparse(NG);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = NG'*M1*NG;
RHS = M0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

E = 4*sort(eig(full(Matrix),full(RHS)));

E(1:min(length(E),20))

end % for N


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%