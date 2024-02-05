%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Multi-element curl-curl problem                                         %
%                                                                         %
% written by Jasper Kreeft (2010)                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns
global xi
global h e
global cc
global globalnr_0 globalnr_1v globalnr_1h
global nr_0 nr_1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

HconvRange = 4;%1:10;

NrCellRange = 4;%1:10;

cc = 0.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start H-convergence loop

error = zeros(10);
er = 0;

for Hconv = HconvRange


numRows     = Hconv;
numColumns  = Hconv;
numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

for N=NrCellRange

disp(['Hconv = ' num2str(Hconv) ', N = ' num2str(N)]);

N2 = N*N;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square');

Mesh = meshgenerator_square('CurlCurl',cc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NGe = normalgrad(N);

NG = zeros(nr_1,nr_0);
M0 = zeros(nr_0);
M1 = zeros(nr_1);
for i=1:numElements

ind0 = globalnr_0(:,i);
ind1 = [ globalnr_1v(:,i) globalnr_1h(:,i) ];

% Normal Gradient
NG(ind1,ind0) = NGe;

% zero-forms
M0e = innerproduct(0,Mesh.J(:,i));
M0(ind0,ind0) = M0(ind0,ind0) + M0e;

% one-forms
M1e = innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)));
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

end

M0 = sparse(M0);
M1 = sparse(M1);
NG = sparse(NG);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and solve system
Neumann_bc
keyboard
L = M1*NG/M0*NG'*M1;

EE = sort(eig(full(L),full(M1)));

E = EE(abs(EE)>.2);

exact = [1 1 2 4 4 5 5 8 9 9]';
nr = min(length(E),10);
er = er+1;
error(1:nr,er) = abs(E(1:nr)-exact(1:nr));


clear L M1 M1e M0 Qinv NG
end % for N
end % for H

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code for convergence plots

% plottenConv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%