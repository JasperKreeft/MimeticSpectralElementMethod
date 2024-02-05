close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N
global xi w
global h e
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square')

[xi,w] = GLLnodes(N);     eta = xi;  % Gauss-Lobotto-Legendre

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

NG = normalgrad(N);

D  = div(N);

[h,e] = MimeticpolyVal(xi,N,1);

M0 = innerproduct(0,Mesh.J);

M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

M2 = innerproduct(2,Mesh.J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

[~,~,boundary,interior] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,[ 0 0 0 0 ]);
boundary = sort(boundary);



Matrix = [ M0 ; M1*NG ];


