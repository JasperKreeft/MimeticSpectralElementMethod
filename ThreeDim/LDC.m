%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% 3D Lid-driven cativy Stokes problem                                     %
%                                                                         %
% written by Jasper Kreeft (2012)                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

if ispc
    path(path,'G:\MSEM\MSEM_codes\ThreeDim\Library');
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
global globalnr_3
global nr_1 nr_2 nr_3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

N = 8;

numElements = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering_LDC;%_single;

[xi,w] = GLLnodes(N);    eta = xi;                 % Gauss-Lobotto-Legendre

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

C = curl_out_3D_assembly();
D = div_out_3D_assembly();

M1 = innerproduct_3D_assembly(1);
M2 = innerproduct_3D_assembly(2);
M3 = innerproduct_3D_assembly(3);

[tangentialVelocity,Wbc,boundary_uvw] = boundaryIntegral_3D();
% [tangentialVelocity,Wbc,boundary_uvw] = boundaryIntegral_3D_single();
interior_uvw = 1:nr_2; interior_uvw(boundary_uvw) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

disp('Matrix assembly')

Matrix = [ spalloc(nr_2,nr_2,0)      M2*C               D'*M3
                C'*M2                 M1          spalloc(nr_1,nr_3,0)
                 M3*D         spalloc(nr_3,nr_1,0) spalloc(nr_3,nr_3,0) ];


RHS = zeros(nr_1+nr_2+nr_3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

% [Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_square(Mesh,FunctionType,[ 0 0 0 0 ]);
% [Pbc,UVbc,boundary_uv,interior_uv] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,[ 0 0 0 0 ]);

% 1
RHS(nr_2+1:nr_2+nr_1) = RHS(nr_2+1:nr_2+nr_1) - Wbc*tangentialVelocity;

% 2
Matrix(:,boundary_uvw) = [];
Matrix(boundary_uvw,:) = [];

RHS(boundary_uvw,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

disp('Solve system')

UVWOP = Matrix\RHS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

disp('Postprocessen')

ind1 = nr_2-length(boundary_uvw);
UVW_in = UVWOP(1:ind1);
ind2 = ind1+nr_1;
O = UVWOP(ind1+1:ind2);
ind3 = ind2+nr_3;
P = UVWOP(ind2+1:ind3);

% O  = zeros(nr_0,1);
UVW = zeros(nr_2,1);
% P  = zeros(nr_2,1);

% O(interior_w)   = O_in;
% O(boundary_w)   = Vorticity_bc;
UVW(interior_uvw) = UVW_in;
% UVW(boundary_uvw) = UVWbc;
% P = P_in;

UU = UVW(globalnr_2x);
VV = UVW(globalnr_2y);
WW = UVW(globalnr_2z);
OOx = O(globalnr_1z);
OOy = O(globalnr_1y);
OOz = O(globalnr_1x);
PP = P(globalnr_3);
if N==1
    PP = PP';
end

divU = D*UVW;
divU = divU(globalnr_3);
% 
% curlW = NG*W;
% curlW1 = curlW(globalnr_1v);
% curlW2 = curlW(globalnr_1h);

post_LDC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

if ispc
    rmpath('G:\MSEM\MSEM_codes\ThreeDim\Library');
else
    rmpath('/media/My Passport/MSEM/MSEM_codes/ThreeDim/Library');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%