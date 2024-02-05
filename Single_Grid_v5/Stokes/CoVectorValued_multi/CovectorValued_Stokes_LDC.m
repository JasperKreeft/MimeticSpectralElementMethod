%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Multi-element Covector-Valued Stokes problem                            %
%                                                                         %
% written by Jasper Kreeft (2012)                                         %
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
global xi w e
global nr_1 nr_2 nr_11
global globalnr_1h globalnr_1v globalnr_2
global globalnr_11_xx globalnr_11_yx globalnr_11_xy globalnr_11_yy
global nn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'LidDrivenCavity';
Domain       = 'SinDeformGrid';
DomInfo      = 0.0;

plot_figures  = 1;

N = 16;
H = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numRows    = H;
numColumns = H;
numElements = numRows*numColumns;

[xi,w] = GLLnodes(N);

Mesh = meshgenerator_square(Domain,DomInfo);

numbering_CoVect('square')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[~,e] = MimeticpolyVal(xi,N,1);

D = divergence_assembly();

CD = Covect_div_assembly();

T = Tensorinnerproduct_assembly();

% Tr = Trace_assembly();
Tr = Trace();

M1 = innerproduct_assembly(1,Mesh);

M2 = innerproduct_assembly(2,Mesh);

[Wbc,U_lid] = toplid();

I = Inclusion(Mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = [          T                   CD'*Tr'*M1          spalloc(nr_11,nr_2,0)
                M1*Tr*CD            spalloc(nr_1,nr_1,0)          D'*M2
           spalloc(nr_2,nr_11,0)           M2*D            spalloc(nr_2,nr_2,0)    ];

% Matrix = [          T                   CD'*Tr'*M1          spalloc(nr_11,nr_2,0)
%                 M1*Tr*CD            eps*speye(nr_1)          D'*M2
%            spalloc(nr_2,nr_11,0)           M2*D            eps*speye(nr_2)    ];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

RHS = [   Wbc*U_lid
         zeros(nr_1,1)
         zeros(nr_2,1) ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

[Pbc,UVbc,boundary,interior] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,[0 0 0 0]);
boundary = boundary + nr_11;

Matrix(:,boundary) = [];
Matrix(boundary,:) = [];

RHS(boundary,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

tau_u_p = Matrix\RHS;

Tau = tau_u_p(1:nr_11);
UVin = tau_u_p(nr_11+(1:nr_1-length(boundary)));
UV = zeros(nr_1,1);
UV(interior) = UVin;
UV(boundary-nr_11) = UVbc;
P = tau_u_p(nr_11+nr_1-length(boundary)+(1:nr_2));

UU = UV(globalnr_1v);
VV = UV(globalnr_1h);
PP = P(globalnr_2);

Txx = Tau(globalnr_11_xx);
Tyx = Tau(globalnr_11_yx);
Txy = Tau(globalnr_11_xy);
Tyy = Tau(globalnr_11_yy);

Pi = I*P;
Pi = Pi/Mesh.J(1,1);  % Include uniform linear scaling

Pixx = Pi(1:N*(N+2));
Piyx = Pi(N*(N+2)+(1:(N+1)^2));
Piyy = Pi(N*(N+2)+(N+1)^2+(1:N*(N+2)));
Pixy = Pi(2*N*(N+2)+(N+1)^2+(1:(N+1)^2));

M = CD*(-Tau+Pi);
Mx = M(1:N*(N+1));
My = M(N*(N+1)+(1:N*(N+1)));
TrM = Tr*M;
TrMx = TrM(globalnr_1v);
TrMy = TrM(globalnr_1h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hglp,eglp] = postproces_grid_square(Domain,DomInfo);
[hegp,eegp] = MimeticpolyVal(Meshp.xip,N,3);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hglp,eglp,Meshp);
pp           = reconstruct(2,PP,eglp,Meshp); %pp = pp-mean(mean(pp));
[txx,tyx,txy,tyy] = reconstruct(11,Txx,Tyx,Txy,Tyy,hglp,eglp,hegp,eegp,Meshp);
[pxx,pyx,pxy,pyy] = reconstruct(11,Pixx,Piyx,Pixy,Piyy,hglp,eglp,hegp,eegp,Meshp);
[tmx,tmy,mm] = reconstruct(1,TrMx,TrMy,hglp,eglp,Meshp);
[mmx,mmy] = reconstruct(21,Mx,My,eglp,eegp,Meshp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotten

if plot_figures
    disp('creation colorplots')
%     meshplot
    plotten_me
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%