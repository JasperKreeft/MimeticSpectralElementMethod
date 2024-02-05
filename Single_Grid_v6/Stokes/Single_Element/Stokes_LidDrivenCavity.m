%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Single-element Stokes solver for Lid-Driven Cavity problem              %
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

global N nn
global xi
global h e
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'LidDrivenCavity';
Domain       = 'SinDeformGrid'; %'MuST'; %'CosDeformGrid'; % 
DomInfo      = 0.0;

N = 14;

plot_figures  = 1;
casenr        = FunctionType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square')

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

NG = normalgrad(N);

D  = div(N);

[h,e] = MimeticpolyVal(xi,N,1);

M0 = innerproduct(0,Mesh.J);

M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

M2 = innerproduct(2,Mesh.J);

Wbc0 = boundaryIntegral([0 1]);

Wbc1 = boundaryIntegral([1 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

[Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_square(Mesh,FunctionType,[ 0 0 0 0 ]);

[Pbc,UVbc,boundary_uv,interior_uv] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,[ 0 0 0 0 ]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = [ spalloc(nr_1,nr_1,0)    M1*NG               D'*M2
               NG'*M1                M0          spalloc(nr_0,nr_2,0)
                M2*D        spalloc(nr_2,nr_0,0) spalloc(nr_2,nr_2,0) ];


RHS = zeros(nr_1+nr_0+nr_2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

% 0
RHS(nr_1+1:nr_1+nr_0) = RHS(nr_1+1:nr_1+nr_0) - Wbc0*TangentialVelocity_bc;
if ~isempty(boundary_w)
RHS = RHS - Matrix(:,boundary_w)*Vorticity_bc;
end

% 1
RHS(1:nr_1) = RHS(1:nr_1) + Wbc1*Pbc;
if ~isempty(boundary_uv)
RHS = RHS - Matrix(:,boundary_uv)*UVbc;
end

ind = sort([boundary_uv boundary_w+nr_1]);

Matrix(:,ind) = [];
Matrix(ind,:) = [];

RHS(ind) = [];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

UVWP = Matrix\RHS;

ind1 = nr_1-length(boundary_uv);
UV_in = UVWP(1:ind1);
ind2 = ind1+nr_0-length(boundary_w);
W_in = UVWP(ind1+1:ind2);
ind3 = ind2+nr_2;
P_in = UVWP(ind2+1:ind3);

UV = zeros(nr_1,1);
W  = zeros(nr_0,1);
P  = zeros(nr_2,1);

UV(interior_uv) = UV_in;
UV(boundary_uv) = UVbc;
W(interior_w)   = W_in;
W(boundary_w)   = Vorticity_bc;
P = P_in;

UU = UV(globalnr_1v);
VV = UV(globalnr_1h);
WW = W(globalnr_0);
PP = P(globalnr_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hGLp,eGLp] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hGLp,eGLp,Meshp);
ww           = reconstruct(0,WW,hGLp);
pp           = reconstruct(2,PP,eGLp,Meshp); pp = pp-mean(mean(pp));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotten

if plot_figures

meshplot
plotten

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hodge decomposition

% dw = NG*W;
% dwu = dw(globalnr_1v);
% dwv = dw(globalnr_1h);
% [dwwu,dwwv,dwwvelo] = reconstruct_oneforms(dwu,dwv,hGLp,eGLp,Meshp);
% figure
% pcolor(Xp,Yp,reshape(dwwvelo,nn,nn))
% shading interp
% hold on
% nn = size(Xp,1);
% dwwu = reshape(dwwu,nn,nn);
% dwwv = reshape(dwwv,nn,nn);
% quiver(Xp(1:4:nn-1,1:4:nn-1),Yp(1:4:nn-1,1:4:nn-1),dwwu(1:4:nn-1,1:4:nn-1),dwwv(1:4:nn-1,1:4:nn-1),'w')
% title('Lid-Driven Cavity Stokes: Hodge Decomposition, dw')
% axis equal
% axis(XYlim)
% colorbar
% set(gca,'clim',[0 10])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stream function

PSI1 = zeros(N+1);
for i=1:N+1
    for j=1:N
        PSI1(i,j+1) = PSI1(i,j) - UU(i+(j-1)*(N+1));
    end
end

PSI2 = zeros(N+1);
for j=1:N+1
    for i=1:N
        PSI2(i+1,j) = PSI2(i,j) + VV(j+(i-1)*(N+1));
    end
end

PSI = (PSI1+PSI2)/2;

psi = reconstruct(0,reshape(PSI,[],1),hGLp);

ppsi = reshape(psi,nn,nn);

if plot_figures
    figure
    pcolor(Xp,Yp,ppsi)
    shading interp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%