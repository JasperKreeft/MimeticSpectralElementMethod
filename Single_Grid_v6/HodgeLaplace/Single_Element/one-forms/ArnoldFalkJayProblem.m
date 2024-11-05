%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Single-element Arnold-Falk-Jay problem                                  %
%                                                                         %
% written by Jasper Kreeft (2012)                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
% run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N w
global h e
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

BC = 2;

Domain = 'SinDeformGrid_01';
if BC==1
    FunctionType = 'arnoldfalkjay1';
elseif BC==2
    FunctionType = 'arnoldfalkjay2';
end
DomInfo = 0.0;


NrCellRange = 4:2:12;

plot_figures  = 1;
error_figures = 1;
casenr        = FunctionType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

for N=NrCellRange

disp(['N = ' num2str(N)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

[xi,w] = GLLnodes(N);     eta = xi;                % Gauss-Lobotto-Legendre

Mesh = meshgenerator_square(Domain,DomInfo);

numbering('square')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NG = normalgrad(N);

D = div(N);

M0 = innerproduct(0,Mesh.J);

M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

M2 = innerproduct(2,Mesh.J);

[fx fy] = forcefunction(1,1,1,FunctionType,Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and solve system

Matrix = -[   M0     NG'*M1
             M1*NG -D'*M2*D ];

RHS = [ zeros(nr_0,1)
        M1*[ fx ; fy ] ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions
if BC==1
    boundary = [];
elseif BC==2
    ind = [ 1:N+1:N*(N+1) N+1:N+1:N*(N+1) ];
    boundary = sort([ ind N*(N+1)+ind ]);
end
interior = 1:nr_1; interior(boundary) = [];

Matrix(:,nr_0+boundary) = [];
Matrix(nr_0+boundary,:) = [];

RHS(nr_0+boundary) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

SigmaUV = Matrix\RHS;

Sigma = SigmaUV(1:nr_0);
UV_in = SigmaUV(nr_0+1:end);

UV    = zeros(nr_1,1);
UV(interior) = UV_in;
UV(boundary) = zeros(length(boundary),1);

SSigma = Sigma(globalnr_0);
UU = UV(globalnr_1v);
VV = UV(globalnr_1h);

DivUV = D*UV;
NG_sigma = NG*Sigma;
NG_sigma_x = NG_sigma(globalnr_1h);
NG_sigma_y = NG_sigma(globalnr_1v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hGLp,eGLp] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
ssigma        = reconstruct(0,SSigma,hGLp);
[uu,vv,velo]  = reconstruct(1,UU,VV,hGLp,eGLp,Meshp);
[ffx,ffy]     = reconstruct(1,fx,fy,hGLp,eGLp,Meshp);
[ng_sx,ng_sy] = reconstruct(1,NG_sigma_x,NG_sigma_y,hGLp,eGLp,Meshp);
div_uv        = reconstruct(2,DivUV,eGLp,Meshp);


% Interpolated solution

% Exact Solution
sigma_ex          = exact_solution(Meshp.X,Meshp.Y,FunctionType,'zero');
[u_ex,v_ex]       = exact_solution(Meshp.X,Meshp.Y,FunctionType,'one');
[fx_ex,fy_ex]     = exact_solution(Meshp.X,Meshp.Y,FunctionType,'force');
[ng_x_ex,ng_y_ex] = exact_solution(Meshp.X,Meshp.Y,FunctionType,'d_zero');
div_ex            = exact_solution(Meshp.X,Meshp.Y,FunctionType,'d_one');

if plot_figures
    plotten_AFG
end

if error_figures
    fout
end

end % for N

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Convergence plots

if error_figures
    foutplot
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
% GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%