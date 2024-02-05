clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N nr_tau xi w e
global nr_1 nr_2
global globalnr_1h globalnr_1v globalnr_2
global nn

% FunctionType = 'AFG';
Domain       = 'SinDeformGrid';
DomInfo      = 0.0;

plot_figures  = 1;

N = 16;
Re = 50;

nr_tau = 2*( N*(N+2) + (N+1)^2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

[xi,w] = GLLnodes(N);

numbering('square')

Mesh = meshgenerator_square(Domain,DomInfo);

% Grid, basis-functions and weights for post-processen
[Meshp,hGLp,eGLp] = postproces_grid_square(Domain,DomInfo);
Xp = reshape(Meshp.X,nn,nn);
Yp = reshape(Meshp.Y,nn,nn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices
disp('vector divergence')
D  = div(N);

disp('tensor divergence')
CD = covect_div(N);

[~,e] = MimeticpolyVal(xi,N,1);

disp('tensor product')
T = tensorinnerproduct();

disp('trace')
Tr = Trace();

disp('inner product one-forms')
M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

disp('inner product two-forms')
M2 = innerproduct(2,Mesh.J);

disp('boundary condition')
[Wbc,UVbc] = toplid();

Conv_Weights = Convection_matrices();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Iterations
if sum(N==[12 13 14 16])
    disp('load Stokes solution')
    load(['Velocity_Stokes_E1_N' num2str(N) '.mat'])
else
    UV = zeros(2*N*(N+1),1);
end

error = 1; z = 0;
while error>1e-4
z = z+1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside
disp('convection')
Conv = Convection(UV,Conv_Weights);

UV_old = UV;

disp('Matrix assembly')
% Matrix = [          T             1/Re*CD'*Tr'*M1  + Conv     spalloc(nr_tau,nr_2,0)
%                 1/Re*M1*Tr*CD           spalloc(nr_1,nr_1,0)          D'*M2
%            spalloc(nr_2,nr_tau,0)           M2*D            spalloc(nr_2,nr_2,0)    ];

Matrix = [        Re*T             CD'*Tr'*M1+Re*Conv   spalloc(nr_tau,nr_2,0)
                M1*Tr*CD            eps*speye(nr_1)          D'*M2
           spalloc(nr_2,nr_tau,0)           M2*D        spalloc(nr_2,nr_2,0)    ]; %eps*speye(nr_2)


%  [fx fy] = forcefunction(1,1,1,FunctionType,Domain,DomInfo);

RHS = [    Wbc*UVbc
         zeros(nr_1,1)
         zeros(nr_2,1) ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

ind = 1:N+1:N*(N+1);
boundary = nr_tau + sort([ ind ind+N ind+N*(N+1) ind+N*(N+2) ]);
interior = 1:nr_1; interior(boundary-nr_tau) = [];

Matrix(:,boundary) = [];
Matrix(boundary,:) = [];

RHS(boundary,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data
disp('Solve matrix system')
tau_u_p = Matrix\RHS;

UV_in = tau_u_p(nr_tau+(1:nr_1-length(boundary)));
UV = zeros(nr_1,1);
UV(interior) = UV_in;

error = max(max(abs(diff(UV-UV_old))))

UU = UV(globalnr_1v);
VV = UV(globalnr_1h);

[uu,vv,velo] = reconstruct(1,UU,VV,hGLp,eGLp,Meshp);
Vp = reshape(velo,nn,nn);
clf
contourf(Xp,Yp,Vp,20)
colorbar
title('Velocity magnitude')
axis equal
pause(2)

end

Tau = tau_u_p(1:nr_tau);
UV_in = tau_u_p(nr_tau+(1:nr_1-length(boundary)));
UV = zeros(nr_1,1);
UV(interior) = UV_in;
P = tau_u_p(nr_tau+nr_1-length(boundary)+(1:nr_2));

UU = UV(globalnr_1v);
VV = UV(globalnr_1h);
PP = P(globalnr_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Reconstruction
pp           = reconstruct(2,PP,eGLp,Meshp); pp = pp-mean(mean(pp));
%  [ffx,ffy]    = reconstruct(1,fx,fy,hGLp,eGLp,Meshp);



if plot_figures
    
Xp = reshape(Meshp.X,nn,nn);
Yp = reshape(Meshp.Y,nn,nn);
up = reshape(uu,nn,nn);
vp = reshape(vv,nn,nn);
Vp = reshape(velo,nn,nn);
Pp = reshape(pp,nn,nn);
%  fxp = reshape(ffx,nn,nn);
%  fyp = reshape(ffy,nn,nn);

XYlim = [min(min(Xp)) max(max(Xp)) min(min(Yp)) max(max(Yp))];

figure
subplot(2,3,1)
contourf(Xp,Yp,up,20)
colorbar
title('u')
axis equal
axis(XYlim)
subplot(2,3,2)
contourf(Xp,Yp,vp,20)
colorbar
title('v')
axis equal
axis(XYlim)
subplot(2,3,3)

subplot(2,3,4)
contourf(Xp,Yp,Pp,20)
colorbar
title('p')
axis equal
axis(XYlim)
%  subplot(2,3,5)
%  contourf(Xp,Yp,fxp,20)
%  colorbar
%  title('f_x')
%  axis equal
%  axis(XYlim)
%  subplot(2,3,6)
%  contourf(Xp,Yp,fyp,20)
%  colorbar
%  title('f_y')
%  axis equal
%  axis(XYlim)

end

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

figure
contours = -[1.5e-3 1e-3 5e-4 2.5e-4 1e-4 5e-5 1e-5 1e-6 0 1e-10 -1e-5 -1e-4 -1e-2 -3e-2 -5e-2 -7e-2 -9e-2 -0.1 -0.11 -0.115 -0.1175];
% contour(Xp,Yp,ppsi,contours)
contour(Xp,Yp,ppsi,20)
axis([-1 1 -1 1])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%