clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N xi w e
global nr_1 nr_2
global globalnr_1h globalnr_1v globalnr_2
global nn

FunctionType = 'AFG';
Domain       = 'SinDeformGrid_01';
DomInfo      = 0.0;

plot_figures  = 1;

N = 12;

nr_tau = 2*( N*(N+2) + (N+1)^2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

[xi,w] = GLLnodes(N);

Mesh = meshgenerator_square(Domain,DomInfo);

numbering('square')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

D  = div(N);

CD = covect_div(N);

[~,e] = MimeticpolyVal(xi,N,1);

T = tensorinnerproduct();

Tr = Trace();

In = Inclusion();

M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

M2 = innerproduct(2,Mesh.J);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

% Matrix = [          T                   CD'*Tr'*M1          spalloc(nr_tau,nr_2,0)
%                 M1*Tr*CD            spalloc(nr_1,nr_1,0)          D'*M2
%            spalloc(nr_2,nr_tau,0)           M2*D            spalloc(nr_2,nr_2,0)    ];

Matrix = [          T                   CD'*Tr'*M1          spalloc(nr_tau,nr_2,0)
                M1*Tr*CD            eps*speye(nr_1)          D'*M2
           spalloc(nr_2,nr_tau,0)           M2*D            eps*speye(nr_2)    ];


[fx fy] = forcefunction(1,1,1,FunctionType,Domain,DomInfo);

RHS = [ zeros(nr_tau,1)
        -M1*[ fx ; fy ]
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

tau_u_p = Matrix\RHS;

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

% Grid, basis-functions and weights for post-processen
[Meshp,hGLp,eGLp] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hGLp,eGLp,Meshp);
pp           = reconstruct(2,PP,eGLp,Meshp); pp = pp-mean(mean(pp));
[ffx,ffy]    = reconstruct(1,fx,fy,hGLp,eGLp,Meshp);



if plot_figures
    
Xp = reshape(Meshp.X,nn,nn);
Yp = reshape(Meshp.Y,nn,nn);
up = reshape(uu,nn,nn);
vp = reshape(vv,nn,nn);
Vp = reshape(velo,nn,nn);
Pp = reshape(pp,nn,nn);
fxp = reshape(ffx,nn,nn);
fyp = reshape(ffy,nn,nn);

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
subplot(2,3,5)
contourf(Xp,Yp,fxp,20)
colorbar
title('f_x')
axis equal
axis(XYlim)
subplot(2,3,6)
contourf(Xp,Yp,fyp,20)
colorbar
title('f_y')
axis equal
axis(XYlim)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%