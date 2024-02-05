%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Multi-element vector Poisson problem                                    %
%                                                                         %
% written by Jasper Kreeft (2011)                                         %
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

global N numRows numColumns
global xi
global h e
global globalnr_0 globalnr_1v globalnr_1h globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

N = 4;

numRows    = 1;
numColumns = 4;
RC = numRows*numColumns;

DomInfo = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering_in('multi');

Mesh = meshgenerator('Annulus',DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

Ge = grad_in(N);
Ce = curl_in(N);

G  = zeros(nr_1,nr_0);
C  = zeros(nr_2,nr_1);
M0 = zeros(nr_0);
M1 = zeros(nr_1);
M2 = zeros(nr_2);
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;

% ind0 = reshape(globalnr_0((c-1)*N+(1:N+1),(r-1)*N+(1:N+1)),1,[]);
ind0 = globalnr_0((c-1)*N+(1:N+1),(r-1)*N+(1:N+1));
% ind1 = [ reshape(globalnr_1h((c-1)*N+(1:N),(r-1)*N+(1:N+1)),1,N*(N+1)) ...
%          reshape(globalnr_1v((c-1)*N+(1:N+1),(r-1)*N+(1:N))',1,N*(N+1))];
ind1 = [ globalnr_1h((c-1)*N+(1:N),(r-1)*N+(1:N+1)) ...
         globalnr_1v((c-1)*N+(1:N+1),(r-1)*N+(1:N))'];
ind2 = globalnr_2((c-1)*N+(1:N),(r-1)*N+(1:N));

% Normal Gradient
G(ind1,ind0) = Ge;

% Divergence operator
C(ind2,ind1) = Ce;

% zero-forms
M0e = innerproduct(0,Mesh.J(:,rc));
M0(ind0,ind0) = M0(ind0,ind0) + M0e;

% one-forms
M1e = innerproduct(1,Mesh.J(:,rc),Mesh.Qinv(:,3*(rc-1)+(1:3)));
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

% two-forms
M2e = innerproduct(2,Mesh.J(:,rc));
M2(ind2,ind2) = M2e;

    end
end

G = sparse(G);
C = sparse(C);
M0 = sparse(M0);
M1 = sparse(M1);
M2 = sparse(M2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and solve system

F = force_annulus(Mesh.X,Mesh.Y);
[f1 f2] = force_annulus(Mesh.X,Mesh.Y);

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid('Annulus',DomInfo);

% Reconstruction
[ffx,ffy,ff] = reconstruct_oneforms_in(f1,f2,hp,ep,Meshp);

Xp = Meshp.X;
Yp = Meshp.Y;
XYlim = [-1 1 -1 1];

figure
for rc=1:numRows*numColumns
subplot(2,2,1)
pcolor(Xp(:,:,rc),Yp(:,:,rc),ffx(:,:,rc)/(4*pi))
hold on
shading interp
colorbar
title('fx')
axis equal
axis(XYlim)
subplot(2,2,2)
pcolor(Xp(:,:,rc),Yp(:,:,rc),ffy(:,:,rc)/(4*pi))
hold on
shading interp
colorbar
title('fy')
axis equal
axis(XYlim)
subplot(2,2,3)
pcolor(Xp(:,:,rc),Yp(:,:,rc),ff(:,:,rc)/(4*pi))    %/(4*pi)
hold on
shading interp
colorbar
title('ff')
axis equal
axis(XYlim)
subplot(2,2,4)
quiver(Xp(:,:,rc),Yp(:,:,rc),ffx(:,:,rc),ffy(:,:,rc))
hold on
axis equal
axis(XYlim)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

Matrix = [ C'*M2*C M1*G
           G'*M1    M0  ];

RHS = [      M1*F
        zeros(nr_0,1) ];


Matrix(:,globalnr_1v(1,:)) = Matrix(:,globalnr_1v(1,:)) + Matrix(:,globalnr_1v(end,:));
% Matrix(:,globalnr_0(1,:)+nr_1) = Matrix(:,globalnr_0(1,:)+nr_1) + Matrix(:,globalnr_0(end,:)+nr_1);

Matrix(:,globalnr_0(end,:)+nr_1) = [];
Matrix(:,globalnr_1v(end,:))     = [];
Matrix(globalnr_0(end,:)+nr_1,:) = [];
Matrix(globalnr_1v(end,:),:)     = [];
RHS(globalnr_0(end,:)+nr_1)      = [];
RHS(globalnr_1v(end,:))          = [];

us = Matrix\RHS;

u = us(1:nr_1-N);



UU = u(globalnr_1h);
VV = zeros(size(globalnr_1v));
VV(1:4*N,:) = u(globalnr_1v(1:4*N,:));
VV(4*N+1,:) = u(globalnr_1v(1,:));

[uu,vv,velo] = reconstruct_oneforms_in2(UU,VV,hp,ep,Meshp);

figure
for rc=1:numRows*numColumns
subplot(2,2,1)
pcolor(Xp(:,:,rc),Yp(:,:,rc),uu(:,:,rc))
hold on
shading interp
colorbar
title('uu')
axis equal
axis(XYlim)
subplot(2,2,2)
pcolor(Xp(:,:,rc),Yp(:,:,rc),vv(:,:,rc))
hold on
shading interp
colorbar
title('vv')
axis equal
axis(XYlim)
subplot(2,2,3)
pcolor(Xp(:,:,rc),Yp(:,:,rc),velo(:,:,rc))
hold on
shading interp
colorbar
title('velo')
axis equal
axis(XYlim)
subplot(2,2,4)
quiver(Xp(:,:,rc),Yp(:,:,rc),uu(:,:,rc),vv(:,:,rc))
hold on
axis equal
axis(XYlim)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

% in = 'finish';
% GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%