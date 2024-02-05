%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Multi-element vector Poisson problem                                    %
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

global N numRows numColumns numElements nn
global xi
global h e
global globalnr_0 globalnr_1v globalnr_1h globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

N = 8;

numRows    = 1;
numColumns = 4;
numElements = numRows*numColumns;

FunctionType = 'AFW_annulus';
Domain = 'AFW_annulus';
DomInfo = 0.5;

plot_figures = 1;
Tecplot = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering_annulus();

Mesh = meshgenerator_annulus(DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NG = normalgradient_assembly();

D = divergence_assembly();

M0 = innerproduct_assembly(0,Mesh);

M1 = innerproduct_assembly(1,Mesh);

M2 = innerproduct_assembly(2,Mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and solve system

% F = force_annulus([1 DomInfo],[-pi/2 0 pi/2 pi 3/2*pi]);
% H = harmonicform_annulus([1 DomInfo],[-pi/2 0 pi/2 pi 3/2*pi]);

F = force_annulus([1 DomInfo],[-pi/4 pi/4 3*pi/4 5*pi/4 7/4*pi]);
H = harmonicform_annulus([1 DomInfo],[-pi/4 pi/4 3*pi/4 5*pi/4 7/4*pi]);

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_annulus(DomInfo);

fv = F(globalnr_1v);
fh = F(globalnr_1h);
% % Reconstruction
[ffx,ffy,ff] = reconstruct(1,fv,fh,hp,ep,Meshp);

h_v = H(globalnr_1v);
h_h = H(globalnr_1h);
% % Reconstruction
[hhx,hhy,hh] = reconstruct(1,h_v,h_h,hp,ep,Meshp);

Xp = Meshp.X;
Yp = Meshp.Y;
XYlim = [-1 1 -1 1];

if isempty(nn)
    nn = sqrt(size(Meshp.X,1));
end

if plot_figures
    figure
    for i=1:numElements

        Xp = reshape(Meshp.X(:,i),nn,nn);
        Yp = reshape(Meshp.Y(:,i),nn,nn);

    subplot(2,2,1)
    pcolor(Xp,Yp,reshape(ffx(:,i),nn,nn))
    hold on
    shading interp
    colorbar
    title('fx')
    axis equal
    axis(XYlim)
    subplot(2,2,2)
    pcolor(Xp,Yp,reshape(ffy(:,i),nn,nn))
    hold on
    shading interp
    colorbar
    title('fy')
    axis equal
    axis(XYlim)
    subplot(2,2,3)
    pcolor(Xp,Yp,reshape(ff(:,i),nn,nn))
    hold on
    shading interp
    colorbar
    title('ff')
    axis equal
    axis(XYlim)
    subplot(2,2,4)
    quiver(Xp,Yp,reshape(ffx(:,i),nn,nn),reshape(ffy(:,i),nn,nn))
    hold on
    axis equal
    axis(XYlim)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

disp('Assemble and solve system')

Matrix = [   M0              NG'*M1   spalloc(nr_0,1,0)
            M1*NG           -D'*M2*D   M1*H
            spalloc(1,nr_0,0) H'*M1     0               ];

RHS = [ zeros(nr_0,1)
        -M1*F 
          0    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Boundary conditions

ind_bc_s = unique(globalnr_0([1:N+1 N*(N+1)+1:(N+1)^2],:));

ind_bc_u = unique(globalnr_1h([1:N+1:N*(N+1) N+1:N+1:N*(N+1)],:));

ind_bc = [ ind_bc_s ; ind_bc_u+nr_0 ];

ind_in_s = 1:nr_0; ind_in_s(ind_bc_s) = [];
ind_in_u = 1:nr_1; ind_in_u(ind_bc_u) = [];

Matrix(ind_bc,:) = [];
Matrix(:,ind_bc) = [];
RHS(ind_bc) = [];

su = Matrix\RHS;

S = zeros(nr_0,1);
u = zeros(nr_1,1);

S(ind_in_s) = su(1:nr_0-length(ind_bc_s));

ind = nr_0-length(ind_bc_s)+(1:nr_1-length(ind_bc_u));
u(ind_in_u) = su(ind);

C = su(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessing

disp('Postprocessing')

UU = u(globalnr_1v);
VV = u(globalnr_1h);

[uu,vv,velo] = reconstruct(1,UU,VV,hp,ep,Meshp);

if plot_figures
    figure
    for i=1:numElements

        Xp = reshape(Meshp.X(:,i),nn,nn);
        Yp = reshape(Meshp.Y(:,i),nn,nn);

    subplot(2,2,1)
    pcolor(Xp,Yp,reshape(uu(:,i),nn,nn))
    hold on
    shading interp
    colorbar
    title('uu')
    axis equal
    axis(XYlim)
    subplot(2,2,2)
    pcolor(Xp,Yp,reshape(vv(:,i),nn,nn))
    hold on
    shading interp
    colorbar
    title('vv')
    axis equal
    axis(XYlim)
    subplot(2,2,3)
    pcolor(Xp,Yp,reshape(velo(:,i),nn,nn))
    hold on
    shading interp
    colorbar
    title('velo')
    axis equal
    axis(XYlim)
    subplot(2,2,4)
    quiver(Xp,Yp,reshape(uu(:,i),nn,nn),reshape(vv(:,i),nn,nn))
    hold on
    axis equal
    axis(XYlim)
    end
end

hodgedecomposition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matlab to Tecplot

if Tecplot

disp('writing to TecPlot')

for i=1:numElements

name = ['AFW_annulus_'  num2str(N)];

if i<10
    filename = strcat(name,'_0',num2str(i));
else
    filename = strcat(name,'_',num2str(i));
end

title = filename;
variablenames = '"X" "Y" "U" "V" "Abs V" "Uz" "Vz" "Uperp" "Vperp" "Hx" "Hy"';
meshsize = [nn nn];
data = [ Meshp.X(:,i) Meshp.Y(:,i) uu(:,i) vv(:,i) velo(:,i) uuz(:,i) vvz(:,i) uuperp(:,i) vvperp(:,i) hhx(:,i) hhy(:,i) ];

MatlabToTecplot('IJK',filename,title,variablenames,meshsize,data,2)

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%