%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Multi-element vector Poisson problem                                    %
%                                                                         %
% written by Jasper Kreeft (2012)                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resetpath
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
global xi w
global h e
global globalnr_0 globalnr_1v globalnr_1h globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

N = 16;

numRows    = 1;
numColumns = 1;
numElements = numRows*numColumns;

FunctionType = 'HD';
Domain = 'SinDeformGrid';
DomInfo = 0.0;

plot_figures = 1;
Tecplot = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square');

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NG = normalgradient_assembly();

D = divergence_assembly();

M0 = innerproduct_assembly(0,Mesh);

M1 = innerproduct_assembly(1,Mesh);

M2 = innerproduct_assembly(2,Mesh);

boundaryintegral_HD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and solve system

f  = zeros(nr_1,1); Fx = 0*globalnr_1v; Fy = 0*globalnr_1h;
for i=1:numElements
    r = ceil(i/numColumns); c = i-(r-1)*numColumns;
    [fx fy] = forcefunction(1,r,c,FunctionType,Domain,DomInfo);
    Fx(:,i) = fx; Fy(:,i) = fy;
end
f(globalnr_1v) = Fx; f(globalnr_1h) = Fy;

% F = force_annulus([1 DomInfo],[-pi/4 pi/4 3*pi/4 5*pi/4 7/4*pi]);
% H = harmonicform_annulus([1 DomInfo],[-pi/4 pi/4 3*pi/4 5*pi/4 7/4*pi]);

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

[ffx,ffy,ff] = reconstruct(1,Fx,Fy,hp,ep,Meshp);

% h_v = H(globalnr_1v);
% h_h = H(globalnr_1h);
% % % Reconstruction
% [hhx,hhy,hh] = reconstruct(1,h_v,h_h,hp,ep,Meshp);

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
    quiver(Xp,Yp,reshape(ffx(:,i),nn,nn),reshape(-ffy(:,i),nn,nn))
    hold on
    axis equal
    axis(XYlim)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessing

disp('Postprocessing')

u = f;

% hodgedecomposition
hodgedecomposition_2

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