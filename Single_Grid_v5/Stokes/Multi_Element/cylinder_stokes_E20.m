%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Multi-element Stokes solver for Lid-Driven Cavity problem               %
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

global N numElements nn
global xi
global h e
global globalnr_0 globalnr_1v globalnr_1h globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

N = 6;

numElements = 20;
% 
plot_figures = 1;
Tecplot = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('cylinder_E20')

R = 0.5;
Mesh = meshgenerator_cylinder_E20(R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NGe = normalgrad(N);

De = div(N);

NG = zeros(nr_1,nr_0);
D  = zeros(nr_2,nr_1);
M0 = zeros(nr_0);        % spalloc(nr_0,nr_0,nr_0);
M1 = zeros(nr_1);
M2 = zeros(nr_2);
Wbc0 = zeros(nr_0);
Wbc1 = zeros(nr_1);

for i=1:numElements

ind0 = globalnr_0(:,i);
ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
ind2 = globalnr_2(:,i);

% Normal Gradient operator
NG(ind1,ind0) = NGe;

% Divergence operator
if i==7 % Orientation correction element 7
    De7 = De;
    for j=1:N
        De7(1+(j-1)*N,1+(j-1)*(N+1)) = 1;
    end
    D(ind2,ind1) = De7;
elseif i==15 % Orientation correction element 15
    De15 = De;
    for j=1:N
        De15(j+N*(N-1),j*(N+1)+N*(N+1)) = -1;
    end
    D(ind2,ind1) = De15;
else
    D(ind2,ind1) = De;
end
        
% zero-forms
M0e = innerproduct(0,Mesh.J(:,i));
M0(ind0,ind0) = M0(ind0,ind0) + M0e;

% one-forms
M1e = innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)));
if i==7 % Orientation correction element 7
    ind = 1:N+1:N*(N+1);
    M1e(:,ind) = -M1e(:,ind);
elseif i==15 % Orientation correction element 15
    ind = N*(N+1)+(N+1:N+1:N*(N+1));
    M1e(:,ind) = -M1e(:,ind);
end
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

% two-forms
M2e = innerproduct(2,Mesh.J(:,i));
M2(ind2,ind2) = M2e;

% boundary matrix zero form
Wbc0(ind0,ind0) = Wbc0(ind0,ind0) + boundaryIntegral([0 1]);

% boundary matrix one form
Wbc1(ind1,ind1) = Wbc1(ind1,ind1) + boundaryIntegral([1 2]);

end

D    = sparse(D);
NG   = sparse(NG);
M0   = sparse(M0);
M1   = sparse(M1);
M2   = sparse(M2);
Wbc0 = sparse(Wbc0);
Wbc1 = sparse(Wbc1);

[Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_stokes_cylinder_E20(Mesh);

[Pbc,UVbc,boundary_uv,interior_uv] = boundaryconditions_1_stokes_cylinder_E20(Mesh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = [       M0           NG'*M1        spalloc(nr_0,nr_2,0)
               M1*NG   spalloc(nr_1,nr_1,0)      D'*M2
         spalloc(nr_2,nr_0,0)  M2*D         spalloc(nr_2,nr_2,0) ];

RHS = zeros(nr_0+nr_1+nr_2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove Boundary Conditions

% 0
RHS(1:nr_0) = RHS(1:nr_0) - Wbc0*TangentialVelocity_bc;
if ~isempty(boundary_w)
RHS = RHS - Matrix(:,boundary_w)*Vorticity_bc;
end

% 1
RHS(nr_0+1:nr_0+nr_1) = RHS(nr_0+1:nr_0+nr_1) + Wbc1*Pbc;
if ~isempty(boundary_uv)
RHS = RHS - Matrix(:,boundary_uv)*UVbc;
end

ind = sort([boundary_w ; boundary_uv+nr_0]);

Matrix(:,ind) = [];
Matrix(ind,:) = [];

RHS(ind,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

WUVP = Matrix\RHS;

ind1 = nr_0-length(boundary_w);
W_in = WUVP(1:ind1);
ind2 = ind1+nr_1-length(boundary_uv);
UV_in = WUVP(ind1+1:ind2);
ind3 = ind2+nr_2;
P_in = WUVP(ind2+1:ind3);

UV = zeros(nr_1,1);
W  = zeros(nr_0,1);
P  = zeros(nr_2,1);

UV(interior_uv) = UV_in;
W(interior_w)   = W_in;
P = P_in;
UV(boundary_uv) = UVbc;
W(boundary_w)   = Vorticity_bc;

UU = UV(globalnr_1v);
VV = UV(globalnr_1h);
WW = W(globalnr_0);
PP = P(globalnr_2);

% Orientation correction element 7 & 15
UU(1:N+1:N*(N+1),7)    = -UU(1:N+1:N*(N+1),7);
VV(N+1:N+1:N*(N+1),15) = -VV(N+1:N+1:N*(N+1),15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_cylinder_E20(R);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hp,ep,Meshp);
ww           = reconstruct(0,WW,hp);
pp           = reconstruct(2,PP,ep,Meshp); pp = pp-mean(mean(mean(pp)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotten

if plot_figures
    meshplot
    plotten_me
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matlab to Tecplot

if Tecplot
    
for i=1:numElements

name = ['liddrivencylinder_'  num2str(N)];

if i<10
    filename = strcat(name,'_0',num2str(i));
else
    filename = strcat(name,'_',num2str(i));
end

title = filename;
variablenames = '"X" "Y" "U" "V" "Abs V" "P" "W"';
meshsize = [nn nn];
data = [ Meshp.X(:,i) Meshp.Y(:,i) uu(:,i) vv(:,i) velo(:,i) pp(:,i) ww(:,i) ];

MatlabToTecplot(filename,title,variablenames,meshsize,data,2)

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

% in = 'finish';
% GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%