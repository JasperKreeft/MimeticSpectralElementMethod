clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_TwoForms/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements
global xi
global h e
global nr_1 nr_2
global globalnr_1v globalnr_1h globalnr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'cylinder_E16';
Domain       = 'naca0012_E16';
DomInfo      = 0.0;

N = 10;

numElements = 16;

plot_figures  = 1;
Tecplot = 0;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Potential flow for two-forms ')
disp(['Testfunction is ' FunctionType])
disp(['Computational Domain is ' Domain ' with settings ' num2str(DomInfo)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('cylinder_E16')

Mesh = meshgenerator_naca0012_E16(DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

M1 = innerproduct_assembly(1,Mesh);

M2 = innerproduct_assembly(2,Mesh);

D = divergence_assembly();

Wbce = boundaryIntegral([1 2]);

Wbc = zeros(nr_1);

for i=1:numElements

ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];

% boundary matrix
Wbc(ind1,ind1) = Wbc(ind1,ind1) + Wbce;

end

Wbc = sparse(Wbc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Forcing function and boundary conditions

[PHIbc,Qbc,boundary_flux,interior_flux] = boundaryconditions_1_cylinder_E16(Mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = [   M1        D'*M2
           M2*D spalloc(nr_2,nr_2,0) ];

RHS = zeros(nr_1+nr_2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Remove boundary Conditions

% RHS(1:nr_1) = RHS(1:nr_1) + Wbc*PHIbc;
if ~isempty(boundary_flux)
RHS = RHS - Matrix(:,boundary_flux)*Qbc;

Matrix(:,boundary_flux) = [];
Matrix(boundary_flux,:) = [];

RHS(boundary_flux) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Solve system and add boundary conditions to data

QPHI = Matrix\RHS;

ind = nr_1-length(boundary_flux);
Q_in = QPHI(1:ind);
PHI  = QPHI(ind+1:end);

Q = zeros(nr_1,1);
Q(interior_flux) = Q_in;
Q(boundary_flux) = Qbc;

PHI = PHI(globalnr_2);
Qxi  = Q(globalnr_1v);
Qeta = Q(globalnr_1h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_naca0012_E16(DomInfo);

% Reconstruction
phi          = reconstruct(2,PHI,ep,Meshp);
[qx,qy,qMag] = reconstruct(1,Qxi,Qeta,hp,ep,Meshp);

% Plotten
if plot_figures
    meshplot
    plotten
end

% figure
% subplot(3,1,1)
% for i=1:numElements
%     Xp = reshape(Meshp.X(:,i),nn,nn);
%     Yp = reshape(Meshp.Y(:,i),nn,nn);
%     Qp = reshape(qx(:,i),nn,nn);
%     pcolor(Xp,Yp,Qp)
%     axis equal
%     axis(XYaxis)
%     hold on
% end
% shading interp
% colorbar
% 
% subplot(3,1,2)
% for i=1:numElements
%     Xp = reshape(Meshp.X(:,i),nn,nn);
%     Yp = reshape(Meshp.Y(:,i),nn,nn);
%     Qp = reshape(qy(:,i),nn,nn);
%     pcolor(Xp,Yp,Qp)
%     axis equal
%     axis(XYaxis)
%     hold on
% end
% shading interp
% colorbar
% 
% subplot(3,1,3)
% for i=1:numElements
%     Xp = reshape(Meshp.X(:,i),nn,nn);
%     Yp = reshape(Meshp.Y(:,i),nn,nn);
%     Qp = reshape(qMag(:,i),nn,nn);
%     pcolor(Xp,Yp,Qp)
%     axis equal
%     axis(XYaxis)
%     hold on
% end
% shading interp
% colorbar

derivplot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matlab to Tecplot

if Tecplot
    
for i=1:numElements

name = ['naca0012_'  num2str(N)];

if i<10
    filename = strcat(name,'_0',num2str(i));
else
    filename = strcat(name,'_',num2str(i));
end

title = filename;
variablenames = '"X" "Y" "phi" "U" "V" "Abs V"';
meshsize = [nn nn];
data = [ Meshp.X(:,i) Meshp.Y(:,i) phi(:,i) qx(:,i) qy(:,i) qMag(:,i) ];

MatlabToTecplot(filename,title,variablenames,meshsize,data,2)

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%