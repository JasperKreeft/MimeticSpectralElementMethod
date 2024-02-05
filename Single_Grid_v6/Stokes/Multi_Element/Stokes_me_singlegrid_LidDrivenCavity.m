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

% s = matlabpool('size'); if s==0; eval(['matlabpool open ' num2str(feature('numCores'))]); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns nn
global xi
global h e
global globalnr_0 globalnr_1v globalnr_1h globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'LidDrivenCavity';
Domain       = 'SinDeformGrid';%'MuST';%'CosDeformGrid';%
DomInfo      = 0.0;

N = 3;
H = 5;

numRows    = H;
numColumns = H;
numElements = H*H;

plot_figures  = 1;
Tecplot = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square')

Mesh = meshgenerator_square(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Topology and Metric matrices

[h,e] = MimeticpolyVal(xi,N,1);

NGe = normalgrad(N);

De = div(N);

NG = zeros(nr_1,nr_0);
D  = zeros(nr_2,nr_1);
M0 = spalloc(nr_0,nr_0,nr_0);
M1 = zeros(nr_1);
M2 = zeros(nr_2);
f  = zeros(nr_1,1);
Wbc0 = zeros(nr_0);
Wbc1 = zeros(nr_1);

tic
for i=1:numElements

ind0 = globalnr_0(:,i);
ind1 = [ globalnr_1v(:,i) ; globalnr_1h(:,i) ];
ind2 = globalnr_2(:,i);

% Normal Gradient operator
NG(ind1,ind0) = NGe;

% Divergence operator
D(ind2,ind1) = De;
        
% zero-forms
M0e = innerproduct(0,Mesh.J(:,i));
M0(ind0,ind0) = M0(ind0,ind0) + M0e;

% one-forms
M1e = innerproduct(1,Mesh.J(:,i),Mesh.Qinv(:,3*(i-1)+(1:3)));
M1(ind1,ind1) = M1(ind1,ind1) + M1e;

% two-forms
M2e = innerproduct(2,Mesh.J(:,i));
M2(ind2,ind2) = M2e;

% boundary matrix zero form
Wbc0(ind0,ind0) = Wbc0(ind0,ind0) + boundaryIntegral([0 1]);

% boundary matrix one form
Wbc1(ind1,ind1) = Wbc1(ind1,ind1) + boundaryIntegral([1 2]);

end
toc

D    = sparse(D);
NG   = sparse(NG);
M0   = sparse(M0);
M1   = sparse(M1);
M2   = sparse(M2);
Wbc0 = sparse(Wbc0);
Wbc1 = sparse(Wbc1);

[Vorticity_bc,TangentialVelocity_bc,boundary_w,interior_w] = boundaryconditions_0_square(Mesh,FunctionType,[ 0 0 0 0 ]);
ind = globalnr_0(N*(N+1)+1,numColumns*(numRows-1)+1);
TangentialVelocity_bc(ind,1)  = TangentialVelocity_bc(ind,1)/2;
TangentialVelocity_bc(nr_0,1) = TangentialVelocity_bc(nr_0,1)/2;


[Pbc,UVbc,boundary_uv,interior_uv] = boundaryconditions_1_square(FunctionType,Domain,DomInfo,[ 0 0 0 0 ]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct system matrix and righthandside

Matrix = [ spalloc(nr_1,nr_1,0)      M1*NG               D'*M2
                NG'*M1                 M0          spalloc(nr_0,nr_2,0)
                 M2*D         spalloc(nr_2,nr_0,0) spalloc(nr_2,nr_2,0) ];

RHS = [ zeros(nr_1,1)
        zeros(nr_0,1)
        zeros(nr_2,1) ];

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

ind = sort([boundary_uv ; boundary_w+nr_1]);

Matrix(:,ind) = [];
Matrix(ind,:) = [];

RHS(ind,:) = [];

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

DivU = D*UV;
DivU = DivU(globalnr_2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hp,ep,Meshp);
ww           = reconstruct(0,WW,hp);
pp           = reconstruct(2,PP,ep,Meshp); pp = pp-mean(mean(mean(pp)));
divu         = reconstruct(2,DivU,ep,Meshp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Stream function

PSI1 = zeros(size(globalnr_0));
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
    
for i=1:N+1
    for j=1:N
        ij = i+(j-1)*(N+1);
        
        PSI1(ij+N+1,rc) = PSI1(ij,rc) + UU(i+(j-1)*(N+1),rc);
        
    end
end

if r<numRows
    PSI1(1:N+1,rc+numColumns) = PSI1(N*(N+1)+1:(N+1)^2,rc);
end

    end
end



PSI2 = zeros(size(globalnr_0));
for r=1:numRows
    for c=1:numColumns
        rc = c+(r-1)*numColumns;
    
for j=1:N+1
    for i=1:N
        ij = i+(j-1)*(N+1);
        
        PSI2(ij+1,rc) = PSI2(ij,rc) - VV(j+(i-1)*(N+1),rc);
        
    end
end

if c<numColumns
    PSI2(1:N+1:(N+1)^2,rc+1) = PSI2(N+1:N+1:(N+1)^2,rc);
end

    end
end

PSI = (PSI1+PSI2)/2;

psi  = reconstruct(0,PSI,hp);

psi = psi/2; % domain [-1,1] to [0,1]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotten

if plot_figures
    meshplot
    plotten_me

figure
for i=1:numElements
    
    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
    up = reshape(uu(:,i),nn,nn);
    vp = reshape(vv(:,i),nn,nn);
    Vp = reshape(velo(:,i),nn,nn);
    
pcolor(Xp,Yp,Vp)
hold on
quiver(Xp(1:4:nn-1,1:4:nn-1),Yp(1:4:nn-1,1:4:nn-1),up(1:4:nn-1,1:4:nn-1),vp(1:4:nn-1,1:4:nn-1),'w')
end
shading interp
axis equal
axis(XYlim)
axis off
set(gca,'clim',[0 1])
colorbar
title('Lid-Driven Cavity Stokes')

figure
for i=1:numElements
    
    Xp   = reshape(Meshp.X(:,i),nn,nn);
    Yp   = reshape(Meshp.Y(:,i),nn,nn);
    psip = reshape(psi(:,i),nn,nn);
    
contourf(Xp,Yp,psip,[-0.11 -0.09 -0.07 -0.05 -0.03 -0.01 -0.001 -0.0001 -0.00001 0 0.00001 0.0001 0.001 0.01 ])
hold on
end
% shading interp
axis equal
axis(XYlim)
axis off
% set(gca,'clim',[0 1])
colorbar
title('Lid-Driven Cavity Stokes: Stream function')

figure
for i=1:numElements

    Xp   = reshape(Meshp.X(:,i),nn,nn);
    Yp   = reshape(Meshp.Y(:,i),nn,nn);
    divup = reshape(divu(:,i),nn,nn);

    pcolor(Xp,Yp,divup)
    hold on
end
shading interp
axis equal
axis(XYlim)
axis off
% set(gca,'clim',[0 1])
colorbar
title('Lid-Driven Cavity Stokes: divergence velocity')


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matlab to Tecplot

if Tecplot
    
for i=1:numElements

name = ['liddrivensquare2_'  num2str(N)];

if i<10
    filename = strcat(name,'_0',num2str(i));
else
    filename = strcat(name,'_',num2str(i));
end

title = filename;
% variablenames = '"X" "Y" "U" "V" "Abs V" "P" "W" "PSI"';
variablenames = '"X" "Y" "U" "V" "Abs V" "P" "W" "DivU"';
meshsize = [nn nn];
% data = [ Meshp.X(:,i) Meshp.Y(:,i) uu(:,i) vv(:,i) velo(:,i) pp(:,i) ww(:,i) psi(:,i) ];
data = [ Meshp.X(:,i) Meshp.Y(:,i) uu(:,i) vv(:,i) velo(:,i) pp(:,i) ww(:,i) 1e15*divu(:,i)];

MatlabToTecplot('IJK',filename,title,variablenames,meshsize,data,2)

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hodge decomposition
% % keyboard
% dw = -NG*W;     % WAAROM HIER EEN MINTEKEN ????
% dwu = dw(globalnr_1v);
% dwv = dw(globalnr_1h);
% [dwwu,dwwv,dwwvelo] = reconstruct_oneforms(dwu,dwv,hp,ep,Meshp);
% nn = size(Xp,1);
% figure
% subplot(2,4,1)
% for rc = 1:numRows*numColumns
% pcolor(Xp(:,:,rc),Yp(:,:,rc),dwwvelo(:,:,rc))
% hold on
% % quiver(Xp(1:4:nn-1,1:4:nn-1,rc),Yp(1:4:nn-1,1:4:nn-1,rc),dwwu(1:4:nn-1,1:4:nn-1,rc),dwwv(1:4:nn-1,1:4:nn-1,rc),'w')
% end
% shading interp
% title('Hodge Decomposition, dw')
% axis equal
% axis(XYlim)
% % colorbar
% set(gca,'clim',[0 10])
% 
% subplot(2,4,5)
% for rc = 1:numRows*numColumns
% quiver(Xp(1:4:nn-1,1:4:nn-1,rc),Yp(1:4:nn-1,1:4:nn-1,rc),dwwu(1:4:nn-1,1:4:nn-1,rc),dwwv(1:4:nn-1,1:4:nn-1,rc),'k')
% hold on
% end
% title('Hodge Decomposition, dw')
% axis equal
% axis(XYlim)
% 
% 
% % FIRST WE NEED TO FIND Pbc, USING INTERPOLATION OF P TOWARDS THE BOUNDARY !!!
% dstarp = inv(M1)*(D'*M2*P-Wbc1*Pbc);
% dstarpu = dstarp(globalnr_1v);
% dstarpv = dstarp(globalnr_1h);
% [dstarppu,dstarppv,dstarppvelo] = reconstruct_oneforms(dstarpu,dstarpv,hp,ep,Meshp);
% nn = size(Xp,1);
% subplot(2,4,2)
% for rc = 1:numRows*numColumns
% pcolor(Xp(:,:,rc),Yp(:,:,rc),dstarppvelo(:,:,rc))
% hold on
% % quiver(Xp(1:4:nn-1,1:4:nn-1,rc),Yp(1:4:nn-1,1:4:nn-1,rc),dstarppu(1:4:nn-1,1:4:nn-1,rc),dstarppv(1:4:nn-1,1:4:nn-1,rc),'w')
% end
% shading interp
% title('Hodge Decomposition, dstarp')
% axis equal
% axis(XYlim)
% % colorbar
% set(gca,'clim',[0 10])
% 
% subplot(2,4,6)
% for rc = 1:numRows*numColumns
% quiver(Xp(1:4:nn-1,1:4:nn-1,rc),Yp(1:4:nn-1,1:4:nn-1,rc),dstarppu(1:4:nn-1,1:4:nn-1,rc),dstarppv(1:4:nn-1,1:4:nn-1,rc),'k')
% hold on
% end
% title('Hodge Decomposition, dstarp')
% axis equal
% axis(XYlim)
% 
% 
% 
% HDu = dw + dstarp;
% HDuu = HDu(globalnr_1v);
% HDuv = HDu(globalnr_1h);
% [HDuuu,HDuvv,HDuvelo] = reconstruct_oneforms(HDuu,HDuv,hp,ep,Meshp);
% nn = size(Xp,1);
% subplot(2,4,3)
% for rc = 1:numRows*numColumns
% pcolor(Xp(:,:,rc),Yp(:,:,rc),HDuvelo(:,:,rc))
% hold on
% % quiver(Xp(1:4:nn-1,1:4:nn-1,rc),Yp(1:4:nn-1,1:4:nn-1,rc),HDuuu(1:4:nn-1,1:4:nn-1,rc),HDuvv(1:4:nn-1,1:4:nn-1,rc),'w')
% end
% shading interp
% title('Hodge Decomposition, u_{HD}')
% axis equal
% axis(XYlim)
% % colorbar
% set(gca,'clim',[0 10])
% 
% subplot(2,4,7)
% for rc = 1:numRows*numColumns
% quiver(Xp(1:4:nn-1,1:4:nn-1,rc),Yp(1:4:nn-1,1:4:nn-1,rc),HDuuu(1:4:nn-1,1:4:nn-1,rc),HDuvv(1:4:nn-1,1:4:nn-1,rc),'k')
% hold on
% end
% title('Hodge Decomposition, u_{HD}')
% axis equal
% axis(XYlim)
% 
% subplot(2,4,4)
% for rc = 1:numRows*numColumns
% pcolor(Xp(:,:,rc),Yp(:,:,rc),velo(:,:,rc))
% hold on
% % quiver(Xp(1:4:nn-1,1:4:nn-1,rc),Yp(1:4:nn-1,1:4:nn-1,rc),HDuuu(1:4:nn-1,1:4:nn-1,rc),HDuvv(1:4:nn-1,1:4:nn-1,rc),'w')
% end
% shading interp
% title('Hodge Decomposition, u_{num}')
% axis equal
% axis(XYlim)
% % colorbar
% set(gca,'clim',[0 10])
% 
% subplot(2,4,8)
% for rc = 1:numRows*numColumns
% quiver(Xp(1:4:nn-1,1:4:nn-1,rc),Yp(1:4:nn-1,1:4:nn-1,rc),uu(1:4:nn-1,1:4:nn-1,rc),vv(1:4:nn-1,1:4:nn-1,rc),'k')
% hold on
% end
% title('Hodge Decomposition, u_{num}')
% axis equal
% axis(XYlim)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%