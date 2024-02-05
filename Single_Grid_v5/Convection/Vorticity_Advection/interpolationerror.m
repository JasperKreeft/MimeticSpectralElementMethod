% Zero-forms, Outer-oriented, constant vector-field

clear all
% close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library_Convection/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N nn numElements numRows numColumns
global globalnr_0 globalnr_1h globalnr_1v globalnr_2
global nr_0 nr_1 nr_2
global w h e
global Mesh
global Re

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'SinDeformGrid_01';
DomInfo = 0.2;
FunctionType = 'reguralizedLDC';

plot_fig  = 0;
avimatlab = 0;
Tecplot   = 0;

NrCellRange = 4;
HconvRange = [ 2 4 8 16 ];
T = 2000;
dt = 1
Re = 100;

filename = ['Zero_Vort_ME_N' num2str(N) '_dt' num2str(dt)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start H-convergence loop

conv_i = 0;
L2_u = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_v = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_w = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_p = zeros(1,max(length(NrCellRange),length(HconvRange)));

for Hconv = HconvRange

numRows    = Hconv;
numColumns = Hconv;
numElements = numRows*numColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start P-convergence loop

for N=NrCellRange

disp('  ')
disp(['N = ' num2str(N) ', H = ' num2str(Hconv)])
disp('  ')

conv_i = conv_i+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

numbering('square');

Mesh = meshgenerator_square(Domain,DomInfo);

[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

[psi_ex,u_ex,v_ex,w_ex,p_ex,Fy_ex] = ReguralizedLDC(Re,Meshp.X,Meshp.Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize figure and plot Mesh

W  = zeros(nr_0,1);

for i=1:numElements

ind0 = globalnr_0(:,i);

w = interpolationfunction(0,Mesh.X(:,i),Mesh.Y(:,i),FunctionType);
W(ind0) = w;

end

disp('creation force function')
UV  = zeros(nr_1,1);
u = 0*globalnr_1v;
v = 0*globalnr_1h;
for i=1:numElements
    r = ceil(i/numColumns);
    c = i-(r-1)*numColumns;
    [U V] = interpolationfunction(1,r,c,FunctionType,Domain,DomInfo);
    u(:,i) = U;
    v(:,i) = V;
end

UV(globalnr_1v) = u;
UV(globalnr_1h) = v;


P = zeros(nr_2,1);
for r=1:numRows
    for c=1:numColumns
        p = interpolationfunction(2,r,c,FunctionType,Domain,DomInfo);
        i = c+(r-1)*numColumns;
        P(globalnr_2(:,i)) = p;
    end
end

WW = W(globalnr_0);
UU = UV(globalnr_1v);
VV = UV(globalnr_1h);
PP = P(globalnr_2);

Wp = reconstruct(0,WW,hp);
[Up,Vp,Velo] = reconstruct(1,UU,VV,hp,ep,Meshp);
Pp = reconstruct(2,PP,ep,Meshp);

disp('calculation of errors')
fouten_interp;

% Clim = [min(min(Wp)) max(max(Wp))];
% posten_vorticity
% pause
end % for N
end % for H

for i=1:length(HconvRange)-1
   ax = 1./(HconvRange);
   rate_w(i) = (log(L2_w(i+1))-log(L2_w(i))) / (log(ax(i+1))-log(ax(i)));
   rate_U(i) = (log(L2_U(i+1))-log(L2_U(i))) / (log(ax(i+1))-log(ax(i)));
   rate_p(i) = (log(L2_p(i+1))-log(L2_p(i))) / (log(ax(i+1))-log(ax(i)));
end


ax = 1./(HconvRange);
axisXY = [ 1e-2 1 1e-10 1];
str2 = 'c';
str3 = ['--v' str2];
figure(11)
loglog(ax,L2_U,str3,'markerface','w')
hold on
axis(axisXY)
title('u error')
figure(12)
loglog(ax,L2_w,str3,'markerface','w')
hold on
axis(axisXY)
title('w error')
figure(13)
loglog(ax,L2_p,str3,'markerface','w')
hold on
axis(axisXY)
title('p error')



% filename = [ 'LDC_Hconv_N' num2str(N) '_c' num2str(10*DomInfo) ];
% save(filename,'NrCellRange','HconvRange','L2_w','L2_u','L2_v','L2_U','L2_du','Energy','Energy_exact','Enstrophy','Enstrophy_exact')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
run Library_Convection/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%