clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run Library/GetLibrary.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N numElements numRows numColumns
global xi
global h e
global globalnr_0 globalnr_1v globalnr_1h globalnr_2
global nr_0 nr_1 nr_2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

FunctionType = 'AFG';
Domain       = 'SinDeformGrid_01';
DomInfo      = 0.0;

NrCellRange = 3;
HconvRange  = 2.^(1:6);

plot_figures  = 1;
error_figures = 0;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('Stokes Problem ')
disp('The Arnold-Falk-Gopalakrishnan (2012) testcase')
disp(['Testfunction is ' FunctionType])
disp(['Computational Domain is ' Domain ' with settings ' num2str(DomInfo)])
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start H-convergence loop

er = 0;
L2_u = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_v = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_p = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_mx = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_my = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_txx = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_tyx = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_txy = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_tyy = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_U = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_M = zeros(1,max(length(NrCellRange),length(HconvRange)));
L2_T = zeros(1,max(length(NrCellRange),length(HconvRange)));

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
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh generation

numbering('square');

Mesh = meshgenerator_square(Domain,DomInfo);
[Meshx,Meshy] = meshgenerator_square_covector(Domain,DomInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

PP = zeros(N^2,Hconv^2);
for hy=1:Hconv
    for hx=1:Hconv
        hh = hx+(hy-1)*Hconv;
        for j=1:N
            for i=1:N
                IJ = i+(j-1)*N;
                ij = i+(j-1)*(N+1);
                ji = j+(i-1)*(N+1);
PP(IJ,hh) = 1/6 * ( (Mesh.X(ij+1,hh)-0.5)^6*Mesh.Y(ij+N+1,hh) + (Mesh.Y(ij+N+1,hh)-0.5)^6*Mesh.X(ij+1,hh) ...
                  -(Mesh.X(ij,hh)-0.5)^6*Mesh.Y(ij+N+1,hh)   - (Mesh.Y(ij+N+1,hh)-0.5)^6*Mesh.X(ij,hh) ...
                  -(Mesh.X(ij+1,hh)-0.5)^6*Mesh.Y(ij,hh)     - (Mesh.Y(ij,hh)-0.5)^6*Mesh.X(ij+1,hh) ...
                  +(Mesh.X(ij,hh)-0.5)^6*Mesh.Y(ij,hh)       + (Mesh.Y(ij,hh)-0.5)^6*Mesh.X(ij,hh) );
            end
        end
    end
end


UU = zeros(N*(N+1),Hconv^2);
for hy=1:Hconv
    for hx=1:Hconv
        hh = hx+(hy-1)*Hconv;
        for j=1:N
            for i=1:N+1
                ij = i+(j-1)*(N+1);
UU(ij,hh) = -(2*Mesh.X(ij,hh)^4-4*Mesh.X(ij,hh)^3+2*Mesh.X(ij,hh)^2) * ( (1/2*Mesh.Y(ij+N+1,hh)^4-Mesh.Y(ij+N+1,hh)^3+1/2*Mesh.Y(ij+N+1,hh)^2) - (1/2*Mesh.Y(ij,hh)^4-Mesh.Y(ij,hh)^3+1/2*Mesh.Y(ij,hh)^2) );
            end
        end
    end
end

VV = zeros(N*(N+1),Hconv^2);
for hy=1:Hconv
    for hx=1:Hconv
        hh = hx+(hy-1)*Hconv;
        for i=1:N
            for j=1:N+1
                ji = j+(i-1)*(N+1);
                ij = i+(j-1)*(N+1);
VV(ji,hh) = (2*Mesh.Y(ij,hh)^4-4*Mesh.Y(ij,hh)^3+2*Mesh.Y(ij,hh)^2) * ( (1/2*Mesh.X(ij+1,hh)^4-Mesh.X(ij+1,hh)^3+1/2*Mesh.X(ij+1,hh)^2) - (1/2*Mesh.X(ij,hh)^4-Mesh.X(ij,hh)^3+1/2*Mesh.X(ij,hh)^2) );
            end
        end
    end
end

MMx = zeros(N*(N+1),Hconv^2);
for hy=1:Hconv
    for hx=1:Hconv
        hh = hx+(hy-1)*Hconv;
        for j=1:N
            for i=1:N+1
                IJ = i+(j-1)*(N+1);
                ij = i+(j-1)*(N+2);
MMx(IJ,hh) = -((2/5*Meshx.X(ij+1,hh)^5-Meshx.X(ij+1,hh)^4+2/3*Meshx.X(ij+1,hh)^3) - (2/5*Meshx.X(ij,hh)^5-Meshx.X(ij,hh)^4+2/3*Meshx.X(ij,hh)^3) ) * ( (1/2*Meshx.Y(ij+N+2,hh)^4-Meshx.Y(ij+N+2,hh)^3+1/2*Meshx.Y(ij+N+2,hh)^2) - (1/2*Meshx.Y(ij,hh)^4-Meshx.Y(ij,hh)^3+1/2*Meshx.Y(ij,hh)^2) );
            end
        end
    end
end

MMy = zeros(N*(N+1),Hconv^2);
for hy=1:Hconv
    for hx=1:Hconv
        hh = hx+(hy-1)*Hconv;
        for i=1:N
            for j=1:N+1
                ji = j+(i-1)*(N+1);
                ij = i+(j-1)*(N+1);
MMy(ji,hh) = ( (2/5*Meshy.Y(ij+N+1,hh)^5-Meshy.Y(ij+N+1,hh)^4+2/3*Meshy.Y(ij+N+1,hh)^3) - (2/5*Meshy.Y(ij,hh)^5-Meshy.Y(ij,hh)^4+2/3*Meshy.Y(ij,hh)^3) ) * ( (1/2*Meshy.X(ij+1,hh)^4-Meshy.X(ij+1,hh)^3+1/2*Meshy.X(ij+1,hh)^2) - (1/2*Meshy.X(ij,hh)^4-Meshy.X(ij,hh)^3+1/2*Meshy.X(ij,hh)^2) );
            end
        end
    end
end

Txx = zeros(N*(N+2),Hconv^2);
for hy=1:Hconv
    for hx=1:Hconv
        hh = hx+(hy-1)*Hconv;
        for j=1:N
            for i=1:N+2
                ij = i+(j-1)*(N+2);
Txx(ij,hh) = -2*( 8*Meshx.X(ij,hh)^3-12*Meshx.X(ij,hh)^2+4*Meshx.X(ij,hh) ) * ( (1/2*Meshx.Y(ij+N+2,hh)^4-Meshx.Y(ij+N+2,hh)^3+1/2*Meshx.Y(ij+N+2,hh)^2)-(1/2*Meshx.Y(ij,hh)^4-Meshx.Y(ij,hh)^3+1/2*Meshx.Y(ij,hh)^2) );
            end
        end
    end
end

Tyx = zeros((N+1)^2,Hconv^2);
for hy=1:Hconv
    for hx=1:Hconv
        hh = hx+(hy-1)*Hconv;
        for j=1:N+1
            for i=1:N+1
                IJ = i+(j-1)*(N+1);
                ij = i+(j-1)*(N+2);
Tyx(IJ,hh) = -( 6*Meshx.Y(ij,hh)^2-6*Meshx.Y(ij,hh)+1 ) * ( (2/5*Meshx.X(ij+1,hh)^5-Meshx.X(ij+1,hh)^4+2/3*Meshx.X(ij+1,hh)^3)-(2/5*Meshx.X(ij,hh)^5-Meshx.X(ij,hh)^4+2/3*Meshx.X(ij,hh)^3) ) ...
             +(2*Meshx.Y(ij,hh)^4-4*Meshx.Y(ij,hh)^3+2*Meshx.Y(ij,hh)^2) * ( (2*Meshx.X(ij+1,hh)^3-3*Meshx.X(ij+1,hh)^2+Meshx.X(ij+1,hh)) - (2*Meshx.X(ij,hh)^3-3*Meshx.X(ij,hh)^2+Meshx.X(ij,hh)) );
            end
        end
    end
end


Txy = zeros((N+1)^2,Hconv^2);
for hy=1:Hconv
    for hx=1:Hconv
        hh = hx+(hy-1)*Hconv;
        for j=1:N+1
            for i=1:N+1
                ij = i+(j-1)*(N+1);
                ji = j+(i-1)*(N+1);
Txy(ji,hh) = ( 6*Meshy.X(ij,hh)^2-6*Meshy.X(ij,hh)+1 ) * ( (2/5*Meshy.Y(ij+N+1,hh)^5-Meshy.Y(ij+N+1,hh)^4+2/3*Meshy.Y(ij+N+1,hh)^3)-(2/5*Meshy.Y(ij,hh)^5-Meshy.Y(ij,hh)^4+2/3*Meshy.Y(ij,hh)^3) ) ...
             -(2*Meshy.X(ij,hh)^4-4*Meshy.X(ij,hh)^3+2*Meshy.X(ij,hh)^2) * ( (2*Meshy.Y(ij+N+1,hh)^3-3*Meshy.Y(ij+N+1,hh)^2+Meshy.Y(ij+N+1,hh)) - (2*Meshy.Y(ij,hh)^3-3*Meshy.Y(ij,hh)^2+Meshy.Y(ij,hh)) );
            end
        end
    end
end




Tyy = zeros(N*(N+2),Hconv^2);
for hy=1:Hconv
    for hx=1:Hconv
        hh = hx+(hy-1)*Hconv;
        for i=1:N
            for j=1:N+2
                ji = j+(i-1)*(N+2);
                ij = i+(j-1)*(N+1);
Tyy(ji,hh) = 2*(8*Meshy.Y(ij,hh)^3-12*Meshy.Y(ij,hh)^2+4*Meshy.Y(ij,hh)) * ( (1/2*Meshy.X(ij+1,hh)^4-Meshy.X(ij+1,hh)^3+1/2*Meshy.X(ij+1,hh)^2) - (1/2*Meshy.X(ij,hh)^4-Meshy.X(ij,hh)^3+1/2*Meshy.X(ij,hh)^2) );
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Postprocessen

disp('Postprocessen')

% Grid, basis-functions and weights for post-processen
[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);
[hegp,eegp] = MimeticpolyVal(Meshp.xip,N,3);

% Reconstruction
[uu,vv,velo] = reconstruct(1,UU,VV,hp,ep,Meshp);
[mmx,mmy] = reconstruct(21,MMx,MMy,ep,eegp,Meshp);
pp           = reconstruct(2,PP,ep,Meshp);
[txx,tyx,txy,tyy] = reconstruct(11,Txx,Tyx,Txy,Tyy,hp,ep,hegp,eegp,Meshp);


x = Meshp.X; y = Meshp.Y;
% curlw2_ex =  8*y.^2.*(y-1).^2.*(x-1)+4*y.^2.*(y-1).^2.*(2*x-1)+8*y.^2.*(y-1).^2.*x+...
%     4*x.*(x-1).^2.*(2*y-1).*(y-1)+4*x.^2.*(x-1).*(2*y-1).*(y-1)+8*x.*(x-1).^2.*y.*(y-1)+...
%     8*x.^2.*(x-1).*y.*(y-1)+4*x.*(x-1).^2.*y.*(2*y-1)+4*x.^2.*(x-1).*y.*(2*y-1);
% curlw1_ex = -4*y.*(y-1).^2.*(2*x-1).*(x-1)-4*y.^2.*(y-1).*(2*x-1).*(x-1)-...
%     8*y.*(y-1).^2.*x.*(x-1)-8*y.^2.*(y-1).*x.*(x-1)-4*y.*(y-1).^2.*x.*(2*x-1)-4*y.^2.*(y-1).*x.*(2*x-1)-...
%     8*x.^2.*(x-1).^2.*(y-1)-4*x.^2.*(x-1).^2.*(2*y-1)-8*x.^2.*(x-1).^2.*y;
% divu_ex = -4*x.*(x-1).^2.*y.*(2*y-1).*(y-1)-4*x.^2.*(x-1).*y.*(2*y-1).*(y-1)+4*y.*(y-1).^2.*x.*(2*x-1).*(x-1)+4*y.^2.*(y-1).*x.*(2*x-1).*(x-1);

u_ex = -(2*x.^4-4*x.^3+2*x.^2).*(2*y.^3-3*y.^2+y);
v_ex =  (2*x.^3-3*x.^2+x).*(2*y.^4-4*y.^3+2*y.^2);
p_ex = (x-1/2).^5+(y-1/2).^5;
mx_ex = u_ex;
my_ex = v_ex;
txx_ex = -2*(8*x.^3-12*x.^2+4*x).*(2*y.^3-3*y.^2+y);
tyx_ex = -(2*x.^4-4*x.^3+2*x.^2).*(6*y.^2-6*y+1)+(2*y.^4-4*y.^3+2*y.^2).*(6*x.^2-6*x+1);
txy_ex = tyx_ex;
tyy_ex = 2*(8*y.^3-12*y.^2+4*y).*(2*x.^3-3*x.^2+x);




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% plotten
% 
% if plot_figures
%     disp('creation colorplots')
%     meshplot
%     plotten_me
% end
% 
% %% Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if error_figures
%     disp('calculation of errors')
%     fout_Stokes;
% end
% 
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



er = er+1;

if ~exist('numElements','var')
    numElements = 1;
end

for i=1:numElements

JW = Meshp.J(:,i).*Meshp.W;

%% Potential error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error is calculated on standard element
L2_u(er) = L2_u(er) + sum((uu(:,i)-u_ex(:,i)).^2.*JW);
L2_v(er) = L2_v(er) + sum((vv(:,i)-v_ex(:,i)).^2.*JW);
L2_p(er) = L2_p(er) + sum((pp(:,i)-p_ex(:,i)).^2.*JW);
L2_mx(er) = L2_mx(er) + sum((mmx(:,i)-mx_ex(:,i)).^2.*JW);
L2_my(er) = L2_my(er) + sum((mmy(:,i)-my_ex(:,i)).^2.*JW);
L2_txx(er) = L2_txx(er) + sum((txx(:,i)-txx_ex(:,i)).^2.*JW);
L2_tyx(er) = L2_tyx(er) + sum((tyx(:,i)-tyx_ex(:,i)).^2.*JW);
L2_txy(er) = L2_txy(er) + sum((txy(:,i)-txy_ex(:,i)).^2.*JW);
L2_tyy(er) = L2_tyy(er) + sum((tyy(:,i)-tyy_ex(:,i)).^2.*JW);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L2_U(er) = sqrt( L2_u(er) + L2_v(er) );
L2_p(er) = sqrt(L2_p(er));
L2_M(er) = sqrt( L2_mx(er) + L2_my(er) );
L2_T(er) = sqrt( L2_txx(er) + L2_tyx(er) + L2_txy(er) + L2_tyy(er) );


end % for N
end % for H


if 0 %figplot
    
global nn

if isempty(nn)
    nn = sqrt(size(Meshp.X,1));
end

XYlim = [min(min(Meshp.X)) max(max(Meshp.X)) min(min(Meshp.Y)) max(max(Meshp.Y))];

h1 = figure('visible','off');
h5 = figure('visible','off');

for i=1:numElements
% 
    Xp = reshape(Meshp.X(:,i),nn,nn);
    Yp = reshape(Meshp.Y(:,i),nn,nn);
%     
    up = reshape(uu(:,i),nn,nn);
    vp = reshape(vv(:,i),nn,nn);
    Vp = reshape(velo(:,i),nn,nn);
    mxp = reshape(mmx(:,i),nn,nn);
    myp = reshape(mmy(:,i),nn,nn);
    Pp = reshape(pp(:,i),nn,nn);
    txxp = reshape(txx(:,i),nn,nn);
    tyxp = reshape(tyx(:,i),nn,nn);
    txyp = reshape(txy(:,i),nn,nn);
    tyyp = reshape(tyy(:,i),nn,nn);
    
    figure(h1)
    set(h1,'visible','off')
    subplot(2,3,1)
    pcolor(Xp,Yp,up)
    hold on
    shading interp
    % colorbar
    title('u')
    axis equal
    axis(XYlim)
    subplot(2,3,2)
    pcolor(Xp,Yp,vp)
    hold on
    shading interp
    % colorbar
    title('v')
    axis equal
    axis(XYlim)
    subplot(2,3,3)
    pcolor(Xp,Yp,Vp)
    hold on
    shading interp
    % colorbar
    title('abs(velo)')
    axis equal
    axis(XYlim)
    subplot(2,3,4)
    pcolor(Xp,Yp,Pp)
    hold on
    % set(gca,'clim',[-5 5])
    shading interp
    % colorbar
    title('p')
    axis equal
    axis(XYlim)
    subplot(2,3,5)
    pcolor(Xp,Yp,mxp)
    hold on
    shading interp
    % colorbar
    title('m_x')
    axis equal
    axis(XYlim)
    subplot(2,3,6)
    pcolor(Xp,Yp,myp)
    hold on
    shading interp
    % colorbar
    title('m_y')
    axis equal
    axis(XYlim)
%     mesh(reshape(Mesh.X(:,i),N+1,N+1),reshape(Mesh.Y(:,i),N+1,N+1),zeros(N+1),'EdgeColor','black')
%     hold on
%     view([0 0 1])
%     axis equal
%     axis(XYlim)

    figure(h5)
    set(h5,'visible','off')
    subplot(2,2,1)
    pcolor(Xp,Yp,txxp)
    hold on
    shading interp
    colorbar
    title('\tau_x_x')
    axis equal
    axis(XYlim)
    subplot(2,2,2)
    pcolor(Xp,Yp,tyxp)
    hold on
    shading interp
    colorbar
    title('\tau_y_x')
    axis equal
    axis(XYlim)
    subplot(2,2,3)
    pcolor(Xp,Yp,txyp)
    hold on
    shading interp
    colorbar
    title('\tau_x_y')
    axis equal
    axis(XYlim)
    subplot(2,2,4)
    pcolor(Xp,Yp,tyyp)
    hold on
    shading interp
    colorbar
    title('\tau_y_y')
    axis equal
    axis(XYlim)


end

set(h1,'visible','on')
set(h5,'visible','on')

end



ax = 1./(HconvRange);
axisXY = [ 0.05 3 1e-10 1];
str1 = '-sg';
str2 = 'g';
figure(10)
loglog(ax,L2_U,str1,'markerface',str2)
hold on
% axis(axisXY)
title('u error Stokes')
figure(20)
loglog(ax,L2_p,str1,'markerface',str2)
hold on
% axis(axisXY)
title('p error Stokes')
figure(30)
loglog(ax,L2_M,str1,'markerface',str2)
hold on
% axis(axisXY)
title('M error Stokes')
figure(40)
loglog(ax,L2_T,str1,'markerface',str2)
hold on
% axis(axisXY)
title('T error Stokes')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
GetLibrary

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%