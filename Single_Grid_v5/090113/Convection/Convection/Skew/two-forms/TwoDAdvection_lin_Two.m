clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call global variables

global N w e
global nn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Settings

Domain  = 'SinDeformGrid';
DomInfo = 0.0;

plot_fig  = 0;
avimatlab = 0;
Tecplot   = 0;

N = 20;
T = 10;
dt = 0.01;

filename = ['Two_lin_SE_N' num2str(N) '_dt' num2str(dt)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

Mesh = meshgenerator_square(Domain,DomInfo);
x = xi; y = x;

[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

Xp = reshape(Meshp.X,nn,nn);
Yp = reshape(Meshp.Y,nn,nn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize figure and plot Mesh

if ispc && (plot_fig || avimatlab)
    fig = figure('windowstyle','docked');
elseif  plot_fig || avimatlab
    fig = figure('Units','normalized','Position',[0, 0, 0.5, 1]);
end
if plot_fig
    set(gcf,'visible','on')
else
    set(gcf,'visible','off')
end


if plot_fig || avimatlab
    subplot(1,2,1)
    meshplot
    axis([min(Mesh.X) max(Mesh.X) min(Mesh.Y) max(Mesh.Y)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

U0 = initialforms(2,1,1,'gaussian',Domain,DomInfo);

Up = reconstruct(2,U0,ep,Meshp);

if plot_fig || avimatlab
subplot(1,2,2)
% surf(reshape(Meshp.X,nn,nn),reshape(Meshp.Y,nn,nn),reshape(Up,nn,nn))
pcolor(reshape(Meshp.X,nn,nn),reshape(Meshp.Y,nn,nn),reshape(Up,nn,nn))
shading interp
% colorbar('EastOutside')
axis equal
% axis([-1 1 -1 1 0 1])
axis([-1 1 -1 1])
% set(gca,'Ztick',[0 1])
set(gca,'clim',[0 1])
pause(0.05)
end

if Tecplot
    dataXY = [ Meshp.X Meshp.Y ];
    data = [ dataXY reshape(UU0,[],1) ];
    name = strcat(filename,'_00');
    MatlabToTecplot(name,name,'"X" "Y" "U"',[nn nn],data,2);
end

x0 = -1/2; y0 = 0; lambda = 1/8; % MOET ZELFDE ZIJN ALS IN initialforms
Up_ex = exp(-((Meshp.X-x0).^2+(Meshp.Y-y0).^2)/(2*lambda^2));

L2error(1) = sqrt(sum(Meshp.W.*(Up-Up_ex).^2));
L2norm = sqrt(sum(Meshp.W.*Up_ex.^2));
Energy(1) = 1/2*sum(Up.*Meshp.W.*Up);

Ax = 1;
Ay = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization

Ex = kron(eye(N),e');

Ey = zeros(N*(N+1),N*N);
for i=1:N
    Ey(1:N*(N+1),(i-1)*N+(1:N)) = kron(eye(N),e(i,:)');
end

vE = [ Ax*Ex; Ay*Ey ];

D = div(N);

% Periodic BC
for i=2*N:-1:1
    ind1 = (i-1)*(N+1)+1;
    ind2 = i*(N+1);
    
    vE(ind1,:) = (vE(ind1,:)+vE(ind2,:))/2;  % central discretization
    vE(ind2,:) = [];

    D(:,ind1) = D(:,ind1)+D(:,ind2);
    D(:,ind2) = [];
    
end

M2 = innerproduct(2,Mesh.J);

A = M2*D*vE;

K1 = (A-A')/2;  % Skew-symmetric
% K1 = A;         % Standard
K2 = M2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = U0;

t = 0; j=1;
Xh = Meshp.X-x0; Yh = Meshp.Y-y0;
while t<T
    j = j+1;
    t = t+dt
    
    U = timemarching(K2,K1,dt,U,'ERK3');
    
    Up = reconstruct(2,U,ep,Meshp);

    Xh = Xh - Ax*dt + 2*((Xh<-1)-(Xh>1));
    Yh = Yh - Ay*dt + 2*((Yh<-1)-(Yh>1));
    Up_ex = exp(-(Xh.^2+Yh.^2)/(2*lambda^2));

    L2error(j) = sqrt(sum(Meshp.W.*(Up-Up_ex).^2));

    Energy(j) = 1/2*sum(Up.*Meshp.W.*Up);
    
    posten

end
close all

taxis = linspace(0,t,j);
figure
subplot(1,2,1)
plot(taxis,Energy(1:j),'-o')
title('Energy')
subplot(1,2,2)
plot(taxis,Energy(1:j)-Energy(1),'-o')

figure
plot(taxis,L2error(1:j)/L2norm,'-o')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%