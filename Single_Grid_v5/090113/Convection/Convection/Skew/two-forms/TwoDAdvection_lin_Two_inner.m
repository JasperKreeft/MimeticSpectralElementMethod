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

plot_fig  = 1;
avimatlab = 0;
Tecplot   = 0;

N = 13;
T = 4;
dt = 0.05;

filename = ['Two_lin_SE_N_inner' num2str(N) '_dt' num2str(dt)];

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
    for i=1:N+1
        plot([ Mesh.X(i,1) Mesh.X(N*(N+1)+i,1) ],[ Mesh.Y(i,1) Mesh.Y(N*(N+1)+i,1) ],'linewidth',1)
        hold on
        plot([ Mesh.X((i-1)*(N+1)+1,1) Mesh.X(i*(N+1),1) ],[Mesh.Y((i-1)*(N+1)+1,1) Mesh.Y(i*(N+1),1)],'linewidth',1)
    end
    axis equal
    axis([min(Mesh.X) max(Mesh.X) min(Mesh.Y) max(Mesh.Y)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

lambda = 1/8;
x0 = -1/2;
y0 = -1/2;

% Gaussian function
% u0 = exp(-((X-x0).^2+(Y-y0).^2)/(2*lambda^2));

U0 = zeros(N*N,1);
for i=1:N
    for j=1:N
        ij = i+(j-1)*N;
        U0(ij,1) = pi/2*lambda^2*( ...
( erf( (x(i+1)-x0)/(sqrt(2)*lambda) ) - erf( (x(i)-x0)/(sqrt(2)*lambda) ) ) * ...
( erf( (y(j+1)-y0)/(sqrt(2)*lambda) ) - erf( (y(j)-y0)/(sqrt(2)*lambda) ) ) );
    end
end

Up0 = ep'*reshape(U0,N,N)*ep;

if plot_fig || avimatlab
    subplot(1,2,2)
    surf(Xp,Yp,Up0)
    % pcolor(Xp,Yp,Up0)
    shading interp
    % colorbar
    axis equal
    axis([-1 1 -1 1 -0.2 1.2])
    % axis([-1 1 -1 1])
    set(gca,'Ztick',[0 1])
    set(gca,'clim',[0 1])
    % pause(0.05)
end

if Tecplot
    dataXY = [ Meshp.X Meshp.Y ];
    data = [ dataXY reshape(UU0,[],1) ];
    name = strcat(filename,'_00');
    MatlabToTecplot(name,name,'"X" "Y" "U"',[nn nn],data,2);
end


Ax = 1;
Ay = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization


Ex = kron(e',eye(N));

Ey = zeros(N*(N+1),N*N);
for i=1:N+1
    Ey((i-1)*N+(1:N),:) = kron(eye(N),e(:,i)');
end

vE = [ -Ay*Ex; Ax*Ey ];

C = curl_in(N);

% Periodic BC
for i=[ 2*N*(N+1):-1:2*N*(N+1)-N+1 N*(N+1):-1:N*N+1 ]
    ind1 = i-N*N;
    ind2 = i;
    
    vE(ind1,:) = (vE(ind1,:)+vE(ind2,:))/2;  % central discretization
    vE(ind2,:) = [];

    C(:,ind1) = C(:,ind1)+C(:,ind2);
    C(:,ind2) = [];
    
end

M2 = innerproduct(2,Mesh.J);

A = M2*C*vE;

K1 = (A-A')/2;  % Skew-symmetric
% K1 = A;         % Standard
K2 = M2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


U = U0;

t = 0; j=0;
while t<T
    j = j+1
    t = t+dt;
    
    U = timemarching(K2,K1,dt,U,'ESDIRK5');
    
    Up = ep'*reshape(U,N,N)*ep;

    posten

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%