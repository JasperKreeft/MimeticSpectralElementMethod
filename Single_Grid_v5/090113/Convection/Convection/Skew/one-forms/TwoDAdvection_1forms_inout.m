clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load libraries

in = 'start';                                                   %#ok<NASGU>
if ispc
    
elseif isunix
run '/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Convection/Skew/Library_Convection/GetLibrary.m'
end

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

N = 24;
T = 4;
dt = 0.05;

filename = ['One_lin_SE_N' num2str(N) '_dt' num2str(dt)];

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

[U0x,U0y] = initialforms(1,1,1,'gaussian',Domain,DomInfo);

UU0y = ep'*reshape(U0y,N+1,N)'*hp;
UU0x = hp'*reshape(U0x,N+1,N)*ep;

if plot_fig || avimatlab
    subplot(1,2,1)
    surf(Xp,Yp,UU0x)
    shading interp
    axis equal
    set(gca,'clim',[-1 1])
    
    subplot(1,2,2)
    surf(Xp,Yp,UU0y)
    shading interp
    axis equal
    set(gca,'clim',[-1 1])
    pause(0.15)
end

if Tecplot
    dataXY = [ Mesh.X Mesh.Y ];
    data = [ dataXY reshape(UU0x+UU0y,[],1) ];
    name = strcat(filename,'_00');
    MatlabToTecplot(name,name,'"X" "Y" "U"',[nn nn],data,2);
end


Ax = 0;
Ay = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization

NG = normalgrad(N);
D  = div(N);
M1 = innerproduct(1,Mesh.J,Mesh.Qinv);

[vE01 vE12] = interiorproduct(1,Ax,Ay);

A = M1*(vE12*D+NG*vE01);
A = -A';


% K1 = (A-A')/2; % Skew-symmetric
K1 = A;      % Standard
K2 = M1;


U = [ U0x ; U0y ];

t = 0; j=0;
while t<T
    j = j+1
    t = t+dt;
    
    U = timemarching(K2,K1,dt,U,'gauss4');
    
    Ux = reshape(U(1:N*(N+1)),N+1,N);
    Uy = reshape(U(N*(N+1)+1:2*N*(N+1)),N+1,N)';
    
    UUx = hp'*Ux*ep;
    UUy = ep'*Uy*hp;
    
    if plot_fig || avimatlab
        subplot(1,2,1)
    surf(Xp,Yp,UUx)
    shading interp
    axis equal
    axis([-1 1 -1 1 -1.1 0.1])
    set(gca,'clim',[-1 1])
    subplot(1,2,2)
        surf(Xp,Yp,UUy)
    shading interp
    axis equal
    axis([-1 1 -1 1 -0.1 1.1])
    set(gca,'clim',[-1 1])
    pause(0.05)
    end

    if avimatlab
        avimovie(filename,fig,t==dt,t>=T);
    end

    if Tecplot
        data = [ dataXY reshape(UUx+UUy,[],1) ];
        if j<10
            name = strcat(filename,'_0',num2str(j));
        else
            name = strcat(filename,'_',num2str(j));
        end
        MatlabToTecplot(name,name,'"X" "Y" "U"',[nn nn],data,2);
    end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';
if ispc
    
elseif isunix
run '/media/My Passport/MSEM/MSEM_codes/Single_Grid_V5/Convection/Skew/Library_Convection/GetLibrary.m'
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%