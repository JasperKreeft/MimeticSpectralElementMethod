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

N = 20;
T = 2;
dt = 0.05;

filename = ['Two_circ_SE_N' num2str(N) '_dt' num2str(dt)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N); eta = xi;

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

Mesh = meshgenerator_square(Domain,DomInfo);
x = xi; y = x;

[Meshp,hh,ee] = postproces_grid_square(Domain,DomInfo);

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
y0 = 0;


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

% load T1

Up = ee'*reshape(U0,N,N)*ee;

if plot_fig || avimatlab
    subplot(1,2,2)
    surf(Xp,Yp,Up)
    % pcolor(Xp,Yp,Up)
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
    dataXY = [ Mesh.X Mesh.Y ];
    data = [ dataXY reshape(UU0,[],1) ];
    name = strcat(filename,'_00');
    MatlabToTecplot(name,name,'"X" "Y" "U"',[nn nn],data,2);
end

Up = reshape(Up,[],1);
Up_ex = exp(-((Meshp.X-x0).^2+(Meshp.Y-y0).^2)/(2*lambda^2));

L2error(1) = sqrt(sum(Meshp.W.*(Up-Up_ex).^2));
L2norm = sqrt(sum(Meshp.W.*Up_ex.^2));
Energy(1) = 1/2*sum(Up.*Meshp.W.*Up);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization

Ex = kron(ones(N),e');
Ey = zeros(N*(N+1),N*N);
for i=1:N
    Ey(1:N*(N+1),(i-1)*N+(1:N)) = kron(ones(N),e(i,:)');
end

WPe = zeros(N);
for i=1:N
    xp = (xi(i+1)+xi(i))/2+(xi(i+1)-xi(i))/2*xi;
    [~,ep] = MimeticpolyVal(xp,N,1);
    u = -2*pi*xp;
    for k=1:N
        WPe(i,k) = sum(w.*u.*ep(k,:))*(xi(i+1)-xi(i))/2;
    end
end
WPee = kron(WPe,ones(N+1,1));
WP = repmat(WPee,1,N);

WQe = zeros(N);
for j=1:N
    yq = (eta(j+1)+eta(j))/2+(eta(j+1)-eta(j))/2*eta;
    [~,eq] = MimeticpolyVal(yq,N,1);
    v = 2*pi*yq;
    for l=1:N
        WQe(j,l) = sum(w.*v.*eq(l,:))*(eta(j+1)-eta(j))/2;
    end
end
WQ = kron(WQe,ones(N+1,N));

V = [ Ex.*WQ ; Ey.*WP ];


D = div(N);

% Periodic BC
for i=2*N:-1:1
    ind1 = (i-1)*(N+1)+1;
    ind2 = i*(N+1);

    V(ind1,:) = (V(ind1,:)+V(ind2,:))/2;  % central discretization
    V(ind2,:) = [];

    D(:,ind1) = D(:,ind1)+D(:,ind2);
    D(:,ind2) = [];

end

M2 = innerproduct(2,Mesh.J);

A = M2*D*V;

K1 = (A-A')/2;  % Skew-symmetric
% K1 = A;         % Standard
K2 = M2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = U0;

t = 0; j = 1;
while roundn(t,-14)<T
    j = j+1
    t = t+dt;
    
    U = timemarching(K2,K1,dt,U,'ESDIRK5');
    
    Up = ee'*reshape(U,N,N)*ee;

    posten
    
    Xh = Meshp.X-x0*cos(2*pi*t)-y0*sin(2*pi*t);
    Yh = Meshp.Y+x0*sin(2*pi*t)-y0*cos(2*pi*t);
    Up_ex = exp(-(Xh.^2+Yh.^2)/(2*lambda^2));
    Up = reshape(Up,[],1);

    L2error(j) = sqrt(sum(Meshp.W.*(Up-Up_ex).^2));

    Energy(j) = 1/2*sum(Up.*Meshp.W.*Up);

end

taxis = linspace(0,t,j);
figure
subplot(1,2,1)
plot(taxis,Energy(1:j),'-o')
subplot(1,2,2)
plot(taxis,Energy(1:j)-Energy(1),'-o')

figure
plot(taxis,L2error(1:j)/L2norm,'-o')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%