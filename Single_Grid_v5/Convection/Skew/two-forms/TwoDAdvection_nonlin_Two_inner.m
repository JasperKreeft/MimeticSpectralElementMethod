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

Domain  = 'CurlCurl';
DomInfo = 0.0;

plot_fig  = 1;
avimatlab = 0;
Tecplot   = 0;

N = 30;
T = 4;
dt = 0.05;

filename = ['Two_lin_SE_N_inner' num2str(N) '_dt' num2str(dt)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N); eta = xi;

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

Mesh = meshgenerator_square(Domain,DomInfo);
x = Mesh.X(1:N+1); y = x;

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

% lambda = 1/8;
% x0 = -1/2;
% y0 = -1/2;
% 
% % Gaussian function
% % w0 = exp(-((X-x0).^2+(Y-y0).^2)/(2*lambda^2));
% 
% W0 = zeros(N*N,1);
% for i=1:N
%     for j=1:N
%         ij = i+(j-1)*N;
%         W0(ij,1) = pi/2*lambda^2*( ...
% ( erf( (x(i+1)-x0)/(sqrt(2)*lambda) ) - erf( (x(i)-x0)/(sqrt(2)*lambda) ) ) * ...
% ( erf( (y(j+1)-y0)/(sqrt(2)*lambda) ) - erf( (y(j)-y0)/(sqrt(2)*lambda) ) ) );
%     end
% end

a = 0.3;
U = 1;
x0 = -0.4;
y0 = 0;

W0 = zeros(N*N,1);
for i=1:N
    for j=1:N
        ij = i+(j-1)*N;
        W0(ij,1) = U*sqrt(pi*exp(1)/2)*( ...
( (x(i+1)-x0)*exp(-(x(i+1)-x0)^2/(2*a^2))-(x(i)-x0)*exp(-(x(i)-x0)^2/(2*a^2)) ) * ( erf( (y(j+1)-y0)/(sqrt(2)*a) ) - erf( (y(j)-y0)/(sqrt(2)*a) ) ) + ...
( (y(j+1)-y0)*exp(-(y(j+1)-y0)^2/(2*a^2))-(y(j)-y0)*exp(-(y(j)-y0)^2/(2*a^2)) ) * ( erf( (x(i+1)-x0)/(sqrt(2)*a) ) - erf( (x(i)-x0)/(sqrt(2)*a) ) ) );
    end
end

x0 = 0.4;
for i=1:N
    for j=1:N
        ij = i+(j-1)*N;
        W0(ij,1) = W0(ij,1) + U*sqrt(pi*exp(1)/2)*( ...
( (x(i+1)-x0)*exp(-(x(i+1)-x0)^2/(2*a^2))-(x(i)-x0)*exp(-(x(i)-x0)^2/(2*a^2)) ) * ( erf( (y(j+1)-y0)/(sqrt(2)*a) ) - erf( (y(j)-y0)/(sqrt(2)*a) ) ) + ...
( (y(j+1)-y0)*exp(-(y(j+1)-y0)^2/(2*a^2))-(y(j)-y0)*exp(-(y(j)-y0)^2/(2*a^2)) ) * ( erf( (x(i+1)-x0)/(sqrt(2)*a) ) - erf( (x(i)-x0)/(sqrt(2)*a) ) ) );
    end
end


Wp0 = reconstruct(2,W0,ep,Meshp);

if plot_fig || avimatlab
    subplot(1,2,2)
    pcolor(Xp,Yp,reshape(Wp0,nn,nn))
    % pcolor(Xp,Yp,Wp0)
    shading interp
    colorbar
    axis equal
%     axis([-1 1 -1 1 -0.2 1.2])
    axis([-pi pi -pi pi])
    axis off
%     set(gca,'Ztick',[0 1])
%     set(gca,'clim',[0 1])
    % pause(0.05)
end

if Tecplot
    dataXY = [ Meshp.X Meshp.Y ];
    data = [ dataXY reshape(Wp0,[],1) ];
    name = strcat(filename,'_00');
    MatlabToTecplot(name,name,'"X" "Y" "W"',[nn nn],data,2);
end

% % Linear velocity field
% U = kron(ones(N+1,1),diff(x)');
% V = kron(ones(N+1,1),diff(y)');

% Circular velocity field
U = zeros(N*(N+1),1); V = zeros(N*(N+1),1);
for i1=1:N
    for i2=1:N+1
        ij = i1+(i2-1)*N;
        U(ij) =  2*pi*(x(i1+1)-x(i1))*y(i2);
        V(ij) = -2*pi*x(i2)*(y(i1+1)-y(i1));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M1 = innerproduct(1,Mesh.J,Mesh.Qinv);
M2 = innerproduct(2,Mesh.J);

C = curl_in(N);

Matrix = [ -M1 C'*M2 ; M2*C zeros(N*N) ];
Rhs = [ zeros(2*N*(N+1),1) ; M2*W0 ];

SigmaPhi = Matrix\Rhs;
Phi = SigmaPhi(2*N*(N+1)+(1:N*N));

Phip = reconstruct(2,Phi,ep,Meshp);

figure
subplot(1,2,1)
pcolor(Xp,Yp,reshape(Phip,nn,nn))
% pcolor(Xp,Yp,Phip)
shading interp
colorbar
axis equal
% axis([-1 1 -1 1 -0.2 1.2])
axis([-pi pi -pi pi])
axis off
% set(gca,'Ztick',[0 1])
% set(gca,'clim',[0 1])
% pause(0.05)

UV = (M1\C'*M2)*Phi;
U = UV(1:N*(N+1));
V = UV(N*(N+1)+(1:N*(N+1)));

[ux,uy,uMag] = reconstruct_oneforms_in2(U,V,hp,ep,Meshp);

subplot(1,2,2)
quiver(Meshp.X,Meshp.Y,ux,uy)

break

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization


Ex = kron(e',ones(N));

Ey = zeros(N*(N+1),N*N);
for i=1:N+1
    Ey((i-1)*N+(1:N),:) = kron(ones(N),e(:,i)');
end

I = zeros(N+1,N,N+1);
for m=1:N
    x_p = (xi(m+1)+xi(m))/2+(xi(m+1)-xi(m))/2*xi;
    [h_p,e_p] = MimeticpolyVal(x_p,N,1);
    for k=1:N
        for i=1:N+1
            I(i,k,m) = sum(w.*h_p(i,:).*e_p(k,:))*(xi(m+1)-xi(m))/2;
        end
    end
end

VI = zeros(N*(N+1),N*N);
for n=1:N+1
    Vy = sum(reshape(V,N,N+1).*kron(e(:,n),ones(1,N+1)),1);
    for m=1:N
        mn = m+(n-1)*N;
        for k=1:N
            for l=1:N
                kl = k+(l-1)*N;
                VI(mn,kl) = Vy*I(:,k,m);
            end
        end
    end
end

UI = zeros(N*(N+1),N*N);
for m=1:N+1
    Ux = sum(reshape(U,N,N+1).*kron(e(:,m),ones(1,N+1)),1);
    for n=1:N
        mn = n+(m-1)*N;
        for k=1:N
            for l=1:N
                kl = k+(l-1)*N;
                UI(mn,kl) = Ux*I(:,l,n);
            end
        end
    end
end

vE = [ -VI.*Ex ; UI.*Ey ];

% Periodic BC
for i=[ 2*N*(N+1):-1:2*N*(N+1)-N+1 N*(N+1):-1:N*N+1 ]
    ind1 = i-N*N;
    ind2 = i;
    
    vE(ind1,:) = (vE(ind1,:)+vE(ind2,:))/2;  % central discretization
    vE(ind2,:) = [];

    C(:,ind1) = C(:,ind1)+C(:,ind2);
    C(:,ind2) = [];
    
end


A = M2*C*vE;

K1 = (A-A')/2;  % Skew-symmetric
% K1 = A;         % Standard
K2 = M2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


W = W0;

t = 0; j=0;
while t<T
    j = j+1
    t = t+dt;
    
    W = timemarching(K2,K1,dt,W,'ESDIRK5');
    
    Up = ep'*reshape(W,N,N)*ep; % Up ipv Wp ivm posten

    posten

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Close libraries

in = 'finish';                                                  %#ok<NASGU>
run 'C:\Users\Jasper\Documents\MATLAB\MSEM_codes\Single_Grid_V5\Convection\Skew/Library_Convection\GetLibrary.m'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%