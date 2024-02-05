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

N = 20;
T = 10;
dt = 0.05;

filename = ['Zero_circ_SE_N' num2str(N) '_dt' num2str(dt)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N); eta = xi;

[h,e] = MimeticpolyVal(xi,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mesh

Mesh = meshgenerator_square(Domain,DomInfo);
x = xi; y = x;

[Meshp,hp,ep] = postproces_grid_square(Domain,DomInfo);

Xp = reshape(Meshp.X,nn,nn);
Yp = reshape(Meshp.Y,nn,nn);

XY = [ Mesh.X Mesh.Y ];

subplot(1,2,1)
for i=1:N+1
    plot([ Mesh.X(i,1) Mesh.X(N*(N+1)+i,1) ],[ Mesh.Y(i,1) Mesh.Y(N*(N+1)+i,1) ],'linewidth',1)
    hold on
    plot([ Mesh.X((i-1)*(N+1)+1,1) Mesh.X(i*(N+1),1) ],[Mesh.Y((i-1)*(N+1)+1,1) Mesh.Y(i*(N+1),1)],'linewidth',1)
end
axis equal
axis([min(Mesh.X) max(Mesh.X) min(Mesh.Y) max(Mesh.Y)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial solution

lambda = 1/8;
% x0 = -1/2;
% y0 = 1/2;

x0 = 0;
y0 = -1/2;


% Gaussian function
U0 = exp(-((Mesh.X-x0).^2+(Mesh.Y-y0).^2)/(2*lambda^2));

Up = hp'*reshape(U0,N+1,N+1)*hp;

subplot(1,2,2)
% surf(Xp,Yp,Up)
pcolor(Xp,Yp,Up)
shading interp
% colorbar
axis equal
% axis([-1 1 -1 1 0 1])
axis([-1 1 -1 1])
% set(gca,'Ztick',[0 1])
set(gca,'clim',[0 1])
pause(0.05)

Up = reshape(Up,[],1);
Up_ex = exp(-((Meshp.X-x0).^2+(Meshp.Y-y0).^2)/(2*lambda^2));

L2error(1) = sqrt(sum(Meshp.W.*(Up-Up_ex).^2));
L2norm = sqrt(sum(Meshp.W.*Up_ex.^2));
Energy(1) = 1/2*sum(Up.*Meshp.W.*Up);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Discretization

Ex = kron(eye(N+1),e');
Ey = zeros((N+1)^2,N*(N+1));
for i=1:N+1
    Ey((i-1)*(N+1)+(1:N+1),1:N*(N+1)) = kron(eye(N+1),e(:,i)');
end

Ax =  2*pi*XY(:,2);
Ay = -2*pi*XY(:,1);

VE = [ diag(Ax)*Ex diag(Ay)*Ey ];

G = grad_in(N);

M0 = spdiags(kron(w,w)',0,(N+1)^2,(N+1)^2);

A = M0*VE*G;

K1 = (A-A')/2;  % Skew-symmetric
% K1 = A;         % Standard
K2 = M0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = U0;

t = 0; j=1;
while t<T

    t = t+dt;
    j = j+1;

    U = timemarching(K2,K1,dt,U,'gauss4');

    Up = hp'*reshape(U,N+1,N+1)*hp;

%     surf(Xp,Yp,Up)
    pcolor(Xp,Yp,Up)
    shading interp
    colorbar
    axis equal
%     axis([-1 1 -1 1 -0.2 1.2])
    axis([-1 1 -1 1])
%     set(gca,'Ztick',[0 1])
    set(gca,'clim',[0 1])
    pause(0.05)

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