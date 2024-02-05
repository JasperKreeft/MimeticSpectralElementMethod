close all
clear all
clc

global N w

%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotten = 1;
error   = 1;
casenr  = 3;

NrCellRange = 9%9:13;

c = 0.0;

for N=NrCellRange

disp(['N = ' num2str(N)])

%% Build grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[xi,w] = GLLnodes(N);     eta = xi;  % Gauss-Lobotto-Legendre

Xi = xi'*ones(1,N+1); Eta = Xi';

%% Topology %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NGp_h = zeros(N*(N+1),(N+1)^2);
for i=1:N
    unit = sparse([ zeros(1,i-1) -1 1 zeros(1,N-i) ]);
    NGp_h((i-1)*(N+1)+(1:N+1),:) = kron(speye(N+1),unit);
end

NGp_v = spdiags([-ones(N*(N+1),1) ones(N*(N+1),1)],[0 N+1],N*(N+1),(N+1)^2);

NGp = [ NGp_v ; -NGp_h ];


nr_cells = N*N;

% Construct the divergence matrix for the internal points.
Dp = zeros(N*N,N*(N+1)+N*(N+1));
for i = 1:N
    for j = 1:N
        Dp(i + (j-1)*N,i+(j-1)*(N+1))   = -1 ;
        Dp(i + (j-1)*N,i+1+(j-1)*(N+1)) =  1 ;
        Dp(i + (j-1)*N,N*(N+1) + j+(i-1)*(N+1))     = -1 ;
        Dp(i + (j-1)*N,N*(N+1) + j+1+(i-1)*(N+1))   =  1 ;
    end
end
Dp = sparse(Dp);

%% Basis functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h,e] = MimeticpolyVal(xi,N,1);

%% Grid generator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Uniform mapping [-1,1]^2 --> [0,pi]^2
% dXdXi  = pi/2*ones(size(Xi));
% dXdEta = zeros(size(Xi));
% dYdXi  = zeros(size(Xi));
% dYdEta = pi/2*ones(size(Xi));

% SinSin mapping
dXdXi  = 1+pi*c*cos(pi*Xi).*sin(pi*Eta);
dXdEta = pi*c*sin(pi*Xi).*cos(pi*Eta);
dYdXi  = pi*c*cos(pi*Xi).*sin(pi*Eta);
dYdEta = 1+pi*c*sin(pi*Xi).*cos(pi*Eta);

JGLLGLL = dXdXi.*dYdEta-dXdEta.*dYdXi;

J = reshape(JGLLGLL,1,(N+1)^2);

Jacobian = spdiags(kron(J,ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

qinv11 = kron(reshape(( dXdXi./JGLLGLL),1,(N+1)^2),[1 0])';
qinv22 = kron(reshape((dYdEta./JGLLGLL),1,(N+1)^2),[0 1])';
qinv12 = kron(reshape((dXdEta./JGLLGLL),1,(N+1)^2),[0 1])';
qinv21 = kron(reshape(( dYdXi./JGLLGLL),1,(N+1)^2),[1 0])';
Qinv = spdiags([qinv21 qinv11+qinv22 qinv12],-1:1,2*(N+1)^2,2*(N+1)^2);


%% Inner product 0-forms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M0 = innerproduct_zeroforms(J);

%% Inner product 1-forms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M1 = innerproduct_oneforms(e,J,Qinv);

%% Inner product 2-forms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M2 = innerproduct_twoforms(e,J);

%% Matrix compilation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Matrix = [ zeros(2*N*(N+1))        -M1*NGp             Dp'*M2
               -NGp'*M1               M0          zeros((N+1)^2,N^2)
                M2*Dp        zeros(N^2,(N+1)^2)    zeros(N^2)       ];

%% Forcing function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch casenr

    case 2
% case 2
fx = -sin(pi/2*xi)'*(sin(pi/2*eta(2:N+1))-sin(pi/2*eta(1:N)))...
     -2*pi*sin(pi*xi)'*(sin(pi*eta(2:N+1))-sin(pi*eta(1:N)));
fy = -(sin(pi/2*xi(2:N+1))-sin(pi/2*xi(1:N)))'*sin(pi/2*eta)...
     +2*pi*(sin(pi*xi(2:N+1))-sin(pi*xi(1:N)))'*sin(pi*eta);

    case 3
% Case 3
fx = 2*pi*cos(pi/2*xi)'*(cos(pi/2*eta(2:N+1))-cos(pi/2*eta(1:N)));
fy = zeros(N,N+1);

end
%% Righthandside %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = [ -M1*[reshape(fx,N*(N+1),1) ; reshape(fy',N*(N+1),1)]
      zeros((N+1)^2,1)
      zeros(N^2,1)     ];

%% Remove for Boundary Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

indu = sort([ 1:N+1:N*(N+1) N+1:N+1:N*(N+1) ]);
indv = indu + N*(N+1);
indw = 2*N*(N+1)+...
    sort([ 1:N+1 N+2:N+1:N*(N+1) 2*(N+1):N+1:N*(N+1) (N*(N+1)+1):(N+1)^2]);

ind = [ indu indv indw ];

Matrix(:,ind) = [];
Matrix(ind,:) = [];

F(ind,:) = [];
% keyboard
%% Solve system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uvwp = Matrix\F;

ind1 = N*(N+1)-length(indu);
ind2 = ind1+N*(N+1)-length(indv);
ind3 = ind2+(N+1)^2-length(indw);
ind4 = ind3+N^2;
U = uvwp(1:ind1);
V = uvwp(ind1+1:ind2);
W = uvwp(ind2+1:ind3);
P = uvwp(ind3+1:ind4);

%% Add Boundary Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UU = [ zeros(1,N) ; reshape(U,N-1,N) ; zeros(1,N) ] ;
VV = [ zeros(1,N) ; reshape(V,N-1,N) ; zeros(1,N) ]';
WW = [ zeros(N+1,1) [ zeros(1,N-1) ; reshape(W,N-1,N-1) ; zeros(1,N-1) ] zeros(N+1,1) ];
PP = reshape(P,N,N);

%% Postprocessen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plotten
nn = 100;
xx = linspace(-1,1,nn);
XX = xx'*ones(1,nn);
YY = XX';
[hh,dhhdxx] = LagrangeVal(xx,N,1);
ee = EdgeVal(dhhdxx);

uu = hh'*UU*ee;
vv = ee'*VV*hh;
ww = hh'*WW*hh;
pp = ee'*PP*ee;    pp = pp-mean(mean(pp));
ffx = hh'*fx*ee;
ffy = ee'*fy*hh;

figure(1)
subplot(2,3,1)
contourf(XX,YY,uu,20)
colorbar
title('u')
axis equal
axis([-1 1 -1 1])
subplot(2,3,2)
contourf(XX,YY,vv,20)
colorbar
title('v')
axis equal
axis([-1 1 -1 1])
subplot(2,3,3)
contourf(XX,YY,ww,20)
colorbar
title('w')
axis equal
axis([-1 1 -1 1])
subplot(2,3,4)
contourf(XX,YY,pp,20)
colorbar
title('p')
axis equal
axis([-1 1 -1 1])
subplot(2,3,5)
contourf(XX,YY,ffx,20)
colorbar
title('f_x')
axis equal
axis([-1 1 -1 1])
subplot(2,3,6)
contourf(XX,YY,ffy,20)
colorbar
title('f_y')
axis equal
axis([-1 1 -1 1])
pause(1)
end

%% Error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if error

nn = 100;
[xig,wg] = Gnodes(nn); etag = xig;
Xigg = xig'*ones(1,nn);  Etagg = Xigg';
Wgg = wg'*wg;

[hh ee ] = MimeticpolyVal(xig,N,1);

uu = hh'*UU*ee;
vv = ee'*VV*hh;
ww = hh'*WW*hh;
pp = ee'*PP*ee;   pp = pp-mean(mean(pp));
ffx = hh'*fx*ee;
ffy = ee'*fy*hh;

%% Exact solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch casenr

    case 2

    u_ex  =        -sin(pi*Xigg).*cos(pi*Etagg);
    v_ex  =         cos(pi*Xigg).*sin(pi*Etagg);
    w_ex  =   -2*pi*sin(pi*Xigg).*sin(pi*Etagg);
    p_ex  =       cos(pi/2*Xigg).*cos(pi/2*Etagg); p_ex = p_ex-mean(mean(p_ex));
    fx_ex = -pi/2*sin(pi/2*Xigg).*cos(pi/2*Etagg)-2*pi^2*sin(pi*Xigg).*cos(pi*Etagg);
    fy_ex = -pi/2*cos(pi/2*Xigg).*sin(pi/2*Etagg)+2*pi^2*cos(pi*Xigg).*sin(pi*Etagg);

    case 3

    u_ex  =      -cos(pi/2*Xigg).*sin(pi/2*Etagg);
    v_ex  =       sin(pi/2*Xigg).*cos(pi/2*Etagg);
    w_ex  =    pi*cos(pi/2*Xigg).*cos(pi/2*Etagg);
    p_ex  =   -pi*sin(pi/2*Xigg).*sin(pi/2*Etagg); p_ex = p_ex-mean(mean(p_ex));
    fx_ex = -pi^2*cos(pi/2*Xigg).*sin(pi/2*Etagg);
    fy_ex = 0*Xi;

end
%% L2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L2_u(N) = sqrt( sum(sum( (uu-u_ex).^2.*Wgg )) );
L2_v(N) = sqrt( sum(sum( (vv-v_ex).^2.*Wgg )) );
L2_w(N) = sqrt( sum(sum( (ww-w_ex).^2.*Wgg )) );
L2_p(N) = sqrt( sum(sum( (pp-p_ex).^2.*Wgg )) );

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

if error
figure
semilogy(L2_u,'-o')
hold on
semilogy(L2_v,'-^g')
semilogy(L2_w,'-xr')
semilogy(L2_p,'-sc')
legend('u','v','\omega','p')
end