% Code steady diffusion
%
% d^2u/dx^2 = f(x)
% with x=0: q0 = -(du/dx)|_{x=0}
%      x=L: u(x=L) = u_L

clear all
clf
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input section                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length of domain
L = 1.0;

% number of elements
Ne = 3;

% element polynomial order
Pe = [2 3 5];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generation of grid (h) and       %
% polynomials (p)                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ng,xe,xg,xen,xi,xj,h,c,Le,Le_j,w,Lo_j] = hp(Ne,Pe,L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mass matrix                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass = lumpedmassmatrix(Ne,Ng,Pe,w,c,h);
Mass = massmatrix(Ne,Ng,Pe,w,c,h,Le_j,Lo_j,xi,xj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diffusion matrix                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Diff = diffusionmatrix(Ne,Ng,Pe,w,c,h,Le_j,xj);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% force function                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [f,u_exact,u_exact_L2,bc,bc_l,bc_r] = exact_steadydiffusion(L,xg);
f = xg'; u_exact = [xg; -11/30*xg+1/10+1/6*xg.^3]; bc = 1; bc_l = .1; bc_r = -.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vector                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = zeros(Ng,1);

if bc == 2
    % Neumann boundary condition right
    b(Ng) = bc_r;
elseif bc == 3
    % Neumann boundary condition left
    b(1) = bc_l;
end

for j=1:Ng  
    b(j) = b(j)+Mass(j,:)*f;
end

% Left Dirichlet boundary condition
for kk = 1:Pe(1)+1
    cc = c(1,kk);
    b(cc) = b(cc)+Diff(cc,1)*bc_l;
end

% Right Dirichlet boundary condition
for kk = 1:Pe(Ne)+1
    cc = c(Ne,kk);
    b(cc) = b(cc)+Diff(cc,Ng)*bc_r;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solving system                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bc == 1
    % Left and right boundary is Dirichlet bc, so
    Diff = Diff(2:Ng-1,2:Ng-1);
    b = b(2:Ng-1);
elseif bc == 2
    % Left boundary is Dirichlet bc, so
    Diff = Diff(2:Ng,2:Ng);
    b = b(2:Ng);
elseif bc == 3
    % Right boundary is Dirichlet bc, so
    Diff = Diff(1:Ng-1,1:Ng-1);
    b = b(1:Ng-1);
end

% linear solver
u = -Diff\b;

% Adding Dirichlet bc to the solution
if bc == 1
    u = [bc_l; u ; bc_r];
elseif bc == 2
    u = [bc_l ; u];
elseif bc == 3
    u = [u ; bc_r];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocessen                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PostProcessen(Ne,xg,xe,xen,Pe,c,u,u_exact)