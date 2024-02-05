% Code steady diffusion
% 
% d^2u/dx^2 = f(x)
% with x=0: q0 = -(du/dx)|_{x=0}
%      x=L: u(x=L) = u_L

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input section                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% length of domain
L = 1.0;

% number of elements
Ne = 3;

% element polynomial order
Pe = [5 5 5];

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
Diff_string = '-D2u+u';
[f,u_exact,u_exact_L2,bc,bc_l,bc_r] = exact_solution(Diff_string,L,xg);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solving system                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bc == 1
    % Left and right boundary is Dirichlet bc, so
    Mass = Mass(2:Ng-1,2:Ng-1);
    Diff = Diff(2:Ng-1,2:Ng-1);
    b = b(2:Ng-1);
elseif bc == 2
    % Left boundary is Dirichlet bc, so
    Mass = Mass(2:Ng,2:Ng);
    Diff = Diff(2:Ng,2:Ng);
    b = b(2:Ng);
elseif bc == 3
    % Right boundary is Dirichlet bc, so
    Diff = Diff(1:Ng-1,1:Ng-1);
    Mass = Mass(1:Ng-1,1:Ng-1);
    b = b(1:Ng-1);
end

% linear solver
u = (Diff+Mass)\b;

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