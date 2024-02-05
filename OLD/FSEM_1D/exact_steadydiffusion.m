function [rhs,u_exact,u_exact_L2,bc,bc_l,bc_r] = exact_steadydiffusion(L,xg)

% clear all
% close all
% clc

% exact
% 
% L = 1; xg = linspace(0,L,3);

syms x a b

f = input('Enter forcing function: f(x) = ');

rhs = subs(f,x,xg)';

% boundary conditions
bc = 0;
disp('  ')
disp('Choose type of boundary conditions:')
disp('  ')
disp('  1. Left and right Dirichlet bc`s')
disp('  2. Left Dirichlet & right Neumann bc`s')
disp('  3. Left Neumann & right Dirichlet bc`s')
disp('  ')  
while bc<1 || bc>3 || rem(bc,1)~=0
    bc = input('enter the temperal method you want to use:  ');
    disp('  ')
end


u_hom = a*x+b;

F = int(f,x);

u_part = int(F,x);

u = u_hom+u_part;

if bc == 1
    u_l = input('enter value of left bc: u_l = ');
    disp('  ')
    u_r = input('enter value of right bc: u_r = ');
    disp('  ')
    
    bb = u_l-subs(u_part,x,0);
    aa = (u_r-bb-subs(u_part,x,L))/L;
    
    bc_l = u_l;
    bc_r = u_r;
    
elseif bc == 2
    u_l = input('enter value of left bc: u_l = ');
    disp('  ')
    du_r = input('enter value of right bc: du_r = ');
    disp('  ')
    
    aa = du_r-subs(F,x,L);
    bb = u_l-subs(u_part,x,0);
    
    bc_l = u_l;
    bc_r = du_r;

elseif bc == 3
    du_l = input('enter value of left bc: du_l = ');
    disp('  ')
    u_r = input('enter value of right bc: u_r = ');
    disp('  ')
    
    aa = du_l-subs(F,x,0);
    bb = u_r-subs(u_part,x,l)-aa*l;

    bc_l = du_l;
    bc_r = u_r;
    
end


u = subs(u,{a b},[aa bb])

% uu = dsolve(['2u=x','u(0)=0','u(L)=0','x')

xx = linspace(0,L,201);
u_exact = [xx; subs(u,x,xx)];
% plot(xx,u_exact(2,:))

u_exact_L2 = [];