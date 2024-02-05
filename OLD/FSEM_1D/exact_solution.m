function [rhs,u_exact,u_exact_L2,bc,bc_l,bc_r] = exact_solution(Diff_string,L,xg)

% clear all
% close all
% clc
% 
% % exact
% % 
% L = 1; xg = linspace(0,L,3);

syms x

% Diff_string = input('Enter differential equation (default D2u+u):   L[u(x)] = ','s');
% if isempty(Diff_string)
%     Diff_string = '-D2u+u';
% end

f_string = input('Enter forcing function: f(x) = ','s');
f = sym(f_string);

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

if bc == 1
    u_l = input('enter value of left bc: u_l = ');
    disp('  ')
    u_r = input('enter value of right bc: u_r = ');
    disp('  ')
    
    bc_l = u_l;
    bc_r = u_r;
    
    bc_l_string = ['u(',num2str(xg(1)),')=',num2str(u_l)];
    bc_r_string = ['u(',num2str(xg(end)),')=',num2str(u_r)];
    
elseif bc == 2
    u_l = input('enter value of left bc: u_l = ');
    disp('  ')
    du_r = input('enter value of right bc: du_r = ');
    disp('  ')
    
    bc_l = u_l;
    bc_r = du_r;

    bc_l_string = ['u(',num2str(xg(1)),')=',num2str(u_l)];
    bc_r_string = ['Du(',num2str(xg(end)),')=',num2str(du_r)];
    
elseif bc == 3
    du_l = input('enter value of left bc: du_l = ');
    disp('  ')
    u_r = input('enter value of right bc: u_r = ');
    disp('  ')

    bc_l = du_l;
    bc_r = u_r;

    bc_l_string = ['Du(',num2str(xg(1)),')=',num2str(Du_l)];
    bc_r_string = ['u(',num2str(xg(end)),')=',num2str(u_r)];

end

u = dsolve([Diff_string,'=',f_string],bc_l_string,bc_r_string,'x');

xx = linspace(0,L,201);
u_exact = [xx; subs(u,x,xx)];
plot(xx,u_exact(2,:))

u_exact_L2 = [];