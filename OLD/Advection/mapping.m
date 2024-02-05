clear all
clf
clc


n = 5;

% x_i = -1 + 2.*sort(rand(1,n));
x_i = -1.5 + 3.*sort(rand(1,n));
% x_i = Gnodes(n);

m = 100;

xi_full = linspace(-1.5,1.5,m);

for i=1:n
    x_p     = x_i;
    x_p(i)  = [];
    c_i  = prod(x_i(i)-x_p);
    for j=1:m
        hb_i_full(i,j) = prod(xi_full(j)-x_p)/c_i;
    end
end

subplot(1,2,1)
plot(xi_full,hb_i_full,'b')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 4;

xi = GLLnodes(m);

for i=1:n
    x_p     = x_i;
    x_p(i)  = [];
    c_i  = prod(x_i(i)-x_p);
    for j=1:m+1
        hb_i(i,j) = prod(xi(j)-x_p)/c_i;
    end
end

hold on
plot(xi,hb_i,'sk')
plot(x_i,ones(1,n),'dg')
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_j = LagrangeVal(xi,m,1);
plot(xi,diag(h_j),'^r')
h_j_full = LagrangeVal(xi_full,m,1);
plot(xi_full,h_j_full,'r')


ub_i = rand(1,n);

u1 = ub_i*hb_i_full;

subplot(1,2,2)
plot(xi_full,u1)
hold on
plot(x_i,ub_i,'og')

u_j = ub_i*hb_i;

u2 = u_j*h_j_full;

plot(xi_full,u2,'--r')
plot(xi,u_j,'xk')