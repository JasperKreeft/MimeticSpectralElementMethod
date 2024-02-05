clear all
clf %ose all
clc

N = 2;
M = N;

xp = linspace(0,1,N+1);
yp = linspace(0,1,N+1);

xd = [0 1/4 3/4 1];
yd = [0 1/4 3/4 1];

dxd = diff(xd);
dyd = diff(yd);

F = zeros((N+1)*(N+1),1);
phi_exact = zeros((N+1)*(N+1),1);
for i=1:N+1
    for j=1:N+1
        k = i+(j-1)*(N+1);
        F(k) = 2*(cos(xd(i))-cos(xd(i)+dxd(i)))*(cos(yd(j))-cos(yd(j)+dyd(j)));%2*sin(xp(i))*sin(yp(j));%
        phi_exact(k) = sin(xp(i))*sin(yp(j));
    end
end
phi_bc([1:4 6:9],1) = phi_exact([1:4 6:9]);


hold on
for i=1:N+1
    for j=1:N+1
        k = i+(j-1)*(N+1);
        surf([xd(i) xd(i+1)],[yd(j) yd(j+1)],[F(k) F(k); F(k) F(k)])
    end
end


[Gp,Cp]=topology_old(N,M);
  
Dd = Gp';  

Grad = Gp*phi_exact;

Curl = Cp*Gp*phi_exact;


H = eye(length(Gp));

A = Dd*H*Gp;

AA = A(5,5);

FF = F(5)-A(5,[1:4 6:9])*phi_bc([1:4 6:9]);

phi_in    = zeros((N+1)*(N+1),1);
phi_in(5) = AA\FF;

phi = phi_in+phi_bc;

PHI = zeros(N+1);
for i=1:N+1
    for j=1:N+1
        k = i+(j-1)*(N+1);
        PHI(i,j) = phi(k);
    end
end

figure
surf(xp,yp,PHI')

%% Exact

xx = linspace(0,1,100);
yy = linspace(0,1,100);

phi_ex = zeros(100);
for i=1:100
    for j=1:100
        phi_ex(i,j) = sin(xx(i))*sin(yy(j));
    end
end

hold on

surf(xx,yy,phi_ex')