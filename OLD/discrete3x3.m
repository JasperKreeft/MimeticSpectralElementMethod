clear all
clf %ose all
clc

N = 3;
M = N;

beta = 1;

xp = linspace(0,1,N+1);
yp = linspace(0,1,N+1);

dxp = diff(xp);
dyp = diff(yp);

xd = [0 linspace(1/(2*N),1-1/(2*N),3) 1];
yd = [0 linspace(1/(2*N),1-1/(2*N),3) 1];

dxd = diff(xd);
dyd = diff(yd);

F = zeros((N+1)*(N+1),1);
phi_exact = zeros((N+1)*(N+1),1);
for i=1:N+1
    for j=1:N+1
        k = i+(j-1)*(N+1);
        F(k) = 2*(cos(xd(i))-cos(xd(i)+dxd(i)))*(cos(yd(j))-cos(yd(j)+dyd(j)));
        phi_exact(k) = sin(xp(i))*sin(yp(j));
    end
end
phi_bc([1:5 8:9 12:16],1) = phi_exact([1:5 8:9 12:16]);


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

A = Dd*beta*Gp;

AA = A([6 7 10 11],[6 7 10 11]);

FF = F([6 7 10 11])-A([6 7 10 11],[1:5 8:9 12:16])*phi_bc([1:5 8:9 12:16]);

phi_in = zeros((N+1)*(N+1),1);
phi_in([6 7 10 11]) = AA\FF;

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

% Exact

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