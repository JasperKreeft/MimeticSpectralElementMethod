clear all
close all
clc

m = 1;

for N=1:2:20

xp = linspace(-1,1,N+1);
yp = linspace(-1,1,N+1);

dxp = diff(xp);
dyp = diff(yp);

xd = linspace((xp(1)+xp(2))/2,(xp(N)+xp(N+1))/2,N);
yd = linspace((yp(1)+yp(2))/2,(yp(N)+yp(N+1))/2,N);

xd_ex = [-1 xd 1];
yd_ex = [-1 yd 1];

F = zeros(N*N,1);
for j=1:N
    for i=1:N
        k = i+(j-1)*N;
        F(k) = -2*(cos(m*pi*xp(i+1))-cos(m*pi*xp(i)))*(cos(m*pi*yp(j+1))-cos(m*pi*yp(j)));
    end
end
phi_exact = zeros((N+2)*(N+2),1);
for j=1:N+2
    for i=1:N+2
        k = i+(j-1)*(N+2);
        phi_exact(k) = sin(m*pi*xd_ex(i))*sin(m*pi*yd_ex(j));
    end
end

PHI_EXACT = reshape(phi_exact,N+2,N+2);

[Dp,Gd] = topology(N);

diag_dx = 2./([dxp 0]+[0 dxp]);
diag_dy = 2./([dyp 0]+[0 dyp]);

diag_v = zeros(1,N*(N+1));
for j=1:N+1
    for i=1:N
        k=i+(j-1)*N;
        diag_v(k) = dxp(i)*diag_dy(j);
    end
end

H = [ kron(diag(dyp),diag(diag_dx)) zeros(N*(N+1));
      zeros((N+1)*N)                diag(diag_v)  ];


A = Dp*H*Gd;


% Additing boundary conditions
for i=1:N
    for j=1:N
        k = i+(j-1)*N;
        if i==1
            F(k,1) = F(k,1) + 2*sin(m*pi*xd_ex(1))*sin(m*pi*yd_ex(j+1));
        elseif i==N
            F(k,1) = F(k,1) + 2*sin(m*pi*xd_ex(N+2))*sin(m*pi*yd_ex(j+1));
        end
        if j==1
            F(k,1) = F(k,1) + 2*sin(m*pi*xd_ex(i+1))*sin(m*pi*yd_ex(1));
        elseif j==N
            F(k,1) = F(k,1) + 2*sin(m*pi*xd_ex(i+1))*sin(m*pi*yd_ex(N+2));
        end
    end
end

phi_in = A\F;

PHI = zeros(N+2); l=0;
for i=1:N+2
    for j=1:N+2
        k = i+(j-1)*(N+2);
        if i==1 || i==N+2
            PHI(i,j) = phi_exact(k);
        elseif j==1 || j==N+2
            PHI(i,j) = phi_exact(k);
        else            
        l=l+1;
        PHI(i,j) = phi_in(l);
        end
    end
end

%% Exact

xx = linspace(-1,1,100);
yy = linspace(0-1,1,100);

phi_ex = zeros(100);
for i=1:100
    for j=1:100
        phi_ex(i,j) = sin(m*pi*xx(i))*sin(m*pi*yy(j));
    end
end

% surf(xx,yy,phi_ex')
contour(xx,yy,phi_ex','k')

hold on

% surf(xd_ex,yd_ex,PHI)
contour(xd_ex,yd_ex,PHI)
colorbar; set(gca,'clim',[-1 1])
axis('square')
hold off
pause

end