clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

mmax = 1;

for m=1:mmax

Z = 25;%1:25;
errorL1 = zeros(size(Z));
errorL2 = zeros(size(Z));

for N=Z
disp(['N = ' num2str(N)])

xp = linspace(-1,1,N+1);
yp = linspace(-1,1,N+1);

dxp = diff(xp);
dyp = diff(yp);

xd = linspace((xp(1)+xp(2))/2,(xp(N)+xp(N+1))/2,N);
yd = linspace((yp(1)+yp(2))/2,(yp(N)+yp(N+1))/2,N);

xd_ex = [-1 xd 1];
yd_ex = [-1 yd 1];

dxd_ex = diff(xd_ex); dyd_ex = diff(yd_ex);

F = zeros(N*N,1);
for j=1:N
    for i=1:N
        k = i+(j-1)*N;
        F(k) = -2*(cos(m*pi*xp(i+1))-cos(m*pi*xp(i)))*(cos(m*pi*yp(j+1))-cos(m*pi*yp(j)));
    end
end

[Dp,Gd] = topology_old(N);

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
% None

phi_in = A\F;

% Exact %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_exact = zeros((N+2)*(N+2),1);
for j=1:N+2
    for i=1:N+2
        k = i+(j-1)*(N+2);
        phi_exact(k) = sin(m*pi*xd_ex(i))*sin(m*pi*yd_ex(j));
    end
end

PHI_EXACT = reshape(phi_exact,N+2,N+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% surf(xd_ex,yd_ex,PHI)
contour(xd_ex,yd_ex,PHI)
colorbar; set(gca,'clim',[-1 1])
axis('square')
hold off
pause

postprocessen_v1

c(N) = cond(A);

end

% error plots

end

toc