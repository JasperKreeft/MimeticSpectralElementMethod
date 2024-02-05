clear all
close all
clc

N = 2;

xp = linspace(0,1,N+1);
yp = linspace(0,1,N+1);

dxp = diff(xp);
dyp = diff(yp);

xd = linspace(xp(2)/2,1-xp(2)/2,N);
yd = linspace(yp(2)/2,1-yp(2)/2,N);

xd_ex = [0 linspace(xp(2)/2,1-xp(2)/2,N) 1];
yd_ex = [0 linspace(yp(2)/2,1-yp(2)/2,N) 1];

F = zeros(N*N,1);

for j=1:N
    for i=1:N
        k = i+(j-1)*N;
        F(k) = 2*(cos(xp(i))-cos(xp(i)+dxp(i)))*(cos(yp(j))-cos(yp(j)+dyp(j)));
    end
end
phi_exact = zeros((N+2)*(N+2),1);
for j=1:N+2
    for i=1:N+2
        k = i+(j-1)*(N+2);
        phi_exact(k) = sin(xd_ex(i))*sin(yd_ex(j));
    end
end
phi_corner = phi_exact([1 N+2 1+(N+1)*(N+2) (N+2)*(N+2)]);
phi_exact_2 = phi_exact;
phi_exact_2([1 N+2 1+(N+1)*(N+2) (N+2)*(N+2)]) = [];

% u_bc = zeros(N*(N+1)+(N+1)*N,1);
% for j=1:N
%     k = 1+(j-1)*(N+1);
%     l = 1+j*(N+2);
%     u_bc(k)   =  phi_exact(l);
%     u_bc(k+N) = -phi_exact(l+N+1);
% end
% u_bc((N+1)*N+(1:N)) = phi_exact(2:N+1);
% u_bc(N*(N+1)+(N+1)*N-(N-1:-1:0)) = -phi_exact((N+2)*(N+2)-(N:-1:1));

q_bc = zeros(N*(N+1)+(N+1)*N,1);
for j=1:N
    k = 1+(j-1)*(N+1);
    q_bc(k)   = cos(yp(j))-cos(yp(j+1));
    q_bc(k+N) = cos(1)*(cos(yp(j+1))-cos(yp(j)));
end
q_bc((N+1)*N+(1:N)) = cos(xp(1:N))-cos(xp(2:N+1));
q_bc(N*(N+1)+(N+1)*N-(N-1:-1:0)) = cos(1)*(cos(xp(2:N+1))-cos(xp(1:N)));


Dp = topology(N);

Gd = Dp';

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
% H = eye(length(Gd));

phi_in = (Dp*H*Gd)\(F-Dp*q_bc);

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

PHI_EXACT = zeros(N+2);
for i=1:N+2
    for j=1:N+2
        k = i+(j-1)*(N+2);
        PHI_EXACT(i,j) = phi_exact(k);
    end
end

% figure
% hold on
% for i=1:N
%     for j=1:N
%         k = i+(j-1)*N;
%         surf([xp(i) xp(i+1)],[yp(j) yp(j+1)],[F(k) F(k); F(k) F(k)])
%     end
% end

% figure
surf(xd_ex,yd_ex,PHI)
% surf(xd_ex,yd_ex,abs(PHI-PHI_EXACT)')

%% Exact

% xx = linspace(0,1,100);
% yy = linspace(0,1,100);
% 
% phi_ex = zeros(100);
% for i=1:100
%     for j=1:100
%         phi_ex(i,j) = sin(xx(i))*sin(yy(j));
%     end
% end
% 
% hold on
% 
% surf(xx,yy,phi_ex')