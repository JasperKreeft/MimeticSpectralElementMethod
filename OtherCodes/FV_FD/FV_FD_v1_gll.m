clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic

mmax = 1;

for m=1:mmax

Z = 40;
errorL1 = zeros(size(Z));
errorL2 = zeros(size(Z));

for N=Z
disp(['N = ' num2str(N)])

xp = GLLnodes(N);
yp = GLLnodes(N);

dxp = diff(xp);
dyp = diff(yp);

xd = Gnodes(N);
yd = Gnodes(N);

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

% phi_in = A\F;

% postprocessen_v1
% 
% c(N) = cond(A);

end

% error plots



end

toc