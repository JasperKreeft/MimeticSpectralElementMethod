clear all
% close all
clc

mmax = 1;

nn = 50;
[xp,w] = Gnodes(nn); yp=xp;
Xp = xp'*ones(1,nn); Yp = Xp';
ww = w'*w;

color = 'brgcmky';
leg = zeros(mmax,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m=1:mmax

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exact_v2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z = 3;
errorL1 = zeros(size(Z));
errorL2 = zeros(size(Z));
c = zeros(size(Z));
for N=Z
disp(['N = ' num2str(N)])

[xGLL,wGLL] = GLLnodes(N);
yGLL = xGLL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h   ,dhdxi  ] = LagrangeVal(xGLL,N,1);

e    = LineVal(dhdxi);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A11 = zeros(N*(N+1)); B11 = zeros(N*(N+1));
for l=1:N
    for p=1:N+1
        pl = p+(l-1)*(N+1);
        % A11
        for q=1:N+1
            for i=1:N
                iq  = i+(q-1)*N;
                A11(pl,iq) = e(i,p)*e(l,q)*wGLL(p)*wGLL(q);
            end
        end
        % B11
        for j=1:N
            pl = p+(l-1)*(N+1);
            pj = p+(j-1)*(N+1);
            B11(pl,pj) = wGLL(p)*sum(wGLL.*e(j,:).*e(l,:));
        end
    end
end

A22 = zeros(N*(N+1)); B22 = zeros(N*(N+1));
for q=1:N+1
    for k=1:N
        kq = k+(q-1)*N;
        % A22
        for j=1:N
            for p=1:N+1
                ij  = p+(j-1)*(N+1);
                A22(kq,pj) = e(k,p)*e(j,q)*wGLL(p)*wGLL(q);
            end
        end
        % B22
        for i=1:N
            qk = k+(q-1)*N;
            qi = i+(q-1)*N;
            B22(qk,qi) = wGLL(q)*sum(wGLL.*e(i,:).*e(k,:));
        end
    end
end

H = [inv(B11)*A11 zeros(N*(N+1)); zeros(N*(N+1)) inv(B22)*A22];

[D,G] = topology_singlegrid(N);

C = D*H*G;

disp('boundary elements are removed because of zero boundary condition')

% top row
C = C(:,1:N*(N+1)-1);

for i=N-2:-1:1
    C(:,(i+1)*(N+1)+1) = [];
    C(:,(i+1)*(N+1)  ) = [];
end

% lower row
C = C(:,N+3:N*(N-1)+3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = zeros(N*N,1);
for j=1:N
    for i=1:N
        k = i+(j-1)*N;
        F(k) = -2*(cos(m*pi*xGLL(i+1))-cos(m*pi*xGLL(i)))*(cos(m*pi*yGLL(j+1))-cos(m*pi*yGLL(j)));
%         F(k) = -1/2*(sin(m*pi*xGLL(i+1))-sin(m*pi*xGLL(i)))*(sin(m*pi*yGLL(j+1))-sin(m*pi*yGLL(j)))+...
%                -m*pi/4*( (yGLL(j+1)-yGLL(j))*(sin(m*pi*xGLL(i+1))-sin(m*pi*xGLL(i)))+(xGLL(i+1)-xGLL(i))*(sin(m*pi*yGLL(j+1))-sin(m*pi*yGLL(j))) );
%         F(k) = 2/3*(xGLL(i+1)^3-xGLL(i)^3)*(yGLL(j+1)-yGLL(j))+2/3*(xGLL(i+1)-xGLL(i))*(yGLL(j+1)^3-yGLL(j)^3)-4*(xGLL(i+1)-xGLL(i))*(yGLL(j+1)-yGLL(j));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HIERRRR C is not a square
phi_in = C\F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

postprocessen_v2_GLL

% c(N) = cond(C);
end

errorplot

% figure
% semilogy(Z,c(Z),'-')

end