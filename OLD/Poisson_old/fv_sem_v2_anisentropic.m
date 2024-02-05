clear all
close all
clc

nn = 50;
[xp,w] = Gnodes(nn); yp=xp;
Xp = xp'*ones(1,nn); Yp = Xp';
ww = w'*w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mmax = 1;
% for m=1:mmax
m = mmax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exact_v2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Z = 1:20;
errorL1 = zeros(size(Z));
errorL2 = zeros(size(Z));
c = zeros(size(Z));
for N=Z
disp(['N = ' num2str(N)])

gridchoice = 2; % 1 uniform, 2 non-uniform

if gridchoice == 1
    xGLL = linspace(-1,1,N+1);
    xG   = linspace((xGLL(1)+xGLL(2))/2,(xGLL(N)+xGLL(N+1))/2,N);
    xEG  = [-1 xG 1];
    [LaPoly,dLaPoly,wGLL] = LagrangePoly(xGLL); % primal polynomials
    [LaPolyw,dLaPolyw,wEG] = LagrangePoly(xEG); % dual polynomials (w = wiggle)
    wG = wEG(2:N+1);
elseif gridchoice == 2
    [xGLL,wGLL] = GLLnodes(N);
    [xG,wG] = Gnodes(N);
    xEG = [-1 xG 1];
end

yGLL = xGLL; yG = xG; yEG = xEG;
dxG = diff(xEG); dyG = dxG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h   ,dhdxi  ] = LagrangeVal(xGLL,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xEG,N,1);
[hw_w,dhwdxiw] = LagrangeVal(xEG,N,3);

e    = LineVal(dhdxi);
ew   = LineVal(dhwdxi);
e_w  = LineVal(dhdxiw);
ew_w = LineVal(dhwdxiw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A11 = zeros(N*(N+1));
A12 = zeros(N*(N+1));
for l=1:N
    for k=1:N+1
        kl = k+(l-1)*(N+1);
        % A11
        for j=1:N
            for i=1:N+1
                ij  = i+(j-1)*(N+1);
                A11(kl,ij) = ew(i,k)*e_w(l,j+1)*wGLL(k)*wG(j);
            end
        end
        % A12
        for j=1:N+1
            for i=1:N
                ij = i+(j-1)*N;
                A12(kl,ij) = h_w(k,i+1)*wG(i)*sum(ew(j,:).*e(l,:).*wGLL);
            end
        end
    end
end
% B11
B11 = zeros(N*(N+1));
for p=1:N+1
    for l=1:N
        for j=1:N
            pl = p+(l-1)*(N+1);
            pj = p+(j-1)*(N+1);
            B11(pl,pj) = wGLL(p)*sum(wG.*e_w(j,2:N+1).*e_w(l,2:N+1));
        end
    end
end

A21 = zeros(N*(N+1));
A22 = zeros(N*(N+1));
for l=1:N+1
    for k=1:N
        kl = k+(l-1)*N;
        % A21
        for j=1:N
            for i=1:N+1
                ij  = i+(j-1)*(N+1);
                A21(kl,ij) = h_w(l,j+1)*wG(j)*sum(ew(i,:).*e(k,:).*wGLL);
            end
        end
        % A22
        for j=1:N+1
            for i=1:N
                ij  = i+(j-1)*N;
                A22(kl,ij) = e_w(k,i+1)*ew(j,l)*wG(i)*wGLL(l);
            end
        end
    end
end
% B22
B22 = zeros(N*(N+1));
for q=1:N+1
    for k=1:N
        for i=1:N
            qk = k+(q-1)*N;
            qi = i+(q-1)*N;
            B22(qk,qi) = wGLL(q)*sum(wG.*e_w(i,2:N+1).*e_w(k,2:N+1));
        end
    end
end

k11 = 10;   A11 = k11*A11;
k12 = 3;    A12 = k12*A12;
k21 = k12;  A21 = k21*A21;
k22 = 1;    A22 = k22*A22;


B11i = inv(B11); B22i = inv(B22);
H = [B11i*A11 B11i*A12; B22i*A21 B22i*A22];

[Dp,Gd] = topology(N);

C = Dp*H*Gd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = zeros(N*N,1);
for j=1:N
    for i=1:N
        k = i+(j-1)*N;
        F(k) = -(k11+k22)*(cos(m*pi*xGLL(i+1))-cos(m*pi*xGLL(i)))*(cos(m*pi*yGLL(j+1))-cos(m*pi*yGLL(j)))+...
               +(k12+k21)*(sin(m*pi*xGLL(i+1))-sin(m*pi*xGLL(i)))*(sin(m*pi*yGLL(j+1))-sin(m*pi*yGLL(j)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_in = C\F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% postprocessen_v2_G
postprocessen_v2_EG

% c(N) = cond(C);
end

errorplot

% figure
% semilogy(Z,c(Z),'-')

% end