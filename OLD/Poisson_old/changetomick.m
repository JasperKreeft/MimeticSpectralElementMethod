clear all
% close all
clc

mmax = 3;

nn = 100;
[xp,w] = Gnodes(nn); yp=xp;
Xp = xp'*ones(1,nn); Yp = Xp';
Wgg = w'*w;

color = 'brgcmky';
leg = zeros(mmax,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m=1:mmax

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exact_v2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NrCellRange = 10;
errorL1 = zeros(size(NrCellRange));
errorL2 = zeros(size(NrCellRange));
c = zeros(size(NrCellRange));
for N=NrCellRange
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

e    = EdgeVal(dhdxi);
ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);
ew_w = EdgeVal(dhwdxiw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A11 = zeros(N*(N+1)); B11 = zeros(N*(N+1));
for l=1:N
    for k=1:N+1
        kl = k+(l-1)*(N+1);
        % A11 & B11
        for j=1:N
            for i=1:N+1
                ij  = i+(j-1)*(N+1);
                A11(kl,ij) = ew(i,k)*e_w(l,j+1)*wGLL(k)*wG(j);
            end
        end
    end
end
for p=1:N+1
    for l=1:N
        for j=1:N
            pl = p+(l-1)*(N+1);
            pj = p+(j-1)*(N+1);
            B11(pl,pj) = wGLL(p)*sum(wG.*e_w(j,2:N+1).*e_w(l,2:N+1));
        end
    end
end

A22 = zeros(N*(N+1)); B22 = zeros(N*(N+1));
for l=1:N+1
    for k=1:N
        kl = k+(l-1)*N;
        % A22 & B22
        for j=1:N+1
            for i=1:N
                ij  = i+(j-1)*N;
                A22(kl,ij) = e_w(k,i+1)*ew(j,l)*wG(i)*wGLL(l);
            end
        end
    end
end
for q=1:N+1
    for k=1:N
        for i=1:N
            qk = k+(q-1)*N;
            qi = i+(q-1)*N;
            B22(qk,qi) = wGLL(q)*sum(wG.*e_w(i,2:N+1).*e_w(k,2:N+1));
        end
    end
end


H = [inv(B11)*A11 zeros(N*(N+1)); zeros(N*(N+1)) inv(B22)*A22];

[Dp,Gd] = topology(N);

C = Dp*H*Gd;

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

phi_in = C\F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% postprocessen_v2_G
postprocessen_v2_EG

c(N) = cond(C);
end

errorplot

figure
semilogy(NrCellRange,c(NrCellRange),'-')

end