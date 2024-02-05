clear all
close all
clc

nn = 100;
[xp,w] = Gnodes(nn); yp=xp;
Xp = xp'*ones(1,nn); Yp = Xp';
Wgg = w'*w;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 1;
exact_v2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NrCellRange = 8;
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

k11 = ones(N+1);%xGLL'*ones(1,N+1);
k12 = ones(N+1);
k21 = ones(N+1);
k22 = ones(N+1);

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
                A11(kl,ij) = k11(k,j)*ew(i,k)*e_w(l,j+1)*wGLL(k)*wG(j);
            end
        end
        % A12
        for j=1:N+1
            for i=1:N
                ij = i+(j-1)*N;
                A12(kl,ij) = h_w(k,i+1)*wG(i)*sum(k12(k,:).*ew(j,:).*e(l,:).*wGLL);
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
                A21(kl,ij) = h_w(l,j+1)*wG(j)*sum(k21(:,l)'.*ew(i,:).*e(k,:).*wGLL);
            end
        end
        % A22
        for j=1:N+1
            for i=1:N
                ij  = i+(j-1)*N;
                A22(kl,ij) = k22(i,l)*e_w(k,i+1)*ew(j,l)*wG(i)*wGLL(l);
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

% k11 = 1; A11 = k11*A11;
% k12 = 0; A12 = k12*A12;
% k21 = 0; A21 = k21*A21;
% k22 = 1; A22 = k22*A22;



B11i = inv(B11); B22i = inv(B22);
H = [B11i*A11 B11i*A12; B22i*A21 B22i*A22];

[Dp,Gd] = topology(N);

C = Dp*H*Gd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% F = zeros(N*N,1);
% for j=1:N
%     for i=1:N
%         k = i+(j-1)*N;
%         F(k) = ((xGLL(i+1)+1)*cos(pi*xGLL(i+1))-(xGLL(i)+1)*cos(pi*xGLL(i)))*(cos(pi*yGLL(j))-cos(pi*yGLL(j+1)));
%     end
% end
F = reshape(force_int(N,xGLL,yGLL),N*N,1);

% x = linspace(-1,1,100)'*ones(1,100); y = x';
% f = -pi^2*(x+1).*sin(pi*x).*sin(pi*y)+pi*cos(pi*x).*sin(pi*y);
% surf(x,y,f); shading interp; title('f'); colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_in = C\F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% postprocessen_v2_G
postprocessen_v2_EG

end

% errorplot