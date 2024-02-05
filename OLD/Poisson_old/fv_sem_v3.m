clear all
close all
clc

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=0.3;

mmax = 1;

nn = 100;
[xp,Gw] = Gnodes(nn); yp = xp;
Xp = xp'*ones(1,nn); Yp = Xp';

color = 'brgcmky';
leg = zeros(mmax,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for m=1:mmax

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% exact_v3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NrCellRange = 2:2:20;
errorL1 = zeros(size(NrCellRange));
errorL2 = zeros(size(NrCellRange));
for N=NrCellRange
disp(['N = ' num2str(N)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

curvedgrid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h   ,dhdxi  ] = LagrangeVal(xiGLL,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xiGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);
[hw_w,dhwdxiw] = LagrangeVal(xiEG,N,3);

e    = EdgeVal(dhdxi);
ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);
ew_w = EdgeVal(dhwdxiw);

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
                A11(kl,ij) = S12J(k,j)*ew(i,k)*e_w(l,j+1)*wGLL(k)*wG(j);
            end
        end
        % A12
        for j=1:N+1
            for i=1:N
                ij = i+(j-1)*N;
                A12(kl,ij) = h_w(k,i+1)*wG(i)*sum(S22J(i,:).*ew(j,:).*e(l,:).*wGLL);
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
                A21(kl,ij) = -1* h_w(l,j+1)*wG(j)*sum(S11J(:,j)'.*ew(i,:).*e(k,:).*wGLL);
            end
        end
        % A22
        for j=1:N+1
            for i=1:N
                ij  = i+(j-1)*N;
                A22(kl,ij) = -1* S21J(i,l)*e_w(k,i+1)*ew(j,l)*wG(i)*wGLL(l);
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

B11i = inv(B11); B22i = inv(B22);
H = [B11i*A11 B11i*A12; B22i*A21 B22i*A22];

[Dp,Gd] = topology(N);

C = Dp*H*Gd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[FF] = force(N,xiGLL,etaGLL,c);
F = reshape(FF,N*N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi_in = C\F;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% postprocessen_v3_EG
fout
end

errorplot

end

toc
% semilogy(NrCellRange,errorL2(NrCellRange))
% hold on
% semilogy(NrCellRange,errorL2_interp(NrCellRange),'--r')