clear all
close all
clc

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=0.2;

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

NrCellRange = 5:10;
errorL1 = zeros(size(NrCellRange));
errorL2 = zeros(size(NrCellRange));
for N=NrCellRange
disp(['N = ' num2str(N)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

curvedgrid2;

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
A21 = zeros(N*(N+1));
A22 = zeros(N*(N+1));

% A11
for p=1:N+1
    for l=1:N
        for j=1:N
            pl = p+(l-1)*(N+1);
            pj = p+(j-1)*(N+1);
            A11(pl,pj) = wGLL(p)*sum(wG.*S11J(p,:).*e_w(j,2:N+1).*e_w(l,2:N+1));
        end
    end
end
% A12
for l=1:N
    for k=1:N+1
        kl = k+(l-1)*(N+1);
        for j=1:N+1
            for i=1:N
                ij = i+(j-1)*N;
                A12(kl,ij) = wGLL(k)*wGLL(j)*S12J(k,j)*e(i,k)*e(l,j);
            end
        end
    end
end

% A21
for l=1:N+1
    for k=1:N
        kl = k+(l-1)*N;
        for j=1:N
            for i=1:N+1
                ij  = i+(j-1)*(N+1);
                A21(kl,ij) = wGLL(i)*wGLL(l)*S21J(i,l)*e(k,i)*e(j,l);
            end
        end
    end
end

% A22
for q=1:N+1
    for k=1:N
        for i=1:N
            qk = k+(q-1)*N;
            qi = i+(q-1)*N;
            A22(qk,qi) = wGLL(q)*sum(wG.*S22J(:,q)'.*e_w(i,2:N+1).*e_w(k,2:N+1));
        end
    end
end

B = zeros(N*N);
for l=1:N
    for k=1:N
        kl = k+(l-1)*N;
        for q=1:N
            for p=1:N
                pq = p+(q-1)*N;
                B(pq,kl) = wG(p)*wG(q)*e_w(k,p+1)*e_w(l,q+1);
            end
        end
    end
end

A = [ A11 A12 ; A21 A22 ];

[Dp,Gd] = topology(N);

C = B*Dp*inv(A)*Gd*B';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[FF] = force(N,xiGLL,etaGLL,c);
F = reshape(FF,N*N,1);

F = B*F;

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