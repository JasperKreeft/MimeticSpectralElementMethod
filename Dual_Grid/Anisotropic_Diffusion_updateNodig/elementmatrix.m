function [Dp,H,Gd] = elementmatrix(N,xiGLL,xiEG,wG,wGLL,numRows,numColumns,k11,k12,k21,k22)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[h   ,dhdxi  ] = LagrangeVal(xiGLL,N,1);
[hw  ,dhwdxi ] = LagrangeVal(xiGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);
% [hw_w,dhwdxiw] = LagrangeVal(xiEG,N,3);

e    = LineVal(dhdxi);
ew   = LineVal(dhwdxi);
e_w  = LineVal(dhdxiw);
% ew_w = LineVal(dhwdxiw);

J = 1/(numRows*numColumns);

Su = (1/numRows)^2;
Sv = (1/numColumns)^2;
% Suv = Svu = 1/(numRows*numColumns);

SuJ = Su/J;
SvJ = Sv/J;
SuvJ = 1; SvuJ = 1;

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

B11i = inv(B11); B22i = inv(B22);
H = [B11i*(SuJ*A11) B11i*(SuvJ*A12); B22i*(SvuJ*A21) B22i*(SvJ*A22)]; % HIERRRR

[Dp,Gd] = topology(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%