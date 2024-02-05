function [H] = elementmatrix(N,xiGLL,xiEG,wG,wGLL,numRows,numColumns)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hw  ,dhwdxi ] = LagrangeVal(xiGLL,N,3);
[h_w ,dhdxiw ] = LagrangeVal(xiEG,N,1);

ew   = EdgeVal(dhwdxi);
e_w  = EdgeVal(dhdxiw);

J = 1/(numRows*numColumns);

Su = (1/numRows)^2;
Sv = (1/numColumns)^2;

SuJ = Su/J;
SvJ = Sv/J;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Au = zeros(N*(N+1)); Bu = zeros(N*(N+1));
for l=1:N
    for k=1:N+1
        kl = k+(l-1)*(N+1);
        % Au & Bu
        for j=1:N
            for i=1:N+1
                ij  = i+(j-1)*(N+1);
                Au(kl,ij) = ew(i,k)*e_w(l,j+1)*wGLL(k)*wG(j);
            end
        end
    end
end
for p=1:N+1
    for l=1:N
        for j=1:N
            pl = p+(l-1)*(N+1);
            pj = p+(j-1)*(N+1);
            Bu(pl,pj) = wGLL(p)*sum(wG.*e_w(j,2:N+1).*e_w(l,2:N+1));
        end
    end
end

Av = zeros(N*(N+1)); Bv = zeros(N*(N+1));
for l=1:N+1
    for k=1:N
        kl = k+(l-1)*N;
        % Av & Bv
        for j=1:N+1
            for i=1:N
                ij  = i+(j-1)*N;
                Av(kl,ij) = e_w(k,i+1)*ew(j,l)*wG(i)*wGLL(l);
            end
        end
    end
end
for q=1:N+1
    for k=1:N
        for i=1:N
            qk = k+(q-1)*N;
            qi = i+(q-1)*N;
            Bv(qk,qi) = wGLL(q)*sum(wG.*e_w(i,2:N+1).*e_w(k,2:N+1));
        end
    end
end


H = [inv(Bu)*(SuJ*Au) zeros(N*(N+1)); zeros(N*(N+1)) inv(Bv)*(SvJ*Av)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%