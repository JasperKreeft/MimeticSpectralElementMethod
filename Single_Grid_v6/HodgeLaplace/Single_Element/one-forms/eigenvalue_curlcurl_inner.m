%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Single-element curl-curl problem on inner-oriented grid                 %
% WERKT ALLEEN VOOR UNIFORME GRIDS                                        %
% written by Jasper Kreeft (2011)                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

path(path,'C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Library_SingleGrid')

N = 12;

C = curl_in(N);

[x,w] = GLLnodes(N);
[h,e] = MimeticpolyVal(x,N,1);

M2 = zeros(N*N);
for i=1:N
    for j=1:N
        ij = i+(j-1)*N;
        for k=1:N
            for l=1:N
                kl = k+(l-1)*N;
                M2(kl,ij) = sum(w.*e(i,:).*e(k,:))*sum(w.*e(j,:).*e(l,:));
            end
        end
    end
end

M1xi = zeros(N*(N+1));
for k=1:N
    for q=1:N+1
        kq = k+(q-1)*N;
        for i=1:N
            iq = i+(q-1)*N;
            M1xi(kq,iq) = w(q)*sum(w.*e(i,:).*e(k,:));
        end
    end
end
M1eta = zeros(N*(N+1));
for l=1:N
    for p=1:N+1
        pl = l+(p-1)*N;
        for j=1:N
            pj = j+(p-1)*N;
            M1eta(pl,pj) = w(p)*sum(w.*e(j,:).*e(l,:));
        end
    end
end
M1 = [ M1xi       zeros(N*(N+1))
       zeros(N*(N+1)) M1eta      ];

L = C'*M2*C;

E = sort(eig(full(L),full(M1)));

E(abs(E)<.2)=[];

E = E*2/E(1);

E(1:min(20,length(E)))

rmpath('C:\Users\Zolder\Documents\MATLAB\MSEM\Single_Grid_v6\Library_SingleGrid')