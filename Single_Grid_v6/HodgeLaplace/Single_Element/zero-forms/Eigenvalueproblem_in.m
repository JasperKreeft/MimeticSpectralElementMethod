%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Single-element grad-grad problem on inner-oriented grid                 %
% WERKT ALLEEN VOOR UNIFORME GRIDS                                        %
% written by Jasper Kreeft (2011)                                         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

N = 12;

G = grad_in(N);

[x,w] = GLLnodes(N);
[h,e] = MimeticpolyVal(x,N,1);

J = 1; %!!!!!!!!!!
M0 = spdiags(J.*kron(w,w)',0,(N+1)^2,(N+1)^2);

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

Matrix = G'*M1*G;
RHS = M0;

E = sort(eig(full(Matrix),full(RHS)));

% E(abs(E)<.2)=[];

E = E/E(2);

E(1:min(20,length(E)))