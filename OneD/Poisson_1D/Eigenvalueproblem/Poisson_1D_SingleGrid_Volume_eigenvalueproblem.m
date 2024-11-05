clear all
close all
clc

global N w e

NrCellRange = 2:2:20;

Jac = pi/2;

error = zeros(10,length(NrCellRange)); er = 0;

n=0;
for N=NrCellRange
disp(['N = ' num2str(N)])
    
[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

D = full(topology1D(N));

M0 = diag(w*Jac);

M1 = zeros(N);
for i=1:N
    for j=1:N
        M1(i,j) = sum(w.*e(i,:).*e(j,:)/Jac);
    end
end

L = M1*D/M0*D'*M1;

E = sort(eig(full(L),full(M1)));

E(abs(E)<.2)=[];

exact = ((1:length(E)).^2)';
nr = min(length(E),10);
er = er+1;
error(1:nr,er) = abs(E(1:nr)-exact(1:nr));

end

plot_convergence_eigenvalues("N",NrCellRange,error,exact)
