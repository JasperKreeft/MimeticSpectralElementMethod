clear all
close all
clc

global N w e

NrElementRange = 2.^(2:1:10);
NrCellRange = 2;%:20;

error = zeros(10,length(NrCellRange)*length(NrElementRange)); er = 0;

n=0;
for M=NrElementRange
    
Jac = pi/(2*M);

for N=NrCellRange
disp(['M = ' num2str(M) ', N = ' num2str(N)])

[xi,w] = GLLnodes(N);

[h,e] = MimeticpolyVal(xi,N,1);

G = full(topology1D(N*M));

A1 = zeros(N);
for i=1:N
    for j=1:N
        A1(i,j) = sum(w.*e(i,:).*e(j,:)/Jac);
    end
end

A = kron(eye(M),A1);

L = G'*A*G;

% L = L*((M/Jac)^2); % Jacobiaan

A0 = zeros(1,N*M+1);
for m=1:M
    ind = N*(m-1)+(1:N+1);
    A0(ind) = A0(ind) + w*Jac;
end
A0 = diag(A0);


E = sort(eig(full(L),full(A0)));

E(abs(E)<.2)=[];

exact = ((1:length(E)).^2)';
nr = min(length(E),10);
er = er+1;
error(1:nr,er) = abs(E(1:nr)-exact(1:nr));


end
end

%%
if length(NrCellRange)>=9
    plot_convergence_eigenvalues("N",NrCellRange,error,exact)
elseif length(NrElementRange)>=9
    plot_convergence_eigenvalues("H",NrElementRange,error,exact)
end