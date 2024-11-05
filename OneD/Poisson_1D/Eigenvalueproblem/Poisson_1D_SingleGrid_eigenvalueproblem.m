clear
close all
clc

global N w e

NrCellRange = 2:2:20;

exact = ((1:NrCellRange(end)).^2)';

error = zeros(10,length(NrCellRange));

Jac = pi/2;

k = 0;
for N=NrCellRange
    disp("N = "+N)
    
    k=k+1;

    [xgl,w] = GLLnodes(N);
    % ygl = (1+xgl)*Jac;

    [~,e]=MimeticpolyVal(xgl,N,1);

    NG = topology1D(N);

    M0 = innerproduct_1D(0,Jac);

    M1 = innerproduct_1D(1,Jac);

    L = NG'*M1*NG;

    E = sort(eig(full(L),full(M0)));
    E(abs(E)<.2)=[];
    
    nr = min(length(E),10);
    error(1:nr,k) = abs(E(1:nr)-exact(1:nr));

end

plot_convergence_eigenvalues("N",NrCellRange,error,exact)
