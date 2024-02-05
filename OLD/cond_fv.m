clear all

Z = 1:40;

for N=Z

A = 4*diag(ones(1,N*N))+diag(-ones(1,N*N-1),+1)+diag(-ones(1,N*N-1),-1)+...
    diag(-ones(1,N*(N-1)),+N)+diag(-ones(1,N*(N-1)),-N);

c_fv(N) = cond(A);

end