function [Conv] = ElementConvectionMatrix(V,T)

global N N2
global ew_w

a = V(1);
b = V(2);

IT = kron(speye(N),ew_w');

Conv = [ a*IT spalloc(N2+2*N,N*(N+1),0)
         spalloc(N2+2*N,N*(N+1),0) b*IT ];

Conv = T*Conv;