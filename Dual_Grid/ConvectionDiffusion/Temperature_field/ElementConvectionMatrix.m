function [Conv] = ElementConvectionMatrix(V,T)

global N N2
global xiEG
global ew_w

etaEG = xiEG;

a = -1/2*(1-etaEG.^2);
b = 0;

Conv = [ kron(speye(N),(ones(N+1,1)*a.*ew_w)') spalloc(N2+2*N,N*(N+1),0)
            spalloc(N2+2*N,N*(N+1),0)          b*kron(speye(N),ew_w') ];

Conv = T*Conv;