function [D,G] = topology1D(N)

D = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);

G = -D';

