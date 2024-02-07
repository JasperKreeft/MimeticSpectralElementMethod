function [x,w]=EGnodes(N)

[x,w] = Gnodes(N);
x = [ -1 x 1 ];
w = [ 0 w 0 ];