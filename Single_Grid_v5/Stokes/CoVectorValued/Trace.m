function Tr = Trace()

global N xi

[~,ew] = MimeticpolyVal(xi,N,3);

Tr = kron(sparse([1 0 ; 0 -1]),kron(speye(N),ew'));