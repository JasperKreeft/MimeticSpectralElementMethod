function In = Inclusion()

global N

xiw = EGnodes(N);

[~,e_w] = MimeticpolyVal(xiw,N,1);

xx = kron(speye(N),e_w');

yy = spalloc(N*(N+2),N^2,2*N*(N+2));
for k=1:N
    ind = (k-1)*N+(1:N);
    yy(:,ind) = kron(speye(N),e_w(k,:)');
end

xy = spalloc((N+1)^2,N^2,0);
yx = xy;

In = [ xx
       yx
       yy
       xy ];
