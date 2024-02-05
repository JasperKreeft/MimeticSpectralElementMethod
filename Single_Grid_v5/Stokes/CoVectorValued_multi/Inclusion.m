function I = Inclusion(Mesh)

global N

xig = Gnodes(N);
[~,e_w] = MimeticpolyVal([-1 xig 1],N,1);

Ixx = kron(speye(N),e_w');

Iyy = zeros(N*(N+2),N*N);
for i=1:N
    ind = (1:N)+(i-1)*N;
    Iyy(:,ind) = kron(speye(N),e_w(i,:)');
end
Iyy = sparse(Iyy);

Iyx = spalloc((N+1)^2,N^2,0);
Ixy = Iyx;

I = [ Ixx ; Iyx ; Iyy ; Ixy ];