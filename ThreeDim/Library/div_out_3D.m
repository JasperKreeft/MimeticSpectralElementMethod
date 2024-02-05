function D = div_out_3D(N)

D1 = kron(speye(N*N),spdiags([-ones(N,1) ones(N,1)],0:1,N,N+1));

part = zeros(N*N,N*(N+1));
for i=1:N
    part((1:N)+(i-1)*N,:) = kron(speye(N),[zeros(1,i-1) -1 1 zeros(1,N-i) ]);
end
D2 = kron(speye(N),sparse(part));

D3 = zeros(N*N*N,N*N*(N+1));
for i=1:N
    D3((1:N*N)+(i-1)*N*N,:) = kron(speye(N*N),[zeros(1,i-1) -1 1 zeros(1,N-i) ]);
end


D = [ D1 D2 D3 ];