function NG = normalgrad(N)

NG_h = zeros(N*(N+1),(N+1)^2);
for i=1:N
    unit = sparse([ zeros(1,i-1) -1 1 zeros(1,N-i) ]);
    NG_h((i-1)*(N+1)+(1:N+1),:) = kron(speye(N+1),unit);
end
    
NG_v = spdiags([ones(N*(N+1),1) -ones(N*(N+1),1)],[0 N+1],N*(N+1),(N+1)^2);

NG = [ NG_v ; NG_h ];
