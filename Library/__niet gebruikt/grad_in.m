function Gp = grad_in(N)

unit = [diag(-ones(N,1)) zeros(N,1)] + ...
         [zeros(N,1) diag(ones(N,1))];
     
Gp_u = kron(speye(N+1),unit);

Gp_v = zeros(N*(N+1),(N+1)*(N+1));
for i=1:N
    for j=1:N+1
        Gp_v((j-1)*N+i,j+(i-1)*(N+1)) = -1;
        Gp_v((j-1)*N+i,j+i*(N+1))     = +1;
    end
end

Gp = sparse([Gp_u ; Gp_v]);