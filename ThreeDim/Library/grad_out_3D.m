function G = grad_out_3D(N)

unit = spdiags([-ones(N,1) ones(N,1)],[0 1],N,N+1);

G_u = kron(speye((N+1)^2),unit);

unit = kron(speye(N+1),[-1;zeros(N-1,1)]);
for i=2:N
    unit_part = kron(speye(N+1),[zeros(i-2,1) ; 1 ; -1 ; zeros(N-i,1)]);
    unit = [ unit unit_part ];
end
unit_part = kron(speye(N+1),[zeros(N-1,1) ; 1]);
unit = [unit unit_part];

G_v = kron(speye(N+1),unit);


G_w = kron(speye((N+1)^2),[-1;zeros(N-1,1)]);
for i=2:N
    unit_part = kron(speye((N+1)^2),[zeros(i-2,1) ; 1 ; -1 ; zeros(N-i,1)]);
    G_w = [ G_w unit_part ];
end
unit_part = kron(speye((N+1)^2),[zeros(N-1,1) ; 1]);
G_w = [G_w unit_part];

G = [ G_w
      G_v
      G_u ];