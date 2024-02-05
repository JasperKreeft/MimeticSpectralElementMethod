function W20 = wedgeproduct_2P0Dform(e_w)

global N wG

I2_GLLG_eta  = kron(speye(N),e_w(:,2:N+1));
I2_GLLG_xi   = kron(e_w(:,2:N+1),speye(N));

W2 = spdiags(kron(wG,wG)',0,N^2,N^2);

W20 = full(I2_GLLG_xi*I2_GLLG_eta*W2);