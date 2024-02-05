function W11 = wedgeproduct_1P1Dform(ew,e_w)

global N wGL wG


I2_GLLG = kron(speye(2),kron(e_w(:,2:N+1),speye(N+1)));

W2 = spdiags([kron(wG,wGL) kron(wG,wGL)]',0,2*N*(N+1),2*N*(N+1));

I2_GGLL = kron(speye(2),kron(speye(N),ew));

W11 = I2_GLLG*W2*I2_GGLL';