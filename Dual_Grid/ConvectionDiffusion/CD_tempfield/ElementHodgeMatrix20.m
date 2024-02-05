function H20 = ElementHodgeMatrix20(J)
% Transformation is not yet in this function

global N N2
global e_w

Jinv = spdiags(1./J,0,N2,N2);

I_GLLG_xi  = kron(speye(N),e_w(:,2:N+1)');
I_GLLG_eta = kron(e_w(:,2:N+1)',speye(N));
% keyboard
H02 = full(Jinv*I_GLLG_eta*I_GLLG_xi);