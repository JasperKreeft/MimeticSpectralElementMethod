function H = ElementHodgeMatrix(J,str)
% Transformation is not yet in this function

global N N2
global e_w

I_GLLG_xi  = kron(speye(N),e_w(:,2:N+1)');
I_GLLG_eta = kron(e_w(:,2:N+1)',speye(N));

if strcmp(str,'02')
    Jinv = spdiags(1./J,0,N2,N2);
    H = full(Jinv*I_GLLG_eta*I_GLLG_xi);
elseif strcmp(str,'20')
    J = spdiags(J,0,N2,N2);
    H = full(J*I_GLLG_eta*I_GLLG_xi);
end