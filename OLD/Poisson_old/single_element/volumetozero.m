function [C]=volumetozero(N,xiG)
% Transformation is not yet in this function

[h_w,dhdxw] = LagrangeVal(xiG,N,1);
e_w = EdgeVal(dhdxw);

% C = IGLLGT

I_GLLG_xi  = kron(speye(N),e_w');
I_GLLG_eta = kron(e_w',speye(N));

C = (I_GLLG_eta*I_GLLG_xi)';