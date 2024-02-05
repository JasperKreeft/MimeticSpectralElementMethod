function M0 = innerproduct_zeroforms(J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Inner-product matrix for 0-forms (or (n-2)-forms with n=2)              %
% for a single spectral element                                           %
%                                                                         %
% written by Jasper Kreeft                                                %
% date: 28-10-2010                                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N w

if size(J) == [ N+1 N+1 ]                                      %#ok<*BDSCA>
    J = reshape(J,1,(N+1)^2)';
elseif size(J) == [ 1 (N+1)^2 ]
    J = J';
end

Jacobian = spdiags(J,0,(N+1)^2,(N+1)^2);
W = spdiags(kron(w,w)',0,(N+1)^2,(N+1)^2);

M0 = Jacobian.*W;