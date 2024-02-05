function M1 = innerproduct_oneforms(e,J,Qinv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Inner-product matrix for 1-forms (or (n-1)-forms with n=2)              %
% for a single spectral element                                           %
%                                                                         %
% written by Jasper Kreeft                                                %
% date: 28-10-2010                                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N w

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W1 = spdiags(kron(kron(w,w),ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

if size(J) == [ N+1 N+1 ]                                      %#ok<*BDSCA>
    J = reshape(J,1,(N+1)^2);
elseif size(J) == [ (N+1)^2 1 ]
    J = J';
end

Jacobian = spdiags(kron(J,ones(1,2))',0,2*(N+1)^2,2*(N+1)^2);

I1 = spalloc(2*(N+1)^2,2*N*(N+1),2*N*(N+1)^2);
I1(1:2:2*(N+1)^2,1:N*(N+1)) = kron(e',speye(N+1));
for i=1:N
    I1(2:2:2*(N+1)^2,(N+i-1)*(N+1)+(1:N+1)) = kron(speye(N+1),e(i,:)');
end

% inner-product matrix for 1-forms ( or (n-1)-forms with n=2 )
M1 = I1'*Qinv'*(W1.*Jacobian)*Qinv*I1;
