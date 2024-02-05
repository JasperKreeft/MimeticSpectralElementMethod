function M2 = innerproduct_twoforms(e,J)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% Inner-product matrix for 2-forms (or n-forms with n=2)                  %
% for a single spectral element                                           %
%                                                                         %
% written by Jasper Kreeft                                                %
% date: 28-10-2010                                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global N w

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1
    J = ones(N+1);
end


if size(J) == [ (N+1)^2 1 ]                                    %#ok<*BDSCA>
    J = reshape(J,N+1,N+1);
elseif size(J) == [ 1 (N+1)^2 ]
    J = reshape(J,N+1,N+1);
end

% inner-product matrix for 2-forms ( or n-forms with n=2 )
M2 = zeros(N^2);
for k=1:N
    for l=1:N
        kl = k+(l-1)*N;
        for i=1:N
            for j=1:N
                ij = i+(j-1)*N;
                M2(kl,ij) = M2(kl,ij) + (w.*e(i,:).*e(k,:))*J*(w.*e(j,:).*e(l,:))';
            end
        end
    end
end