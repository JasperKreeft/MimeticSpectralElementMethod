function D = divergence(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% D = divergence(N)     D:  divergence matrix on primary grid             %
%                                                                         %
% Written by Jasper Kreeft - 2010                                         %
% Contact: j.j.kreeft@tudelft.nl                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct the divergence matrix for the internal points.

D_h = kron(speye(N),spdiags([-ones(N,1) ones(N,1)],0:1,N,N+1));

D_v = zeros(N*N,N*(N+1));
for i = 1:N
    for j = 1:N
        ind1 = i + (j-1)*N;
        ind2 = j+(i-1)*(N+1);
        D_v(ind1,ind2)   = -1;
        D_v(ind1,ind2+1) =  1;
    end
end
D_v = sparse(D_v);

D = [ D_h D_v ];
