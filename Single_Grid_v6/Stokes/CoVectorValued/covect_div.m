function CD = covect_div(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% CD = covectdiv(N)     CD:  divergence matrix on covectorvalued grid     %
%                                                                         %
% Written by Jasper Kreeft - 2012                                         %
% Contact: j.j.kreeft@tudelft.nl                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Construct the divergence matrix for the internal points.

CD_h = kron(speye(N),spdiags([-ones(N+1,1) ones(N+1,1)],0:1,N+1,N+2));

CD_v = zeros(N*(N+1),(N+1)^2);
for i=1:N
    CD_v((1:N+1)+(i-1)*(N+1),:) = kron(speye(N+1),[zeros(1,i-1) -1 1 zeros(1,N-i) ]);
end
CD_v = sparse(CD_v);

CD_x = [ CD_h CD_v ];

CD = kron(speye(2),CD_x);